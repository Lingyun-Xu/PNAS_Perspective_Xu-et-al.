#!/usr/bin/env python3
"""
Fetch protein sequences from RCSB PDB, annotate with NIH taxonomy and superkingdom,
and write unique sequences to JSON.

Output JSON records contain:
  - sequence: unique protein sequence (one-letter code)
  - pdb_ids: list of PDB entry IDs sharing that sequence
  - source_organism: biological source species (gene source / natural source),
                     not the expression/production host organism
  - species: NCBI taxonomy group for source_organism (from NIH taxonomy lookup)
  - superkingdom: broad category assigned from species (Fig2_Pie_chart.ipynb logic)

Designed for long runs: retries, throttling, periodic checkpointing, resume.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote_plus, urlencode

import requests

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
NIH_TAXONOMY_URL = "https://www.ncbi.nlm.nih.gov/taxonomy/"

DEFAULT_CHECKPOINT = "pdb_fetch_checkpoint.json"
DEFAULT_OUTPUT = "unique_protein_sequences.json"
DEFAULT_SPECIES_CACHE = "nih_species_cache.json"
REPORT_DIV_NOT_FOUND = "Report div not found."

REPLACEMENT_DICT = {
    "MSE": "M",
    "CME": "C",
    "MLY": "K",
    "SNN": "K",
    "SEP": "S",
    "TPO": "T",
    "LLP": "K",
    "UNK": "X",
    "OCS": "C",
    "KCR": "K",
    "CSO": "C",
    "PCA": "Q",
    "HYP": "P",
    "FME": "M",
    "CSD": "C",
    "M3L": "K",
    "M2L": "K",
    "MLZ": "K",
    "KCX": "K",
    "SME": "M",
}

NUCLEOTIDE_PATTERN = re.compile(r"^[ACGTU]+$")
SPECIAL_NUCLEOTIDE_PATTERN = re.compile(r"\((DA|DC|DT|DG|1MG|PSU|PED|C49|5CM|GTP)\)")
NONCANONICAL_PATTERN = re.compile(r"\((.*?)\)")

# NIH species -> superkingdom mappings from Fig2_Pie_chart.ipynb (df_final_meltome loop order).
VIRUSES_SPECIES = frozenset(['species, viruses', 'isolate, viruses', 'viruses', 'virus', 'Report div not found.', 'serotype, viruses', '(Epstein-Barr virus(strain B95-8)), isolate, viruses'])

BACTERIA_SPECIES = frozenset(['(E. coli), species, enterobacteria', 'species, g-proteobacteria', 'species, a-proteobacteria', 'species, bacteria', 'species, high G+C Gram-positive bacteria', 'strain, enterobacteria', 'species, cyanobacteria', 'species, GNS bacteria', '(Bacteria stick insect), genus, walking sticks', 'enterobacteria', 'species, enterobacteria', 'strain, b-proteobacteria', 'species, proteobacteria', 'g-proteobacteria', 'strain, e-proteobacteria', 'strain, g-proteobacteria', 'strain, a-proteobacteria', 'strain, green sulfur bacteria', 'strain, bacteria', 'species, green sulfur bacteria', 'strain, high G+C Gram-positive bacteria', 'species, b-proteobacteria', 'strain, CFB group bacteria', 'genus, CFB group bacteria', 'species, CFB group bacteria', 'subspecies, high G+C Gram-positive bacteria', 'strain, fusobacteria', 'species, actinobacteria', 'species, fusobacteria', 'bacteria', 'strain, cyanobacteria', 'species, e-proteobacteria', 'strain, GNS bacteria', 'species, d-proteobacteria', 'strain, proteobacteria', 'genus, a-proteobacteria', 'isolate, bacteria', 'strain, actinobacteria', 'isolate, b-proteobacteria', 'isolate, cyanobacteria', 'forma specialis, enterobacteria', 'forma specialis, g-proteobacteria', 'isolate, CFB group bacteria', 'CFB group bacteria', 'subspecies, a-proteobacteria', 'genus, high G+C Gram-positive bacteria', 'a-proteobacteria', 'species, firmicutes', 'strain, firmicutes', 'genus, firmicutes', 'subspecies, firmicutes', 'firmicutes', 'isolate, firmicutes', 'species, aquificales', 'strain, aquificales', 'species, chlamydias', 'strain, chlamydias', 'species, planctomycetes', 'species, thermotogales', 'strain, thermotogales', 'species, spirochetes', 'strain, spirochetes', 'species, verrucomicrobia', 'strain, verrucomicrobia', 'strain, mycoplasmas', 'species, mycoplasmas', 'cyanobacteria', 'b-proteobacteria', 'subspecies, g-proteobacteria'])

FUNGI_SPECIES = frozenset(['(corn smut), species, smut fungi & allies', 'strain, budding yeasts & allies', 'species, ascomycete fungi', '(fission yeast), species, ascomycete fungi', 'species, smut fungi', 'species, mucoromycotan fungi', 'strain, ascomycete fungi', 'strain, basidiomycete fungi', 'species, basidiomycete fungi', 'species, blastocladiomycotan fungi', 'species, zoopagomycotan fungi', 'strain, mucoromycotan fungi', '(maitake), species, basidiomycete fungi', '(shiitake mushroom), species, basidiomycete fungi', 'strain, zoopagomycotan fungi', '(rose black spot fungus), species, ascomycete fungi', 'varietas, basidiomycete fungi', "(witches' butter), species, basidiomycete fungi", '(Monilia pod rot of cocoa), species, basidiomycete fungi', '(rice blast fungus), species, ascomycete fungi', '(peach leaf curl fungus), species, ascomycete fungi', 'strain, fungi', 'species, rust fungi', 'strain, blastocladiomycotan fungi', '(Karnal bunt), species, basidiomycete fungi', 'forma specialis, rust fungi', '(black scurf of potato), species, basidiomycete fungi', 'strain, rust fungi', 'varietas, ascomycete fungi', '(white Piedmont truffle), species, ascomycete fungi', '(whitish truffle), species, ascomycete fungi', 'species, fungi', 'fungi', "(brewer's yeast), species, budding yeasts", 'species, budding yeasts', 'strain, budding yeasts', 'species, chytrids & allies', 'strain, chytrids & allies', 'strain, glomeromycetes', 'species, glomeromycetes', 'species, monoblepharidomycetes', 'strain, powdery mildews', 'species, powdery mildews', 'species, microsporidians'])

INSECTA_SPECIES = frozenset(['(fruit fly), species, flies', 'species, flies', 'subspecies, flies', 'species, butterflies', '(Australian sheep blowfly), species, flies', '(redheaded pine sawfly), species, hymenopterans', '(emerald ash borer), species, beetles', '(red flour beetle), species, beetles', '(jewel wasp), species, wasps, ants & bees', 'species, beetles', '(mountain pine beetle), species, beetles', '(honey bee), species, bees', 'species, wasps, ants & bees', 'species, termites', '(Formosan subterranean termite), species, termites', '(migratory locust), species, grasshoppers', '(Asian citrus psyllid), species, psyllids', '(winter moth), species, moths', '(tobacco budworm), species, moths', '(domestic silkworm), species, moths', 'species, moths', 'species, bugs', '(small brown planthopper), species, bugs', '(pea aphid), species, aphids', '(rose-grain aphid), species, aphids', '(cowpea aphid), species, aphids', '(soybean aphid), species, aphids', 'species, ants', '(African malaria mosquito), species, mosquitos', '(yellow fever mosquito), species, mosquitos', '(American malaria mosquito), species, mosquitos', '(Asian malaria mosquito), species, mosquitos', '(western flower thrips), species, thrips', 'species, thrips', '(German cockroach), species, roaches', 'flies', '(oriental rat flea), species, fleas', 'species, spiders', '(fruit fly), species, flies', 'species, flies', 'subspecies, flies', 'species, butterflies', '(Australian sheep blowfly), species, flies', '(redheaded pine sawfly), species, hymenopterans', '(emerald ash borer), species, beetles', '(red flour beetle), species, beetles', '(jewel wasp), species, wasps, ants & bees', 'species, beetles', '(mountain pine beetle), species, beetles', '(honey bee), species, bees', 'species, wasps, ants & bees', 'species, termites', '(Formosan subterranean termite), species, termites', '(migratory locust), species, grasshoppers', '(Asian citrus psyllid), species, psyllids', '(winter moth), species, moths', '(tobacco budworm), species, moths', '(domestic silkworm), species, moths', 'species, moths', 'species, bugs', '(small brown planthopper), species, bugs', '(pea aphid), species, aphids', '(rose-grain aphid), species, aphids', '(cowpea aphid), species, aphids', '(soybean aphid), species, aphids', 'species, ants', '(African malaria mosquito), species, mosquitos', '(yellow fever mosquito), species, mosquitos', '(American malaria mosquito), species, mosquitos', '(Asian malaria mosquito), species, mosquitos', '(western flower thrips), species, thrips', 'species, thrips', '(German cockroach), species, roaches'])

OTHER_MAMMALIA_SPECIES = frozenset(['(gray short-tailed opossum), species, marsupials', '(rabbit), species, rabbits & hares', '(cattle), species, even-toed ungulates & whales', '(pig), species, even-toed ungulates & whales', '(goat), species, even-toed ungulates & whales', '(Wild Bactrian camel), species, even-toed ungulates & whales', "(Reeves' muntjac), species, even-toed ungulates & whales", '(domestic cattle), species, even-toed ungulates & whales', '(horse), species, odd-toed ungulates', '(house mouse), species, rodents', '(Norway rat), species, rodents', '(black rat), species, rodents', 'rodents', '(Chinese hamster), species, rodents', '(prairie deer mouse), subspecies, rodents', '(naked mole-rat), species, rodents', '(degu), species, rodents', '(golden hamster), species, rodents', "(Ord's kangaroo rat), species, rodents", '(Damara mole-rat), species, rodents', '(sperm whale), species, whales & dolphins', 'subspecies, whales & dolphins', '(Yangtze finless porpoise), subspecies, whales & dolphins', '(Indo-pacific humpbacked dolphin), species, whales & dolphins', '(Sumatran orangutan), species, primates', '(small-eared galago), species, primates', '(white-tufted-ear marmoset), species, primates', '(northern white-cheeked gibbon), species, primates', '(crab-eating macaque), species, primates', '(dog), subspecies, carnivores', 'subspecies, carnivores', '(puma), species, carnivores', '(meerkat), species, carnivores', '(cattle), species, even-toed ungulates & whales', '(sperm whale), species, whales & dolphins', '(pig), species, even-toed ungulates & whales', '(goat), species, even-toed ungulates & whales', '(Wild Bactrian camel), species, even-toed ungulates & whales', 'subspecies, whales & dolphins', '(Yangtze finless porpoise), subspecies, whales & dolphins', '(Indo-pacific humpbacked dolphin), species, whales & dolphins', "(Reeves' muntjac), species, even-toed ungulates & whales", '(rabbit), species, rabbits & hares', '(house mouse), species, rodents', '(Norway rat), species, rodents', '(black rat), species, rodents', '(Chinese hamster), species, rodents', '(prairie deer mouse), subspecies, rodents', '(naked mole-rat), species, rodents', '(degu), species, rodents', '(golden hamster), species, rodents', "(Ord's kangaroo rat), species, rodents", '(Damara mole-rat), species, rodents', '(black flying fox), species, bats', '(Chinese tree shrew), species, placentals', '(African savanna elephant), species, placentals', '(Cape golden mole), species, placentals', '(pale spear-nosed bat), species, bats', '(common wombat), species, marsupials', 'species, bats', '(platypus), species, monotremes', '(common wombat), species, marsupials', '(koala), species, marsupials'])

HOMO_SPECIES = frozenset(['human', '(human), species, primates', 'Homo Sapiens'])

NEMATODA_SPECIES = frozenset(['species, nematodes', '(nodule worm), species, nematodes', '(common roundworm), species, nematodes', '(American rat lungworm), species, nematodes', '(barber pole worm), species, nematodes', '(New World hookworm), species, nematodes', '(agent of lymphatic filariasis), species, nematodes', '(herring worm), species, nematodes', '(bovine lungworm), species, nematodes', '(Guinea worm), species, nematodes', '(pine wood nematode), species, nematodes', '(dog roundworm), species, nematodes', '(human whipworm), species, nematodes', '(pig roundworm), species, nematodes', '(African eye worm), species, nematodes', '(human pinworm), species, nematodes', 'species, nematodes', '(nodule worm), species, nematodes', '(common roundworm), species, nematodes', '(American rat lungworm), species, nematodes', '(barber pole worm), species, nematodes', '(New World hookworm), species, nematodes', '(agent of lymphatic filariasis), species, nematodes', '(herring worm), species, nematodes', '(bovine lungworm), species, nematodes', '(Guinea worm), species, nematodes', '(pine wood nematode), species, nematodes', '(dog roundworm), species, nematodes', '(human whipworm), species, nematodes', '(pig roundworm), species, nematodes', '(African eye worm), species, nematodes', 'nematodes'])

VIRIDIPLANTAE_SPECIES = frozenset(['(flower-of-an-hour), species, eudicots', '(winter squash), species, eudicots', '(snapdragon), species, eudicots', '(thale cress), species, eudicots', '(common sunflower), species, eudicots', '(garden lettuce), species, eudicots', '(rotklee), species, eudicots', 'genus, eudicots', '(soybean), species, eudicots', 'species, eudicots', '(rape), species, eudicots', '(spinach), species, eudicots', '(common jujube), species, eudicots', '(cassava), species, eudicots', '(tree cotton), species, eudicots', '(common tobacco), species, eudicots', '(mile-a-minute), species, eudicots', '(sea-island cotton), species, eudicots', '(spotted monkey flower), species, eudicots', '(peanut), species, eudicots', '(white lupine), species, eudicots', 'subspecies, eudicots', '(durian), species, eudicots', '(apricot), species, eudicots', '(tomato), species, eudicots', '(English walnut), species, eudicots', '(sweet wormwood), species, eudicots', '(pomegranate), species, eudicots', '(globe artichoke), varietas, eudicots', 'varietas, eudicots', '(barrel medic), species, eudicots', '(Rose-of-Sharon), species, eudicots', '(apple), species, eudicots', '(sesame), species, eudicots', '(narrow-leaved blue lupine), species, eudicots', '(muskmelon), species, eudicots', '(castor bean), species, eudicots', '(radish), species, eudicots', '(saltwater cress), species, eudicots', '(wild soybean), species, eudicots', '(pigeon pea), species, eudicots', '(wine grape), species, eudicots', '(biennial turnip rape), subspecies, eudicots', '(clementine), species, eudicots', '(rocotillo), species, eudicots', '(sweet cherry), species, eudicots', '(kidney bean), species, eudicots', '(black cottonwood), species, eudicots', '(table beet), subspecies, eudicots', '(China rose), species, eudicots', '(wood tobacco), species, eudicots', '(scarlet sage), species, eudicots', '(rubber tree), species, eudicots', '(white poplar), species, eudicots', '(coffee), species, eudicots', '(shan jing zi), species, eudicots', '(cucumber), species, eudicots', '(rose gum), species, eudicots', '(barley), species, monocots', '(pineapple), species, monocots', '(domesticated barley), subspecies, monocots', 'species, monocots', '(Japanese rice), monocots', '(maize), species, monocots', '(garden asparagus), species, monocots', '(African oil palm), species, monocots', '(stiff brome), species, monocots', '(Balbis banana), species, monocots', 'subspecies, monocots', '(Australian wild rice), species, monocots', '(weeping love grass), species, monocots', '(wild Malaysian banana), subspecies, monocots', '(proso millet), species, monocots', '(sorghum), species, monocots', '(date palm), species, monocots', 'varietas, monocots', '(African wild rice), species, monocots', '(malo sina), species, monocots', 'species, flowering plants', '(Rocky Mountain columbine), species, flowering plants', '(sacred lotus), species, flowering plants', '(opium poppy), species, flowering plants', 'forma, flowering plants', 'species, green plants', 'species, green algae', 'strain, green algae', 'forma, green algae', 'species, mosses', 'isolate, green algae', 'species, club-mosses', 'subspecies, liverworts', '(liverwort), species, liverworts', '(belladonna), species, eudicots'])

OTHER_VERTEBRATA_SPECIES = frozenset(['(Mozambique tilapia), species, bony fishes', '(Japanese quail), species, birds', '(guppy), species, bony fishes', '(Okinawa rail), species, birds', '(skipjack tuna), species, bony fishes', '(ocean pout), species, bony fishes', '(Antarctic eel pout), species, bony fishes', '(Japanese medaka), species, bony fishes', 'species, bony fishes', '(viviparous blenny), species, bony fishes', '(rainbow trout), species, bony fishes', '(ocean sunfish), species, bony fishes', '(zebrafish), species, bony fishes', '(large yellow croaker), species, bony fishes', '(Nile tilapia), species, bony fishes', '(spotted green pufferfish), species, bony fishes', '(spotted gar), species, bony fishes', '(mangrove rivulus), species, bony fishes', '(northern pike), species, bony fishes', '(Atlantic herring), species, bony fishes', '(sterlet), species, bony fishes', '(three-spined stickleback), species, bony fishes', '(Siamese fighting fish), species, bony fishes', '(javanese ricefish), species, bony fishes', '(southern platyfish), species, bony fishes', '(orbiculate cardinalfish), species, bony fishes', '(yellow perch), species, bony fishes', '(climbing perch), species, bony fishes', '(ballan wrasse), species, bony fishes', '(torafugu), species, bony fishes', '(striped catfish), species, bony fishes', '(Asian bonytongue), species, bony fishes', '(sheepshead minnow), species, bony fishes', '(huchen), species, bony fishes', '(goldfish), species, bony fishes', '(pinecone soldierfish), species, bony fishes', '(western mosquitofish), species, bony fishes', '(Mexican tetra), species, bony fishes', '(rohu), species, bony fishes', '(milkfish), species, bony fishes', '(Atlantic salmon), species, bony fishes', '(Amazon molly), species, bony fishes', '(electric eel), species, bony fishes', '(gilthead seabream), species, bony fishes', '(river trout), species, bony fishes', '(chicken), species, birds', '(scaled quail), species, birds', '(speckled mousebird), species, birds', '(chimney swift), species, birds', '(Downy woodpecker), species, birds', '(cuckoo roller), species, birds', '(white-throated tinamou), species, birds', '(Dalmatian pelican), species, birds', '(Chinese bamboo-partridge), species, birds', '(Kakapo), species, birds', '(rifleman), species, birds', '(common cuckoo), species, birds', '(Collared flycatcher), species, birds', 'subspecies, birds', '(northern leopard frog), species, frogs & toads', '(African clawed frog), species, frogs & toads', '(tropical clawed frog), species, frogs & toads', '(American bullfrog), species, frogs & toads', '(Chinese alligator), species, vertebrates', '(Gila monster), species, lizards & snakes', '(king cobra), species, snakes', '(tropical rattlesnake), subspecies, snakes', '(Common wall lizard), species, lizards & snakes', 'species, snakes', '(central bearded dragon), species, lizards & snakes', '(green anole), species, lizards & snakes', '(American alligator), species, vertebrates', '(Green sea turtle), species, turtles', "(Agassiz's desert tortoise), species, turtles", '(elephant shark), species, chimaeras', '(cloudy catshark), species, sharks & rays', '(brownbanded bambooshark), species, sharks & rays', '(coelacanth), species, coelacanths', '(Gabon caecilian), species, caecilians', 'species, caecilians', 'bony fishes', '(sea lamprey), species, lampreys'])

ARCHAEA_SPECIES = frozenset(['species, archaea', 'species, euryarchaeotes', 'strain, euryarchaeotes', 'species, crenarchaeotes', 'strain, crenarchaeotes'])

OTHERS_SPECIES = frozenset(['species, metagenomes', 'species, cellular slime molds', 'strain, cellular slime molds', '(North American leech), species, segmented worms', '(medicinal leech), species, segmented worms', 'species, segmented worms', '(common earthworm), species, segmented worms', 'species, mites & ticks', '(black-legged tick), species, mites & ticks', '(two-spotted spider mite), species, mites & ticks', '(European house dust mite), species, mites & ticks', '(scrub typhus mite), species, mites & ticks', 'subspecies, kinetoplastids', 'species, kinetoplastids', 'strain, kinetoplastids', 'species, myxozoans', '(Florida lancelet), species, lancelets', "(Belcher's lancelet), species, lancelets", '(Pacific oyster), species, bivalves', 'species, bivalves', '(Yesso scallop), species, bivalves', 'species, stony corals', '(cauliflower coral), species, stony corals', 'species, downy mildews', '(potato late blight agent), species, downy mildews', '(black shank of tobacco agent), species, downy mildews', 'strain, downy mildews', 'species, eukaryotes', 'strain, eukaryotes', 'isolate, eukaryotes', 'species, flatworms', '(bird nasal schistosome), species, flatworms', '(oriental liver fluke), species, flatworms', '(Southeast Asian liver fluke), species, flatworms', 'species, ciliates', 'strain, ciliates', '(Japanese sea cucumber), species, sea cucumbers', 'species, oomycetes', 'strain, oomycetes', 'species, animals', 'species, tardigrades', '(California two-spot octopus), species, cephalopods', 'species, springtails', 'species, sponges', 'species, brachiopods', '(bloodfluke planorb), species, gastropods', '(owl limpet), species, gastropods', 'species, gastropods', '(cloth-of-gold cone), species, gastropods', '(eastern emerald elysia), species, gastropods', 'species, choanoflagellates', '(purple sea urchin), species, sea urchins', '(starlet sea anemone), species, sea anemones', '(Australian red waratah sea anemone), species, sea anemones', '(Egyptian scorpion), species, scorpions', 'species, tunicates', '(vase tunicate), species, tunicates', '(common water flea), species, crustaceans', '(swimming crab), species, crustaceans', 'species, crustaceans', 'species, pelagophytes', 'species, centipedes', 'species, diatoms', 'strain, diatoms', 'species, haptophytes', 'species, dinoflagellates', 'species, apicomplexans', 'isolate, apicomplexans', 'strain, apicomplexans', 'species, brown algae', 'strain, cryptomonads', 'species, apicomplexans', 'strain, apicomplexans', 'strain, red algae', 'species, red algae', 'isolate, trichomonads', 'species, rotifers', 'species, placozoans', 'species, forams', 'species, hydrozoans', '(malaria parasite P. falciparum), species, apicomplexans'])

SYNTHETIC_SPECIES = frozenset(['synthetic', 'species, other sequences', 'species', 'synthesis'])

SUPERKINGDOM_RULES: tuple[tuple[str, frozenset[str]], ...] = (
    ("Bacteria", BACTERIA_SPECIES),
    ("Fungi", FUNGI_SPECIES),
    ("Others", OTHERS_SPECIES),
    ("Other Mammalia", OTHER_MAMMALIA_SPECIES),
    ("Homo", HOMO_SPECIES),
    ("Nematoda", NEMATODA_SPECIES),
    ("Insecta", INSECTA_SPECIES),
    ("Viridiplantae", VIRIDIPLANTAE_SPECIES),
    ("Other Vertebrata", OTHER_VERTEBRATA_SPECIES),
    ("Archaea", ARCHAEA_SPECIES),
    ("Viruses", VIRUSES_SPECIES),
    ("Synthetic", SYNTHETIC_SPECIES),
)


def species_to_superkingdom(species: str | None) -> str | None:
    """Map NIH species text to superkingdom label (Fig2_Pie_chart.ipynb logic)."""
    if species is None:
        return None
    value = str(species).strip()
    if not value:
        return None
    for label, members in SUPERKINGDOM_RULES:
        if value in members:
            return label
    # NIH often returns "(common name), rank, group" — match taxonomy suffix too.
    for label, members in SUPERKINGDOM_RULES:
        for member in members:
            if "," in member and len(member) >= 12 and value.endswith(member):
                return label
    return None


def fetch_protein_pdb_ids(batch_size: int = 1000, max_ids: int | None = None) -> list[str]:
    """Fetch PDB entry IDs where polymer entity type is Protein."""
    all_pdb_ids: list[str] = []
    start = 0

    while True:
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": "Protein",
                },
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": batch_size},
                "results_content_type": ["experimental"],
            },
        }

        data = post_json(SEARCH_URL, query)
        pdb_ids = [item["identifier"] for item in data.get("result_set", [])]
        all_pdb_ids.extend(pdb_ids)

        if max_ids is not None and len(all_pdb_ids) >= max_ids:
            return all_pdb_ids[:max_ids]

        if len(pdb_ids) < batch_size:
            break

        start += batch_size
        time.sleep(0.2)

    return all_pdb_ids


def post_json(
    url: str,
    payload: dict[str, Any],
    *,
    retries: int = 5,
    timeout: int = 30,
    base_sleep: float = 1.0,
) -> dict[str, Any]:
    """POST JSON with exponential backoff on transient failures."""
    headers = {"Content-Type": "application/json"}
    last_error: Exception | None = None

    for attempt in range(retries):
        try:
            response = requests.post(url, json=payload, headers=headers, timeout=timeout)
            if response.status_code in {429, 500, 502, 503, 504}:
                raise requests.HTTPError(
                    f"HTTP {response.status_code}: {response.text[:200]}", response=response
                )
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt < retries - 1:
                sleep_for = base_sleep * (2**attempt)
                print(f"Request failed ({exc}); retrying in {sleep_for:.1f}s...", file=sys.stderr)
                time.sleep(sleep_for)
            else:
                break

    raise RuntimeError(f"Request failed after {retries} attempts: {last_error}") from last_error


def fetch_pdb_data(query: str, variables: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """Execute GraphQL query and return the `data` block."""
    payload: dict[str, Any] = {"query": query}
    if variables:
        payload["variables"] = variables

    try:
        result = post_json(GRAPHQL_URL, payload)
    except RuntimeError as exc:
        print(f"GraphQL request failed: {exc}", file=sys.stderr)
        return None

    if "errors" in result:
        print(f"GraphQL errors: {result['errors']}", file=sys.stderr)
        return None

    return result.get("data")


def get_source_organism(entity: dict[str, Any]) -> str | None:
    """
    Return biological source species, not expression/production host.

    Priority (matches RCSB polymer_entities API):
      1) rcsb_entity_source_organism.scientific_name
      2) rcsb_entity_source_organism.ncbi_scientific_name
      3) entity_src_gen.pdbx_gene_src_scientific_name
      4) entity_src_nat.pdbx_organism_scientific
      5) uniprots.rcsb_uniprot_protein.source_organism.scientific_name
    """
    invalid = {"", "?", "unknown", "Unknown", "UNKNOWN"}

    for org in entity.get("rcsb_entity_source_organism") or []:
        for field in ("scientific_name", "ncbi_scientific_name"):
            name = (org or {}).get(field)
            if name and name.strip() not in invalid:
                return name.strip()

    for src in entity.get("entity_src_gen") or []:
        name = (src or {}).get("pdbx_gene_src_scientific_name")
        if name and name.strip() not in invalid:
            return name.strip()

    for src in entity.get("entity_src_nat") or []:
        name = (src or {}).get("pdbx_organism_scientific")
        if name and name.strip() not in invalid:
            return name.strip()

    for uni in entity.get("uniprots") or []:
        protein = (uni or {}).get("rcsb_uniprot_protein") or {}
        name = ((protein.get("source_organism") or {}).get("scientific_name"))
        if name and name.strip() not in invalid:
            return name.strip()

    return None


def is_dna_rna(sequence: str) -> bool:
    if not sequence:
        return True
    return bool(NUCLEOTIDE_PATTERN.match(sequence)) or bool(
        SPECIAL_NUCLEOTIDE_PATTERN.search(sequence)
    )


def replace_non_canonical(sequence: str) -> str:
    for match in NONCANONICAL_PATTERN.findall(sequence):
        replacement = REPLACEMENT_DICT.get(match)
        if replacement is not None:
            sequence = sequence.replace(f"({match})", replacement)
    return sequence


def get_entity_sequence(entity: dict[str, Any]) -> str | None:
    """Return protein one-letter sequence from entity_poly."""
    poly = entity.get("entity_poly") or {}
    return poly.get("pdbx_seq_one_letter_code") or poly.get("pdbx_seq_one_letter_code_can")


def is_protein_entity(entity: dict[str, Any]) -> bool:
    poly = entity.get("entity_poly") or {}
    polymer_type = poly.get("rcsb_entity_polymer_type")
    if polymer_type:
        return polymer_type == "Protein"
    seq_type = poly.get("type")
    return seq_type == "polypeptide(L)" if seq_type else True


def process_pdb_id(pdb_id: str) -> list[dict[str, Any]]:
    """Fetch protein entities for one PDB entry."""
    entity_query = f"""
    {{
      entry(entry_id: "{pdb_id}") {{
        polymer_entities {{
          rcsb_id
          entity_poly {{
            pdbx_seq_one_letter_code
            pdbx_seq_one_letter_code_can
            rcsb_entity_polymer_type
            type
          }}
          rcsb_entity_source_organism {{
            scientific_name
            ncbi_scientific_name
            source_type
          }}
          entity_src_gen {{
            pdbx_gene_src_scientific_name
            pdbx_host_org_scientific_name
          }}
          entity_src_nat {{
            pdbx_organism_scientific
          }}
          uniprots {{
            rcsb_uniprot_protein {{
              source_organism {{
                scientific_name
              }}
            }}
          }}
        }}
      }}
    }}
    """

    entry_data = fetch_pdb_data(entity_query)
    if not entry_data or not entry_data.get("entry"):
        return []

    entities = entry_data["entry"].get("polymer_entities") or []
    rows: list[dict[str, Any]] = []

    for entity in entities:
        if not is_protein_entity(entity):
            continue

        sequence = get_entity_sequence(entity)
        if not sequence:
            continue

        if is_dna_rna(sequence):
            continue

        cleaned_sequence = replace_non_canonical(sequence)
        if "(" in cleaned_sequence or ")" in cleaned_sequence:
            continue

        rows.append(
            {
                "pdb_id": pdb_id.upper(),
                "entity_id": entity.get("rcsb_id"),
                "sequence": cleaned_sequence,
                "source_organism": get_source_organism(entity),
            }
        )

    return rows


def load_checkpoint(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"processed_pdb_ids": [], "records_by_sequence": {}}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_checkpoint(path: Path, processed_pdb_ids: list[str], records_by_sequence: dict[str, Any]) -> None:
    payload = {
        "processed_pdb_ids": processed_pdb_ids,
        "records_by_sequence": records_by_sequence,
    }
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    tmp_path.replace(path)


def merge_record(
    records_by_sequence: dict[str, dict[str, Any]],
    row: dict[str, Any],
) -> None:
    sequence = row["sequence"]
    pdb_id = row["pdb_id"]
    source_organism = row.get("source_organism")

    if sequence not in records_by_sequence:
        records_by_sequence[sequence] = {
            "sequence": sequence,
            "pdb_ids": [pdb_id],
            "source_organism": source_organism,
        }
        return

    record = records_by_sequence[sequence]
    if pdb_id not in record["pdb_ids"]:
        record["pdb_ids"].append(pdb_id)

    if record.get("source_organism") is None and source_organism is not None:
        record["source_organism"] = source_organism


def format_output_record(record: dict[str, Any]) -> dict[str, Any]:
    """Ensure each JSON record has the required fields in a stable order."""
    return {
        "sequence": record.get("sequence"),
        "pdb_ids": sorted(set(record.get("pdb_ids") or [])),
        "source_organism": record.get("source_organism"),
        "species": record.get("species"),
        "superkingdom": record.get("superkingdom"),
    }


def records_to_output(records_by_sequence: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    output = [format_output_record(record) for record in records_by_sequence.values()]
    output.sort(key=lambda item: item["sequence"])
    return output


def write_output_json(path: Path, records_by_sequence: dict[str, dict[str, Any]]) -> None:
    output = records_to_output(records_by_sequence)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)


def fetch_species_from_nih(organ_name: str, *, timeout: int = 20) -> str:
    """
    Query NCBI taxonomy for an organism name and return the NIH description text.

    Mirrors Fig2_Pie_chart.ipynb `find_species`: parses div.rprt > div.supp > p.desc.
    Returns REPORT_DIV_NOT_FOUND when the report block is missing.
    """
    try:
        from bs4 import BeautifulSoup
    except ImportError as exc:
        raise ImportError(
            "beautifulsoup4 is required for NIH species annotation. "
            "Install with: pip install beautifulsoup4"
        ) from exc

    params = urlencode({"term": organ_name.strip()}, quote_via=quote_plus)
    url = f"{NIH_TAXONOMY_URL}?{params}"

    try:
        response = requests.get(
            url,
            headers={"User-Agent": "Mozilla/5.0"},
            timeout=timeout,
        )
        response.raise_for_status()
        html = response.text
    except requests.RequestException as exc:
        print(f"NIH request failed for {organ_name}: {exc}", file=sys.stderr)
        return f"Error fetching page: {exc}"

    soup = BeautifulSoup(html, "html.parser")
    try:
        rprt_div = soup.find("div", class_="rprt")
        if not rprt_div:
            return REPORT_DIV_NOT_FOUND

        supp_div = rprt_div.find("div", class_="supp")
        if not supp_div:
            return "Supplementary div not found."

        desc_div = supp_div.find("p", class_="desc")
        if not desc_div:
            return "Description paragraph not found."

        return desc_div.get_text(strip=True)
    except Exception as exc:
        print(f"NIH parse failed for {organ_name}: {exc}", file=sys.stderr)
        return f"Error processing page: {exc}"


def load_species_cache(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return {str(key): str(value) for key, value in data.items()}


def save_species_cache(path: Path, cache: dict[str, str]) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
    tmp_path.replace(path)


def unique_source_organisms(records: list[dict[str, Any]]) -> list[str]:
    organisms: list[str] = []
    seen: set[str] = set()
    for record in records:
        organism = record.get("source_organism")
        if not organism or not str(organism).strip():
            continue
        name = str(organism).strip()
        if name not in seen:
            seen.add(name)
            organisms.append(name)
    return organisms


def build_organism_species_map(
    organisms: list[str],
    *,
    cache_path: Path,
    sleep_seconds: float = 0.35,
) -> dict[str, str]:
    """Look up NIH taxonomy species text for each unique source organism."""
    cache = load_species_cache(cache_path)
    total = len(organisms)

    for index, organism in enumerate(organisms, start=1):
        if organism in cache:
            continue

        print(f"[{index}/{total}] NIH lookup: {organism}...", file=sys.stderr)
        cache[organism] = fetch_species_from_nih(organism)
        save_species_cache(cache_path, cache)
        time.sleep(sleep_seconds)

    return cache


def annotate_records_with_species(
    records: list[dict[str, Any]],
    *,
    cache_path: Path,
    sleep_seconds: float = 0.35,
) -> list[dict[str, Any]]:
    """Add a `species` field to each record from NIH taxonomy lookups."""
    organisms = unique_source_organisms(records)
    if not organisms:
        print("No source organisms to annotate.", file=sys.stderr)
        for record in records:
            record["species"] = None
        return records

    print(f"Looking up NIH species for {len(organisms)} unique organisms...", file=sys.stderr)
    organism_species_map = build_organism_species_map(
        organisms,
        cache_path=cache_path,
        sleep_seconds=sleep_seconds,
    )

    for record in records:
        organism = record.get("source_organism")
        if organism and str(organism).strip():
            record["species"] = organism_species_map.get(str(organism).strip())
        else:
            record["species"] = None

    return records


def annotate_records_with_superkingdom(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Add a `superkingdom` field to each record from its NIH `species` text."""
    mapped = 0
    unmapped = 0

    for record in records:
        superkingdom = species_to_superkingdom(record.get("species"))
        record["superkingdom"] = superkingdom
        if superkingdom is None:
            unmapped += 1
        else:
            mapped += 1

    print(
        f"Assigned superkingdom to {mapped} records ({unmapped} unmapped).",
        file=sys.stderr,
    )
    return records


def annotate_records_with_taxonomy(
    records: list[dict[str, Any]],
    *,
    cache_path: Path,
    sleep_seconds: float = 0.35,
    annotate_superkingdom: bool = True,
) -> list[dict[str, Any]]:
    """Add NIH species and optional superkingdom annotations."""
    records = annotate_records_with_species(
        records,
        cache_path=cache_path,
        sleep_seconds=sleep_seconds,
    )
    if annotate_superkingdom:
        records = annotate_records_with_superkingdom(records)
    return records


def load_output_records(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON list in {path}")
    return data


def run_fetch(
    pdb_ids: list[str],
    *,
    checkpoint_path: Path,
    output_path: Path,
    sleep_seconds: float = 0.25,
    checkpoint_every: int = 25,
) -> list[dict[str, Any]]:
    checkpoint = load_checkpoint(checkpoint_path)
    processed = set(checkpoint.get("processed_pdb_ids", []))
    records_by_sequence: dict[str, dict[str, Any]] = checkpoint.get("records_by_sequence", {})

    total = len(pdb_ids)
    processed_count = 0
    failed_ids: list[str] = []

    for index, pdb_id in enumerate(pdb_ids, start=1):
        if pdb_id in processed:
            continue

        print(f"[{index}/{total}] Fetching {pdb_id}...", file=sys.stderr)

        try:
            rows = process_pdb_id(pdb_id)
            for row in rows:
                merge_record(records_by_sequence, row)
            processed.add(pdb_id)
            processed_count += 1
        except Exception as exc:
            print(f"Failed on {pdb_id}: {exc}", file=sys.stderr)
            failed_ids.append(pdb_id)

        if processed_count > 0 and processed_count % checkpoint_every == 0:
            save_checkpoint(checkpoint_path, sorted(processed), records_by_sequence)
            write_output_json(output_path, records_by_sequence)
            print(
                f"Checkpoint saved: {len(processed)} PDB IDs, "
                f"{len(records_by_sequence)} unique sequences",
                file=sys.stderr,
            )

        time.sleep(sleep_seconds)

    save_checkpoint(checkpoint_path, sorted(processed), records_by_sequence)
    output = records_to_output(records_by_sequence)
    write_output_json(output_path, records_by_sequence)

    if failed_ids:
        failed_path = output_path.with_name(output_path.stem + "_failed_ids.json")
        with failed_path.open("w", encoding="utf-8") as handle:
            json.dump(sorted(set(failed_ids)), handle, indent=2)
        print(f"Wrote failed PDB IDs to {failed_path}", file=sys.stderr)

    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch unique protein sequences and source organisms from RCSB PDB."
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Output JSON path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--checkpoint",
        default=DEFAULT_CHECKPOINT,
        help=f"Checkpoint JSON path for resume (default: {DEFAULT_CHECKPOINT})",
    )
    parser.add_argument(
        "--pdb-id-list",
        default=None,
        help="Optional JSON file containing a list of PDB IDs to process",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=None,
        help="Limit number of PDB entries fetched/processed (useful for testing)",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.25,
        help="Delay between PDB entry requests (default: 0.25)",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="Save checkpoint every N newly processed entries (default: 25)",
    )
    parser.add_argument(
        "--skip-fetch-ids",
        action="store_true",
        help="Require --pdb-id-list and do not query the search API",
    )
    parser.add_argument(
        "--skip-species-annotation",
        action="store_true",
        help="Skip NIH taxonomy species lookup after PDB fetch",
    )
    parser.add_argument(
        "--skip-superkingdom-annotation",
        action="store_true",
        help="Skip superkingdom assignment from NIH species text",
    )
    parser.add_argument(
        "--species-only",
        action="store_true",
        help="Only annotate an existing output JSON with species/superkingdom (no PDB fetch)",
    )
    parser.add_argument(
        "--input",
        default=None,
        help="Input JSON for --species-only (defaults to --output path)",
    )
    parser.add_argument(
        "--species-cache",
        default=DEFAULT_SPECIES_CACHE,
        help=f"Cache file for NIH organism->species lookups (default: {DEFAULT_SPECIES_CACHE})",
    )
    parser.add_argument(
        "--species-sleep-seconds",
        type=float,
        default=0.35,
        help="Delay between NIH taxonomy requests (default: 0.35)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_path = Path(args.output)
    checkpoint_path = Path(args.checkpoint)
    species_cache_path = Path(args.species_cache)
    annotate_species = not args.skip_species_annotation
    annotate_superkingdom = not args.skip_superkingdom_annotation

    if args.species_only:
        input_path = Path(args.input) if args.input else output_path
        if not input_path.exists():
            print(f"Error: input file not found: {input_path}", file=sys.stderr)
            return 2

        records = load_output_records(input_path)
        if annotate_species:
            records = annotate_records_with_species(
                records,
                cache_path=species_cache_path,
                sleep_seconds=args.species_sleep_seconds,
            )
        elif annotate_superkingdom and not any("species" in record for record in records):
            print(
                "Error: records lack a species column; run species annotation first "
                "or omit --skip-species-annotation.",
                file=sys.stderr,
            )
            return 2

        if annotate_superkingdom:
            records = annotate_records_with_superkingdom(records)

        records = [format_output_record(record) for record in records]
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(records, handle, indent=2)
        print(f"Wrote {len(records)} annotated records to {output_path}", file=sys.stderr)
        return 0

    if args.pdb_id_list:
        with Path(args.pdb_id_list).open("r", encoding="utf-8") as handle:
            pdb_ids = json.load(handle)
        pdb_ids = [str(item).upper() for item in pdb_ids if item]
    elif args.skip_fetch_ids:
        print("Error: --skip-fetch-ids requires --pdb-id-list", file=sys.stderr)
        return 2
    else:
        print("Fetching protein PDB IDs from RCSB search API...", file=sys.stderr)
        pdb_ids = fetch_protein_pdb_ids(max_ids=args.max_entries)

    if args.max_entries is not None and args.pdb_id_list is None:
        pdb_ids = pdb_ids[: args.max_entries]

    if not pdb_ids:
        print("No PDB IDs to process.", file=sys.stderr)
        return 1

    print(f"Processing {len(pdb_ids)} PDB entries...", file=sys.stderr)
    output = run_fetch(
        pdb_ids,
        checkpoint_path=checkpoint_path,
        output_path=output_path,
        sleep_seconds=args.sleep_seconds,
        checkpoint_every=args.checkpoint_every,
    )

    print(f"Wrote {len(output)} unique sequences to {output_path}", file=sys.stderr)

    if annotate_species:
        output = annotate_records_with_species(
            output,
            cache_path=species_cache_path,
            sleep_seconds=args.species_sleep_seconds,
        )
        print(f"Added NIH species annotations to {output_path}", file=sys.stderr)

    if annotate_superkingdom:
        output = annotate_records_with_superkingdom(output)
        print(f"Added superkingdom annotations to {output_path}", file=sys.stderr)

    if annotate_species or annotate_superkingdom:
        output = [format_output_record(record) for record in output]
        with output_path.open("w", encoding="utf-8") as handle:
            json.dump(output, handle, indent=2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

