# Data and code: PNAS manuscript 2026-25352

Datasets, Python scripts and statistics tables of the Supporting Information of PNAS manuscript 2026-25352: the datasets
of the SI (Datasets S1-S17), the Python scripts that built them, drew the SI figures and computed the SI statistics,
and the resulting statistics tables.

## Contents

```
datasets/     Datasets S1-S17, numbered as in the SI; datasets/README.txt describes every file and its columns,
              Dataset_S16_* and Dataset_S17_* have their own README.txt (source data of the main-text and SI figures)
scripts/      analysis/ (figures and statistics), datasets/ (building the datasets), checks/ (verifications);
              scripts/README.md lists every script
statistics/   statistics tables behind the SI (tab-separated)
```

## Large files

The 7 dataset files of 50 MB or more (about 1.5 GB) are in the Zenodo record https://doi.org/10.5281/zenodo.23925266, not here;
`datasets/LARGE_FILES_ON_ZENODO.md` lists them with their MD5 sums and where they belong.

## Environment

Python 3.10; package versions are in `requirements.txt`. File paths in the scripts are placeholders (`/path/to/...`)
for the folders of the original analysis; some scripts need modules of the original pipelines that are not included
(see `scripts/README.md`).

## Data sources

The datasets are derived from public resources, cited in the SI: AlphaFold Protein Structure Database (versions 4 and
6), RCSB Protein Data Bank, UniProtKB, OpenProteinSet, ESM Metagenomic Atlas, published GenEra gene ages, FireProtDB,
the mega-scale stability dataset, Meltome Atlas, MaveDB and mdCATH. Please also follow the terms of these resources.

## Licence

Code: MIT (`../LICENSE`). Data and tables in `datasets/` and `statistics/`: CC BY 4.0 (`../LICENSE-DATA.md`).

## Citation

Citation details will be added on publication.
