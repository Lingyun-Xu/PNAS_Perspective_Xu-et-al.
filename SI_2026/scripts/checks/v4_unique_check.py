import pandas as pd
d = pd.read_parquet("/path/to/project/SI_2026/datasets/AFDB_v4_extracted/AF2_48_model_oganisms.parquet", columns=["AAseq_one_letter_code", "organism_scientific"])
d["seq"] = d.AAseq_one_letter_code.str.lstrip(";").str.replace(r"\s+", "", regex=True)
for name, pat in (("E. coli", "Escherichia coli"), ("yeast", "Saccharomyces cerevisiae"), ("human", "Homo sapiens"), ("plants", "Arabidopsis|Glycine max|Oryza sativa|Zea mays")):
    s = d[d.organism_scientific.str.contains(pat)]
    print(f"{name:8s} models {len(s):>7,} | unique sequences {s.seq.nunique():>7,}")
