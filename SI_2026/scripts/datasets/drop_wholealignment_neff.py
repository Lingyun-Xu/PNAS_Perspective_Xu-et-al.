#!/usr/bin/env python3
"""Remove the whole-alignment Neff (the sum over the alignment of 1/n, n = sequences with >= 80% identity; columns
neff, neff_per_len, neff_uniclust30, log10_neff and their summaries) from the SI datasets (author request,
13 Nov 2026). The SI reports only Neff as defined in the AlphaFold2 paper (column AF2-style_Neff) and Neff/length
(AF2-style_Neff / query_length).

Run after the scripts that write these columns: build_si_datasets.py, add_af2style_neff_other_datasets.py,
figS8_published_ages_9species.py, table_s19_human.py and plot_neff_sets.py. The source data of the first version of
main-text Fig. 1 (Dataset S16, folder Fig1b_Fig1h-k) is left as it was. Scripts that analysed the whole-alignment
Neff (the checks of af2style_neff_stats.py and af2style_neff_length_stats.py, fig1_stats.py, fig1hk_v6.py,
jsd_grid.py, neff_two_peaks.py) need the datasets as they were before this change (repository history).

Every file is written to <out root> under the same relative path; every column that is kept is checked to be
unchanged (text files: byte-identical fields; parquet and xlsx: equal values).
Usage: python drop_wholealignment_neff.py <SI_2026/datasets> <out root>"""
import sys
from pathlib import Path

import openpyxl
import pandas as pd
import pyarrow.parquet as pq

DROP = {
    "Dataset_S10_alignment_depth/merged_neff.parquet": ["neff", "neff_per_len"],
    "Dataset_S13_homolog_class/panel_strata.parquet": ["neff"],
    "Dataset_S13_homolog_class/panel_bin_summary.csv": ["neff_median", "neff_q25", "neff_q75"],
    "Dataset_S14_published_gene_ages/genera_balanced_sample_ages.csv": ["neff"],
    "Dataset_S14_published_gene_ages/genera_balanced_sample_summary.csv": ["neff_median"],
    "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv": ["neff"],
    "Dataset_S14_published_gene_ages/suite_bin_summary.csv": ["neff_median"],
    "Dataset_S17_source_data_SI_figures/FigS19/neff_sets.tsv": ["neff", "log10_neff"],
    "Dataset_S17_source_data_SI_figures/FigS19/neff_sets_tests.tsv": [
        "n_with_neff", "median_neff", "share_neff_lt_10", "cliffs_delta_neff_vs_pdb3825", "p_neff_vs_pdb3825"],
    "Dataset_S17_source_data_SI_figures/FigS20/sel3800_genera_ages.csv": ["neff"],
    "Dataset_S17_source_data_SI_figures/FigS20/sel690_genera_ages.csv": ["neff"],
    "Dataset_S2_sel3800.xlsx": ["neff_uniclust30"],
    "Dataset_S3_sel690.xlsx": ["neff_uniclust30"],
}
README_OLD, README_NEW = "Neff and AF2-style Neff (OpenProteinSet", "AF2-style Neff (OpenProteinSet"


def text_file(src, dst, cols):
    sep = "\t" if src.suffix == ".tsv" else ","
    raw = src.read_bytes().decode("utf-8")
    assert '"' not in raw, f"{src}: quoted fields, not handled"
    nl = "\r\n" if "\r\n" in raw else "\n"
    end = raw.endswith(nl)
    lines = raw[:-len(nl)].split(nl) if end else raw.split(nl)
    head = lines[0].split(sep)
    drop = [head.index(c) for c in cols]
    keep = [i for i in range(len(head)) if i not in drop]
    out = []
    for k, line in enumerate(lines):
        f = line.split(sep)
        assert len(f) == len(head), (src, k, len(f), len(head))
        out.append(sep.join(f[i] for i in keep))
    dst.write_bytes((nl.join(out) + (nl if end else "")).encode("utf-8"))
    a = pd.read_csv(src, sep=sep, dtype=str, keep_default_na=False)
    b = pd.read_csv(dst, sep=sep, dtype=str, keep_default_na=False)
    assert list(b.columns) == [c for c in a.columns if c not in cols] and a.drop(columns=cols).equals(b)
    return len(a), len(b.columns)


def parquet_file(src, dst, cols):
    codec = pq.ParquetFile(src).metadata.row_group(0).column(0).compression.lower()
    a = pd.read_parquet(src)
    b = a.drop(columns=cols)
    b.to_parquet(dst, index=False, compression=codec)
    c = pd.read_parquet(dst)
    assert c.equals(b.reset_index(drop=True)) and list(c.columns) == [x for x in a.columns if x not in cols]
    return len(c), len(c.columns), codec


def xlsx_file(src, dst, cols):
    wb = openpyxl.load_workbook(src)
    ws = wb["data"]
    for col in cols:
        head = [c.value for c in ws[1]]
        ws.delete_cols(head.index(col) + 1)
    n = 0
    for row in wb["README"].iter_rows():
        for cell in row:
            if isinstance(cell.value, str) and README_OLD in cell.value:
                cell.value = cell.value.replace(README_OLD, README_NEW)
                n += 1
    assert n == 1, (src, "README cells changed", n)
    wb.save(dst)
    a = pd.read_excel(src, sheet_name="data")
    b = pd.read_excel(dst, sheet_name="data")
    assert a.drop(columns=cols).equals(b), src
    return len(b), len(b.columns)


def main():
    src_root, dst_root = Path(sys.argv[1]), Path(sys.argv[2])
    for rel, cols in DROP.items():
        src, dst = src_root / rel, dst_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.suffix == ".parquet":
            info = parquet_file(src, dst, cols)
        elif src.suffix == ".xlsx":
            info = xlsx_file(src, dst, cols)
        else:
            info = text_file(src, dst, cols)
        print(f"{rel}: dropped {cols}; rows, columns kept{', codec' if len(info) == 3 else ''}: {info}", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
