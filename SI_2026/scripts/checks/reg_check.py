"""AIUPred regressions of Dataset S12 (all proteins): coefficients for the SI table, and the counts behind them."""
import pandas as pd

D = "/path/to/project/SI_2026/datasets/Dataset_S12_disorder/"
r = pd.read_csv(D + "regressions.csv")
r = r[(r.disorder_source == "aiupred") & (r.data == "all proteins")].copy()
r["term"] = r.term.str.replace(r"C\(superkingdom, Treatment\('Viridiplantae'\)\)\[T\.(\w+)\]", r"group:\1", regex=True)
pd.set_option("display.width", 200)
print(r[["model", "term", "coef", "se_cluster_proteome", "p", "n", "r2"]].round(4).to_string(index=False))

m = pd.read_parquet(D + "merged_per_protein.parquet",
                    columns=["proteome", "is_plant", "msa_raw_depth", "plddt_ordered", "superkingdom"])
use = m.dropna(subset=["msa_raw_depth"])
print("n with depth", len(use), "| proteomes", use.proteome.nunique(),
      "| plant proteomes", use.loc[use.is_plant == 1, "proteome"].nunique(),
      "| with ordered pLDDT", int(use.plddt_ordered.notna().sum()))
print("groups:", use.superkingdom.value_counts().to_dict())
print("median log10 depth, plants:", round(float(use.loc[use.is_plant == 1, "msa_raw_depth"].clip(lower=1).pipe(
    lambda s: s.apply(lambda v: __import__('math').log10(v))).median()), 3))
