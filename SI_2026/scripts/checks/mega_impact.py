"""READ-ONLY estimate: how the Fig. 2c thermostability composition changes if all 1,841,285 mega-scale rows are counted
(the plotted table kept the first 1,048,575, the Excel row limit). Group of each mega-scale entry name taken from the
rows that were mapped (firemega_complete_species.csv); names only in the dropped rows stay unmapped."""
import pandas as pd

H = "/path/to/home"
full = pd.read_csv(f"{H}/megadata_thermalstability/mega_entry_names.csv", low_memory=False)
fm = pd.read_csv(f"{H}/GEO_dataset/firemega_complete_species.csv", low_memory=False)
me = pd.read_csv(f"{H}/GEO_dataset/meltome_altas_with_organisms_species_superkindoms.csv", low_memory=False)
print("mega_entry_names columns", list(full.columns), len(full), "| firemega columns", list(fm.columns), len(fm))
name_col = full.columns[0]
kept = full.iloc[:1048575]
dropped = full.iloc[1048575:]
print(f"kept rows {len(kept):,} ({kept[name_col].nunique():,} names) | dropped rows {len(dropped):,} ({dropped[name_col].nunique():,} names, "
      f"{(~dropped[name_col].isin(set(kept[name_col]))).sum():,} rows whose name is not among the kept rows)")
mega_rows = fm.iloc[-1048575:] if len(fm) > 1048575 else fm   # FireProtDB rows come first (53,445)
key = [c for c in fm.columns if c.lower() in ("name", "names")]
if key:
    g = dict(zip(mega_rows[key[0]], mega_rows["superkingdom"]))
    dropped_groups = dropped[name_col].map(g)
    print("dropped rows with a known group:", int(dropped_groups.notna().sum()))
    plotted = pd.concat([fm.superkingdom, me.superkingdom]).replace({"Nematoda": "Others"})
    allrows = pd.concat([plotted, dropped_groups.replace({"Nematoda": "Others"})])
    tot_p, tot_a = len(plotted), int(allrows.notna().sum())
    for grp in ["Viridiplantae", "Bacteria", "Synthetic", "Homo"]:
        print(f"  {grp:14s} plotted {100 * (plotted == grp).sum() / tot_p:5.1f}%  with all mapped mega rows {100 * (allrows == grp).sum() / tot_a:5.1f}%")
