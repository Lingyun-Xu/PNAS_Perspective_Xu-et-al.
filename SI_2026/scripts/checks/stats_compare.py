import pandas as pd, numpy as np
old_dir = "/path/to/project/SI_2026/statistics"
new_dir = "/path/to/scratch/stats_new"
for name, keys, vals in (("dtemplate_group_tests.tsv", ["measure", "comparison"], ["n", "statistic", "effect", "p_holm", "median_plants", "median_other"]),
                         ("dtemplate_among_plant_dunn.tsv", ["measure", "a", "b"], ["dunn_z", "cliffs_delta_a_vs_b", "p_holm"]),
                         ("dtemplate_correlation_tests.tsv", ["measure", "subset"], ["n", "rho", "p_holm", "statistic"])):
    o = pd.read_csv(f"{old_dir}/{name}", sep="\t"); n = pd.read_csv(f"{new_dir}/{name}", sep="\t")
    m = o.merge(n, on=keys, suffixes=("_old", "_new"))
    worst = 0.0
    for v in vals:
        if f"{v}_old" in m:
            a, b = m[f"{v}_old"].astype(float), m[f"{v}_new"].astype(float)
            ok = (a.isna() & b.isna()) | np.isclose(a, b, rtol=1e-6, atol=1e-12)
            worst = max(worst, float((~ok).sum()))
    print(f"{name}: old rows {len(o)}, matched {len(m)}, new rows {len(n)} (added {len(n) - len(m)}); cells differing: {int(worst)}")
