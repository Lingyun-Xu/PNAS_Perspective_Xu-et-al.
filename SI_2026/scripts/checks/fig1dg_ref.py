import numpy as np, torch
raw = "/path/to/scratch/s16_raw"
d12 = torch.load("/path/to/home/AF2_clustering/min_pairwise_dist_Plants_021125.pt", map_location="cpu", weights_only=False).numpy().ravel()
for name in ("fig3_pdb_all", "fig3_pdb_excl"):
    r0 = torch.load(f"{raw}/{name}/round_0/min_pairwise_dist.pt", map_location="cpu", weights_only=False).numpy().ravel()
    same = len(r0) == len(d12) and np.allclose(r0, d12)
    print(f"{name} round 0 vs Fig. 1d-g plant tensor: identical={same}; medians {np.median(r0):.3f} vs {np.median(d12):.3f}; zeros {int((r0 < 1e-6).sum())} vs {int((d12 < 1e-6).sum())}")
