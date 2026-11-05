"""Are the two peaks made by the HHblits diversity filter (-diff 1000)? Read-only on the alignments.
1) uniclust30_unfiltered: rows and per-column coverage of a random sample of alignments.
2) OpenProteinSet pdb folders of the original analysis (three alignments per chain):
   rows of each alignment for 5,000 chains; AF2-style Neff and old Neff of each alignment for 600 of them.
Usage: python two_peaks_server.py <chains.txt> <outdir>"""
import random
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, "/path/to/server_data/neff_v6_af2style/scripts")
from fetch_af2style_neff import GAP, LUT, af2style, parse_a3m  # noqa: E402

STORE = Path("/path/to/server_data/neff_v6/a3m")
PDB = Path("/path/to/server_data/work/openproteinset/pdb")
FILES = ("bfd_uniclust_hits", "uniref90_hits", "mgnify_hits")
chains_file, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)


def coverage(seqs):
    """rows kept, query length, and column coverage: min, 10th percentile, median, 90th percentile, share of columns with >= 1,000 rows"""
    length = len(seqs[0])
    seqs = [s for s in seqs if len(s) == length]
    arr = np.frombuffer(b"".join(seqs), dtype=np.uint8).reshape(len(seqs), length)
    cov = (LUT[arr] < GAP).sum(0)
    return len(seqs), length, int(cov.min()), float(np.percentile(cov, 10)), float(np.median(cov)), float(np.percentile(cov, 90)), float((cov >= 1000).mean())


def one_cluster(cid):
    try:
        text = zlib.decompress((STORE / f"{cid % 1000:03d}" / f"{cid}.deflate").read_bytes(), -15)
        return (cid,) + coverage(parse_a3m(text))
    except Exception as exc:
        print("skip", cid, str(exc)[:60], flush=True)
        return None


def rows_of(chain):
    row = [chain]
    for name in FILES:
        p = PDB / chain / "a3m" / f"{name}.a3m"
        try:
            data = p.read_bytes()
            row.append(data.count(b"\n>") + data.startswith(b">"))
        except OSError:
            row.append(-1)
    return row


if __name__ == "__main__":
    rng = random.Random(20261104)
    # ---- 1
    ids = [int(l.split("\t", 1)[0]) for l in open("/path/to/server_data/neff_v6_af2style/af2style_neff_by_cluster.tsv").read().splitlines()[1:]]
    sample = rng.sample(ids, 8000)
    with Pool(12) as pool, open(out / "uniclust_sample_coverage.tsv", "w") as fh:
        fh.write("rep_id\trows\tquery_length\tcov_min\tcov_p10\tcov_median\tcov_p90\tshare_columns_ge_1000\n")
        for k, r in enumerate(pool.imap_unordered(one_cluster, sample, chunksize=16)):
            if r:
                fh.write("\t".join(map(str, r)) + "\n")
            if k % 2000 == 0:
                print("uniclust sample", k, flush=True)
    print("part 1 done", flush=True)

    # ---- 2, rows
    chains = open(chains_file).read().split()
    print("chains", len(chains), flush=True)
    rows_sample = rng.sample(chains, min(5000, len(chains)))
    with ThreadPoolExecutor(8) as pool, open(out / "pdb_rows.tsv", "w") as fh:
        fh.write("chain\t" + "\t".join(FILES) + "\n")
        for k, r in enumerate(pool.map(rows_of, rows_sample)):
            fh.write("\t".join(map(str, r)) + "\n")
            if k % 1000 == 0:
                print("rows", k, flush=True)
    print("part 2 done", flush=True)

    # ---- 3, depth of each alignment of the same chain
    with open(out / "pdb_af2style.tsv", "w") as fh:
        fh.write("chain\talignment\taf2style_neff\tneff\tn_rows\tquery_length\tcov_median\n")
        for k, chain in enumerate(rows_sample[:600]):
            for name in FILES:
                p = PDB / chain / "a3m" / f"{name}.a3m"
                try:
                    seqs = parse_a3m(p.read_bytes())
                    med, mean, whole, old, n, length, dev = af2style(seqs)
                    cov = coverage(seqs)
                    fh.write(f"{chain}\t{name}\t{med:.4f}\t{old:.4f}\t{n}\t{length}\t{cov[4]}\n")
                except Exception as exc:
                    print("skip", chain, name, str(exc)[:80], flush=True)
            fh.flush()
            if k % 50 == 0:
                print("depth", k, flush=True)
    print("ALL DONE", flush=True)
