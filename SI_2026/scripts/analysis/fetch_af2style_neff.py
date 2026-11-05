"""AF2-style Neff of the OpenProteinSet uniclust30 alignments used for Neff (fetch_neff.py), as defined in Jumper et al.
2021 (Methods, test set of recent PDB sequences; Fig. 5 legend): per-position effective sequence counts, median across
positions.

The alignments are read as in fetch_neff.py: a3m lowercase insertions removed, 20 amino acids plus gap, any other letter
read as a gap, rows whose length differs from the query left out. Then, on the same rows:
  identity of two rows = identical residues / columns that are non-gap in either row
  n_i = number of rows with identity >= 0.8 to row i, itself included; weight w_i = 1 / n_i
  Neff at a column = sum of w_i over the rows with a residue in that column
  AF2-style Neff = median of the column values over all columns of the alignment (the query positions)

Also written, from the same counts: af2style_mean (mean over columns), neff_nongap_id (sum of w_i, the whole-row value
with this identity) and neff_recomputed (the fetch_neff.py definition: identity = columns in the same state / all
columns, a shared gap counting as a match, sum of 1/n_i), which checks the run against neff_by_cluster.tsv.

Raw deflate members are kept under --store (zlib.decompress(data, -15) restores the a3m). --also-read names a second
store that is only read, such as the October copies of fetch_neff.py; an alignment that it does not have or cannot
return is taken from --store or downloaded into it.
Resumable: ids already in --out are skipped.

Usage: python fetch_af2style_neff.py --ids ids.txt --out af2style_neff_by_cluster.tsv --store DIR [--also-read DIR]
"""
import argparse
import os
import queue
import string
import struct
import threading
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import torch

DIR = "/path/to/server_data/neff_v6/"
ALPHABET = "ARNDCQEGHILKMFPSTWYV-"
GAP = 20
LUT = np.full(256, GAP, dtype=np.int64)
for k, a in enumerate(ALPHABET):
    LUT[ord(a)] = k
RM_LC = bytes.maketrans(b"", b"")
LOWER = string.ascii_lowercase.encode()
COL_BLOCK = 2048
ROW_BLOCK = 4096


def fetch(url, a, b, tries=6):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception:
            if t == tries - 1:
                raise
            time.sleep(2 ** t)


def parse_a3m(text: bytes):
    seqs, buf = [], []
    for line in text.split(b"\n"):
        if not line or line.startswith(b"#"):
            continue
        if line.startswith(b">"):
            if buf:
                seqs.append(b"".join(buf))
                buf = []
            continue
        buf.append(line.translate(RM_LC, LOWER))
    if buf:
        seqs.append(b"".join(buf))
    return seqs


@torch.no_grad()
def _counts(idx, length, seq_id, device, row_block, col_block):
    """n_i of every row under both identities. fp16 matmuls are exact for counts over column blocks of <= 2048"""
    n = idx.shape[0]
    dtype = torch.float16 if device == "cuda" else torch.float32
    same_blocks, res_blocks = [], []
    residues = torch.zeros(n, dtype=torch.float32, device=device)
    for c0 in range(0, length, col_block):
        part = idx[:, c0:c0 + col_block]
        same_blocks.append(torch.nn.functional.one_hot(part, GAP + 1)[..., :GAP].to(dtype).reshape(n, -1))
        res = (part < GAP)
        residues += res.sum(1).float()
        res_blocks.append(res.to(dtype))
    n_af2 = torch.zeros(n, dtype=torch.float32, device=device)
    n_old = torch.zeros(n, dtype=torch.float32, device=device)
    for r0 in range(0, n, row_block):
        r1 = min(n, r0 + row_block)
        same = torch.zeros(r1 - r0, n, dtype=torch.float32, device=device)
        both = torch.zeros(r1 - r0, n, dtype=torch.float32, device=device)
        for oh, res in zip(same_blocks, res_blocks):
            same += (oh[r0:r1] @ oh.T).float()
            both += (res[r0:r1] @ res.T).float()
        either = residues[r0:r1, None] + residues[None, :] - both
        n_af2[r0:r1] = (same >= seq_id * either - 1e-3).sum(1).float()
        n_old[r0:r1] = (same + (length - either) >= seq_id * length - 1e-3).sum(1).float()
    return n_af2.double().cpu().numpy(), n_old.double().cpu().numpy()


def af2style(seqs, seq_id=0.8):
    """Returns (af2style_neff, af2style_mean, neff_nongap_id, neff_recomputed, n_rows, query_length, device_used)."""
    length = len(seqs[0])
    seqs = [s for s in seqs if len(s) == length]
    arr = np.frombuffer(b"".join(seqs), dtype=np.uint8).reshape(len(seqs), length)
    idx_np = LUT[arr]
    idx = torch.from_numpy(idx_np)
    n = idx.shape[0]
    n_af2 = None
    for row_block, col_block in ((ROW_BLOCK, COL_BLOCK), (512, 512)):
        try:
            n_af2, n_old = _counts(idx.cuda(), length, seq_id, "cuda", row_block, col_block)
            dev = "cuda"
            break
        except RuntimeError as exc:
            if "out of memory" not in str(exc):
                raise
            torch.cuda.empty_cache()
    if n_af2 is None:
        n_af2, n_old = _counts(idx, length, seq_id, "cpu", 1024, 256)
        dev = "cpu"
    w = 1.0 / n_af2
    has_residue = idx_np < GAP
    per_column = np.zeros(length)
    for r0 in range(0, n, 4096):
        per_column += w[r0:r0 + 4096] @ has_residue[r0:r0 + 4096].astype(np.float64)
    return (float(np.median(per_column)), float(per_column.mean()), float(w.sum()), float((1.0 / n_old).sum()), n,
            length, dev)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--store", required=True)
    ap.add_argument("--also-read", default=None)
    ap.add_argument("--threads", type=int, default=48)
    args = ap.parse_args()

    ids = np.load(DIR + "msa_ids.npy")
    zips = np.load(DIR + "msa_zip.npy")
    meth = np.load(DIR + "msa_method.npy")
    comp = np.load(DIR + "msa_compsize.npy")
    offs = np.load(DIR + "msa_offset.npy")
    keys = open(DIR + "zip_keys.txt").read().splitlines()

    want = sorted({int(x) for x in open(args.ids).read().split()})
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as fh:
            next(fh, None)
            done = {int(l.split("\t", 1)[0]) for l in fh if l.rstrip("\n").endswith(("\tok", "\tno_msa", "\tok_cpu"))}
    todo = [w for w in want if w not in done]
    print(f"ids {len(want):,} done {len(done):,} todo {len(todo):,}", flush=True)
    os.makedirs(args.store, exist_ok=True)

    def get(cid):
        i = int(np.searchsorted(ids, np.uint64(cid)))
        if i >= len(ids) or int(ids[i]) != cid:
            return cid, None, "no_msa"
        path = os.path.join(args.store, f"{cid % 1000:03d}", f"{cid}.deflate")
        if args.also_read:
            try:
                with open(os.path.join(args.also_read, f"{cid % 1000:03d}", f"{cid}.deflate"), "rb") as fh:
                    return cid, zlib.decompress(fh.read(), -15), "ok"
            except (OSError, zlib.error):
                pass
        if os.path.exists(path):
            raw = open(path, "rb").read()
        else:
            url = "https://openfold.s3.amazonaws.com/" + keys[zips[i]]
            off, size = int(offs[i]), int(comp[i])
            hdr = fetch(url, off, off + 511)
            if hdr[:4] != b"PK\x03\x04":
                return cid, None, "bad_header"
            n_len, e_len = struct.unpack_from("<HH", hdr, 26)
            start = off + 30 + n_len + e_len
            raw = fetch(url, start, start + size - 1) if size else b""
            if meth[i] == 0:
                raw = zlib.compress(raw, 6)[2:-4]
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = f"{path}.{os.getpid()}.{threading.get_ident()}.tmp"
            with open(tmp, "wb") as fh:
                fh.write(raw)
            os.replace(tmp, path)
        return cid, zlib.decompress(raw, -15), "ok"

    def get_safe(cid):
        try:
            return get(cid)
        except Exception as exc:
            return cid, None, f"fetch_error:{str(exc)[:60]}"

    q = queue.Queue(maxsize=64)

    def producer():
        with ThreadPoolExecutor(args.threads) as pool:
            for start in range(0, len(todo), 4 * args.threads):
                for res in pool.map(get_safe, todo[start:start + 4 * args.threads]):
                    q.put(res)
        q.put(None)

    threading.Thread(target=producer, daemon=True).start()
    new = not os.path.exists(args.out)
    t0, k = time.time(), 0
    with open(args.out, "a") as out:
        if new:
            out.write("rep_id\taf2style_neff\taf2style_mean\tneff_nongap_id\tneff_recomputed\tn_rows\tquery_length\tstatus\n")
        while True:
            item = q.get()
            if item is None:
                break
            cid, text, status = item
            row = [cid, "", "", "", "", "", "", status]
            if text is not None:
                try:
                    med, mean, whole, old, n, length, dev = af2style(parse_a3m(text))
                    row = [cid, f"{med:.4f}", f"{mean:.4f}", f"{whole:.4f}", f"{old:.4f}", n, length,
                           "ok" if dev == "cuda" else "ok_cpu"]
                except Exception as exc:
                    row[-1] = f"neff_error:{str(exc)[:60]}"
            out.write("\t".join(map(str, row)) + "\n")
            k += 1
            if k % 2000 == 0:
                out.flush()
                rate = k / (time.time() - t0)
                print(f"{k:,}/{len(todo):,} {rate:.1f}/s eta {(len(todo) - k) / rate / 3600:.1f} h", flush=True)
    print("DONE", k, flush=True)


if __name__ == "__main__":
    main()
