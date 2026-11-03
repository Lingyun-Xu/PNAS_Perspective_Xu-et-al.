"""Re-export the figures of an existing plotting script as label-free panels (style of main-text Fig. 1a).

Usage (any machine):
    NOLABEL_OUT=<dir> NOLABEL_NAMES="orig_stem=new_stem,..." python nolabel_hook.py <script.py> [script args]

The script runs unchanged, except that
  * Figure.savefig no longer writes the original file; instead each panel (axes) of a figure whose file stem is
    listed in NOLABEL_NAMES is saved to NOLABEL_OUT/<new_stem>_<a|b|c...>.png, and each colorbar to
    <new_stem>_<letter>_colorbar.png (figures not listed are skipped);
  * all other file writes (open for writing, DataFrame.to_parquet) outside NOLABEL_OUT and temp/cache folders go
    to /dev/null, so rerunning a script cannot change its original outputs.
Panel style (from plddt_vs_species.ipynb, cell 42): 12 x 8 in, 700 dpi; no titles, axis labels, tick labels,
legends or text; ticks kept (width 3, length 6); axis lines width 3; top/right axis lines removed. Each panel is
enlarged as a whole, so line widths and marker sizes are scaled with it and the panel keeps its original look.
"""
import builtins
import io
import os
import pickle
import runpy
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.figure as mfig  # noqa: E402
from matplotlib.collections import LineCollection, PathCollection, PolyCollection  # noqa: E402

OUT = Path(os.environ["NOLABEL_OUT"]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
NAMES = dict(kv.split("=", 1) for kv in os.environ.get("NOLABEL_NAMES", "").split(",") if "=" in kv)
W, H = float(os.environ.get("NOLABEL_W", 12)), float(os.environ.get("NOLABEL_H", 8))
DPI = int(os.environ.get("NOLABEL_DPI", 700))
TICK_W, TICK_L, SPINE_W = 3, 6, 3
MIN_LW = 2.0
_done = set()

# ---------------------------------------------------------------- block writes outside the output folder
_ALLOWED = [str(OUT), tempfile.gettempdir(), "/tmp", "/dev/", str(Path.home() / ".cache"), str(Path.home() / ".config"),
            str(Path.home() / ".matplotlib")]
_open = builtins.open


def _blocked(file, mode):
    if not isinstance(file, (str, bytes, os.PathLike)) or not any(m in mode for m in "wax+"):
        return False
    p = os.path.abspath(os.fsdecode(file))
    return not any(p.startswith(a) for a in _ALLOWED)


def _guarded_open(file, mode="r", *args, **kwargs):
    if _blocked(file, mode):
        print(f"[nolabel] write skipped: {file}", file=sys.stderr)
        return _open(os.devnull, mode, *args, **kwargs)
    return _open(file, mode, *args, **kwargs)


builtins.open = _guarded_open
io.open = _guarded_open
try:
    import pandas as pd

    def _no_parquet(self, path=None, *a, **k):
        if path is None or not _blocked(path, "w"):
            return _orig_parquet(self, path, *a, **k)
        print(f"[nolabel] write skipped: {path}", file=sys.stderr)

    _orig_parquet = pd.DataFrame.to_parquet
    pd.DataFrame.to_parquet = _no_parquet
except ImportError:
    pass


# ---------------------------------------------------------------- panel export
def _is_cbar(ax):
    return getattr(ax, "_colorbar", None) is not None or ax.get_label() == "<colorbar>"


def _has_content(ax):
    return bool(ax.lines or ax.collections or ax.images or len(ax.patches) > 0 or ax.containers)


def _strip(ax, f, colorbar=False):
    for loc in ("center", "left", "right"):
        ax.set_title("", loc=loc)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(which="both", labelbottom=False, labelleft=False, labeltop=False, labelright=False)
    ax.tick_params(which="major", width=TICK_W, length=TICK_L)
    ax.tick_params(which="minor", width=TICK_W * 0.67, length=TICK_L * 0.5)
    if ax.get_legend() is not None:
        ax.get_legend().set_visible(False)
    for t in ax.texts:
        t.set_visible(False)
    for child in ax.child_axes:
        child.set_visible(False)
    for name, s in ax.spines.items():
        s.set_linewidth(SPINE_W)
        if not colorbar and name in ("top", "right") and not ax.images:
            s.set_visible(False)
    if colorbar:
        return
    def lw(w):  # enlarge with the panel; drawn lines at least MIN_LW points wide
        return 0 if w <= 0 else max(w * f, MIN_LW)

    for ln in ax.lines:
        ln.set_linewidth(lw(ln.get_linewidth()))
        ln.set_markersize(ln.get_markersize() * f)
        ln.set_markeredgewidth(ln.get_markeredgewidth() * f)
    for c in ax.collections:
        if isinstance(c, PathCollection):
            c.set_sizes(c.get_sizes() * f * f)
            c.set_linewidths([w * f for w in c.get_linewidths()])
        elif isinstance(c, LineCollection):
            c.set_linewidths([lw(w) for w in c.get_linewidths()])
        elif isinstance(c, PolyCollection):
            c.set_linewidths([w * f for w in c.get_linewidths()])
    for p in ax.patches:
        p.set_linewidth(lw(p.get_linewidth()))


def _export(fig, new_stem):
    try:
        fig = pickle.loads(pickle.dumps(fig))
    except Exception as e:  # unpicklable figure: modify in place (the script's own save is suppressed anyway)
        print(f"[nolabel] pickle failed ({e}); exporting in place", file=sys.stderr)
    fig_w, fig_h = fig.get_size_inches()
    axes = [a for a in fig.axes if a.get_visible()]
    for t in list(fig.texts):
        t.set_visible(False)
    for attr in ("_suptitle", "_supxlabel", "_supylabel"):
        if getattr(fig, attr, None) is not None:
            getattr(fig, attr).set_visible(False)
    for leg in fig.legends:
        leg.set_visible(False)
    # group axes that share a position (twin axes) so they are drawn together
    groups = []
    for a in axes:
        if _is_cbar(a):
            continue
        pos = tuple(round(v, 4) for v in a.get_position().bounds)
        for g in groups:
            if g[0] == pos:
                g[1].append(a)
                break
        else:
            groups.append((pos, [a]))
    groups = [(pos, g) for pos, g in groups if any(_has_content(a) for a in g)]
    cbars = []  # (axes, owner panel index, vertical), decided before any axes is moved
    for cb in (a for a in axes if _is_cbar(a)):
        b = cb.get_position().bounds
        if not groups:
            break
        owner = min(range(len(groups)), key=lambda j: abs(groups[j][0][0] + groups[j][0][2] - b[0])
                    + abs(groups[j][0][1] - b[1]))
        cbars.append((cb, owner, b[3] * fig_h >= b[2] * fig_w))
    written = []
    for i, (pos, group) in enumerate(groups):
        letter = chr(ord("a") + i)
        f = min(W * 0.78 / (pos[2] * fig_w), H * 0.78 / (pos[3] * fig_h))
        for a in axes:
            a.set_visible(False)
        fig.set_size_inches(W, H)
        for a in group:
            a.set_visible(True)
            a.set_position([0.12, 0.12, 0.78, 0.78])
            _strip(a, f)
        path = OUT / f"{new_stem}_{letter}.png"
        _orig_savefig(fig, path, dpi=DPI)
        written.append(path.name)
        for a in group:
            a.set_visible(False)
        # colorbars that belong to this panel (closest by position)
        for k, (cb, owner, vertical) in enumerate(cbars):
            if owner != i:
                continue
            fig.set_size_inches((1.6, H) if vertical else (W, 1.6))
            cb.set_visible(True)
            cb.set_position([0.35, 0.12, 0.3, 0.78] if vertical else [0.12, 0.35, 0.78, 0.3])
            _strip(cb, 1, colorbar=True)
            n_own = sum(1 for c in cbars if c[1] == i)
            path = OUT / (f"{new_stem}_{letter}_colorbar.png" if n_own == 1 else f"{new_stem}_{letter}_colorbar{k + 1}.png")
            _orig_savefig(fig, path, dpi=DPI)
            written.append(path.name)
            cb.set_visible(False)
    print(f"[nolabel] {new_stem}: {', '.join(written)}", file=sys.stderr)


_orig_savefig = mfig.Figure.savefig


def _savefig(self, fname, *args, **kwargs):
    if not isinstance(fname, (str, bytes, os.PathLike)):
        return
    p = Path(os.fsdecode(fname))
    key = f"{p.parent.name}/{p.stem}" if f"{p.parent.name}/{p.stem}" in NAMES else p.stem
    if key in NAMES and NAMES[key] not in _done:
        _done.add(NAMES[key])
        _export(self, NAMES[key])


mfig.Figure.savefig = _savefig

if __name__ == "__main__":
    script = sys.argv[1]
    sys.argv = sys.argv[1:]
    sys.path.insert(0, str(Path(script).resolve().parent))
    runpy.run_path(script, run_name="__main__")
    missing = [v for v in NAMES.values() if v not in _done]
    print(f"[nolabel] done; not produced: {missing}" if missing else "[nolabel] done; all requested figures exported",
          file=sys.stderr)
