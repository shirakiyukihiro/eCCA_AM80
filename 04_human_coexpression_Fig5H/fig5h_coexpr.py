#!/usr/bin/env python3
# ============================================================
# Fig 5H — per-cell Meflin/Thy1 co-expression (human dCCA dual-IF)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming
#              strategy in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Cellpose nuclei segmentation + per-cell Otsu positivity; input = Keyence field folders.
# Paths made relative for public release: put inputs in ./data, outputs in ./output.
# ============================================================

"""
fig5h_coexpr.py  (v2)
=====================
Per-cell Meflin / Thy1 co-expression quantification for human dCCA dual-IF (Fig 5H).

v2 change (important): positivity is decided per channel by **Otsu on the
per-cell intensities** (default), NOT by "K x sigma over background". The
background method fails for Thy1 because abundant extracellular/fibrillar Thy1
inflates the background noise and pushes the threshold above almost all cells
(Thy1 was under-called). Otsu finds the natural per-channel split and is
symmetric between channels. The legacy method is still available (--method bgsigma).

Channel mapping is auto-detected from the folder name "Meflin(G)_Thy1(R)"
(Fig 1 was ISLR=R; here Meflin=G, Thy1=R -- never hard-code).

USAGE
    python fig5h_coexpr.py "<field folder OR study root>" --out OUT
Options:
    --method otsu|bgsigma   (default otsu)
    --ring-um 2.0           perinuclear ring width (um)
    --mult 1.0              multiply the Otsu threshold (sensitivity operating point)
    --sweep 0.7,0.85,1.0,1.15,1.3   threshold multipliers reported for robustness
    --k 3.0                 only used by --method bgsigma
    --no-cellpose           force the classical watershed fallback
"""
import os, re, sys, glob, zipfile, struct, argparse
import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, gaussian
from skimage.feature import peak_local_max
from skimage.segmentation import watershed, expand_labels
from skimage.measure import regionprops_table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ----------------------------------------------------------------- BCF calibration
def bcf_um_per_px(bcf_path):
    try:
        z = zipfile.ZipFile(bcf_path)
        xml = z.read("GroupFileProperty/Image/properties.xml").decode("utf-8", "replace")
        m = re.search(r"<Calibration[^>]*>(.*?)</Calibration>", xml)
        if not m:
            return None
        nm = struct.unpack("<d", struct.pack("<q", int(m.group(1))))[0]
        return nm / 1000.0
    except Exception:
        return None

def bcf_channels(bcf_path):
    """Fluorophores actually shot (e.g. 'DAPI+GFP+TRITC'). Target identity is NOT in the BCF."""
    try:
        z = zipfile.ZipFile(bcf_path); names = set(z.namelist()); shot = []
        for c in range(8):
            p = f"GroupFileProperty/Channel{c}/properties.xml"
            if p not in names:
                break
            x = z.read(p).decode("utf-8", "replace")
            def g(t):
                m = re.search(rf"<{t}[^>]*>(.*?)</{t}>", x); return m.group(1) if m else None
            if g("IsShot") == "True" and g("IsOverlay") != "True":
                shot.append(g("Comment") or "(unnamed)")
        return "+".join(shot)
    except Exception:
        return ""

def pixel_coloc(mef, thy, dapi):
    """Segmentation-FREE co-localization within the tissue mask (Manders/Jaccard/Pearson).
    Jaccard/Pearson are swap-invariant; M1/M2 swap if the two markers are swapped."""
    tot = gaussian(mef + thy + dapi, sigma=3)
    if tot.max() <= 0:
        return {}
    tis = tot > threshold_otsu(tot)
    def pmask(ch):
        v = ch[tis]
        if v.size < 10 or not (v > 0).any():
            return np.zeros(ch.shape, bool)
        return (ch > threshold_otsu(v[v > 0])) & tis
    mm, tm = pmask(mef), pmask(thy)
    inter = int((mm & tm).sum()); union = int((mm | tm).sum())
    sm, st = mef[tis].sum(), thy[tis].sum()
    a, b = mef[tis], thy[tis]
    return {
        "coloc_jaccard": round(inter / union, 3) if union else np.nan,
        "manders_M1_MefInThy": round(mef[tm].sum() / sm, 3) if sm > 0 else np.nan,
        "manders_M2_ThyInMef": round(thy[mm].sum() / st, 3) if st > 0 else np.nan,
        "pixel_pearson_tissue": round(float(np.corrcoef(a, b)[0, 1]), 3) if a.size > 2 else np.nan,
    }

# ----------------------------------------------------------------- IO
def channel_map_from_name(folder, swap=False):
    name = os.path.basename(folder.rstrip("/"))
    mef = re.search(r"Meflin\((\w)\)", name, re.I)
    thy = re.search(r"Thy1\((\w)\)", name, re.I)
    mef = mef.group(1).upper() if mef else "G"
    thy = thy.group(1).upper() if thy else "R"
    if swap:
        mef, thy = thy, mef
    dapi = (set("RGB") - {mef, thy}).pop()
    return {"Meflin": mef, "Thy1": thy, "DAPI": dapi}

def patient_id(folder):
    m = re.search(r"(\d+_\d+)", os.path.basename(folder.rstrip("/")))
    return m.group(1) if m else os.path.basename(folder.rstrip("/"))

def load_zstack_mip(folder, ch_num):
    """Max-intensity projection of Img_Z*_CH{ch_num}.tif (full-focus approximation)."""
    files = sorted(glob.glob(os.path.join(folder, f"Img_Z*_CH{ch_num}.tif")))
    if not files:
        return None
    planes = []
    for f in files:
        s = tifffile.imread(f)
        planes.append(s if s.ndim == 2 else s.max(axis=-1))
    return np.max(np.stack(planes), axis=0).astype(np.float32)

def _read_hrff(folder, letter):
    path = os.path.join(folder, f"HR_FF_Img_{letter}.tif")
    if not os.path.exists(path):
        return None
    a = tifffile.imread(path)
    if a.ndim == 3 and a.shape[-1] in (3, 4):
        a = a[..., {"R": 0, "G": 1, "B": 2}[letter]]
    elif a.ndim == 3:
        a = a.max(axis=0)
    return a.astype(np.float32)

def has_hrff(folder):
    return all(os.path.exists(os.path.join(folder, f"HR_FF_Img_{l}.tif")) for l in "RGB")

def learn_ch_mapping(folder):
    """Determine Z-stack CH(1/2/3) -> display R/G/B by correlating each CH's MIP with
    the folder's HR_FF planes. Returns ({'R':ch,'G':ch,'B':ch}, {letter: corr})."""
    chans = {n: load_zstack_mip(folder, n) for n in (1, 2, 3, 4)}
    chans = {n: m for n, m in chans.items() if m is not None}
    if not chans:
        return None, None
    mapping, corrs, used = {}, {}, set()
    for letter in "BGR":                      # assign B(DAPI) first, then G, R
        ref = _read_hrff(folder, letter)
        if ref is None:
            return None, None
        rf = ref.ravel().astype(np.float64)
        best, bestr = None, -2.0
        for n, m in chans.items():
            if n in used or m.shape != ref.shape:
                continue
            r = np.corrcoef(rf, m.ravel().astype(np.float64))[0, 1]
            if r > bestr:
                bestr, best = r, n
        mapping[letter], corrs[letter], _ = best, round(float(bestr), 3), used.add(best)
    return mapping, corrs

def load_plane(folder, letter, ch_map=None, prefer_zstack=False):
    """HR_FF plane if present (unless prefer_zstack); else Z-stack MIP via ch_map."""
    if not prefer_zstack:
        a = _read_hrff(folder, letter)
        if a is not None:
            return a
    if ch_map and ch_map.get(letter):
        m = load_zstack_mip(folder, ch_map[letter])
        if m is not None:
            return m
    a = _read_hrff(folder, letter)
    if a is not None:
        return a
    raise FileNotFoundError(f"No HR_FF or mapped Z-stack channel for '{letter}' in {folder}")

# ----------------------------------------------------------------- segmentation
def segment_nuclei(dapi, diam_px, use_cellpose=True):
    if use_cellpose:
        try:
            from cellpose import models
            diam = float(diam_px) if diam_px else None
            try:                                   # cellpose v3.x
                model = models.Cellpose(gpu=True, model_type="nuclei")
                res = model.eval(dapi, channels=[0, 0], diameter=diam)
            except (TypeError, AttributeError):    # cellpose v4.x
                model = models.CellposeModel(gpu=True)
                res = model.eval(dapi, diameter=diam)
            return np.asarray(res[0]).astype(np.int32), "cellpose"
        except Exception as e:
            sys.stderr.write(f"[cellpose unavailable -> watershed] {e}\n")
    sm = gaussian(dapi, sigma=max(diam_px / 8.0, 1.0))
    if sm.max() <= 0:
        return np.zeros(dapi.shape, np.int32), "empty"
    mask = ndi.binary_opening(ndi.binary_fill_holes(sm > threshold_otsu(sm)), iterations=1)
    dist = ndi.distance_transform_edt(mask)
    coords = peak_local_max(dist, min_distance=max(int(diam_px * 0.5), 3),
                            labels=mask, exclude_border=False)
    markers = np.zeros(dapi.shape, np.int32)
    for i, (r, c) in enumerate(coords, 1):
        markers[r, c] = i
    if markers.max() == 0:
        return ndi.label(mask)[0].astype(np.int32), "watershed"
    return watershed(-dist, markers, mask=mask).astype(np.int32), "watershed"

# ----------------------------------------------------------------- thresholds
def channel_threshold(vals, method, bg, sigma, k=3.0, mult=1.0):
    """Return positivity threshold for one channel's per-cell intensities."""
    v = vals[vals > 0]
    if method == "bgsigma":
        return k * sigma
    # otsu (default): natural split of the cell-intensity distribution
    if v.size < 10 or v.max() <= 0:
        return k * sigma                      # too few cells -> safe fallback
    try:
        return float(threshold_otsu(v)) * mult
    except Exception:
        return k * sigma

# ----------------------------------------------------------------- core
def analyse_field(folder, method="otsu", k=3.0, mult=1.0, ring_um=2.0,
                  use_cellpose=True, swap=False, ch_map=None, prefer_zstack=False,
                  sweep=(0.7, 0.85, 1.0, 1.15, 1.3), outdir="out"):
    cmap = channel_map_from_name(folder, swap=swap)
    pid, fid = patient_id(folder), os.path.basename(folder.rstrip("/"))
    src = "zstack_MIP" if (prefer_zstack or not has_hrff(folder)) else "HR_FF"
    dapi = load_plane(folder, cmap["DAPI"], ch_map, prefer_zstack)
    mef = load_plane(folder, cmap["Meflin"], ch_map, prefer_zstack)
    thy = load_plane(folder, cmap["Thy1"], ch_map, prefer_zstack)

    bcf = os.path.join(folder, "Img.bcf")
    umpx = bcf_um_per_px(bcf) if os.path.exists(bcf) else None
    fluor = bcf_channels(bcf) if os.path.exists(bcf) else ""
    if fluor:
        sys.stderr.write(f"   [{fid}] BCF fluorophores shot: {fluor} | "
                         f"mapping used Meflin={cmap['Meflin']}, Thy1={cmap['Thy1']}, DAPI={cmap['DAPI']}"
                         f"{'  (SWAPPED)' if swap else ''}\n")
    ring_px = max(int(round(ring_um / umpx)), 1) if umpx else 3
    diam_px = max(int(round(7.0 / umpx)), 6) if umpx else 14

    labels, seg = segment_nuclei(dapi, diam_px, use_cellpose=use_cellpose)
    n_nuc = int(labels.max())
    if n_nuc == 0:
        sys.stderr.write(f"[{fid}] no nuclei; skipped\n")
        return None, None

    cells = expand_labels(labels, distance=ring_px)
    noncell = cells == 0
    def bg_sigma(ch):
        x = ch[noncell]
        bg = float(np.median(x))
        return bg, float(1.4826 * np.median(np.abs(x - bg)) + 1e-6)
    bg_m, sg_m = bg_sigma(mef)
    bg_t, sg_t = bg_sigma(thy)

    props = regionprops_table(cells, intensity_image=np.dstack([mef, thy]),
                              properties=("label", "centroid", "area", "intensity_mean"))
    df = pd.DataFrame(props).rename(columns={
        "centroid-0": "y", "centroid-1": "x",
        "intensity_mean-0": "mef_raw", "intensity_mean-1": "thy_raw"})
    df["mef"] = (df["mef_raw"] - bg_m).clip(lower=0)
    df["thy"] = (df["thy_raw"] - bg_t).clip(lower=0)
    df["patient"], df["field"] = pid, fid

    thr_m = channel_threshold(df["mef"].values, method, bg_m, sg_m, k, mult)
    thr_t = channel_threshold(df["thy"].values, method, bg_t, sg_t, k, mult)
    df["mef_pos"] = df["mef"] > thr_m
    df["thy_pos"] = df["thy"] > thr_t
    df["any_pos"] = df["mef_pos"] | df["thy_pos"]
    df["class"] = np.where(df["mef_pos"] & df["thy_pos"], "double",
                  np.where(df["mef_pos"], "Meflin_only",
                  np.where(df["thy_pos"], "Thy1_only", "negative")))

    pos = df[df["any_pos"]]; n_pos = len(pos)
    r = float(np.corrcoef(np.log1p(df["mef"]), np.log1p(df["thy"]))[0, 1])
    summ = {"patient": pid, "field": fid, "segmenter": seg, "method": method,
            "um_per_px": round(umpx, 5) if umpx else "", "n_nuclei": n_nuc,
            "n_pos": n_pos,
            "n_Meflin_only": int((df["class"] == "Meflin_only").sum()),
            "n_Thy1_only": int((df["class"] == "Thy1_only").sum()),
            "n_double": int((df["class"] == "double").sum()),
            "thr_Meflin": round(thr_m, 3), "thr_Thy1": round(thr_t, 3),
            "corr_logMef_logThy": round(r, 3)}
    summ["bcf_fluorophores"] = fluor
    summ["source"] = src
    summ["Meflin_ch"], summ["Thy1_ch"] = cmap["Meflin"], cmap["Thy1"]
    summ.update(pixel_coloc(mef, thy, dapi))
    if n_pos:
        summ["pct_Meflin_only"] = round(100 * summ["n_Meflin_only"] / n_pos, 2)
        summ["pct_Thy1_only"] = round(100 * summ["n_Thy1_only"] / n_pos, 2)
        summ["pct_double"] = round(100 * summ["n_double"] / n_pos, 2)

    # robustness: double%% across threshold multipliers (otsu) or K (bgsigma)
    base_m = thr_m / mult if (method == "otsu" and mult) else thr_m
    base_t = thr_t / mult if (method == "otsu" and mult) else thr_t
    sweep_pct = {}
    for f in sweep:
        if method == "otsu":
            tm, tt = base_m * f, base_t * f
            key = f"pct_double_x{f}"
        else:
            tm, tt = f * sg_m, f * sg_t
            key = f"pct_double_k{f}"
        mp = df["mef"] > tm; tp = df["thy"] > tt; ap = mp | tp
        val = round(100 * (mp & tp).sum() / ap.sum(), 2) if ap.sum() else np.nan
        summ[key] = val; sweep_pct[f] = val

    _qc(folder, outdir, dapi, mef, thy, df, thr_m, thr_t, sweep, sweep_pct,
        method, mult, fid, seg)
    return summ, df

# ----------------------------------------------------------------- QC images
def _qc(folder, outdir, dapi, mef, thy, df, thr_m, thr_t, sweep, sweep_pct,
        method, mult, fid, seg):
    os.makedirs(outdir, exist_ok=True)
    def n(x, p=99.5):
        v = x[x > 0]; hi = np.percentile(v, p) if v.size else 1
        return np.clip(x / max(hi, 1e-6), 0, 1)
    col = {"double": "gold", "Meflin_only": "lime", "Thy1_only": "red",
           "negative": (0.55, 0.55, 0.55)}

    rgb = np.dstack([n(thy), n(mef), n(dapi)])      # Thy1=R, Meflin=G, DAPI=B
    fig, ax = plt.subplots(figsize=(7, 5)); ax.imshow(rgb); ax.axis("off")
    for cls, c in col.items():
        s = df[df["class"] == cls]
        ax.scatter(s["x"], s["y"], s=6, facecolors="none", edgecolors=c,
                   linewidths=0.5, label=f"{cls} ({len(s)})")
    ax.legend(loc="upper right", fontsize=6, framealpha=0.4)
    ax.set_title(f"{fid}  [{seg}/{method}]", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/{fid}__seg_overlay.png", dpi=130); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    for cls, c in col.items():
        s = df[df["class"] == cls]
        ax.scatter(s["mef"] + 1, s["thy"] + 1, s=8, c=[c], alpha=0.6,
                   edgecolors="none", label=f"{cls} ({len(s)})")
    ax.axvline(thr_m, color="gray", ls="--", lw=0.8)
    ax.axhline(thr_t, color="gray", ls="--", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Meflin (bg-sub,+1)"); ax.set_ylabel("Thy1 (bg-sub,+1)")
    ax.set_title(f"{fid}  per-cell ({method}) threshold-free view", fontsize=8)
    ax.legend(fontsize=6, framealpha=0.4)
    fig.tight_layout(); fig.savefig(f"{outdir}/{fid}__scatter.png", dpi=130); plt.close(fig)

    xs = list(sweep); ys = [sweep_pct[f] for f in xs]
    fig, ax = plt.subplots(figsize=(5, 3.2)); ax.plot(xs, ys, "o-")
    ax.axvline(mult if method == "otsu" else 3.0, color="red", ls="--", lw=0.8,
               label=("operating mult" if method == "otsu" else "K"))
    ax.set_xlabel("Otsu threshold multiplier" if method == "otsu"
                  else "K (x sigma over background)")
    ax.set_ylabel("double-positive %% of positive cells")
    ax.set_ylim(0, 100); ax.legend(fontsize=7)
    ax.set_title(f"{fid}  threshold sensitivity", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/{fid}__sensitivity.png", dpi=130); plt.close(fig)

# ----------------------------------------------------------------- driver
def find_field_folders(root, exclude=None):
    pat = re.compile(exclude) if exclude else None
    def ok(d):
        return not (pat and pat.search(os.path.basename(d.rstrip("/"))))
    def is_field(fs):
        return (any(f.startswith("HR_FF_Img_") and f.endswith(".tif") for f in fs)
                or any(re.match(r"Img_Z\d+_CH\d+\.tif$", f) for f in fs))
    if os.path.isdir(root) and is_field(os.listdir(root)):
        return [root] if ok(root) else []
    out = [d for d, _, fs in os.walk(root) if ok(d) and is_field(fs)]
    return sorted(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--out", default="fig5h_out")
    ap.add_argument("--method", choices=["otsu", "bgsigma"], default="otsu")
    ap.add_argument("--k", type=float, default=3.0)
    ap.add_argument("--mult", type=float, default=1.0)
    ap.add_argument("--ring-um", type=float, default=2.0)
    ap.add_argument("--no-cellpose", action="store_true")
    ap.add_argument("--swap", action="store_true",
                    help="swap Meflin<->Thy1 channel assignment (apply uniformly if folder names are mislabelled)")
    ap.add_argument("--sweep", default="0.7,0.85,1.0,1.15,1.3")
    ap.add_argument("--exclude", default=None,
                    help="regex on folder name to skip (e.g. '_x40' to keep 20x only)")
    ap.add_argument("--zstack-all", action="store_true",
                    help="use Z-stack MIP for ALL folders (uniform processing); HR_FF used only to calibrate channels")
    a = ap.parse_args()
    sweep = tuple(float(x) for x in a.sweep.split(","))

    folders = find_field_folders(a.path, exclude=a.exclude)
    if not folders:
        sys.exit(f"No field folders under {a.path}")

    # learn Z-stack CH -> display R/G/B from a folder that has HR_FF (for Z-stack-only folders)
    ch_map = None
    refs = [f for f in folders if has_hrff(f)]
    if refs:
        ch_map, corrs = learn_ch_mapping(refs[0])
        print(f"CH->color calibration (from {os.path.basename(refs[0])}): {ch_map}  corr={corrs}")
    n_zonly = sum(1 for f in folders if not has_hrff(f))
    if n_zonly and not ch_map:
        sys.exit("ERROR: %d folders lack HR_FF and no HR_FF reference exists to calibrate "
                 "Z-stack channels. Include at least one folder with HR_FF." % n_zonly)
    if n_zonly:
        print(f"{n_zonly} folder(s) without HR_FF will use Z-stack MIP via the calibration above.")

    print(f"Found {len(folders)} field folder(s). method={a.method}"
          f"{'  SWAP' if a.swap else ''}{'  ZSTACK-ALL' if a.zstack_all else ''}")
    rows, cells = [], []
    for f in folders:
        print(f"  -> {os.path.basename(f)}")
        s, df = analyse_field(f, method=a.method, k=a.k, mult=a.mult,
                              ring_um=a.ring_um, use_cellpose=not a.no_cellpose,
                              swap=a.swap, ch_map=ch_map, prefer_zstack=a.zstack_all,
                              sweep=sweep, outdir=a.out)
        if s:
            rows.append(s); cells.append(df)
    if not rows:
        sys.exit("No fields produced results.")
    os.makedirs(a.out, exist_ok=True)
    img_df = pd.DataFrame(rows)
    img_df.to_csv(f"{a.out}/per_image_summary.csv", index=False, encoding="cp932")
    pd.concat(cells, ignore_index=True).to_csv(f"{a.out}/per_cell.csv",
                                               index=False, encoding="cp932")
    if {"pct_double", "pct_Meflin_only", "pct_Thy1_only"}.issubset(img_df.columns):
        pat = (img_df.groupby("patient")[["pct_Meflin_only", "pct_Thy1_only", "pct_double"]]
               .agg(["mean", "count"]))
        pat.columns = ["_".join(c) for c in pat.columns]
        pat.to_csv(f"{a.out}/per_patient_summary.csv", encoding="cp932")
        m = pat[["pct_Meflin_only_mean", "pct_Thy1_only_mean", "pct_double_mean"]]
        print(f"\nStudy-level (patient as unit, n={len(pat)}):")
        for c in m.columns:
            print(f"  {c.replace('_mean',''):16s}: {m[c].mean():5.1f}% +/- {m[c].sem():.1f}% s.e.m.")
        fig, ax = plt.subplots(figsize=(5, 3.6))
        cats = list(m.columns); labs = ["Meflin only", "Thy1 only", "double"]
        ax.bar(labs, m.mean()[cats], yerr=m.sem()[cats], capsize=4,
               color=["lime", "red", "gold"], alpha=0.7, edgecolor="grey")
        for pid in pat.index:
            ax.scatter(range(3), pat.loc[pid, cats], color="k", s=12, zorder=3)
        ax.set_ylabel("%% of marker-positive cells"); ax.set_ylim(0, 100)
        ax.set_title(f"Fig 5H co-expression (n={len(pat)}, {a.method})", fontsize=9)
        fig.tight_layout(); fig.savefig(f"{a.out}/per_patient_proportions.png", dpi=140)
        plt.close(fig)
    print(f"\nWrote CSVs + QC images to: {a.out}/")

if __name__ == "__main__":
    main()
