#!/usr/bin/env python3
# ============================================================
# Fig 5H — Keyence .bcf channel / consistency QC (supporting)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming
#              strategy in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Verifies fluorophore-channel structure across fields.
# Paths made relative for public release: put inputs in ./data, outputs in ./output.
# ============================================================

"""
inspect_bcf_channels.py
=======================
Walks a study folder, reads every Img.bcf, and reports per folder:
  - patient id, lens, um/px
  - which FLUOROPHORE channels were shot (DAPI / GFP / TRITC / Cy5 ...)
Then checks whether the channel setup is CONSISTENT across all folders.

IMPORTANT: a Keyence .bcf records the fluorophore/filter per channel, NOT which
antibody (Meflin vs Thy1) is on which channel, and no explicit R/G/B display
color. So this tool verifies the *fluorophore structure* and *consistency*; the
target<->channel identity must come from the staining protocol. (And note the
co-expression result -- the double-positive %% -- is invariant to swapping the
two markers, so only the single-positive labels depend on that identity.)

USAGE
    python inspect_bcf_channels.py "data/HBDC_Meflin_Thy1" --csv bcf_channels.csv
"""
import os, re, sys, glob, zipfile, struct, argparse, collections
import pandas as pd


def _double(int_str):
    return struct.unpack("<d", struct.pack("<q", int(int_str)))[0]


def read_bcf(path):
    z = zipfile.ZipFile(path)
    names = set(z.namelist())

    def rd(p):
        return z.read(p).decode("utf-8", "replace") if p in names else ""

    def get(xml, tag):
        m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", xml, re.S)
        return m.group(1) if m else None

    lens = rd("GroupFileProperty/Lens/properties.xml")
    img = rd("GroupFileProperty/Image/properties.xml")
    cal = get(img, "Calibration")
    umpx = round(_double(cal) / 1000.0, 5) if cal else None

    shot = []
    for c in range(8):
        p = f"GroupFileProperty/Channel{c}/properties.xml"
        if p not in names:
            break
        x = rd(p)
        if get(x, "IsShot") == "True" and get(x, "IsOverlay") != "True":
            shot.append(get(x, "Comment") or "(unnamed)")
    return {
        "lens": get(lens, "LensName"),
        "mag_x": (int(get(lens, "Magnification")) / 100.0 if get(lens, "Magnification") else None),
        "um_per_px": umpx,
        "channels_shot": "+".join(shot),
        "n_channels": len(shot),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    bcfs = ([a.root] if a.root.endswith(".bcf")
            else sorted(glob.glob(os.path.join(a.root, "**", "Img.bcf"), recursive=True)
                        + glob.glob(os.path.join(a.root, "**", "*.bcf"), recursive=True)))
    bcfs = sorted(set(bcfs))
    if not bcfs:
        sys.exit(f"No .bcf found under {a.root}")

    rows = []
    for b in bcfs:
        folder = os.path.basename(os.path.dirname(b))
        pid = (re.search(r"(\d+_\d+)", folder) or [None, folder])[1] \
            if re.search(r"(\d+_\d+)", folder) else folder
        try:
            info = read_bcf(b)
        except Exception as e:
            print(f"[ERROR] {b}: {e}")
            continue
        info.update({"patient": pid, "folder": folder})
        rows.append(info)

    df = pd.DataFrame(rows)[["patient", "folder", "mag_x", "um_per_px",
                             "n_channels", "channels_shot", "lens"]]
    pd.set_option("display.max_rows", None, "display.width", 200)
    print(df.to_string(index=False))

    print("\n================ CONSISTENCY ================")
    for col in ["mag_x", "um_per_px", "channels_shot", "n_channels"]:
        vc = collections.Counter(df[col].astype(str))
        verdict = "OK (uniform)" if len(vc) == 1 else "!! MIXED !!"
        print(f"{col:16s}: {verdict}  {dict(vc)}")
    print("\nPer-patient field counts:")
    print(df.groupby("patient").size().to_string())

    if a.csv:
        df.to_csv(a.csv, index=False, encoding="cp932")
        print(f"\nWrote {a.csv} (CP932)")


if __name__ == "__main__":
    main()
