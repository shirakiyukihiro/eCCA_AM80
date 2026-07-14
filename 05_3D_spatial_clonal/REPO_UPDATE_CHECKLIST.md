# GitHub リポジトリ更新チェックリスト

対象: `shirakiyukihiro/eCCA_AM80`
理由: Data availability でコードを公開している以上、**査読者が走らせて再現できる**必要がある。
現状のままだと、走らせた図が本文の Fig. 4G と食い違う。

---

## ① 図 4G のスクリプトを差し替え 【必須】

```bash
git rm 05_3D_spatial_clonal/figure_kde_panelB_style.py
# 納品した figure_4G_kde_map.py を 05_3D_spatial_clonal/ に置く
git add 05_3D_spatial_clonal/figure_4G_kde_map.py
```

**理由:** 旧スクリプトは (a) DBSCAN 上位10クラスタの番号マーカーと (b) inset 概念図を描くが、
どちらも実際の Fig. 4G には無い。新スクリプトは点 + 2本の等高線のみを描く(=実際の図)。

---

## ② Clark–Evans のバグ修正 【必須】

`05_3D_spatial_clonal/analyze_tumor_clones_DoG_Zcrop.py` の **116–117 行目**:

```python
# 修正前(2次元用の係数 0.5 を3次元の式に使っている)
expected_nn = 0.5 / (n_total / hull_vol_um3)**(1/3)
clark_evans = mean_nn / expected_nn
```

```python
# 修正後
from math import gamma, pi
C_3D = gamma(4/3) / ((4/3) * pi) ** (1/3)      # = 0.553961
rho = n_total / hull_vol_um3
expected_nn = C_3D * rho ** (-1/3)
clark_evans = mean_nn / expected_nn
```

→ **R = 0.719 → 0.649**(集簇はより強く出る。結論は不変)

Monte-Carlo による p 値の追加は `clark_evans_fix.patch.py` を参照(任意だが推奨)。
**修正後に必ず再実行し、R = 0.649 を実測で確認すること。**

---

## ③ README_3D.md を全面差し替え 【必須】

現行の `05_3D_spatial_clonal/README_3D.md` は**実在しないファイルを大量に参照**している:

| README_3D が参照 | 実態 |
|---|---|
| `01_napari_workflow/`(17スクリプト) | **存在しない**(実際は `napari_visualisation_workflow.md`) |
| `02_quantitative_analysis/detect_cells_DoG_v3.py` | 実際は `05_3D_spatial_clonal/detect_cells_DoG.py` |
| `02_quantitative_analysis/analyze_tumor_clones_DoG_Zcrop.py` | 実際は同フォルダ直下(フラット構成) |
| `04_figures/figure_kde_panelB_style.py` | 実際は同フォルダ直下 |
| Section 3 の検証スクリプト5本 | **1本も存在しない** |
| Clark–Evans R = 0.719 | → 0.649 に修正が必要 |
| 「Figure ?, panel B」 | → Fig. 4G |

→ 納品した `README_3D.md` で丸ごと置換。

---

## ④ ルート README.md を差し替え 【必須】

`figure_kde_panelB_style.py` への参照と、「**ten largest clusters numbered**」という
記述(実際の図には無い)を修正済み。納品した `README_root.md` を `README.md` として置換。

---

## ⑤ マスクファイルの場所を修正 【必須 — 現状スクリプトが動かない】

```bash
git mv mask_polygons data/mask_polygons
```

**理由:** マスク JSON はリポジトリ直下の `mask_polygons/` にあるが、
`detect_cells_DoG.py`(33/40/47行)、`_dbscan_fallback.py`、ルート README、
`napari_visualisation_workflow.md` はすべて **`data/mask_polygons/`** を参照している。
このままでは `detect_cells_DoG.py` が入力を見つけられない。

---

## ⑥ publication_figures.py の2箇所を修正 【必須】

**(a) 3行目のヘッダコメント**

```python
# 修正前
# 3D clonal spatial analysis — publication figures (KDE / Ripley / clone stats)
# 修正後
# 3D spatial analysis — internal validation figures (detection / cluster statistics)
```

Ripley は**リポジトリのどこにも実装がない**。名前だけ出ているとレビュアーが探して混乱する。

**(b) 35–36 行目の入力パス**

```python
# 修正前(前段の analyze_... は output/ に書き出すので、data/ では見つからない)
NPZ_PATH  = 'data/BDC_tumor_DoG_Zcrop_clones.npz'
CSV_PATH  = 'data/BDC_tumor_DoG_Zcrop_cluster_features.csv'
# 修正後
NPZ_PATH  = 'output/BDC_tumor_DoG_Zcrop_clones.npz'
CSV_PATH  = 'output/BDC_tumor_DoG_Zcrop_cluster_features.csv'
```

---

## ⑦ napari_visualisation_workflow.md の1行を修正 【必須】

8–9行目:

```
see `detect_cells_DoG.py`, `analyze_tumor_clones_DoG_Zcrop.py`, and
`figure_kde_panelB_style.py`.
```
→ `figure_4G_kde_map.py` に変更。

---

## ⑧ 検証スクリプトを追加 【必須】 — 解決済み

Supplementary が引用している Richardson–Lucy デコンボリューションの検証コードが
リポジトリに無かった問題。公開用に整えた2本を追加してください。

```bash
git add 05_3D_spatial_clonal/test_deconvolution_middle_zone.py
git add 05_3D_spatial_clonal/test_deconvolution_deep_zone.py
```

**元スクリプトからの変更点:**

| 変更 | 理由 |
|---|---|
| 絶対パス `/mnt/SSD_1/ito_imaris/...` → 相対パス `data/` `output/` | 他のスクリプトと統一。そのままでは他者の環境で動かない |
| 日本語コメント → 英語 | 公開リポジトリ用 |
| ファイル名 `test_deconvolution_subvolume.py` → `test_deconvolution_middle_zone.py` | 「深部」との対比が名前で分かるように |
| **深部スクリプトを自己完結化** | 元は `BDC_tumor_DoG_clones.npz`(全ボリューム版のクラスタ結果)を要求していたが、**この npz を生成するスクリプトがリポジトリに存在しない**。フォールバックの `result_..._DoG.npz` は `cluster_labels` が全部 −1 なので、実行すると "No clusters found" で終了する。→ 同じ DBSCAN パラメータ(ε=50 µm, min_samples=3)を**スクリプト内で全ボリュームに実行**するよう変更。検出結果 npz だけあれば動く |
| PSF の µm 換算をコメントに明記 | σ_vox=(1.5,1.0,1.0) × voxel(5.22, 2.35, 2.35 µm)= **σ_z 7.83 µm / σ_xy 2.35 µm** → Supp の記載と一致することを検証済み |

**ε スイープ:** Supplementary から削除済み(新版に反映)。
「ε = 50 µm was chosen **because it approximates the observed mean nearest-neighbour
distance**」というデータに基づく根拠に置き換えました。

---

## 更新後の 05_3D_spatial_clonal/ の中身

```
05_3D_spatial_clonal/
├── README_3D.md                      ← ③ で差し替え
├── detect_cells_DoG.py
├── analyze_tumor_clones_DoG_Zcrop.py ← ② でバグ修正
├── figure_4G_kde_map.py              ← ① で新規追加
├── publication_figures.py            ← ⑥ で2箇所修正
├── _dbscan_fallback.py
├── test_deconvolution_middle_zone.py ← ⑧ で新規追加
├── test_deconvolution_deep_zone.py   ← ⑧ で新規追加
├── napari_visualisation_workflow.md  ← ⑦ で1行修正
├── environment.yml
└── requirements.txt
    (figure_kde_panelB_style.py は削除)
```

---

## 最後に

投稿前に、クリーンな環境で `data/` にマスクと `.ims` を置いた状態から
① → ② → ③ の順に実行し、**Fig. 4G が本文の図と一致すること**、
**R = 0.649 が出ること**を確認してください。
