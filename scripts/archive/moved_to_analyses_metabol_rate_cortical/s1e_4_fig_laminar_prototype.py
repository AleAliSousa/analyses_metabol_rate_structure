#!/usr/bin/env python
"""s1e_4 (PROTOTYPE) figure — laminar strata: receptor classes vs BigBrain cell density (human), and
cross-species laminar GLI.

Inputs (written by scripts/s1e_4_laminar_prototype.R):
  tables/s1e/s1e_laminar_strata_panel.csv   (panel_set == "primary": 6 area keys x 3 strata)
  tables/s1e/s1e_laminar_gli_species.csv    (areas 44/45, 7 species; posterior cingulate 23/31, 17 species)
Output: figs/s1e/laminar_prototype.png / .pdf

Row A: one panel per area; x = stratum (supragranular I-III, granular IV, infragranular V-VI); y = z-score within
        the area across its three strata. Receptor classes (Zilles & Palomero-Gallagher 2017, fmol/mg protein;
        3-brain means) as thin coloured lines; BigBrain cell-body density (Dickscheid et al. 2021, cells/0.1 mm^3,
        one donor; thickness-weighted over layers) as the heavy black line. Area 4 is agranular (BigBrain 4a
        layer IV = 1.2 % of depth; Zilles used a 3 % stripe) and is drawn with an open marker.
Row B: GLI (%) by stratum per species; human heavy, other species thin grey; species means over brains.
Run from the repo root: python scripts/s1e_4_fig_laminar_prototype.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = os.path.abspath(os.getcwd())
while not os.path.exists(os.path.join(d, ".git")) and os.path.dirname(d) != d:
    d = os.path.dirname(d)
os.chdir(d)

if "apply_figure_style" in globals():
    apply_figure_style(sizes=(8, 7, 6))
else:
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                         "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlelocation": "left", "axes.titleweight": "normal", "legend.frameon": False})

panel = pd.read_csv("tables/s1e/s1e_laminar_strata_panel.csv")
gli = pd.read_csv("tables/s1e/s1e_laminar_gli_species.csv")
pp = panel[panel.panel_set == "primary"].copy()
strata = ["supragranular", "granular", "infragranular"]
xt = ["I–III", "IV", "V–VI"]
area_order = ["V1", "V2d+V2v", "V3d", "FG1", "44", "4"]
area_title = {"V1": "V1 (hOc1)", "V2d+V2v": "V2 (hOc2)", "V3d": "V3d (hOc3d)", "FG1": "FG1",
              "44": "44 (Broca)", "4": "4 (4a; agranular)"}
classes = [("class_iGlu_z_within_area", "iGlu (AMPA+NMDA+kainate)", "#1b9e77"),
           ("class_GABA_z_within_area", "GABA (A, A/BZ, B)", "#7570b3"),
           ("class_ACh3_z_within_area", "ACh (M1+M2+M3)", "#d95f02"),
           ("class_monoamineGPCR_z_within_area", "monoamine GPCR (α1, α2, 5-HT1A, 5-HT2, D1)", "#e7298a")]
cell_col = "cell_density_per_0.1mm3_z_within_area"

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.15], hspace=0.75, wspace=0.08, top=0.86, bottom=0.08)
axesA = [fig.add_subplot(gs[0, i]) for i in range(6)]
x = np.arange(3)
for ax, a in zip(axesA, area_order):
    sub = pp[pp.area == a].set_index("stratum").loc[strata]
    for col, lab, c in classes:
        ax.plot(x, sub[col].values, "-", color=c, lw=1.0, alpha=0.9, marker="o", ms=2.5, label=lab)
    mk = "o" if a != "4" else "o"
    ax.plot(x, sub[cell_col].values, "-", color="black", lw=2.0, marker=mk, ms=4.5,
            mfc=("white" if a == "4" else "black"), label="cell-body density (BigBrain)")
    ax.set_title(area_title[a])
    ax.set_xticks(x); ax.set_xticklabels(xt)
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.set_ylim(-1.75, 1.75); ax.set_xlim(-0.35, 2.35)
    if ax is not axesA[0]:
        ax.tick_params(labelleft=False)
axesA[0].set_ylabel("z within area (3 strata)")
fig.text(0.5, 0.975, "A  Receptor classes fall from supra- to infragranular stratum in every area; cell-body density peaks in layer IV in 5 of 6",
         ha="center", va="bottom", fontsize=8)
fig.text(0.5, 0.52, "layer (stratum)", ha="center", va="top", fontsize=8)
h, l = axesA[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=5, frameon=False, fontsize=6.5,
           handlelength=1.6, columnspacing=1.0)

# Row B: GLI by stratum per species
axesB = [fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:6])]
specs = [("PalomeroGallagher_Zilles_2018_TableS3", "44", "Area 44, 7 primate species"),
         ("PalomeroGallagher_Zilles_2018_TableS3", "45", "Area 45, 7 primate species"),
         ("Zilles_etal_1986_Table2", "31", "Posterior cingulate area 31, 17 primates")]
short = {"Homo sapiens": "human", "Pan troglodytes": "chimpanzee", "Pan paniscus": "bonobo", "Gorilla gorilla": "gorilla",
         "Pongo pygmaeus": "orangutan", "Hylobates lar": "gibbon", "Macaca fascicularis": "macaque", "Macaca mulatta": "macaque"}
for ax, (src, area, ttl) in zip(axesB, specs):
    g = gli[(gli.source == src) & (gli.area.astype(str) == area)]
    labels = []
    for _, r in g.iterrows():
        y = r[strata].values.astype(float)
        if r.Species == "Homo sapiens":
            ax.plot(x, y, "-", color="black", lw=2.0, marker="o", ms=4, zorder=3)
            labels.append(("human", y[-1], "black"))
        else:
            ax.plot(x, y, "-", color="0.55", lw=0.8, alpha=0.8, marker="o", ms=2, zorder=2)
    nh = g[g.Species != "Homo sapiens"].copy()
    for _, r in pd.concat([nh.nlargest(1, "infragranular"), nh.nsmallest(1, "infragranular")]).iterrows():
        labels.append((short.get(r.Species, r.Species.split()[0]), float(r["infragranular"]), "0.35"))
    # stagger right-end labels so they do not overlap (minimum spacing in data units)
    labels.sort(key=lambda t: t[1]); yr = g[strata].values.astype(float); step = 0.09 * (yr.max() - yr.min())
    ys = [labels[0][1]]
    for _, yv, _ in labels[1:]:
        ys.append(max(yv, ys[-1] + step))
    for (lab, yv, c), yt in zip(labels, ys):
        ax.annotate(lab, (x[-1], yv), xytext=(x[-1] + 0.12, yt), textcoords="data", va="center", fontsize=6.5, color=c,
                    arrowprops=dict(arrowstyle="-", color="0.7", lw=0.5) if abs(yt - yv) > 1e-9 else None)
    ax.set_title(ttl)
    ax.set_xticks(x); ax.set_xticklabels(xt); ax.set_xlim(-0.3, 3.1)
    ax.margins(y=0.08)
axesB[0].set_ylabel("grey level index (GLI, %)")
for ax in axesB[1:]:
    ax.tick_params(labelleft=False)
ymins = [ax.get_ylim()[0] for ax in axesB[:2]]; ymaxs = [ax.get_ylim()[1] for ax in axesB[:2]]
for ax in axesB[:2]:
    ax.set_ylim(min(ymins), max(ymaxs))
fig.text(0.5, 0.465, "B  Laminar GLI by species: granular stratum densest in most species; species rank orders conserved across strata (44 ρ 0.7–0.9, cingulate ρ ≥ 0.9; 45 only ρ 0.3–0.7)",
         ha="center", va="bottom", fontsize=8)
fig.text(0.5, 0.01, "layer (stratum)", ha="center", va="bottom", fontsize=8)

os.makedirs("figs/s1e", exist_ok=True)
fig.savefig("figs/s1e/laminar_prototype.png", dpi=300, bbox_inches="tight")
fig.savefig("figs/s1e/laminar_prototype.pdf", bbox_inches="tight")
print("wrote figs/s1e/laminar_prototype.png/.pdf")
