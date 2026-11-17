"""
Study s1e_3 (PROTOTYPE) figures — human vs macaque (vs rat) region groups: receptor-class z vs CMRglc z.

Reads  tables/s1e/species_group_panel.csv and tables/s1e/species_group_spearman.csv
       (written by scripts/s1e_3_human_macaque_groups.R; the old human_macaque_* names are kept as copies)
Writes figs/s1e/human_macaque_panels.png/.pdf       (rows human, macaque, rat; 4 receptor classes)
       figs/s1e/human_macaque_concordance.png/.pdf  (human z vs macaque z per class + CMRglc, identity line)
       figs/s1e/cmrglc_concordance_three_species.png/.pdf  (pairwise CMRglc z: human-macaque, human-rat, macaque-rat)
Rat (added 2026-11-17): CMRglc = Sokoloff 1977 Table 3 (2-DG, conscious, 30 min, n = 10; 26 groups); receptors = the 10
EBRAINS rat items (Oc1B/Oc1M -> primary_visual, Oc2L/Oc2M -> visual_association, CA1-3 + DG -> hippocampus_dentate),
i.e. 3 receptor groups, of which visual_association has no rat CMRglc and the hippocampal items lack 5 receptors.

All plotted values are WITHIN-SPECIES z-scores (over the groups of that species that have the variable);
raw receptor units differ between species (fmol/mg protein vs fmol/neuron) and are never plotted together.
Descriptive only; not part of run_all.R / output_manifest.csv.
Run from anywhere inside the repo:  python scripts/s1e_3_fig_human_macaque.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = os.path.abspath(os.getcwd())
while not os.path.exists(os.path.join(d, ".git")) and os.path.dirname(d) != d:
    d = os.path.dirname(d)
os.chdir(d)

# figure-style: use the kernel helper when present, else an equivalent rcParams ladder
if "apply_figure_style" in globals():
    apply_figure_style(sizes=(8, 7, 6))
else:
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                         "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.direction": "out", "ytick.direction": "out", "legend.frameon": False,
                         "savefig.dpi": 300, "pdf.fonttype": 42})

panel = pd.read_csv("tables/s1e/species_group_panel.csv")
sp = pd.read_csv("tables/s1e/species_group_spearman.csv")
os.makedirs("figs/s1e", exist_ok=True)

classes = [("class_iGlu", "iGlu (AMPA+NMDA+kainate)"),
           ("class_GABA", "GABA (A, A/BZ, B)"),
           ("class_ACh", "ACh (M1+M2+M3)"),
           ("class_monoamineGPCR", "Monoamine GPCR (α1, α2, 5-HT1A/2, D1)")]
short = {"primary_visual": "V1", "visual_association": "VisAssoc", "primary_somatosensory_motor": "S1/M1",
         "premotor": "PM", "dorsolateral_prefrontal": "DLPFC", "orbitofrontal_ventral_prefrontal": "OFC",
         "anterior_cingulate": "ACC", "posterior_cingulate_retrosplenial": "PCC", "inferior_parietal": "IPL",
         "superior_parietal": "SPL", "inferior_frontal_broca_homolog": "IFG", "frontal_eye_field_area8": "FEF/8",
         "frontal_pole": "FP", "superior_temporal_association": "STG", "inferior_middle_temporal": "ITG/MTG",
         "temporal_pole_perirhinal_entorhinal": "TP/PRh", "fusiform_ventral_occipitotemporal": "FG",
         "primary_auditory": "A1", "intraparietal_sulcus": "IPS",
         # subcortical / other groups (three-species CMRglc figure)
         "hippocampus_dentate": "Hipp", "amygdala": "Amyg", "caudate": "Caud", "putamen": "Put", "pallidum": "GP",
         "nucleus_accumbens": "NAcc", "basal_forebrain_septum": "Sept/BF", "thalamus_mediodorsal": "MD thal",
         "lateral_geniculate": "LGN", "medial_geniculate": "MGN", "hypothalamus": "Hypoth", "mammillary_body": "Mamm",
         "superior_colliculus": "SC", "inferior_colliculus": "IC", "substantia_nigra": "SN", "subthalamic_nucleus": "STN",
         "brainstem_nuclei": "BS nuclei", "cerebellar_cortex": "Cb ctx", "cerebellar_nuclei": "Cb nuclei",
         "white_matter_cerebral": "WM cereb", "white_matter_brainstem_cerebellar_spinal": "WM bs/cb",
         "frontal_lobe": "Frontal", "parietal_lobe": "Parietal", "temporal_lobe": "Temporal", "occipital_lobe": "Occipital",
         "striatum_undivided": "CPu", "olfactory_cortex": "Olf ctx", "thalamus_lateral_nuclei": "Lat thal",
         "thalamus_ventral_nuclear_group": "Vent thal", "red_nucleus": "RN", "optic_pathway": "Optic"}
lab = lambda g: short.get(g, g)

COL_SHARED, COL_OTHER = "#B2182B", "#6E6E6E"   # focal (shared groups) vs other groups of that species
def sp_lookup(species, var, basis="all_groups_with_data"):
    r = sp[(sp.analysis == "within_species") & (sp.species == species) & (sp.variable == var) & (sp.basis == basis)]
    return (int(r.n.iloc[0]), float(r.rho.iloc[0])) if len(r) and np.isfinite(r.rho.iloc[0]) else (int(r.n.iloc[0]) if len(r) else 0, np.nan)

try:
    from adjustText import adjust_text
except ImportError:  # fallback: fixed offsets
    adjust_text = None

def label_points(ax, x, y, names, color, fs=6):
    """Direct labels; collected per axes and repelled once in finish_labels()."""
    store = ax.__dict__.setdefault("_s1e_texts", [])
    for xi, yi, nm in zip(x, y, names):
        store.append(ax.text(xi, yi, nm, fontsize=fs, color=color, zorder=4))

def finish_labels(ax):
    texts = ax.__dict__.get("_s1e_texts", [])
    if not texts: return
    if adjust_text is not None:
        adjust_text(texts, ax=ax, expand=(1.3, 1.6), force_text=(0.4, 0.8), force_points=(0.3, 0.5),
                    arrowprops=dict(arrowstyle="-", color="0.6", lw=0.4), only_move={"points": "xy", "text": "xy"})
    else:
        for k, t in enumerate(texts):
            t.set_position((t.get_position()[0], t.get_position()[1])); t.set_ha("left"); t.set_va("bottom" if k % 2 else "top")

## ---- Figure 1: 3 x 4 panels (rows = species incl. rat, cols = receptor classes) ----
COL_RAT = "#2166AC"
units = {"human": "FDG PET (Heiss 1985)", "macaque": "2-DG (Kennedy 1978, Palombo 1990)", "rat": "2-DG (Sokoloff 1977, 30 min)"}
fig, axes = plt.subplots(3, 4, figsize=(7.2, 6.6), sharex=False, sharey="row")
for r, species in enumerate(["human", "macaque", "rat"]):
    sub = panel[(panel.species == species) & panel.has_receptors & panel.has_cmrglc].copy()
    rec_only = panel[(panel.species == species) & panel.has_receptors]        # rat: 3 receptor groups, 2 with CMRglc
    for c, (var, title) in enumerate(classes):
        ax = axes[r, c]
        x = sub[var + "_z"].to_numpy(); y = sub["cmrglc_z"].to_numpy(); shared = sub.shared_group_primary.to_numpy()
        names = [lab(g) for g in sub.group_id]
        if species == "rat":
            ok = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[ok], y[ok], s=22, color=COL_RAT, zorder=3)
            label_points(ax, x[ok], y[ok], [n for n, o in zip(names, ok) if o], COL_RAT)
            n_rec = int(len(rec_only)); n_cls = int(rec_only[var].notna().sum()); n_xy = int(ok.sum())
            if var == "class_iGlu":
                txt = (f"rat receptor groups: n = {n_rec}\n(V1, VisAssoc, Hipp)\nwith CMRglc: n = {n_xy}\n"
                       f"(no rat CMRglc for VisAssoc)\nz over the 3 groups; no ρ")
            else:
                txt = (f"rat receptor groups: n = {n_rec}\nclass complete in {n_cls} (V1, VisAssoc);\nhippocampal items lack\n"
                       f"GABAA/BZ, GABAB, M3, D1, α4β2\n→ no z (n < 3), no ρ")
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=5, color=COL_RAT, linespacing=1.1)
            if not ok.any():
                ax.set_xlim(-2, 2); ax.set_ylim(-2, 2)
        else:
            ax.scatter(x[~shared], y[~shared], s=18, facecolor="white", edgecolor=COL_OTHER, linewidth=0.8, zorder=2)
            ax.scatter(x[shared], y[shared], s=22, color=COL_SHARED, zorder=3)
            label_points(ax, x[~shared], y[~shared], [n for n, s in zip(names, shared) if not s], COL_OTHER)
            label_points(ax, x[shared], y[shared], [n for n, s in zip(names, shared) if s], COL_SHARED)
            n, rho = sp_lookup(species, var)
            ax.text(0.03, 0.95, f"ρ = {rho:.2f}, n = {n}", transform=ax.transAxes, ha="left", va="top", fontsize=7)
        ax.axhline(0, color="0.85", lw=0.6, zorder=1); ax.axvline(0, color="0.85", lw=0.6, zorder=1)
        ax.margins(0.15)
        if r == 0:
            ax.set_title(title, fontsize=8, loc="left")
        if c == 0:
            ax.set_ylabel(f"{species.capitalize()} CMRglc, z\n{units[species]}")
        else:
            ax.tick_params(labelleft=False)
    lo, hi = axes[r, 0].get_ylim(); axes[r, 0].set_ylim(lo, hi + (0.22 if species != "rat" else 0.9) * (hi - lo))   # headroom for the text
    for ax in axes[r]: finish_labels(ax)
fig.supxlabel("Receptor-class density (within-species z)", fontsize=8)
fig.text(0.01, 0.99, "Filled red = the 6 region groups with receptor and CMRglc data in human and macaque; open grey = groups with data in one species only;\n"
         "blue = rat (EBRAINS rat receptor items; Sokoloff 1977 2-DG). z-scores are computed within species over all groups having the variable;\n"
         "ρ = Spearman over all groups of that species (descriptive). Rat receptor data cover 3 groups only (Oc1, Oc2, hippocampus).",
         fontsize=7, va="top")
fig.tight_layout(rect=(0, 0.01, 1, 0.94))
fig.savefig("figs/s1e/human_macaque_panels.png"); fig.savefig("figs/s1e/human_macaque_panels.pdf")

## ---- Figure 2: concordance human z vs macaque z (shared primary groups) ----
ph = panel[(panel.species == "human")].set_index("group_id"); pm = panel[(panel.species == "macaque")].set_index("group_id")
shared_ids = list(panel.loc[panel.shared_group_primary, "group_id"].unique())
ext_ids = [g for g in panel.loc[panel.shared_group_extended, "group_id"].unique() if g not in shared_ids]
rec_both = [g for g in panel.loc[panel.shared_receptors_both_species, "group_id"].unique() if g not in shared_ids]
vars2 = classes + [("EI_receptor_index", "E/I receptor index (iGlu / GABA)"), ("cmrglc", "CMRglc")]
fig2, axes2 = plt.subplots(2, 3, figsize=(7.2, 5.2))
for ax, (var, title) in zip(axes2.ravel(), vars2):
    zv = var + "_z"
    def pts(ids):
        return (ph.loc[ids, zv].to_numpy(dtype=float), pm.loc[ids, zv].to_numpy(dtype=float))
    if var != "cmrglc":  # receptor classes exist for 13 groups in both species: show the non-shared ones open grey
        xo, yo = pts(rec_both)
        ok = np.isfinite(xo) & np.isfinite(yo)
        ax.scatter(xo[ok], yo[ok], s=16, facecolor="white", edgecolor=COL_OTHER, linewidth=0.8, zorder=2)
        label_points(ax, xo[ok], yo[ok], [lab(g) for g, o in zip(rec_both, ok) if o], COL_OTHER)
    xs, ys = pts(shared_ids)
    ax.scatter(xs, ys, s=22, color=COL_SHARED, zorder=3)
    label_points(ax, xs, ys, [lab(g) for g in shared_ids], COL_SHARED)
    lim = np.nanmax(np.abs(np.concatenate([xs, ys]))) * 1.35 + 0.2
    ax.plot([-lim, lim], [-lim, lim], color="0.75", lw=0.7, ls="--", zorder=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim * 1.35); ax.set_aspect("equal")
    r = sp[(sp.analysis == "cross_species_concordance") & (sp.basis == "shared_groups_primary") & (sp.variable == zv)]
    rho = float(r.rho.iloc[0]); n = int(r.n.iloc[0])
    if var != "cmrglc":
        r13 = sp[(sp.analysis == "cross_species_concordance") & (sp.basis == "all_groups_with_receptors_in_both_species") & (sp.variable == zv)]
        txt = f"shared: ρ = {rho:.2f} (n = {n})\nall 13 receptor groups: ρ = {float(r13.rho.iloc[0]):.2f}"
    else:
        txt = f"shared: ρ = {rho:.2f} (n = {n})"
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=7)
    finish_labels(ax)
    ax.set_title(title, fontsize=8, loc="left")
    ax.set_xlabel("human (z)"); ax.set_ylabel("macaque (z)")
fig2.text(0.01, 0.99, "Cross-species concordance of within-species z-scores. Filled red = 6 shared groups (receptors + CMRglc in both species);\n"
          "open grey = other groups with receptors in both species (no CMRglc in one species); dashed = identity. ρ = Spearman (descriptive).",
          fontsize=7, va="top")
fig2.tight_layout(rect=(0, 0, 1, 0.94))
fig2.savefig("figs/s1e/human_macaque_concordance.png"); fig2.savefig("figs/s1e/human_macaque_concordance.pdf")
## ---- Figure 3: three-species CMRglc concordance (pairwise within-species z over shared CMRglc groups) ----
pr = panel[(panel.species == "rat")].set_index("group_id")
pairs = [("human", "macaque", ph, pm, "shared_cmrglc_human_macaque"),
         ("human", "rat", ph, pr, "shared_cmrglc_human_rat"),
         ("macaque", "rat", pm, pr, "shared_cmrglc_macaque_rat")]
fig3, axes3 = plt.subplots(1, 3, figsize=(7.2, 3.3))
for ax, (a, b, pa, pb, flag) in zip(axes3, pairs):
    ids = [g for g in panel.loc[panel[flag], "group_id"].unique()]
    x = pa.loc[ids, "cmrglc_z"].to_numpy(dtype=float); y = pb.loc[ids, "cmrglc_z"].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y); x, y = x[ok], y[ok]; ids = [g for g, o in zip(ids, ok) if o]
    three = np.array([bool(panel.loc[(panel.group_id == g) & (panel.species == a), "shared_cmrglc_all_three"].iloc[0]) for g in ids])
    low = np.array([str(panel.loc[(panel.group_id == g) & (panel.species == a), "homology_confidence"].iloc[0]).startswith("low") for g in ids])
    ax.scatter(x[~three], y[~three], s=16, facecolor="white", edgecolor=COL_OTHER, linewidth=0.8, zorder=2)
    ax.scatter(x[three], y[three], s=20, color=COL_SHARED, zorder=3)
    ax.scatter(x[low], y[low], s=60, facecolor="none", edgecolor="0.3", linewidth=0.5, zorder=3, marker="s")
    label_points(ax, x, y, [lab(g) for g in ids], "0.25", fs=4.5)
    lim = np.nanmax(np.abs(np.concatenate([x, y]))) * 1.25 + 0.2
    ax.plot([-lim, lim], [-lim, lim], color="0.75", lw=0.7, ls="--", zorder=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim * 1.3); ax.set_aspect("equal")
    r = sp[(sp.analysis == "cross_species_concordance") & (sp.basis == "shared_cmrglc_groups") & (sp.species == f"{a}_vs_{b}") & (sp.variable == "cmrglc_z")]
    r3 = sp[(sp.analysis == "cross_species_concordance") & (sp.basis == "shared_cmrglc_groups_all_three_species") & (sp.species == f"{a}_vs_{b}") & (sp.variable == "cmrglc_z")]
    rx = sp[(sp.analysis == "cross_species_concordance") & (sp.basis == "shared_cmrglc_groups_excl_low_confidence") & (sp.species == f"{a}_vs_{b}") & (sp.variable == "cmrglc_z")]
    ax.text(0.03, 0.97, f"ρ = {float(r.rho.iloc[0]):.2f}, n = {int(r.n.iloc[0])}\n17 groups in all three: ρ = {float(r3.rho.iloc[0]):.2f}\nexcl. low-confidence: ρ = {float(rx.rho.iloc[0]):.2f} (n = {int(rx.n.iloc[0])})",
            transform=ax.transAxes, ha="left", va="top", fontsize=6)
    finish_labels(ax)
    ax.set_title(f"{a} vs {b}", fontsize=8, loc="left")
    ax.set_xlabel(f"{a} CMRglc (z; {units[a].split(' (')[0]})"); ax.set_ylabel(f"{b} CMRglc (z)")
fig3.text(0.01, 0.99, "Regional CMRglc concordance of within-species z-scores over the region groups each pair shares.\n"
          "Human: Heiss 1985/2004 FDG PET; macaque: Kennedy 1978 + Palombo 1990 2-DG; rat: Sokoloff 1977 2-DG (conscious, 30 min).\n"
          "Filled red = the 17 groups with CMRglc in all three species; open grey = pair-only groups; squares = low-confidence (composite) groups; dashed = identity.",
          fontsize=6, va="top")
fig3.tight_layout(rect=(0, 0.02, 1, 0.88))
fig3.savefig("figs/s1e/cmrglc_concordance_three_species.png"); fig3.savefig("figs/s1e/cmrglc_concordance_three_species.pdf")
print("wrote figs/s1e/human_macaque_panels.{png,pdf}, figs/s1e/human_macaque_concordance.{png,pdf} and figs/s1e/cmrglc_concordance_three_species.{png,pdf}")
