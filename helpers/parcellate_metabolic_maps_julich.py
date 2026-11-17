"""s1e_0 prototype: parcellate cortical metabolic maps by Julich-Brain v3.1 areas.

Inputs (all under data_raw/, already downloaded; checksums in
data_raw/neuromaps_metabolic/neuromaps_download_manifest.csv):
  * Raichle-lab group maps (Vaishnavi et al. 2010 PNAS) distributed by neuromaps
    as fsLR-164k surface maps: cmrglc, cmr02, cbf, cbv (left/right).
  * Castrillon et al. 2023 Sci Adv group CMRglc, MNI152 3 mm volume (neuromaps).
  * Julich-Brain v3.1 maximum-probability maps: fsaverage-164k label GIfTIs
    (labels 1-146 cortical areas + GapMaps, 147 CorpusCallosumMask,
    148 BrainStemMask) and MNI152 (ICBM 2009c asym) 1 mm label volumes.

Processing
  1. fsLR-164k -> fsaverage-164k. If Connectome Workbench (`wb_command`) is on
     PATH, neuromaps.transforms.fslr_to_fsaverage is used (ADAP_BARY_AREA with
     the fsLR medial-wall ROI, fsaverage medial wall masked). Otherwise a
     pure-Python barycentric resampler is used on the same neuromaps template
     spheres (tpl-fsLR_space-fsaverage_den-164k sphere -> tpl-fsaverage_den-164k
     sphere): each fsaverage vertex is located in a triangle of the fsLR sphere
     deformed to fsaverage, data are interpolated with barycentric weights,
     source vertices outside the fsLR no-medial-wall ROI get zero weight, and
     the fsaverage no-medial-wall mask is applied to the output. For equal
     164k densities this is the barycentric step of ADAP_BARY_AREA without the
     vertex-area adaptation. The manifest records which path was used.
  2. Vertex-wise glycolytic-index PROXY: residual of cmrglc regressed on cmr02
     (ordinary least squares over all valid cortical vertices, both hemispheres
     pooled). Vaishnavi et al. 2010 defined GI as the residual of CMRGlc
     regressed on CMRO2 on quantile-normalised images in the volume; the maps
     here are the distributed group surface maps in neuromaps' units, so this
     is a proxy with the same sign convention (positive = more glycolysis than
     expected from oxygen consumption), not the published GI.
  3. Parcellation on fsaverage per hemisphere: mean over vertices with finite,
     non-zero data; n_vertices (all vertices with the label) and the fraction
     of vertices excluded as NaN/zero are reported. Labels 0, 147, 148 excluded.
  4. Castrillon CMRglc: Julich MNI label volumes (left, right, merged) are
     resampled to the 3 mm data grid by nearest neighbour through the affines
     (target voxel -> world mm -> source voxel; no template re-registration);
     mean over voxels with finite non-zero data; n_voxels reported.

Outputs
  data_intermediate/s1e_julich_area_metabolism.csv
  data_intermediate/s1e_julich_area_metabolism_manifest.csv
  data_intermediate/s1e_fsaverage164k_metabolic_maps.npz  (resampled vertex data)

Dependencies: numpy, scipy, pandas, nibabel, neuromaps (templates only).
Environment variables: NEUROMAPS_DATA (template cache; defaults to
data_raw/neuromaps_metabolic/templates), S1E_FORCE_PYTHON_RESAMPLER=1 forces the
pure-Python resampler even when wb_command is available.
"""
import csv
import hashlib
import os
import platform
import shutil
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
from scipy.ndimage import affine_transform
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW_MET = ROOT / 'data_raw/neuromaps_metabolic'
RAW_JB = ROOT / 'data_raw/anatomical_atlases/JulichBrain_v31'
OUT = ROOT / 'data_intermediate'
OUT.mkdir(exist_ok=True, parents=True)
os.environ.setdefault('NEUROMAPS_DATA', str(RAW_MET / 'templates'))
Path(os.environ['NEUROMAPS_DATA']).mkdir(exist_ok=True, parents=True)

import neuromaps  # noqa: E402  (after NEUROMAPS_DATA is set)
from neuromaps.datasets import fetch_atlas  # noqa: E402

MEASURES = ['cmrglc', 'cmr02', 'cbf', 'cbv']
EXCLUDE_LABELS = {0, 147, 148}
HEMIS = {'L': 'left', 'R': 'right'}
HEISS_CORTEX_UMOL = 33.5  # Heiss et al. 2004 global cerebral cortex rCMRglc, umol/100 g/min

inputs = []  # (role, path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load(path, role):
    inputs.append((role, Path(path)))
    return nib.load(str(path))


def gifti_data(img, i=0):
    return np.asarray(img.darrays[i].data)


# --------------------------------------------------------------------------
# 1. fsLR-164k -> fsaverage-164k
# --------------------------------------------------------------------------
fslr = fetch_atlas('fsLR', '164k')
fsavg = fetch_atlas('fsaverage', '164k')
atlas_dir = Path(os.environ['NEUROMAPS_DATA']) / 'atlases'


def sphere_paths(hemi):
    src = atlas_dir / 'fsLR' / f'tpl-fsLR_space-fsaverage_den-164k_hemi-{hemi}_sphere.surf.gii'
    trg = Path(getattr(fsavg['sphere'], hemi))
    src_roi = Path(getattr(fslr['medial'], hemi))
    trg_roi = Path(getattr(fsavg['medial'], hemi))
    return src, trg, src_roi, trg_roi


def load_surf(path, role):
    g = load(path, role)
    return gifti_data(g, 0).astype(np.float64), gifti_data(g, 1).astype(np.int64)


def barycentric_weights(src_sphere, trg_sphere, hemi):
    """Locate each target vertex in a triangle of the source sphere.

    Returns (idx, W, inside): source vertex indices (nT,3), weights (nT,3),
    and a boolean saying whether a containing triangle was found (otherwise
    nearest-neighbour fallback with weight 1).
    """
    S, F = load_surf(src_sphere, f'template sphere fsLR deformed to fsaverage {hemi}')
    T, _ = load_surf(trg_sphere, f'template sphere fsaverage {hemi}')
    S /= np.linalg.norm(S, axis=1, keepdims=True)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    nv = len(S)
    counts = np.zeros(nv, int)
    for col in range(3):
        np.add.at(counts, F[:, col], 1)
    maxdeg = counts.max()
    inc = -np.ones((nv, maxdeg), int)
    fill = np.zeros(nv, int)
    for ti, (a, b, c) in enumerate(F):
        for v in (a, b, c):
            inc[v, fill[v]] = ti
            fill[v] += 1
    # candidate triangles: those incident to the 3 nearest source vertices
    _, nn3 = cKDTree(S).query(T, k=3)
    nT = len(T)
    best_tri = -np.ones(nT, int)
    best_w = np.zeros((nT, 3))
    best_min = np.full(nT, -np.inf)
    for j in range(3):
        nn = nn3[:, j]
        for k in range(maxdeg):
            tri = inc[nn, k]
            ok = tri >= 0
            tri_s = np.where(ok, tri, 0)
            A, B, C = S[F[tri_s, 0]], S[F[tri_s, 1]], S[F[tri_s, 2]]
            n = np.cross(B - A, C - A)
            denom = np.einsum('ij,ij->i', T, n)
            denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
            s = np.einsum('ij,ij->i', A, n) / denom  # ray from origin through T hits plane at s*T
            P = T * s[:, None]
            v0, v1, v2 = B - A, C - A, P - A
            d00 = np.einsum('ij,ij->i', v0, v0)
            d01 = np.einsum('ij,ij->i', v0, v1)
            d11 = np.einsum('ij,ij->i', v1, v1)
            d20 = np.einsum('ij,ij->i', v2, v0)
            d21 = np.einsum('ij,ij->i', v2, v1)
            den = d00 * d11 - d01 * d01
            den = np.where(np.abs(den) < 1e-18, 1e-18, den)
            wb = (d11 * d20 - d01 * d21) / den
            wc = (d00 * d21 - d01 * d20) / den
            wa = 1 - wb - wc
            W = np.stack([wa, wb, wc], 1)
            mn = np.where(ok & (s > 0), W.min(1), -np.inf)
            better = mn > best_min
            best_min = np.where(better, mn, best_min)
            best_tri = np.where(better, tri_s, best_tri)
            best_w = np.where(better[:, None], W, best_w)
    inside = best_min >= -1e-6
    idx = F[best_tri]
    W = np.clip(best_w, 0, None)
    W /= W.sum(1, keepdims=True)
    idx[~inside] = nn3[~inside, 0][:, None]
    W[~inside] = [1.0, 0.0, 0.0]
    return idx, W, inside


def python_resample(data, idx, W, src_roi, trg_roi):
    w = W * src_roi[idx]
    ws = w.sum(1)
    out = np.where(ws > 0, (data[idx] * w).sum(1) / np.where(ws > 0, ws, 1.0), 0.0)
    out[~trg_roi] = 0.0
    return out


use_wb = shutil.which('wb_command') is not None and os.environ.get('S1E_FORCE_PYTHON_RESAMPLER', '0') != '1'
resample_method = ('neuromaps.transforms.fslr_to_fsaverage (wb_command -metric-resample ADAP_BARY_AREA, '
                   'fsLR medial-wall ROI, fsaverage medial wall masked)' if use_wb else
                   'pure-Python barycentric interpolation on neuromaps template spheres '
                   '(tpl-fsLR_space-fsaverage_den-164k -> tpl-fsaverage_den-164k), fsLR medial-wall '
                   'ROI weights, fsaverage medial wall masked; wb_command not available')

surf = {}  # (measure, hemi) -> fsaverage-164k vertex array
resample_notes = {}
for hemi in HEMIS:
    src_sphere, trg_sphere, src_roi_p, trg_roi_p = sphere_paths(hemi)
    src_roi = gifti_data(load(src_roi_p, f'fsLR no-medial-wall ROI {hemi}')).astype(bool)
    trg_roi = gifti_data(load(trg_roi_p, f'fsaverage no-medial-wall ROI {hemi}')).astype(bool)
    if not use_wb:
        idx, W, inside = barycentric_weights(src_sphere, trg_sphere, hemi)
        resample_notes[hemi] = f'{(~inside).sum()} of {len(inside)} fsaverage vertices used nearest-neighbour fallback'
    for m in MEASURES:
        p = RAW_MET / f'source-raichle_desc-{m}_space-fsLR_den-164k_hemi-{hemi}_feature.func.gii'
        g = load(p, f'raichle {m} fsLR-164k {hemi}')
        d = gifti_data(g).astype(np.float64)
        if use_wb:
            from neuromaps import transforms
            res, = transforms.fslr_to_fsaverage(str(p), target_density='164k', hemi=hemi)
            out = gifti_data(res).astype(np.float64)
        else:
            out = python_resample(d, idx, W, src_roi.astype(float), trg_roi)
        surf[(m, hemi)] = out

# --------------------------------------------------------------------------
# 2. Glycolytic-index proxy (vertex-wise residual of cmrglc on cmr02)
# --------------------------------------------------------------------------
glc = np.concatenate([surf[('cmrglc', h)] for h in HEMIS])
oxy = np.concatenate([surf[('cmr02', h)] for h in HEMIS])
valid = np.isfinite(glc) & np.isfinite(oxy) & (glc != 0) & (oxy != 0)
slope, intercept = np.polyfit(oxy[valid], glc[valid], 1)
gi = np.zeros_like(glc)
gi[valid] = glc[valid] - (intercept + slope * oxy[valid])
r_glc_oxy = np.corrcoef(oxy[valid], glc[valid])[0, 1]
nL = len(surf[('cmrglc', 'L')])
surf[('gi_proxy', 'L')] = gi[:nL]
surf[('gi_proxy', 'R')] = gi[nL:]
gi_note = (f'gi_proxy = cmrglc - ({intercept:.4f} + {slope:.6f} * cmr02); OLS over {valid.sum()} valid vertices, '
           f'both hemispheres pooled, r(cmr02, cmrglc) = {r_glc_oxy:.4f}; zero where either input is invalid')

# --------------------------------------------------------------------------
# 3. Surface parcellation by Julich-Brain v3.1 fsaverage labels
# --------------------------------------------------------------------------
label_table = None
labels = {}
for hemi, side in HEMIS.items():
    g = load(RAW_JB / f'JulichBrain_v31_fsaverage_labelled_{side}_hemisphere.label.gii', f'Julich v3.1 fsaverage labels {side}')
    labels[hemi] = gifti_data(g).astype(int)
    # right-hemisphere table writes 'Area MFG5 (MFG,IFG)' without the space; normalise to the left/CSV spelling
    lt = {k: ', '.join(part.strip() for part in v.split(',')) for k, v in g.labeltable.get_labels_as_dict().items()}
    if label_table is None:
        label_table = lt
    else:
        assert lt == label_table, 'left/right GIfTI label tables differ beyond comma spacing'
surf_measures = MEASURES + ['gi_proxy']
rows = []
for hemi in HEMIS:
    lab = labels[hemi]
    anyvalid = np.ones(len(lab), bool)
    for m in surf_measures:
        anyvalid &= np.isfinite(surf[(m, hemi)])
    for lid, name in sorted(label_table.items()):
        if lid in EXCLUDE_LABELS:
            continue
        mask = lab == lid
        n = int(mask.sum())
        row = dict(julich_label=lid, julich_region=name, hemisphere=hemi, n_vertices=n)
        for m in surf_measures:
            d = surf[(m, hemi)]
            ok = mask & np.isfinite(d) & (d != 0)
            row[m] = float(d[ok].mean()) if ok.any() else np.nan
            row[f'{m}_sd'] = float(d[ok].std(ddof=1)) if ok.sum() > 1 else np.nan
            if m == 'cmrglc':
                row['n_vertices_valid'] = int(ok.sum())
                row['frac_vertices_excluded'] = float(1 - ok.sum() / n) if n else np.nan
        rows.append(row)
surf_df = pd.DataFrame(rows)

# --------------------------------------------------------------------------
# 4. Castrillon CMRglc parcellation with the Julich MNI label volumes
# --------------------------------------------------------------------------
cas_img = load(RAW_MET / 'source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz', 'castrillon2023 cmrglc MNI152 3mm')
cas = cas_img.get_fdata(dtype=np.float64)
cas_valid = np.isfinite(cas) & (cas != 0)
voxel_mm3 = float(abs(np.linalg.det(cas_img.affine[:3, :3])))


def reslice_labels(img):
    """Target voxel -> world mm -> source voxel; nearest neighbour; no re-registration."""
    transform = np.linalg.inv(img.affine) @ cas_img.affine
    data = np.asanyarray(img.dataobj)
    return affine_transform(data, transform[:3, :3], transform[:3, 3], output_shape=cas_img.shape,
                            order=0, mode='constant', cval=0, prefilter=False).astype(np.int32)


lut = pd.read_csv(RAW_JB / 'JulichBrain_v31_MNI152_labelled_lut.csv')
inputs.append(('Julich v3.1 MNI152 LUT', RAW_JB / 'JulichBrain_v31_MNI152_labelled_lut.csv'))
lut['base'] = lut['region'].str.replace(r'\s+(left|right)$', '', regex=True)
lut['hemi'] = lut['region'].str.extract(r'(left|right)$')[0].map({'left': 'L', 'right': 'R'})


def mni_to_gifti_name(name):
    # 'Frontal-II (GapMap)' -> 'GapMap-Frontal-II'
    if name.endswith(' (GapMap)'):
        return 'GapMap-' + name[:-len(' (GapMap)')]
    return name


lut['gifti_name'] = lut['base'].map(mni_to_gifti_name)
gifti_names = set(label_table.values())
vol = {}
for hemi, side in HEMIS.items():
    vol[hemi] = reslice_labels(load(RAW_JB / f'JulichBrain_v31_MNI152_labelled_{side}_hemisphere.nii.gz', f'Julich v3.1 MNI152 labels {side}'))
vol['B'] = reslice_labels(load(RAW_JB / 'JulichBrain_v31_MNI152_labelled_merged.nii.gz', 'Julich v3.1 MNI152 labels merged'))
assert not ((vol['L'] > 0) & (vol['R'] > 0)).any()
mni_name_to_label = dict(zip(lut['gifti_name'], lut['label']))
vrows = []
for name in gifti_names:
    if name in ('Unknown', 'CorpusCallosumMask', 'BrainStemMask'):
        continue
    lid = mni_name_to_label.get(name)
    for hemi in ['L', 'R', 'B']:
        if lid is None:
            vrows.append(dict(julich_region=name, hemisphere=hemi, n_voxels=0, n_voxels_valid=0, castrillon_cmrglc=np.nan, castrillon_cmrglc_sd=np.nan))
            continue
        mask = vol[hemi] == lid
        ok = mask & cas_valid
        vrows.append(dict(julich_region=name, hemisphere=hemi, mni_label=lid, n_voxels=int(mask.sum()), n_voxels_valid=int(ok.sum()),
                          castrillon_cmrglc=float(cas[ok].mean()) if ok.any() else np.nan,
                          castrillon_cmrglc_sd=float(cas[ok].std(ddof=1)) if ok.sum() > 1 else np.nan))
vol_df = pd.DataFrame(vrows)

# --------------------------------------------------------------------------
# Combine: per hemisphere rows + bilateral row (vertex-weighted means)
# --------------------------------------------------------------------------
bil_rows = []
for (lid, name), g in surf_df.groupby(['julich_label', 'julich_region']):
    row = dict(julich_label=lid, julich_region=name, hemisphere='B', n_vertices=int(g.n_vertices.sum()),
               n_vertices_valid=int(g.n_vertices_valid.sum()))
    row['frac_vertices_excluded'] = 1 - row['n_vertices_valid'] / row['n_vertices'] if row['n_vertices'] else np.nan
    for m in surf_measures:
        w = g.n_vertices_valid.values.astype(float)
        v = g[m].values
        okw = np.isfinite(v) & (w > 0)
        row[m] = float(np.average(v[okw], weights=w[okw])) if okw.any() else np.nan
        row[f'{m}_sd'] = np.nan  # not pooled; see hemisphere rows
    bil_rows.append(row)
combined = pd.concat([surf_df, pd.DataFrame(bil_rows)], ignore_index=True)
combined = combined.merge(vol_df, on=['julich_region', 'hemisphere'], how='left')
combined['frac_voxels_excluded'] = np.where(combined.n_voxels > 0, 1 - combined.n_voxels_valid / combined.n_voxels.replace(0, np.nan), np.nan)
# relative units: ratio to the vertex-weighted cortical mean of the same measure/hemisphere set (labels 1-146)
for m in ['cmrglc', 'cmr02', 'cbf', 'cbv', 'castrillon_cmrglc']:
    for hemi in ['L', 'R', 'B']:
        sel = (combined.hemisphere == hemi) & (combined.julich_label <= 146)
        w = combined.loc[sel, 'n_vertices_valid' if m != 'castrillon_cmrglc' else 'n_voxels_valid'].values.astype(float)
        v = combined.loc[sel, m].values
        ok = np.isfinite(v) & (w > 0)
        cmean = np.average(v[ok], weights=w[ok])
        combined.loc[sel, f'{m}_rel_cortex'] = combined.loc[sel, m] / cmean
        combined.loc[sel, f'{m}_cortex_mean'] = cmean
combined['hemisphere'] = pd.Categorical(combined.hemisphere, ['L', 'R', 'B'])
combined = combined.sort_values(['julich_label', 'hemisphere']).reset_index(drop=True)
combined['hemisphere'] = combined.hemisphere.astype(str)
front = ['julich_label', 'julich_region', 'hemisphere', 'n_vertices', 'n_vertices_valid', 'frac_vertices_excluded',
         'cmrglc', 'cmr02', 'cbf', 'cbv', 'gi_proxy', 'castrillon_cmrglc', 'n_voxels', 'n_voxels_valid', 'frac_voxels_excluded']
combined = combined[front + [c for c in combined.columns if c not in front]]
combined['resample_method'] = 'wb_command' if use_wb else 'python_barycentric'
out_csv = OUT / 's1e_julich_area_metabolism.csv'
combined.to_csv(out_csv, index=False, float_format='%.6g')
np.savez_compressed(OUT / 's1e_fsaverage164k_metabolic_maps.npz',
                    **{f'{m}_{h}': surf[(m, h)].astype(np.float32) for (m, h) in surf},
                    labels_L=labels['L'], labels_R=labels['R'])

# --------------------------------------------------------------------------
# Unit check and manifest
# --------------------------------------------------------------------------
cort = combined[(combined.hemisphere == 'B') & (combined.julich_label <= 146)]
unit_lines = []
for m in ['cmrglc', 'cmr02', 'cbf', 'cbv', 'castrillon_cmrglc']:
    cm = float(cort[f'{m}_cortex_mean'].iloc[0])
    verdict = ('order of magnitude consistent with absolute umol/100 g/min (Heiss 2004 cortex 33.5; group-map cortical mean is lower, partial-volume and cohort differences expected)' if m.endswith('cmrglc') and 15 < cm < 60
               else 'NOT absolute umol/100 g/min; treat as neuromaps-distributed normalised/arbitrary units' if m.endswith('cmrglc')
               else 'units as distributed by neuromaps (not checked against a reference)')
    unit_lines.append(f'{m}: vertex/voxel-weighted cortical mean (Julich labels 1-146, bilateral) = {cm:.4g}; {verdict}')

manifest = []
for role, p in inputs:
    manifest.append(dict(role=role, path=str(p.relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p),
                         bytes=p.stat().st_size, sha256=sha256(p)))
for role, p in [('output', out_csv), ('output', OUT / 's1e_fsaverage164k_metabolic_maps.npz')]:
    manifest.append(dict(role=role, path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha256(p)))
meta = [
    ('software', 'python', platform.python_version()),
    ('software', 'numpy', np.__version__),
    ('software', 'scipy', __import__('scipy').__version__),
    ('software', 'pandas', pd.__version__),
    ('software', 'nibabel', nib.__version__),
    ('software', 'neuromaps', neuromaps.__version__),
    ('software', 'wb_command', shutil.which('wb_command') or 'not found'),
    ('method', 'fslr_to_fsaverage', resample_method),
    ('method', 'gi_proxy', gi_note),
    ('method', 'surface_parcellation', 'mean over vertices with finite non-zero data per Julich label and hemisphere; labels 0/147/148 excluded; bilateral row = valid-vertex-weighted mean of L and R'),
    ('method', 'volume_parcellation', f'Julich MNI152 (ICBM 2009c asym, 1 mm) labels nearest-neighbour resampled to the Castrillon 3 mm grid via affines (voxel {voxel_mm3:.1f} mm3); mean over voxels with finite non-zero data; no template re-registration'),
    ('run', 'date', date.today().isoformat()),
    ('run', 'platform', platform.platform()),
]
for h, note in resample_notes.items():
    meta.append(('method', f'resampler_fallback_{h}', note))
for line in unit_lines:
    meta.append(('units', line.split(':')[0], line.split(': ', 1)[1]))
for role, key, val in meta:
    manifest.append(dict(role=role, path=key, bytes='', sha256=val))
with (OUT / 's1e_julich_area_metabolism_manifest.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['role', 'path', 'bytes', 'sha256'])
    w.writeheader()
    w.writerows(manifest)

print(f'wrote {out_csv} ({len(combined)} rows)')
print(resample_method)
for note in resample_notes.values():
    print(note)
print(gi_note)
for line in unit_lines:
    print(line)
small = cort[cort.n_vertices < 200]
print(f'{len(small)} bilateral regions with < 200 vertices: ' + '; '.join(f"{r.julich_region} ({r.n_vertices})" for r in small.itertuples()))
sys.exit(0)
