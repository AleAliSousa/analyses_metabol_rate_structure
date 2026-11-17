"""Sample public smoothed Johansen Bmax maps in explicit anatomical masks.

Native affine coordinates are honored. Nearest-neighbor reslicing aligns grids,
NOT anatomical template variants. All estimates retain a template-review flag.
Dependencies: numpy, scipy, nibabel, matplotlib (see requirements-atlas.txt).
"""
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import nibabel as nib
from scipy.ndimage import affine_transform
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'johansen-matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir()) / 'johansen-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data_raw'
OUT = ROOT / 'tables/s1c'
FIG = ROOT / 'figs/s1c'
OUT.mkdir(exist_ok=True, parents=True)
FIG.mkdir(exist_ok=True, parents=True)
source_files = []


def load(path):
    source_files.append(path)
    return nib.load(path)


def write(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


mean_img = load(RAW / 'Johansen_2024_atlas/Bmax.mean.MNI152.sm5.nii.gz')
sd_img = load(RAW / 'Johansen_2024_atlas/Bmax.sd.MNI152.sm5.nii.gz')
assert mean_img.shape == sd_img.shape and np.allclose(mean_img.affine, sd_img.affine)
mean = mean_img.get_fdata(dtype=np.float32)
sd = sd_img.get_fdata(dtype=np.float32)
voxel_mm3 = abs(np.linalg.det(mean_img.affine[:3, :3]))
valid = np.isfinite(mean) & (mean != 0)  # Distributed map uses zero outside its support.


def reslice(img, order=0):
    """Target voxel -> world mm -> source voxel; handles axis flips explicitly."""
    transform = np.linalg.inv(img.affine) @ mean_img.affine
    data = np.asanyarray(img.dataobj)
    result = affine_transform(data, transform[:3, :3], transform[:3, 3],
                              output_shape=mean_img.shape, order=order,
                              mode='constant', cval=0, prefilter=False)
    return result


aal_img = load(RAW / 'anatomical_atlases/AAL3v2/AAL3v1_1mm.nii.gz')
massp_img = load(RAW / 'anatomical_atlases/MASSP2/massp2_maxlabel_MNI2009b.nii.gz')
dent_img = load(RAW / 'anatomical_atlases/Diedrichsen/atl-Anatom_space-MNI_dseg.nii')
cit_img = load(RAW / 'anatomical_atlases/CIT168/STH_MNI2009c.nii.gz')
aal = reslice(aal_img).astype(np.int16)
massp = reslice(massp_img).astype(np.int16)
dent = reslice(dent_img).astype(np.int16)
cit = reslice(cit_img, order=1)
assert np.max(cit) <= 1.0001 and np.min(cit) >= 0
aal_labels = {x.findtext('name'): int(x.findtext('index')) for x in
              ET.parse(RAW / 'anatomical_atlases/AAL3v2/AAL3v1_1mm.xml').iter('label')}
massp_labels = {}
for line in (RAW / 'anatomical_atlases/MASSP2/labels.txt').read_text().splitlines()[1:]:
    parts = [x.strip() for x in line.split(':')]
    massp_labels[f'{parts[2]}_{parts[3]}'] = int(parts[0])
dent_labels = {r['name']: int(r['index']) for r in csv.DictReader(
    (RAW / 'anatomical_atlases/Diedrichsen/atl-Anatom.tsv').open(), delimiter='\t')}
specs = []


def add(region, atlas, names, labels, volume, small=True, role='extension', matched=True, note=''):
    ids = [labels[n] for n in names]
    specs.append(dict(anatomy_id=region, mask_atlas=atlas, mask_labels=';'.join(names),
        mask_values=';'.join(map(str,ids)), small_nucleus=small, role=role,
        named_anatomy_match=matched, note=note, mask=np.isin(volume,ids)))


for region, acronym in [('lateral_geniculate_nucleus','lgn'),('medial_geniculate_nucleus','mgn'),
    ('substantia_nigra','sn'),('red_nucleus','rn'),('subthalamic_nucleus','stn'),
    ('inferior_colliculus','ico'),('superior_colliculus','sco'),('internal_capsule','ic')]:
    add(region,'MASSP2_MNI2009b',[acronym+'_l',acronym+'_r'],massp_labels,massp,
        small=acronym!='ic',note='Bilateral union of named masks. MNI2009b-to-Johansen template correspondence is not independently validated.')
add('dentate_nucleus','Diedrichsen_MNI6',['Left_Dentate','Right_Dentate'],dent_labels,dent,
    note='Dentate nucleus only, not dentate gyrus or other deep cerebellar nuclei.')
add('vermis','AAL3_MNI',[f'Vermis_{x}' for x in ['1_2','3','4_5','6','7','8','9','10']],
    aal_labels,aal,small=False,matched=False,
    note='Lobular vermis union. Tissue extent does not establish identity with the Heiss vermis sample; descriptive only.')
# Alternative atlas definitions quantify mask dependence; they are not new regions.
for region, names in [
    ('lateral_geniculate_nucleus',['Thal_LGN_L','Thal_LGN_R']),
    ('medial_geniculate_nucleus',['Thal_MGN_L','Thal_MGN_R']),
    ('substantia_nigra',['SN_pc_L','SN_pc_R','SN_pr_L','SN_pr_R']),
    ('red_nucleus',['Red_N_L','Red_N_R'])]:
    add(region,'AAL3_MNI',names,aal_labels,aal,role='alternative_mask',
        note='Alternative mask definition; SN combines pars compacta and reticulata by voxel count.')
specs.append(dict(anatomy_id='subthalamic_nucleus',mask_atlas='CIT168_MNI2009c',
    mask_labels='STH probability >= 0.5',mask_values='0.5',small_nucleus=True,
    role='alternative_mask',named_anatomy_match=True,
    note='Trilinearly resliced probabilities thresholded at 0.5; threshold is a sensitivity choice.',mask=cit>=0.5))
# Large structures already in Table 2 provide a scale/registration sanity check.
for region, names in [('caudate_nucleus',['Caudate_L','Caudate_R']),
    ('putamen',['Putamen_L','Putamen_R']),('hippocampus',['Hippocampus_L','Hippocampus_R']),
    ('amygdala',['Amygdala_L','Amygdala_R']),('nucleus_accumbens',['N_Acc_L','N_Acc_R'])]:
    add(region,'AAL3_MNI',names,aal_labels,aal,small=False,role='table_overlap_check',
        note='Different masks and 5-mm smoothing: comparison to nonsmoothed Table 2 is a sanity check, not equality validation.')

with (ROOT / 'data_intermediate/heiss_2004_regions.csv').open() as f:
    rates = {r['anatomy_id']: r for r in csv.DictReader(f)}
result = []
fig, axes = plt.subplots(2, 5, figsize=(15, 6.5))
for spec in specs:
    mask = spec['mask']
    n_mask = int(mask.sum())
    assert n_mask > 0, (spec['anatomy_id'], spec['mask_atlas'])
    good = mask & valid
    vox = np.argwhere(mask)
    world = nib.affines.apply_affine(mean_img.affine, vox)
    x = world[:,0]
    values = mean[good]
    sampled_values = mean[tuple(vox.T)]
    sampled_good = np.isfinite(sampled_values) & (sampled_values != 0)
    left_values = sampled_values[(x<0) & sampled_good]
    right_values = sampled_values[(x>0) & sampled_good]
    row = {k:v for k,v in spec.items() if k!='mask'}
    row.update(n_mask_voxels=n_mask, volume_mm3=n_mask*voxel_mm3,
        n_valid_voxels=int(good.sum()), map_coverage_fraction=float(good.sum()/n_mask),
        bmax_pmol_ml=float(values.mean()) if len(values) else '',
        spatial_sd_pmol_ml=float(values.std(ddof=1)) if len(values)>1 else '',
        mean_voxel_between_subject_sd=float(sd[good].mean()) if len(values) else '',
        left_n_voxels=int(np.sum(x<0)),right_n_voxels=int(np.sum(x>0)),
        left_centroid_x=float(x[x<0].mean()) if np.any(x<0) else '',
        right_centroid_x=float(x[x>0].mean()) if np.any(x>0) else '',
        left_bmax_pmol_ml=float(left_values.mean()) if len(left_values) else '',
        right_bmax_pmol_ml=float(right_values.mean()) if len(right_values) else '',
        smoothing_fwhm_mm=5, template_registration_validated=False,
        primary_included=False, rcMRGlc=float(rates[spec['anatomy_id']]['rcmrglc_mean_umol_100g_min']))
    # Deliberate 1-mm translations test location sensitivity; no interpolation
    # of Bmax itself and no wrap-around across image edges.
    shifts=[]
    for axis in range(3):
        for direction in [-1,1]:
            idx=vox.copy(); idx[:,axis]+=direction
            inside=np.all((idx>=0)&(idx<np.array(mean.shape)),axis=1)
            idx=idx[inside]; vals=mean[tuple(idx.T)]
            vals=vals[np.isfinite(vals)&(vals!=0)]
            shifts.append(float(vals.mean()) if len(vals) else np.nan)
    row.update(shift_1mm_mean_min=float(np.nanmin(shifts)),shift_1mm_mean_max=float(np.nanmax(shifts)),
        shift_1mm_max_abs_change=float(np.nanmax(np.abs(np.array(shifts)-row['bmax_pmol_ml']))))
    result.append(row)
    if spec['role']=='extension':
        ax=axes.flat[len([r for r in result if r['role']=='extension'])-1]
        center=np.rint(vox.mean(axis=0)).astype(int)
        z=center[2]
        ax.imshow(mean[:,:,z].T,origin='lower',cmap='gray',vmin=0,vmax=700)
        ax.contour(mask[:,:,z].T,levels=[0.5],colors=['#E9A335'],linewidths=1)
        ax.set_xlim(center[0]-35,center[0]+35); ax.set_ylim(center[1]-30,center[1]+30)
        ax.set_title(spec['anatomy_id'].replace('_',' '),fontsize=9)
        ax.axis('off')
write('johansen_atlas_samples.csv',result)
fig.suptitle('Anatomical mask overlays on the Johansen 5-mm-smoothed map\nNominal MNI alignment; these overlays do not validate template registration',fontsize=12)
fig.tight_layout()
fig.savefig(FIG/'johansen_atlas_mask_qc.png',dpi=160)
plt.close(fig)

current = {r['anatomy_id'] for r in csv.DictReader((ROOT/'tables/s1c/synaptic_density_rcmrglc_matched_regions.csv').open())}
sampled = {r['anatomy_id']:r for r in result if r['role']=='extension'}
reasons={
    'cerebral_cortex':'Global cortical aggregate: volumetric cortical map is visualization-only; no replacement from arbitrarily averaged surface parcels.',
    'medial_thalamic_nuclei':'Heiss medial-thalamic identity is unresolved. Do not substitute mediodorsal or whole thalamus.',
    'basal_forebrain':'Heiss sampled territory is not specified well enough to substitute a cholinergic forebrain mask.',
    'centrum_semiovale':'No source-matched centrum-semiovale mask; whole white matter and internal capsule are different territories.',
}
audit=[]
for key,r in rates.items():
    if key in current: status='published_table_match_retained'; reason='Existing strict Table 2 match.'
    elif key in sampled: status='sampled_exploratory'; reason=sampled[key]['note']
    else: status='not_sampled_anatomy_unresolved'; reason=reasons[key]
    audit.append(dict(anatomy_id=key,region=r['canonical_label'],status=status,reason=reason,
        small_nucleus=sampled[key]['small_nucleus'] if key in sampled else '',
        primary_atlas_included=False))
write('johansen_atlas_coverage_audit.csv',audit)
geometry=[]
for p in source_files:
    img=nib.load(p)
    geometry.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                         shape=list(img.shape),affine=img.affine.tolist()))
(OUT/'johansen_atlas_geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
print(f'Sampled {len(sampled)} missing territories and {len(result)-len(sampled)} alternative/check masks. Four territories remain anatomically unresolved.')
