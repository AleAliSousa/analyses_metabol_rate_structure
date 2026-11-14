"""Build a source-specific coverage inventory; never infer parent-region rates.

Run from any directory. --evo-root overrides the local source repository.
This inventory describes the five audited papers, not every Heiss publication.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--evo-root', type=Path, default=Path('/Users/crossmodal/Library/CloudStorage/OneDrive-AllenInstitute/Species/Evo-M1-Trait-Data'))
args = parser.parse_args()
EVO = args.evo_root
OUT = ROOT / 'tables/anatomy'
OUT.mkdir(parents=True, exist_ok=True)
sources = []


def read(path):
    sources.append({'source_file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(name, rows, fields=None):
    with (OUT / name).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# Deliberately distinct compartments: a sampled broad thalamus is not an atlas
# whole-thalamus mask; lentiform is not putamen; cerebellar WM is not cortex.
families = {
    'Cortex': ['Cerebral cortex (global)', 'Frontal cortex', 'Frontomesial cortex', 'Frontolateral cortex', 'Orbital cortex', 'Inferior frontal cortex', 'Superior frontal cortex', 'Premotor cortex', 'Sensorimotor cortex', 'Parietal cortex', 'Inferior parietal cortex', 'Superior parietal cortex', 'Temporal cortex', 'Medial temporal cortex', 'Lateral temporal cortex', 'Inferior temporal cortex', 'Middle temporal cortex', 'Superior temporal cortex', 'Occipital cortex', 'Primary visual cortex', 'Visual association cortex', 'Insular cortex', 'Cingulate gyrus (unsplit)', 'Anterior cingulate cortex', 'Posterior cingulate cortex'],
    'Limbic / forebrain': ['Hippocampus', 'Hippocampal structures (broader)', 'Amygdala', 'Inferior limbic structures (unspecified)', 'Basal forebrain'],
    'Basal ganglia': ['Caudate', 'Putamen', 'Nucleus accumbens', 'Striatum (aggregate)', 'Lentiform nucleus (composite)', 'Pallidum (unsplit)', 'External pallidum (GPe)', 'Internal pallidum (GPi)', 'Subthalamic nucleus'],
    'Thalamus': ['Thalamus (whole atlas ROI)', 'Thalamus (broad sampled ROI)', 'Medial thalamic region', 'Lateral geniculate nucleus', 'Medial geniculate nucleus', 'Ventral anterior nucleus (VA)', 'Ventral lateral nucleus (VL)', 'Centromedian nucleus (CM)'],
    'Brainstem': ['Substantia nigra', 'Red nucleus', 'Inferior colliculus', 'Superior colliculus', 'Midbrain (mixed)', 'Brainstem (unsplit)', 'Pons (mixed sampled ROI)', 'Pons (selected WM reference ROI)', 'Cerebral peduncle'],
    'Cerebellum': ['Cerebellar cortex', 'Cerebellar gray matter (unspecified)', 'Cerebellum (atlas ROI; tissue uncertain)', 'Cerebellar vermis (tissue uncertain)', 'Dentate nucleus', 'Central cerebellar white matter'],
    'Cerebral white matter': ['Centrum semiovale / oval center', 'Internal capsule', 'Frontal white matter', 'Temporal white matter', 'Occipital white matter', 'Frontal basal white matter', 'Occipitoparietal white matter'],
    'Whole brain': ['Whole brain (weighted sampled regions)'],
}
region_family = {region: family for family, regions in families.items() for region in regions}
legend = {
    'A': 'Regional mean available; named atlas ROI or cortical aggregate. Boundaries still require pairwise matching.',
    'S': 'Regional mean available; selected spatial sample, not a complete structure mask.',
    'U': 'Regional mean available; anatomical label or tissue boundary needs clarification.',
    'P': 'Regional means published, but native table not yet transcribed into a checked machine-readable input.',
    'T': 'Region reported in statistical tests only; no regional mean in the audited table.',
    '-': 'No separately reported measurement in the audited tables. Not evidence of biological absence.',
}
datasets = [
    dict(dataset_id='Heiss_2004', label='Heiss 2004', cohort_id='Heiss2004_HC9', n=9, population='Healthy awake controls', measure='Absolute rCMRGlc', units='umol/100 g/min', variants='One primary table', independence='Overlap with older Heiss cohorts and Borghammer controls not established.', method='HRRT ~2.2 mm. Cortical gray-matter lobar aggregation; small structures sampled with MRI-guided spheres.', role='Primary metabolic reference', source='Heiss_etal_2004_TABLE1.csv; Heiss_etal_2004.pdf, Methods and Table 1', caution='Older comparison column is Heiss 1984, not another independent dataset. Do not inherit the 2004 nuclear labels for those older values.'),
    dict(dataset_id='Heiss_1984', label='Heiss 1984', cohort_id='Heiss1984_HC7', n=7, population='Healthy awake men, 28–38 years', measure='Absolute CMRGlc, kinetic and autoradiographic', units='umol/100 g/min', variants='Two estimators in the same seven people', independence='Overlap with Heiss 1985 unknown. Re-tabulation in 2004 reuses these data.', method='CT-guided sampled ROIs. ~8 mm FWHM, 11 mm slices.', role='Candidate historical convergence; native table extraction pending', source='Heiss 1984 Table 3, p.220; https://doi.org/10.1038/jcbfm.1984.30', caution='22 bilateral/midline territories from 43 side-specific sites. Broad thalamus, mixed midbrain and lentiform are not individual nuclei. Native rates must remain separate from 2004 re-scaled comparisons.'),
    dict(dataset_id='Heiss_1985', label='Heiss 1985', cohort_id='Heiss1985_HC4', n=4, population='Four men after ~28 hours on duty with little/no sleep', measure='Awake absolute CMRGlc; asleep percentage change', units='umol/100 g/min; % change from awake', variants='Three nondreamers and one dreamer, kept separate; awake/asleep repeated', independence='States are repeated measurements; subgroups are 3 + 1, not two four-person cohorts. Overlap with 1984 unknown.', method='~8 mm FWHM, 11 mm slices. Regional procedure refers to Heiss 1984.', role='State-sensitivity comparison', source='Heiss_etal_1985_TABLEI.csv; heiss_etal_1985.pdf, pp.362–363', caution='Awake is sleep-deprived waking, not a rested normative reference. Whole-brain row is a weighted aggregate of sampled regions.'),
    dict(dataset_id='Borghammer_2012', label='Borghammer 2012', cohort_id='Borghammer2012_HC11', n=11, population='Controls only: 8 Aarhus + 3 Cologne', measure='Normalized regional CMRGlc', units='Dimensionless reference ratios', variants='WM; global mean (GM); reference cluster. Same subjects.', independence='Some participants reused from Eggers 2009; overlap of the control subset requires verification. Overlap with Heiss 2004 unknown.', method='HRRT ~2 mm. Histological atlas warped to individual MRI for eight subcortical structures. Two WM reference VOIs.', role='Separate relative-metabolism convergence', source='Borghammer_etal_2012_Table2.csv; article Methods, Fig.2 and Table 2; https://doi.org/10.1111/j.1600-0404.2011.01556.x', caution='GM means global mean. Cerebellum is central cerebellar WM. Pooled absolute regional rates cannot be recovered from group ratios. PD rows excluded.'),
    dict(dataset_id='Andersen_2022', label='Andersen 2022', cohort_id='Andersen2022_HC26', n=26, population='11 young + 15 older controls, kept separate', measure='Paired FDG and UCB-J regional normalized uptake', units='Dimensionless gray-matter-normalized values in S1', variants='Two tracers in the same subjects; two age groups. S2 contains age-comparison tests.', independence='Paired tracers, normalization variants and S2 are not independent cohort replications. Cross-paper overlap unverified.', method='Bilateral Hammers N30R83-derived atlas ROIs. Geometric-transfer-matrix partial-volume correction in S1.', role='Paired metabolic–synaptic convergence; additional source, not selected by Heiss authorship', source='Andersen_etal_2022_TableS1.csv; TableS2.csv; reference_tables/region_labels and errata; article Methods', caution='S1: 13 regions with means. HIP occurs only in S2 statistical tests. SN expansion is uncertain. Young CBL row has an unresolved published mean/relative-difference inconsistency.'),
]
records = []


def add(dataset, source_label, region, status, source, note='', availability='mean_available', anatomy='named_region'):
    assert region in region_family, region
    records.append(dict(dataset_id=dataset, family=region_family[region], region=region,
                        source_label=source_label, coverage=status, mean_availability=availability,
                        anatomical_scope=anatomy, notes=note, source=source))


h04 = read(ROOT / 'data_intermediate/heiss_2004_regions.csv')
h04_regions = ['Cerebral cortex (global)', 'Frontal cortex', 'Parietal cortex', 'Temporal cortex', 'Occipital cortex', 'Insular cortex', 'Hippocampus', 'Caudate', 'Putamen', 'Nucleus accumbens', 'Pallidum (unsplit)', 'Medial thalamic region', 'Lateral geniculate nucleus', 'Medial geniculate nucleus', 'Subthalamic nucleus', 'Basal forebrain', 'Amygdala', 'Inferior colliculus', 'Superior colliculus', 'Substantia nigra', 'Red nucleus', 'Cerebellar cortex', 'Cerebellar vermis (tissue uncertain)', 'Dentate nucleus', 'Centrum semiovale / oval center', 'Internal capsule']
assert len(h04) == len(h04_regions) == 26
assert [int(r['source_row']) for r in h04] == list(range(1, 27))
for i, (r, region) in enumerate(zip(h04, h04_regions)):
    status = 'A' if i < 6 else 'S'
    note = 'Shared region name does not establish identical boundaries across papers.'
    if region == 'Medial thalamic region':
        note = 'Nucleus medial thalami is not whole thalamus; identity with mediodorsal nucleus has not been established.'
    if region == 'Cerebellar vermis (tissue uncertain)':
        status, note = 'U', 'Vermis reported separately; complete vermis gray/white-matter extent is not established.'
    if region == 'Pallidum (unsplit)':
        note = 'Pallidal sample without separate GPe/GPi estimates; do not assign it to either segment.'
    add('Heiss_2004', r['heiss_region'], region, status, 'Heiss 2004 Table 1', note, anatomy='cortical_aggregate' if i < 6 else 'spatial_sample')

h84 = [
    ('Cerebellar cortex', 'Cerebellar cortex'), ('Dentate nucleus', 'Dentate nucleus'), ('Vermis', 'Cerebellar vermis (tissue uncertain)'),
    ('Pons', 'Pons (mixed sampled ROI)'), ('Midbrain', 'Midbrain (mixed)'), ('Cerebral peduncle', 'Cerebral peduncle'),
    ('Thalamus', 'Thalamus (broad sampled ROI)'), ('Lentiform nucleus', 'Lentiform nucleus (composite)'), ('Caudate nucleus', 'Caudate'),
    ('Frontal white matter', 'Frontal white matter'), ('Temporal white matter', 'Temporal white matter'), ('Occipital white matter', 'Occipital white matter'),
    ('Oval center', 'Centrum semiovale / oval center'), ('Cingulate gyrus', 'Cingulate gyrus (unsplit)'), ('Hippocampal structures', 'Hippocampal structures (broader)'),
    ('Frontomesial cortex', 'Frontomesial cortex'), ('Frontolateral cortex', 'Frontolateral cortex'), ('Insular cortex', 'Insular cortex'),
    ('Temporal cortex', 'Temporal cortex'), ('Parietal cortex', 'Parietal cortex'), ('Occipital cortex', 'Occipital cortex'), ('Primary visual cortex', 'Primary visual cortex'),
]
for label, region in h84:
    note = 'Native Table 3 means await transcription; the 2004 comparison column does not expand coverage.'
    if label == 'Thalamus':
        note += ' Sampled thalamic volume ~4.0 cm3, not a complete whole-thalamus mask.'
    if label == 'Lentiform nucleus':
        note += ' Composite territory; cannot relabel as isolated putamen.'
    add('Heiss_1984', label, region, 'P', 'Heiss 1984 Table 3, p.220; DOI 10.1038/jcbfm.1984.30', note, 'published_not_transcribed', 'spatial_sample')

h85 = read(EVO / 'Heiss_etal_1985/Heiss_etal_1985_TABLEI.csv')
h85_map = dict(zip(
    ['Orbital cortex','Inf. front. cortex','Sup. front. cortex','Premotor cortex','Sensorimot. cortex','Insular cortex','Inf. temp. cortex','Middle temp. cortex','Sup. temp. cortex','Inf. par. cortex','Sup. par. cortex','Vis. assoc. cortex','Prim. vis. cortex','Inf. limbic struct.','Hippocampus','Cingulate gyrus','Front. bas. white matter','Occ. par. white matter','Caudate nucleus','Lentiform nucleus','Thalamus','Brainstem','Cerebellar gray matter','Whole brain (all regions weighted according to size)'],
    ['Orbital cortex','Inferior frontal cortex','Superior frontal cortex','Premotor cortex','Sensorimotor cortex','Insular cortex','Inferior temporal cortex','Middle temporal cortex','Superior temporal cortex','Inferior parietal cortex','Superior parietal cortex','Visual association cortex','Primary visual cortex','Inferior limbic structures (unspecified)','Hippocampus','Cingulate gyrus (unsplit)','Frontal basal white matter','Occipitoparietal white matter','Caudate','Lentiform nucleus (composite)','Thalamus (broad sampled ROI)','Brainstem (unsplit)','Cerebellar gray matter (unspecified)','Whole brain (weighted sampled regions)']))
assert set(h85_map) == {r['structure_as_published'] for r in h85}
for label, region in h85_map.items():
    uncertain = label in ['Inf. limbic struct.', 'Cerebellar gray matter']
    add('Heiss_1985', label, region, 'U' if uncertain else 'S', 'Heiss 1985 Table I', 'Regional boundaries refer to 1984 methods. Retain n=3 nondreamers and n=1 dreamer separately; sleep values are percent changes.', anatomy='boundary_uncertain' if uncertain else 'spatial_sample_or_sample_aggregate')

borg = read(EVO / 'Borghammer_etal_2012/Borghammer_etal_2012_Table2.csv')
borg_map = {'Putamen':'Putamen','Caudate':'Caudate','GPe':'External pallidum (GPe)','GPi':'Internal pallidum (GPi)','VA':'Ventral anterior nucleus (VA)','VL':'Ventral lateral nucleus (VL)','CM':'Centromedian nucleus (CM)','STN':'Subthalamic nucleus','Pons':'Pons (selected WM reference ROI)','Cerebellum':'Central cerebellar white matter'}
assert set(borg_map) == {r['structure_as_published'] for r in borg if r['group_as_published'] == 'Controls'}
for label, region in borg_map.items():
    wm = label in ['Pons','Cerebellum']
    add('Borghammer_2012', label, region, 'S' if wm else 'A', 'Borghammer 2012 Table 2; Methods and Fig.2', 'Selected white-matter reference VOI, not the complete gross structure.' if wm else 'Separate left/right means and three normalizations in the same controls. No complete thalamus or unsplit pallidum rate.', anatomy='selected_white_matter_sample' if wm else 'atlas_roi')

a1 = read(EVO / 'Andersen_etal_2022/Andersen_etal_2022_TableS1.csv')
a2 = read(EVO / 'Andersen_etal_2022/Andersen_etal_2022_TableS2.csv')
read(EVO / 'Andersen_etal_2022/reference_tables/Andersen_etal_2022_region_labels.csv')
read(EVO / 'Andersen_etal_2022/reference_tables/Andersen_etal_2022_errata.csv')
amap = {'FCX':'Frontal cortex','Amg':'Amygdala','MTCx':'Medial temporal cortex','LTCx':'Lateral temporal cortex','PCX':'Parietal cortex','OCX':'Occipital cortex','STR':'Striatum (aggregate)','SN':'Substantia nigra','ICX':'Insular cortex','ACC':'Anterior cingulate cortex','PCCx':'Posterior cingulate cortex','CBL':'Cerebellum (atlas ROI; tissue uncertain)','Thal':'Thalamus (whole atlas ROI)'}
assert set(amap) == {r['region_as_published'] for r in a1}
for label, region in amap.items():
    note = 'Paired tracers use the same bilateral atlas ROI. Cross-paper boundary equivalence requires checking.'
    if label == 'SN':
        note = 'Source label SN has no explicit definition in paper/supplement; substantia nigra is a conventional interpretation only.'
    if label == 'CBL':
        note = 'Atlas cerebellum tissue composition not resolved here; do not equate to cortex or central WM. Young CBL printed means/relative difference have an unresolved inconsistency.'
    add('Andersen_2022', label, region, 'U' if label in ['SN','CBL'] else 'A', 'Andersen 2022 Table S1; region-label reference and errata', note, anatomy='boundary_or_label_uncertain' if label in ['SN','CBL'] else 'atlas_roi')
assert {r['region_as_published'].lower() for r in a2} - {x.lower() for x in amap} == {'hip'}
add('Andersen_2022', 'HIP', 'Hippocampus', 'T', 'Andersen 2022 Table S2', 'Age-comparison p-values only; no HIP regional mean in S1. Do not count as an available metabolic rate.', 'tests_only', 'atlas_roi')

# Complete rectangular product, with explicit absences and original source labels.
assert len({(r['dataset_id'], r['region']) for r in records}) == len(records)
assert all(r['coverage'] in legend for r in records)
regions = [r for r in region_family if any(x['region'] == r for x in records)]
lookup = {(r['dataset_id'],r['region']):r for r in records}
matrix = [dict(dataset_id=d['dataset_id'], **{r: lookup.get((d['dataset_id'],r), {}).get('coverage','-') for r in regions}) for d in datasets]
long = [lookup.get((d['dataset_id'],r), dict(dataset_id=d['dataset_id'],family=region_family[r],region=r,source_label='',coverage='-',mean_availability='not_reported',anatomical_scope='not_reported',notes='No separate row in audited tables; parent/child regions are not substituted.',source=d['source'])) for d in datasets for r in regions]
for d in datasets:
    subset = [r for r in records if r['dataset_id'] == d['dataset_id']]
    d['reported_territories'] = len(subset)
    d['means_available_territories'] = sum(r['mean_availability'] == 'mean_available' for r in subset)
assert [d['reported_territories'] for d in datasets] == [26,22,24,10,14]
assert [d['means_available_territories'] for d in datasets] == [26,0,24,10,13]
assert lookup[('Borghammer_2012','Central cerebellar white matter')]['source_label'] == 'Cerebellum'
assert ('Heiss_1984','Putamen') not in lookup
assert ('Heiss_2004','Thalamus (whole atlas ROI)') not in lookup
write_csv('metabolic_dataset_region_coverage.csv', matrix)
write_csv('metabolic_dataset_region_evidence.csv', long)
write_csv('metabolic_dataset_inventory.csv', datasets)
payload = dict(as_of='2026-11-12', legend=legend, regions=[dict(region=r,family=region_family[r]) for r in regions], datasets=datasets, evidence=records, sources=sources, matrix=matrix)
json_path = ROOT / 'data_intermediate/metabolic_dataset_region_coverage.json'
json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+'\n')
print(json.dumps({'datasets':len(datasets),'regions':len(regions),'reported_cells':len(records),'matrix_cells':len(long),'reported_counts':{d['label']:d['reported_territories'] for d in datasets}}, indent=2))
