"""
Validation Script for E:I Ratio Analysis
Checks data quality, classification rates, and result validity
"""
import pandas as pd
import numpy as np
from pathlib import Path

print('='*80)
print('E:I ANALYSIS - DATA VALIDATION')
print('='*80)

# Paths
PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
DATA_DIR = EI_DIR / "data" / "cortical"

validation_results = []

def validate(check_name, condition, message):
    """Record validation result"""
    status = "PASS" if condition else "FAIL"
    validation_results.append((status, check_name, message))
    symbol = "[+]" if condition else "[-]"
    print(f"  [{status}] {symbol} {check_name}: {message}")

print('\n[1/7] Checking file existence and sizes...')
required_files = {
    'WHB-10Xv3_cell_metadata.csv': (600, 700),  # MB
    'WHB_cluster_annotation_term.csv': (0.5, 1),
    'WHB_cluster_to_cluster_annotation_membership.csv': (0.9, 1.1),
    'ei_ratios_whb_cortical_fine_grained.csv': (0.001, 0.01),
    'ei_ratios_aggregated_lobes.csv': (0.0001, 0.001),
    'ei_metabolism_integrated.csv': (0.0001, 0.001),
    'ei_metabolism_statistics.txt': (0.0001, 0.01),
}

for filename, (min_mb, max_mb) in required_files.items():
    file_path = DATA_DIR.parent if 'WHB' in filename and filename.endswith('_cell_metadata.csv') else DATA_DIR
    if 'WHB' in filename and not filename.endswith('_cell_metadata.csv'):
        file_path = EI_DIR
    elif filename.startswith('ei_'):
        file_path = DATA_DIR

    # Find the file
    if filename == 'WHB-10Xv3_cell_metadata.csv':
        file_path = PROJECT_DIR / "cellxgene_integration" / "data" / "cortical" / filename
    elif filename.startswith('WHB_'):
        file_path = EI_DIR / filename
    else:
        file_path = DATA_DIR / filename

    if file_path.exists():
        size_mb = file_path.stat().st_size / (1024**2)
        validate(f"File: {filename}",
                min_mb <= size_mb <= max_mb,
                f"{size_mb:.2f} MB (expected {min_mb}-{max_mb} MB)")
    else:
        validate(f"File: {filename}", False, f"File not found at {file_path}")

print('\n[2/7] Validating cell metadata...')
cell_df = pd.read_csv(PROJECT_DIR / "cellxgene_integration" / "data" / "cortical" / "WHB-10Xv3_cell_metadata.csv")

total_cells = len(cell_df)
validate("Total cells loaded",
        total_cells > 3000000,
        f"{total_cells:,} cells")

required_columns = ['cluster_alias', 'region_of_interest_label', 'anatomical_division_label']
for col in required_columns:
    validate(f"Column: {col}",
            col in cell_df.columns,
            "Present" if col in cell_df.columns else "Missing")

print('\n[3/7] Checking E/I classification rates...')
annotations = pd.read_csv(EI_DIR / 'WHB_cluster_annotation_term.csv')
membership = pd.read_csv(EI_DIR / 'WHB_cluster_to_cluster_annotation_membership.csv')

nt_annot = annotations[annotations['cluster_annotation_term_set_name'] == 'neurotransmitter']
validate("Neurotransmitter annotations found",
        len(nt_annot) > 0,
        f"{len(nt_annot)} annotations")

# Check classification
def classify_ei(nt):
    if pd.isna(nt): return 'Unknown'
    n = str(nt).upper()
    if 'VGLUT' in n: return 'Excitatory'
    if 'GABA' in n or n == 'GLY': return 'Inhibitory'
    return 'Other'

nt_annot['ei_type'] = nt_annot['name'].apply(classify_ei)
excit_count = (nt_annot['ei_type'] == 'Excitatory').sum()
inhib_count = (nt_annot['ei_type'] == 'Inhibitory').sum()

validate("Excitatory annotations",
        excit_count > 0,
        f"{excit_count} excitatory types identified")

validate("Inhibitory annotations",
        inhib_count > 0,
        f"{inhib_count} inhibitory types identified")

# Build cluster mapping
cluster_to_ei = {}
for _, row in nt_annot.iterrows():
    members = membership[membership['cluster_annotation_term_label'] == row['label']]
    for clust_alias in members['cluster_alias'].values:
        cluster_to_ei[clust_alias] = row['ei_type']

# Map to cells
cell_df['ei_type'] = cell_df['cluster_alias'].map(cluster_to_ei).fillna('Unknown')
ei_classified = (cell_df['ei_type'].isin(['Excitatory', 'Inhibitory'])).sum()
classification_rate = ei_classified / len(cell_df) * 100

validate("Cell classification rate",
        classification_rate > 60,
        f"{classification_rate:.1f}% of cells classified as E or I")

print('\n[4/7] Validating cortical region extraction...')
cortical = cell_df[cell_df['anatomical_division_label'] == 'Cerebral cortex']
cortical_pct = len(cortical) / len(cell_df) * 100

validate("Cortical cells extracted",
        len(cortical) > 1000000,
        f"{len(cortical):,} cells ({cortical_pct:.1f}% of total)")

cortical_ei = cortical[cortical['ei_type'].isin(['Excitatory', 'Inhibitory'])]
cortical_ei_pct = len(cortical_ei) / len(cortical) * 100

validate("Cortical E/I classification",
        cortical_ei_pct > 60,
        f"{cortical_ei_pct:.1f}% of cortical cells classified")

print('\n[5/7] Validating E:I ratios...')
ei_df = pd.read_csv(DATA_DIR / 'ei_ratios_whb_cortical_fine_grained.csv')

validate("Fine-grained regions",
        len(ei_df) >= 20,
        f"{len(ei_df)} cortical regions analyzed")

# Check E:I ratio ranges (typical cortical: 2:1 to 6:1)
valid_ratios = ei_df[ei_df['EI_Ratio'].notna()]
min_ratio = valid_ratios['EI_Ratio'].min()
max_ratio = valid_ratios['EI_Ratio'].max()
mean_ratio = valid_ratios['EI_Ratio'].mean()

validate("E:I ratio range",
        1.0 < min_ratio < 8.0 and 2.0 < max_ratio < 10.0,
        f"Range: {min_ratio:.2f} to {max_ratio:.2f}, Mean: {mean_ratio:.2f}")

validate("E:I ratios in expected range",
        (valid_ratios['EI_Ratio'] > 1.5).sum() / len(valid_ratios) > 0.8,
        f"{(valid_ratios['EI_Ratio'] > 1.5).sum()}/{len(valid_ratios)} ratios > 1.5")

# Check cell counts
validate("Minimum cell count per region",
        ei_df['Total'].min() > 100,
        f"Smallest region: {ei_df['Total'].min():,} cells")

print('\n[6/7] Validating lobe aggregation...')
lobe_df = pd.read_csv(DATA_DIR / 'ei_ratios_aggregated_lobes.csv')

validate("Lobes identified",
        len(lobe_df) >= 5,
        f"{len(lobe_df)} lobes")

# Check aggregation logic: sum E, sum I, then divide
for _, row in lobe_df.iterrows():
    calculated_ratio = row['Excitatory'] / row['Inhibitory'] if row['Inhibitory'] > 0 else np.nan
    validate(f"Lobe aggregation: {row['Lobe']}",
            abs(calculated_ratio - row['EI_Ratio_aggregated']) < 0.01,
            f"E:I = {row['EI_Ratio_aggregated']:.2f} (verified)")

print('\n[7/7] Validating rCMRGlc integration...')
integrated = pd.read_csv(DATA_DIR / 'ei_metabolism_integrated.csv')

validate("Lobes with rCMRGlc data",
        len(integrated) >= 5,
        f"{len(integrated)} lobes integrated")

# Check that rCMRGlc values are reasonable (20-50 umol/100g/min)
valid_met = integrated[integrated['rCMRGlc_Mean'].notna()]
validate("rCMRGlc values in range",
        (valid_met['rCMRGlc_Mean'] > 20).all() and (valid_met['rCMRGlc_Mean'] < 50).all(),
        f"Range: {valid_met['rCMRGlc_Mean'].min():.1f} to {valid_met['rCMRGlc_Mean'].max():.1f} umol/100g/min")

# Check correlation statistics
stats_file = DATA_DIR / 'ei_metabolism_statistics.txt'
if stats_file.exists():
    with open(stats_file) as f:
        stats_content = f.read()
    validate("Statistics file created",
            'Pearson' in stats_content and 'Spearman' in stats_content,
            "Contains correlation statistics")
else:
    validate("Statistics file created", False, "File not found")

# Summary
print('\n' + '='*80)
print('VALIDATION SUMMARY')
print('='*80)

passed = sum(1 for status, _, _ in validation_results if status == "PASS")
failed = sum(1 for status, _, _ in validation_results if status == "FAIL")
total = len(validation_results)

print(f'\nTotal checks: {total}')
print(f'Passed: {passed} ({passed/total*100:.1f}%)')
print(f'Failed: {failed} ({failed/total*100:.1f}%)')

if failed > 0:
    print('\nFailed checks:')
    for status, check, msg in validation_results:
        if status == "FAIL":
            print(f'  [-] {check}: {msg}')

# Save validation report
report_path = DATA_DIR / 'validation_report.txt'
with open(report_path, 'w') as f:
    f.write('E:I ANALYSIS - VALIDATION REPORT\n')
    f.write('='*80 + '\n\n')
    f.write(f'Generated: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
    f.write(f'Total checks: {total}\n')
    f.write(f'Passed: {passed} ({passed/total*100:.1f}%)\n')
    f.write(f'Failed: {failed} ({failed/total*100:.1f}%)\n\n')
    f.write('DETAILED RESULTS:\n')
    f.write('='*80 + '\n\n')
    for status, check, msg in validation_results:
        symbol = "[+]" if status == "PASS" else "[-]"
        f.write(f'[{status}] {symbol} {check}\n    {msg}\n\n')

print(f'\nValidation report saved to: {report_path}')
print('='*80)
