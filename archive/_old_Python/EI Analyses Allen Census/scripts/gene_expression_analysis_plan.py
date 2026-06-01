"""
METABOLIC GENE EXPRESSION ANALYSIS - DATA REQUIREMENTS
Note: Full implementation requires downloading expression matrices
"""
import pandas as pd
import numpy as np
from pathlib import Path

print('='*80)
print('METABOLIC GENE EXPRESSION ANALYSIS - DATA ASSESSMENT')
print('='*80)

PROJECT_DIR = Path(r"C:\Users\michaelproulx\MetabolicBrain")
EI_DIR = PROJECT_DIR / "cellxgene_integration" / "ei_analysis"
OUTPUT_DIR = EI_DIR / "data" / "cortical" / "followup_analyses"

print('\nCurrent data available:')
print('  - Cell metadata (3.3M cells): cluster assignments, regions, anatomical divisions')
print('  - Cluster annotations: neurotransmitter, supercluster, cluster, subcluster')
print('  - Cell type classifications: E/I, MGE/CGE interneurons, glia')

print('\nData NOT currently available:')
print('  - Gene expression matrices (count or normalized)')
print('  - Individual gene expression values per cell')

print('\n' + '='*80)
print('PROPOSED METABOLIC GENE EXPRESSION ANALYSIS')
print('='*80)

# Define metabolic genes of interest
METABOLIC_GENES = {
    'Glycolysis': [
        'HK1',      # Hexokinase 1 (first step of glycolysis)
        'HK2',      # Hexokinase 2
        'PFKM',     # Phosphofructokinase, muscle
        'PFKP',     # Phosphofructokinase, platelet
        'PKM',      # Pyruvate kinase, muscle
        'LDHA',     # Lactate dehydrogenase A
        'LDHB',     # Lactate dehydrogenase B
    ],
    'Oxidative_Phosphorylation': [
        'COX4I1',   # Cytochrome c oxidase subunit 4I1
        'COX5A',    # Cytochrome c oxidase subunit 5A
        'COX6C',    # Cytochrome c oxidase subunit 6C
        'ATP5F1A',  # ATP synthase F1 subunit alpha
        'ATP5F1B',  # ATP synthase F1 subunit beta
        'ATP5F1C',  # ATP synthase F1 subunit gamma
        'NDUFA1',   # NADH dehydrogenase (Complex I) subunit
        'NDUFB1',   # NADH dehydrogenase (Complex I) subunit
    ],
    'Glucose_Transport': [
        'SLC2A1',   # GLUT1 (ubiquitous glucose transporter)
        'SLC2A3',   # GLUT3 (neuronal glucose transporter)
        'SLC2A4',   # GLUT4 (insulin-responsive)
    ],
    'Mitochondrial_Function': [
        'TOMM20',   # Mitochondrial import receptor
        'TIMM50',   # Mitochondrial import inner membrane translocase
        'VDAC1',    # Voltage-dependent anion channel 1
        'VDAC2',    # Voltage-dependent anion channel 2
        'SLC25A4',  # ADP/ATP translocase
        'SLC25A5',  # ADP/ATP translocase
    ],
    'Lactate_Shuttle': [
        'SLC16A1',  # MCT1 (monocarboxylate transporter 1)
        'SLC16A3',  # MCT4
        'SLC16A7',  # MCT2 (neuronal)
    ],
    'Activity_Markers': [
        'FOS',      # Immediate early gene
        'JUN',      # Immediate early gene
        'ARC',      # Activity-regulated cytoskeleton-associated
        'EGR1',     # Early growth response 1
        'NPAS4',    # Neuronal PAS domain protein 4
        'NR4A1',    # Nuclear receptor subfamily 4 group A member 1
    ]
}

print('\nMetabolic gene panel (total genes):')
total_genes = 0
for category, genes in METABOLIC_GENES.items():
    print(f'  {category:30s}: {len(genes):2d} genes')
    total_genes += len(genes)
    print(f'    {", ".join(genes[:5])}{"..." if len(genes) > 5 else ""}')

print(f'\nTotal genes to analyze: {total_genes}')

print('\n' + '='*80)
print('ANALYSIS PLAN')
print('='*80)

analysis_steps = """
STEP 1: Download WHB Expression Data
-------------------------------------
Source: Allen Brain Cell Atlas S3 bucket
URL: https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/

Required files:
  - WHB-10Xv3_gene.csv (gene names/metadata, ~MB)
  - WHB-10Xv3_exon.h5ad or similar (expression matrix, ~10-50 GB)

Alternative: Use Allen Brain Cell Atlas API/SDK to query specific genes
  - abc_atlas_access Python package
  - Can query subset of genes to reduce download size

STEP 2: Extract Metabolic Gene Expression
-----------------------------------------
For each cell (3.3M cells) × each metabolic gene (40 genes):
  - Get normalized expression value
  - Calculate mean expression per cell type
  - Calculate mean expression per region

Output: Gene expression matrix (cells × genes)

STEP 3: Aggregate by Region
----------------------------
For each cortical region:
  - Calculate mean expression of each metabolic gene
  - Separate by cell type (Excitatory vs Interneuron)
  - Calculate "metabolic signature" score (weighted average of key genes)

Metrics:
  - Glycolysis score = mean(HK1, HK2, PFKM, PKM)
  - OxPhos score = mean(COX genes, ATP synthase genes)
  - Glucose uptake = mean(SLC2A1, SLC2A3)
  - Activity score = mean(FOS, ARC, EGR1)

STEP 4: Correlate with rCMRGlc
-------------------------------
Test correlations:
  - Glycolysis score vs rCMRGlc
  - OxPhos score vs rCMRGlc
  - Glucose uptake vs rCMRGlc
  - Activity score vs rCMRGlc

Separate analyses:
  - Excitatory neurons only
  - Interneurons only
  - Glia only

STEP 5: Test E:I-Specific Hypotheses
------------------------------------
Question: Do high E:I regions have different metabolic gene expression?

Compare high E:I regions (>3.0) vs low E:I regions (<2.5):
  - Metabolic gene expression in excitatory neurons
  - Metabolic gene expression in interneurons
  - Ratio of excitatory/inhibitory metabolic demand

STEP 6: Activity vs Metabolism
-------------------------------
Hypothesis: Regions with higher activity markers have higher metabolism

Correlate:
  - Activity gene expression (FOS, ARC) vs rCMRGlc
  - Activity genes in E neurons vs E:I ratio
  - Activity genes vs metabolic genes (are they co-expressed?)
"""

print(analysis_steps)

print('\n' + '='*80)
print('ALTERNATIVE APPROACH (WITHOUT FULL EXPRESSION DATA)')
print('='*80)

alt_approach = """
If downloading full expression matrices is not feasible, alternative options:

OPTION 1: Use Existing Single-Cell Atlases with Web Interface
--------------------------------------------------------------
- Allen Brain Map portal: portal.brain-map.org
- CELLxGENE Discover: cellxgene.cziscience.com
- Can query specific genes through web interface
- Manually extract mean expression per cluster/region

OPTION 2: Use Published Summary Statistics
-------------------------------------------
- Look for published papers with WHB dataset
- Many include supplementary tables with mean gene expression per cluster
- Can use cluster-level data instead of cell-level

OPTION 3: Focus on Cluster-Level Analysis
-----------------------------------------
- WHB has 461 clusters with annotations
- Each cluster has characteristic gene expression profile
- Use cluster annotations to infer metabolic characteristics:
  - "Deep-layer corticothalamic" - known to be high metabolic
  - "MGE interneuron" - includes PV cells (high metabolic)
  - etc.

OPTION 4: Literature-Based Inference
------------------------------------
- Use published knowledge about metabolic gene expression
- Known facts:
  - PV interneurons: High PV, high metabolic genes
  - Astrocytes: High glycolytic genes, lactate shuttle genes
  - Excitatory neurons: Variable metabolic profiles

Apply these rules to our cell composition data
"""

print(alt_approach)

# Create summary document
print('\n' + '='*80)
print('CREATING DOCUMENTATION...')
print('='*80)

with open(OUTPUT_DIR / 'analysis3_gene_expression_plan.txt', 'w') as f:
    f.write('METABOLIC GENE EXPRESSION ANALYSIS - IMPLEMENTATION PLAN\n')
    f.write('='*80 + '\n\n')
    f.write(f'Total metabolic genes to analyze: {total_genes}\n\n')
    f.write('GENE PANEL:\n')
    for category, genes in METABOLIC_GENES.items():
        f.write(f'\n{category}:\n')
        for gene in genes:
            f.write(f'  - {gene}\n')
    f.write('\n' + '='*80 + '\n')
    f.write('ANALYSIS STEPS:\n')
    f.write('='*80 + '\n')
    f.write(analysis_steps)
    f.write('\n' + '='*80 + '\n')
    f.write('ALTERNATIVE APPROACHES:\n')
    f.write('='*80 + '\n')
    f.write(alt_approach)

print('  Saved: analysis3_gene_expression_plan.txt')

# Save gene list
gene_list = []
for category, genes in METABOLIC_GENES.items():
    for gene in genes:
        gene_list.append({'Gene': gene, 'Category': category})

pd.DataFrame(gene_list).to_csv(OUTPUT_DIR / 'metabolic_genes_list.csv', index=False)
print('  Saved: metabolic_genes_list.csv')

print('\n' + '='*80)
print('GENE EXPRESSION ANALYSIS DOCUMENTATION COMPLETE')
print('='*80)
print('\nNOTE: Full implementation requires ~10-50 GB expression matrix download')
print('Recommended: Use Allen Brain Cell Atlas API for targeted gene queries')
print('='*80)
