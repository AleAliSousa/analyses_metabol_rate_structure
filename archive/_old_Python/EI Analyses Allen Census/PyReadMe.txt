================================================================================
PYTHON SCRIPTS GUIDE - E:I RATIO vs rCMRGlc ANALYSIS
================================================================================
Created: 2026-02-10
Analysis Status: COMPLETE
Total Scripts: 26 Python files
Scripts Actually Used: 9 (35%)
Scripts Safe to Remove: 17 (65%)
================================================================================

TABLE OF CONTENTS
================================================================================
1. Executive Summary
2. Scripts Used in Final Analysis (KEEP THESE)
3. Exploratory/Unused Scripts (CAN REMOVE)
4. Phase Folders (Planning Only - CAN REMOVE)
5. Execution Order for Final Analysis
6. What Each Script Does (Alphabetical)
7. Recommendations

================================================================================
1. EXECUTIVE SUMMARY
================================================================================

ANALYSIS GOAL:
  Determine if Excitatory:Inhibitory (E:I) neuron ratios correlate with
  regional cerebral glucose metabolism (rCMRGlc).

FINAL RESULT:
  ✓ Analysis COMPLETE
  ✓ Primary finding: E:I ratio does NOT significantly correlate with metabolism
  ✓ Secondary finding: FUNCTIONAL CLASSIFICATION (sensorimotor vs associative)
    DOES correlate significantly (p=0.0089)

WORKFLOW SUMMARY:
  1. Download WHB data (672 MB) + taxonomy files
  2. Calculate E:I ratios using neurotransmitter annotations
  3. Map to brain lobes and correlate with Heiss rCMRGlc data
  4. Generate visualizations and validation
  5. Perform follow-up analyses (subregion, clustering, cell composition)
  6. Test sensorimotor vs associative hypothesis (SIGNIFICANT!)

================================================================================
2. SCRIPTS USED IN FINAL ANALYSIS (KEEP THESE - 9 FILES)
================================================================================

These scripts were successfully executed and produced the final results:

PHASE 1: DATA DOWNLOAD (3 scripts)
----------------------------------
1. download_whb_cortical.py
   Purpose: Downloads Allen Brain WHB cell metadata (672 MB)
   Runtime: ~2-3 minutes
   Output: data/WHB-10Xv3-2_cell_metadata.csv
   Status: SUCCESSFULLY EXECUTED

2. download_whb_taxonomy.py
   Purpose: Downloads taxonomy files for E/I classification
   Runtime: <1 minute
   Output: WHB_cluster*.csv files (3 files, ~1.6 MB total)
   Status: SUCCESSFULLY EXECUTED

3. verify_whb_data.py
   Purpose: Verifies data integrity and structure
   Runtime: <1 minute
   Output: Console report of data structure
   Status: SUCCESSFULLY EXECUTED


PHASE 2: CORE ANALYSIS (3 scripts)
----------------------------------
4. final_ei_analysis.py
   Purpose: Main analysis - calculates E:I ratios and correlates with rCMRGlc
   Runtime: ~0.1 minutes
   Output: ei_ratios*.csv, ei_metabolism_integrated.csv, statistics.txt
   Status: SUCCESSFULLY EXECUTED - PRIMARY ANALYSIS SCRIPT

5. create_visualizations.py
   Purpose: Generates scatter plots and visualizations
   Runtime: ~0.7 minutes
   Output: 4 PNG visualizations in data/cortical/visualizations/
   Status: SUCCESSFULLY EXECUTED

6. validate_analysis.py
   Purpose: Runs quality control checks on results
   Runtime: ~0.2 minutes
   Output: validation_report.txt (30/31 checks PASSED)
   Status: SUCCESSFULLY EXECUTED


PHASE 3: FOLLOW-UP ANALYSES (3 scripts)
---------------------------------------
7. subregion_level_analysis.py
   Purpose: Analyzes E:I at fine-grained subregion level (34 regions)
   Runtime: ~1 minute
   Output: SUBREGION_ANALYSIS_SUMMARY.txt, SUBREGION_TABLE.md
   Status: SUCCESSFULLY EXECUTED

8. comprehensive_followup_analyses.py
   Purpose: Tests sensorimotor vs associative hypothesis + cell composition
   Runtime: ~5 minutes
   Output: FOLLOWUP_ANALYSES_SUMMARY.txt, analysis1/2 visualizations
   Status: SUCCESSFULLY EXECUTED - MAJOR FINDING (p=0.0089)

9. refined_cell_composition_analysis.py
   Purpose: Analyzes interneuron subtypes (MGE, CGE, LAMP5)
   Runtime: ~2 minutes
   Output: Cell composition analysis with 9-panel visualization
   Status: SUCCESSFULLY EXECUTED

TOTAL RUNTIME: ~12 minutes
TOTAL OUTPUT: 3,369,219 cells analyzed, 20+ output files generated

================================================================================
3. EXPLORATORY/UNUSED SCRIPTS (CAN REMOVE - 10 FILES)
================================================================================

These scripts were created for exploration or planning but NOT used in the
final analysis. They can be safely deleted without affecting results:

EXPLORATORY SCRIPTS (6 files)
-----------------------------
1. assess_data_requirements.py
   Why created: Estimate download sizes from CELLxGENE API
   Why unused: Switched to direct Allen Brain downloads instead
   Safe to remove: YES

2. check_annotations.py
   Why created: Explore available annotation sets
   Why unused: Quick exploration only, not needed for analysis
   Safe to remove: YES

3. check_existing_metadata.py
   Why created: Examine metadata structure
   Why unused: Quick exploration only, not needed for analysis
   Safe to remove: YES

4. download_cellxgene_api.py
   Why created: Download via CELLxGENE REST API
   Why unused: Used direct Allen Brain S3 downloads instead
   Safe to remove: YES

5. scan_full_whb.py
   Why created: Identify all unique brain regions
   Why unused: Used pre-defined region mapping instead
   Safe to remove: YES

6. calculate_ei_hmba.py
   Why created: Test methodology on HMBA basal ganglia data
   Why unused: Proof-of-concept only, cortical analysis is final
   Safe to remove: YES


DEPRECATED VERSIONS (4 files)
-----------------------------
7. complete_ei_analysis.py
   Why created: First attempt at complete pipeline
   Why unused: Had annotation issues, replaced by complete_ei_analysis_fixed.py
   Safe to remove: YES - replaced by final_ei_analysis.py

8. complete_ei_analysis_fixed.py
   Why created: Improved version with neurotransmitter annotations
   Why unused: Evolved into final_ei_analysis.py (cleaner code)
   Safe to remove: YES - replaced by final_ei_analysis.py

9. per_lobe_clustering_analysis.py
   Why created: Analyze clustering patterns within lobes
   Why unused: Findings incorporated into comprehensive_followup_analyses.py
   Safe to remove: MAYBE - has unique clustering analysis, but redundant

10. gene_expression_analysis_plan.py
    Why created: Plan for metabolic gene expression analysis
    Why unused: Requires 10-50 GB additional downloads, not executed
    Safe to remove: NO - keep for future work if you want gene expression


================================================================================
4. PHASE FOLDERS (PLANNING ONLY - CAN REMOVE - 7 FILES)
================================================================================

These folders contain planned workflow scripts that were never executed.
The actual analysis used a different approach (direct WHB downloads).

phase1_cortical_metadata/ (4 files) - CAN REMOVE ENTIRE FOLDER
---------------------------------------------------------------
1. create_cortical_region_mapping.py
   Status: Executed once, created mapping file, no longer needed
   Safe to remove: YES (mapping file already created)

2. download_hca_cortical_metadata.py
   Status: Never executed (used WHB instead of HCA)
   Safe to remove: YES

3. download_abc_cortical_metadata.py
   Status: Never executed (used WHB direct downloads)
   Safe to remove: YES

4. download_abc_direct.py
   Status: Never executed (testing script)
   Safe to remove: YES


phase2_ei_calculation/ (1 file) - CAN REMOVE ENTIRE FOLDER
----------------------------------------------------------
5. calculate_ei_ratios_count_based.py
   Status: Never executed (logic incorporated into final_ei_analysis.py)
   Safe to remove: YES


phase3_integration/ (2 files) - CAN REMOVE ENTIRE FOLDER
--------------------------------------------------------
6. aggregate_ei_ratios_to_lobes.py
   Status: Never executed (logic incorporated into final_ei_analysis.py)
   Safe to remove: YES

7. integrate_ei_metabolism.py
   Status: Never executed (logic incorporated into final_ei_analysis.py)
   Safe to remove: YES

RECOMMENDATION: Delete all 3 phase folders entirely. They were planning
artifacts that became obsolete when the simpler WHB-based approach succeeded.

================================================================================
5. EXECUTION ORDER FOR FINAL ANALYSIS
================================================================================

If you need to reproduce the analysis from scratch, run in this order:

STEP 1: DOWNLOAD DATA
---------------------
python download_whb_cortical.py
python download_whb_taxonomy.py
python verify_whb_data.py

Expected output:
  - WHB-10Xv3-2_cell_metadata.csv (672 MB, 3.37M cells)
  - WHB_cluster.csv (83 KB)
  - WHB_cluster_annotation_term.csv (534 KB)
  - WHB_cluster_to_cluster_annotation_membership.csv (995 KB)

Runtime: ~3-4 minutes total


STEP 2: CORE E:I ANALYSIS
-------------------------
python final_ei_analysis.py
python create_visualizations.py
python validate_analysis.py

Expected output:
  - ei_ratios_whb_cortical_fine_grained.csv (34 regions)
  - ei_ratios_aggregated_lobes.csv (6 lobes)
  - ei_metabolism_integrated.csv (final dataset)
  - ei_metabolism_statistics.txt (correlation results)
  - 4 PNG visualizations
  - validation_report.txt (QC)

Runtime: ~1 minute total

Key finding: r=0.668, p=0.147 (not significant)


STEP 3: FOLLOW-UP ANALYSES
--------------------------
python subregion_level_analysis.py
python comprehensive_followup_analyses.py
python refined_cell_composition_analysis.py

Expected output:
  - SUBREGION_ANALYSIS_SUMMARY.txt
  - FOLLOWUP_ANALYSES_SUMMARY.txt (sensorimotor vs associative: p=0.0089!)
  - Multiple analysis visualizations
  - Cell composition data

Runtime: ~8 minutes total

Major finding: Sensorimotor regions have significantly higher metabolism!


STEP 4: REVIEW RESULTS
----------------------
Read in this order:
  1. READ_THIS_FIRST.txt (quick summary)
  2. FOLLOWUP_ANALYSES_SUMMARY.txt (detailed findings)
  3. View visualizations in data/cortical/followup_analyses/
  4. FINAL_ANALYSIS_REPORT.txt (comprehensive methodology)

TOTAL TIME: ~12 minutes from download to complete analysis

================================================================================
6. WHAT EACH SCRIPT DOES (ALPHABETICAL REFERENCE)
================================================================================

assess_data_requirements.py
  Queries CELLxGENE API to check available cortical datasets and estimate
  download sizes without actually downloading data.
  Status: UNUSED (switched to direct Allen Brain downloads)

calculate_ei_hmba.py
  Calculates E:I ratios from existing HMBA basal ganglia data with complete
  E/I annotations to demonstrate methodology applicable to cortical data.
  Status: UNUSED (proof-of-concept only)

check_annotations.py
  Displays all available annotation sets (neighborhoods, classes, subclasses)
  from the downloaded cluster annotation data.
  Status: UNUSED (exploratory only)

check_existing_metadata.py
  Examines the structure, columns, and region content of the cell_metadata.csv
  file to identify cortical vs. subcortical cells.
  Status: UNUSED (exploratory only)

complete_ei_analysis.py
  Runs the full E:I analysis pipeline including loading metadata, classifying
  neurons as excitatory/inhibitory, calculating E:I ratios by region, mapping
  to Heiss lobes, and correlating with rCMRGlc.
  Status: DEPRECATED (replaced by final_ei_analysis.py)

complete_ei_analysis_fixed.py
  Improved version that uses neurotransmitter-level annotations for more
  accurate E/I classification across the complete analysis pipeline.
  Status: DEPRECATED (evolved into final_ei_analysis.py)

comprehensive_followup_analyses.py
  Performs three follow-up analyses: sensorimotor vs. associative cortex
  classification, metabolic gene expression planning, and interneuron subtype
  composition (PV, SST, VIP).
  Status: USED ✓ - Generated major finding (p=0.0089)

create_visualizations.py
  Generates scatter plots, heatmaps, and visualizations showing the relationship
  between E:I ratios and rCMRGlc metabolic rates across cortical regions.
  Status: USED ✓ - Generated 4 key visualizations

download_cellxgene_api.py
  Uses the CELLxGENE REST API to programmatically search for and download
  cortical brain datasets.
  Status: UNUSED (used direct S3 downloads instead)

download_whb_cortical.py
  Downloads cortical cell metadata from the Allen Brain Whole Human Brain (WHB)
  dataset using S3 URLs.
  Status: USED ✓ - Primary data download script

download_whb_taxonomy.py
  Downloads WHB taxonomy files (cluster annotations and membership data) needed
  for E/I neurotransmitter classification.
  Status: USED ✓ - Required for E/I classification

final_ei_analysis.py
  Comprehensive analysis script that classifies E/I from neurotransmitter
  annotations and calculates complete E:I ratios mapped to lobes with Heiss
  metabolic data.
  Status: USED ✓ - PRIMARY ANALYSIS SCRIPT (main pipeline)

gene_expression_analysis_plan.py
  Assesses available data and proposes a framework for metabolic gene expression
  analysis across glycolysis, oxidative phosphorylation, and other metabolic
  pathways.
  Status: PLANNED (requires 10-50 GB additional data, not executed)

per_lobe_clustering_analysis.py
  Analyzes E:I ratios and metabolic patterns within each cortical lobe
  individually and compares clustering characteristics across lobes.
  Status: MAYBE USED (findings incorporated into comprehensive analysis)

refined_cell_composition_analysis.py
  Uses supercluster-level annotations to analyze interneuron subtypes (MGE, CGE,
  LAMP5) and their composition across cortical subregions.
  Status: USED ✓ - Analyzed cell composition vs metabolism

scan_full_whb.py
  Scans the entire WHB dataset in chunks to identify and count all unique brain
  regions present in the complete dataset.
  Status: UNUSED (used pre-defined region mapping instead)

subregion_level_analysis.py
  Keeps cortical subregions independent while assigning lobe-level rCMRGlc
  values to analyze E:I relationships at fine-grained regional resolution.
  Status: USED ✓ - Generated 34-region fine-grained analysis

validate_analysis.py
  Validates E:I analysis quality by checking file existence, sizes, data
  completeness, classification rates, and result validity.
  Status: USED ✓ - QC validation (30/31 checks passed)

verify_whb_data.py
  Loads and verifies the structure of WHB cortical metadata including columns,
  regions, and cell types to prepare for E:I analysis.
  Status: USED ✓ - Data verification after download

================================================================================
7. RECOMMENDATIONS
================================================================================

KEEP (9 files):
--------------
✓ download_whb_cortical.py
✓ download_whb_taxonomy.py
✓ verify_whb_data.py
✓ final_ei_analysis.py
✓ create_visualizations.py
✓ validate_analysis.py
✓ subregion_level_analysis.py
✓ comprehensive_followup_analyses.py
✓ refined_cell_composition_analysis.py

KEEP FOR FUTURE WORK (1 file):
------------------------------
? gene_expression_analysis_plan.py (if you want to pursue gene expression)

SAFE TO DELETE (17 files):
--------------------------
✗ assess_data_requirements.py
✗ calculate_ei_hmba.py
✗ check_annotations.py
✗ check_existing_metadata.py
✗ complete_ei_analysis.py
✗ complete_ei_analysis_fixed.py
✗ download_cellxgene_api.py
✗ per_lobe_clustering_analysis.py (maybe keep if you like clustering analysis)
✗ scan_full_whb.py

✗ phase1_cortical_metadata/ (entire folder - 4 files)
  - create_cortical_region_mapping.py
  - download_hca_cortical_metadata.py
  - download_abc_cortical_metadata.py
  - download_abc_direct.py

✗ phase2_ei_calculation/ (entire folder - 1 file)
  - calculate_ei_ratios_count_based.py

✗ phase3_integration/ (entire folder - 2 files)
  - aggregate_ei_ratios_to_lobes.py
  - integrate_ei_metabolism.py

DISK SPACE SAVINGS: ~50 KB of code (small, but reduces clutter)

ARCHIVE INSTEAD OF DELETE:
--------------------------
If you're hesitant to delete, create an "archive/" or "unused/" subfolder
and move all the unused scripts there. This keeps them for reference but
makes the main folder cleaner.

mkdir archive
move <unused_files> archive/
move phase1_cortical_metadata archive/
move phase2_ei_calculation archive/
move phase3_integration archive/

================================================================================
FINAL WORKFLOW SUMMARY
================================================================================

SUCCESSFUL PIPELINE (9 scripts, ~12 minutes):
  1. Download WHB data (download_whb_cortical.py + download_whb_taxonomy.py)
  2. Verify data (verify_whb_data.py)
  3. Calculate E:I ratios (final_ei_analysis.py)
  4. Create visualizations (create_visualizations.py)
  5. Validate results (validate_analysis.py)
  6. Analyze subregions (subregion_level_analysis.py)
  7. Test functional hypothesis (comprehensive_followup_analyses.py)
  8. Analyze cell composition (refined_cell_composition_analysis.py)
  9. Review results (READ_THIS_FIRST.txt, FOLLOWUP_ANALYSES_SUMMARY.txt)

KEY FINDINGS:
  Primary: E:I ratio does NOT correlate with metabolism (r=0.668, p=0.147)
  Secondary: Functional classification DOES! (sensorimotor vs associative, p=0.0089)

PUBLICATION POTENTIAL:
  "Functional Classification Explains Regional Metabolic Variation in
   Human Cortex: A Single-Cell Atlas Study"

================================================================================
END OF PYTHON SCRIPTS GUIDE
================================================================================
Questions? See READ_THIS_FIRST.txt or FOLLOWUP_ANALYSES_SUMMARY.txt
