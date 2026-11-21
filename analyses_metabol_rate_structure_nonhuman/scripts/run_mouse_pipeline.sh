#!/usr/bin/env bash
# Run the mouse s1a + s1b pipeline (s1c: needs data_raw/SV2A_Mouse_Synaptome_Atlas_Data.xlsx) in order (from anywhere inside the repo).
set -euo pipefail
cd "$(dirname "$0")"
for s in mouse_s1b_0_sync_mouse_rates mouse_s1b_1_extract_allen_wmb10x mouse_s1b_2_validate_roi_map \
         mouse_s1b_3_n_neuronal_composition mouse_s1b_3_nn_nonneuronal_composition mouse_s1b_5_EI_ratio_cortex \
         mouse_s1a_0_sync_isotropic_fractionator mouse_s1a_1_extract_allen_merfish_ccf mouse_s1a_2_regional_cell_composition_vs_rcmr \
         mouse_s1c_0_prepare_sv2a_synaptome mouse_s1c_synaptic_density_vs_rcmr; do
  echo "== $s"; Rscript "$s.R"
done
