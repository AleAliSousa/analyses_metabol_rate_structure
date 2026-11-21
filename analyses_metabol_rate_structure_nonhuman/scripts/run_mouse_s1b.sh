#!/usr/bin/env bash
# Run the mouse s1b pipeline in order (from anywhere inside the repo).
set -euo pipefail
cd "$(dirname "$0")"
for s in mouse_s1b_0_sync_mouse_rates mouse_s1b_1_extract_allen_wmb10x mouse_s1b_2_validate_roi_map \
         mouse_s1b_3_n_neuronal_composition mouse_s1b_3_nn_nonneuronal_composition mouse_s1b_5_EI_ratio_cortex; do
  echo "== $s"; Rscript "$s.R"
done
