analyses_metabol_rate_structure/
├── README.md
├── .gitignore
├── scripts/
├── R/
├── data_raw/
├── data_working/
├── data_final/
├── metadata/
├── tables/
├── figs/
└── archive/

## Data corrections

See `metadata/CORRECTIONS_LOG.md` for the full log.

- **2026-06-03 — Stephan_primates.csv, `Telencephalon`, _Alouatta seniculus_.**
  Value was `37.39` (≈1000× too small). True value from the source table
  (Stephan et al. 1981) is **37,388**; an OCR step misread the thousands-comma
  as a decimal point (`37,388` → `37.388`), later rounded to `37.39`. Detected
  because the five fundamental brain parts did not sum to `Total_brain_net_volume`
  and the neocortex (`NeoWG` = 31,660) exceeded the stated telencephalon.
  Corrected to `37388` in `data_raw/Stephan_primates.csv` and in all 11 archived
  copies. Also removed a stray empty column from `data_raw/Stephan_primates.csv`
  (now 60 × 69). Re-run derivation scripts so downstream files pick up the fix.