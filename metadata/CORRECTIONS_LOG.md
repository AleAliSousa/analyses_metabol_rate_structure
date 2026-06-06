# Data corrections log

A running record of corrections made to the raw/source data files.

---

## 2026-06-03 — `Telencephalon`, *Alouatta seniculus*: OCR comma→period error

**Affected variable:** `Telencephalon` (Stephan et al. 1981 code 10), units mm³
**Affected species:** *Alouatta seniculus* (code `_0179`)

### The error
The stored value was `37.39` (or `37.388` in earlier copies) — roughly **1000×
too small**. The source table (Stephan et al. 1981) gives the value as
**37,388**. An OCR/digitisation step misread the thousands-separator comma as a
decimal point:

```
37,388  (source PDF table)
  └─ OCR ─►  37.388        (early files, e.g. 2019 insula spreadsheets)
              └─ rounded ─►  37.39   (later files, incl. data_raw)
```

This matches the note already in the metadata Timeline:
*"11 Sept 2020 — Some edits made to Stephan 1981 values which copied over with
'.' in place of ','."*

### How it was detected
1. **Additive check:** the five fundamental brain parts
   (`Medulla_oblongata + Cerebellum + Mesencephalon + Diencephalon + Telencephalon`)
   summed to only 11,622 vs. the stated `Total_brain_net_volume` of 49,009
   (76 % short). With `Telencephalon = 37388` they reconcile.
2. **Containment check:** the neocortex (`NeoWG` = 31,660) was *larger* than the
   stated telencephalon (37.39) — impossible, since neocortex is part of the
   telencephalon.
3. **Format check:** in the oldest copy, every other number in the row used a
   comma thousands-separator (`49,009`, `1,593`, `31,660`) while `Telencephalon`
   alone used a period (`37.388`) — confirming the comma→period swap.

### Corrected value
**`37388`** (written as `37,388` in files that use comma thousands-separators).

### Files corrected
| File | Previous | Now |
|---|---|---|
| `data_raw/Stephan_primates.csv` | `37.39` | `37388` |
| `archive/Stephan_primates copy.csv` | `37.39` | `37388` |
| `archive/_old/original with humans/Stephan_primates.csv` | `37.39` | `37388` |
| `archive/_old/original with humans/old/Stephan_primates_order_errors.csv` | `37.39` | `37388` |
| `archive/_old/original with humans/old/Stephan_primates_order_errors.xls` | `37.39` | `37388` |
| `archive/_old/original with humans/old/Stephan_primates.xlsx` | `37.39` | `37388` |
| `archive/_old/AAPA 20192020 analysis/Stephan_NHprimates.csv` | `37.39` | `37388` |
| `archive/_old/original with humans/old/Stephan_NHprimates 08052020.csv` | `37.388` | `37,388` |
| `archive/_old/original with humans/old/Stephan_primates08052020.csv` | `37.388` | `37,388` |
| `archive/_old/original with humans/old/Stephan_insula_primates.csv` | `37.388` | `37,388` |
| `archive/_old/original with humans/old/Stephan_insula_primates_needbrainbodydata.csv` | `37.388` | `37,388` |
| `archive/_old/original with humans/old/Stephan_insula_NHprimates.csv` | `37.388` | `37388` |

Each edit changed only the single target cell (verified cell-by-cell). No active
analysis script reads the archived copies; the live pipeline reads
`data_raw/Stephan_primates.csv`.

### Also fixed
- Removed a stray empty column (unnamed, between `Frontal.motor` and `order`)
  from `data_raw/Stephan_primates.csv`. The file is now **60 × 69**.

### Action required
Re-run the derivation scripts (e.g. `scripts/0_Heiss_Stephan_and_table1_*.R`,
`scripts/s3_predicValuesPGLS_*.R`) so any downstream/intermediate files pick up
the corrected raw value.
