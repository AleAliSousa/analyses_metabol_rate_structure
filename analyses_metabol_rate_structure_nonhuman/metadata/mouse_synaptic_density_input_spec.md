# Input contract for mouse_s1c (synaptic density vs rate)

`mouse_s1c_synaptic_density_vs_rcmr.R` reads
`data_raw/mouse_synaptic_density_by_region.csv` (override with env var
`MOUSE_SYNAPSE_CSV`). The file is **not** in the repo yet: no mouse regional
synaptic-density table exists in Evo-M1-Trait-Data, and the candidate sources
below could not be reached from the cloud container (hosts blocked).

| column | required | meaning |
|---|---|---|
| `cmrgl_region` | yes | one of the Evo-M1 mouse CMRgl region names (e.g. `Sensorimotor_cortex`, `Hippocampus`, `Thalamus`; see `data_raw/Evo_M1_mouse_CMRgl_by_region.csv`) |
| `synaptic_density` | yes | regional mean, one measure per file |
| `measure` | yes | e.g. `SV2A_puncta_density`, `PSD95_puncta_density`, `SV2A_Bmax` |
| `unit` | yes | single unit for the whole file |
| `source` | yes | citation key |
| `synaptic_density_sd`, `n_animals` | no | spread and sample size |
| `comparability_note` | no | e.g. age/sex/method caveats |

Mirrors the human s1c input (Johansen 2024 SV2A Bmax). Region names must match
the rate table; map finer atlas regions to these names before saving, keeping
the mapping in `metadata/`.

Candidate sources (check, then ingest through Evo-M1-Trait-Data):
- SV2A Mouse Synaptome Atlas, Grant lab, University of Edinburgh
  (https://git-pages.ecdf.ed.ac.uk/grantlab/sv2a-mouse-synaptome-atlas/), regional
  densities across ~140 regions at several ages.
- Dataset "Synaptic expression of SV2A in mouse and human brain", Edinburgh Research Explorer
  (https://www.research.ed.ac.uk/en/datasets/synaptic-expression-of-sv2a-in-mouse-and-human-brain/).
- Mouse SV2A PET (e.g. [11C]UCB-J / [18F]SynVesT-1) regional Bmax, if a per-region table is published.
