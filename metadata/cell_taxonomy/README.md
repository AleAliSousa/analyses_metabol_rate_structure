# s1b cell-taxonomy crosswalk

`s1b_jorstad_to_siletti_cell_class_map.csv` maps the neuronal subclasses used
for the Jorstad et al. E:I comparison to the Siletti whole-human-brain
taxonomy used by s1b.

Jorstad class, neighborhood, subclass, and description fields are transcribed
from `Table_S1_subclass_info.csv`. `Table_S5_MERFISH_snRNA_prop_EI.xlsx`
confirms that these subclasses are the categories used for the reported
excitatory and inhibitory proportions.

Most matches can be made at Siletti's broad `supercluster_term` level. The
exception is Jorstad `L5 ET`: Siletti placed putative layer-5
extratelencephalic neurons within the heterogeneous `Miscellaneous`
supercluster. The crosswalk therefore selects only Siletti clusters 113, 114,
117, and 118, whose official taxonomy labels are `L5ET_113`, `L5ET_114`,
`L5ET_117`, and `L5ET_118`. It must never map the entire `Miscellaneous`
supercluster to excitatory neurons.

The `s1b_5c` cerebral-cortex-class analysis uses a deliberately restricted
Siletti definition. It retains neocortical projection neurons, the selected
L5 ET clusters, cortical interneurons, and the amygdalar and hippocampal
excitatory superclusters required by the project's cerebral-cortex scope. It
excludes thalamic, mammillary, midbrain-derived, rhombic-lip, and cerebellar
superclusters even when rare cells with those labels occur in a mapped
cerebral-cortex ROI. Medium and eccentric medium spiny neurons remain excluded
from the primary definition and enter only the labelled MSN sensitivity.

Running `scripts/s1b_5d_EI_cell_class_x_anatomical_scope_definition_table.R`
creates `data_analysis/s1b_5_siletti_supercluster_EI_inclusion_comparison.csv`,
which reports `E`, `I`, `Neither`, or `Partial E` for every observed Siletti
neuronal supercluster under both the s1b_5b and s1b_5c definitions.

Sources:

- Jorstad et al. 2023, Science 382, eadf6812:
  https://doi.org/10.1126/science.adf6812
- Siletti et al. 2023, Science 382, eadd7046:
  https://doi.org/10.1126/science.add7046
- Allen Brain Cell Atlas whole-human-brain taxonomy:
  https://alleninstitute.github.io/abc_atlas_access/_static/WHB-taxonomy/20240330/subcluster.html
