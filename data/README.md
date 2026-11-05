# data/

Numerical source data for every figure in the manuscript and supplementary
materials. The directory is intentionally **flat**: each file is named
descriptively, deduplicated across figures, and referenced directly by the
figure scripts in `scripts/`.

`06_trees/` contains the per-species StrainPhlAn 4 phylogenetic trees used
by the legacy strain-tree analysis.

Larger files (e.g. `ndb.csv.gz`, `mag_metadata_v3.csv.gz`) are stored
gzipped to keep clone size manageable.
