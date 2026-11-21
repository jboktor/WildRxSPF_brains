# VR-3 MeApd vs MeApv test

- **Blocked:** MapMyCells `output.csv` uses snRNA-style cell IDs (`47_02_41__s1`), while seqFISH metadata uses `roi*_run*_n` IDs — zero join.
- seqFISH MapMyCells is only present for `clusters_roi4_block4` (partial).
- Deviation: MeApd/MeApv animal-level test not computed until full seqFISH Allen labels are joined.
- Gradient / node-fraction / hippocampal-audit outputs above are still valid.

