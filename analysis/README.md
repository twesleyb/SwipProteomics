## `SwipProteomics/analysis`

This directory contains executable scripts which perform the proteomics data analysis.

```
┌── 1_WASH-iBioID
│   │ # Analysis of WASH iBioID proteomics with DEP
│   └── 1_BioID-analysis.R
│
├── 2_SWIP-TMT
│   │ # Analysis of SWIP TMT with linear mixed models (MSstatsTMT and lmerTest)
│   ├── 0_PD-data-preprocess.R
│   ├── 1_MSstatsTMT-analysis.R
│   ├── 2_generate-networks.R
│   ├── 3_leidenalg-clustering.py
│   ├── 4_post-leiden.R
│   ├── 5_module-preservation.R
│   ├── 6_lmerTest-module-analysis.R
│   └── 7_module-gsea.R
│
└── 3_Plotting
      # Generate plots
```

