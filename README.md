# IBD-virome

Code for: Unveiling the viral signature of inflammation in murine models of inflammatory bowel disease

## Contents

### Virome analysis

- `virome_pipeline.sh`: virome sequencing data analyses using the [ViroProfiler](https://github.com/deng-lab/viroprofiler) pipeline.
- `alpha.diversity_beta.diversity.R`: alpha and beta diversity script.
- `ANCOM_NetComi.R`: ANCOM and NetComi script.
- `Lme_cor.R`: Lme and cor script.
- `XGBOOST_SHAP.R`: XGBOOST and SHAP script.
- `remap.sh` and `abundance.sh`: Remap public IBD virome data to reference viral contigs and calculate abundance.
- `samtools.sh` and `vOTUs.R`: vOTUs in IBD patient script.

### Bacterial microbiome analysis (16S)

- `amplicon_pipeline.sh`: 16S rRNA gene sequencing data analyses using the [QIIME2](https://qiime2.org/) and [nf-ampseq](https://github.com/deng-lab/nf-ampseq) pipeline.
