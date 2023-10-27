# Single cell data 전처리 파이프라인

전처리 과정에서 사용한 snakemake pipeline

- scAtlas.yaml
  Conda 환경 구축에 사용.
  ``conda env create --file scAtlas.yaml``
  - Cell Ranger 버전 및 reference  
    !(Cell Ranger v.6.0.2)[https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0]
    Reference set  
    GRCh38-2020-A ``curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz``

