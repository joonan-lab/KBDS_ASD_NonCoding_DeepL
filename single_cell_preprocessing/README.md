# Single cell data 전처리 파이프라인

전처리 과정에서 사용한 pipeline

- scAtlas.yaml
  Conda 환경 구축에 사용.
  ``conda env create --file scAtlas.yaml``
  - Cell Ranger 버전 및 reference  
    [Cell Ranger v.6.0.2 다운로드](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0)  
    GRCh38-2020-A   ``curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz``

본 파이프라인은 [snakemake](https://snakemake.readthedocs.io/en/stable/)를 이용하여 구축함.
- Snakefile  
  워크플로우 파일.  
- config.yaml  
  설정파일.
  
