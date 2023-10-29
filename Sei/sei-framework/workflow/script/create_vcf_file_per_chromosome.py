'''
Description:
    This script creates temporary input vcf files per chromosome for sei-framework 'run_pipeline.sh'. e.g. [VCF_file]_chr1.vcf

'''
import pandas as pd
import os, sys

vcf_file = sys.argv[1]
outdir = sys.argv[2]

## temporary path
tmp_dir = os.path.join(outdir, "tmp/input")

if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)

filename = os.path.basename(vcf_file)
filename = os.path.splitext(filename)[0]

