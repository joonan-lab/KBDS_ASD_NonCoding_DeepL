'''
Description:
    This script is run after 'sequence_class.py'. It annotates as the class of
    maximum sequence class score for the input variants and outputs the results
    as TSV file.

Usage:
    sequence_max_class.py <seqclass-result> <result-dir>
    sequence_max_class.py -h | --help

Options:
    <seqclass-result>       The result TSV file of sequence_class.py
    <result-dir>            Path to result directory.
'''

import os
from docopt import docopt
import pandas as pd
from natsort import natsort_keygen
from tqdm import tqdm

if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')

    seqclass_res = arguments['<seqclass-result>']
    outdir = arguments['<result-dir>']
    output = os.path.join(outdir, "max_sequence_class_score.tsv")

    seqclass_scores = pd.read_csv(seqclass_res, sep="\t")
    max_seqclass = pd.DataFrame(columns=["chrom","pos","id","ref","alt", "max_seqclass", "seqclass_max_absdiff", "seqclass_pred_score"])

    for i in tqdm(seqclass_scores.index):
        var_seqclass = seqclass_scores.loc[i]
        max_score = var_seqclass["seqclass_max_absdiff"]
        seqclass = dict(map(reversed, dict(var_seqclass[9:].abs()).items()))
        max_class = seqclass[max_score]
        seq_score = var_seqclass[max_class]

        max_seqclass = max_seqclass.append({"chrom":var_seqclass["chrom"], "pos":var_seqclass["pos"], "id":var_seqclass["id"], "ref":var_seqclass["ref"],
                                            "alt":var_seqclass["alt"], "max_seqclass":max_class, "seqclass_max_absdiff":max_score, "seqclass_pred_score":seq_score}, ignore_index=True)

    max_seqclass.sort_values(["chrom","pos"], key=natsort_keygen(), inplace=True)
    max_seqclass.to_csv(output, sep="\t", index=False)

    max_seqclass["SampleID"] = max_seqclass.id.str.split(";").str[1].str.split("=").str[1]
    max_seqclass_cnt = pd.pivot_table(data=max_seqclass, index='max_seqclass', columns='SampleID', values='id', aggfunc='count')
    max_seqclass_cnt.fillna(0, inplace=True)
    max_seqclass_cnt = max_seqclass_cnt.applymap(int)
    max_seqclass_cnt.insert(0, "max_seqclass", max_seqclass_cnt.index)
    max_seqclass_cnt.reset_index(drop=True, inplace=True)
    max_seqclass_cnt.sort_values("max_seqclass", key=natsort_keygen(), inplace=True)

    output = os.path.join(outdir, "max_sequence_class_counts_per_sample.tsv")
    max_seqclass_cnt.to_csv(output, sep="\t", index=False)

    print("Done!")