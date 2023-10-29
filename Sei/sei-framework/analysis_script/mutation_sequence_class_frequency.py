import pandas as pd
import os, sys
from natsort import natsort_keygen
import argparse


parser = argparse.ArgumentParser(description="Count the DNM frequency each sequence class for Sei result")
parser.add_argument("--sei_result", help="Sei result file [max_sequence_class_score.tsv]")
parser.add_argument("--adj_file", required=False, help="Adjust factor file name")
parser.add_argument("--threshold", type=float, help="Cutoff of absolute seqeunce class score")
parser.add_argument("--outdir", help="Output file directory")
parser.add_argument("--raw_output", help="Output file name for raw frequencies result")
parser.add_argument("--adj_output", help="Output file name for adjusted frequencies result")

args = parser.parse_args()
sei_res = args.sei_result
adj_file = args.adj_file
cutoff = args.threshold
outdir = args.outdir
raw_output = args.raw_output
adj_output = args.adj_output

if not os.path.exists(outdir):
    os.makedirs(outdir)


adj_factor = pd.read_csv(adj_file, sep="\t")


cnames = pd.read_csv("/home/sonic/sei/sei-manuscript/resources/cnames.tsv", sep="\t")
all_seqclass = pd.DataFrame(cnames.loc[:, ["ID","name"]].apply(lambda x: " ".join(list(x)), axis=1), columns=["seqclass"])

seq_score = pd.read_csv(sei_res, sep="\t")
seq_score = seq_score.loc[seq_score.seqclass_max_absdiff > cutoff]
seq_score["sampleID"] = seq_score.id.str.split(";").str[1].str.split("=").str[1]

## raw frequency
seq_raw_freq = pd.pivot_table(data=seq_score, index="max_seqclass", columns="sampleID", values="id", aggfunc="count")
seq_raw_freq.insert(0, "seqclass", seq_raw_freq.index)

#if len(seq_raw_freq.seqclass.unique()) != len(all_seqclass):
#    seq_raw_freq = pd.merge(all_seqclass, seq_raw_freq, on="seqclass", how="left")

seq_raw_freq.insert(1, "label", seq_raw_freq.seqclass.str.split(" ").str[0])
seq_raw_freq.insert(2, "group_label", seq_raw_freq.label.str.replace("\d", "", regex=True))
seq_raw_freq.sort_values("label", key=natsort_keygen(), inplace=True)
seq_raw_freq.reset_index(drop=True, inplace=True)
seq_raw_freq.iloc[:, 3:] = seq_raw_freq.iloc[:, 3:].fillna(0).applymap(int)

## adjusted frequency
seq_adj_freq = seq_raw_freq.copy()

samples = seq_score.sampleID.unique().tolist()
for sample in samples:
    adj_f = float(adj_factor.loc[adj_factor.SAMPLE==sample, "AdjustFactor"])
    seq_adj_freq[sample] = seq_adj_freq[sample] * adj_f

seq_raw_freq.to_csv(os.path.join(outdir, raw_output), sep="\t", index=False)
seq_adj_freq.to_csv(os.path.join(outdir, adj_output), sep="\t", index=False)