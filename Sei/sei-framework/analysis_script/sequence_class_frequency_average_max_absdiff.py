import pandas as pd
import numpy as np
from natsort import natsort_keygen
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

##
cnames = pd.read_csv("/home/sonic/sei/sei-manuscript/resources/cnames.tsv", sep="\t")
cnames_dict = {}
for i in cnames.index:
    cnames_id = cnames.loc[i, "ID"]
    cnames_group = cnames.loc[i, "group"]
    cnames_dict[cnames_id] = cnames_group


##
adj_factor = pd.read_csv("/home/sonic/ASD_samples/Korean/adjustFile_patAgeOnly_Korean_ASD_WGS_v3_final_autosome_DNV_874_samples.20220816.txt", sep="\t")

case_sei = pd.read_csv("/data2/sei-framework/output/Korean_ASD_WGS/hg38/cases/max_sequence_class_score.tsv", sep="\t")
ctrl_sei = pd.read_csv("/data2/sei-framework/output/Korean_ASD_WGS/hg38/controls/max_sequence_class_score.tsv", sep="\t")


##
color_dict = {"E":"#984ef3",
              "CTCF":"#fdb462",
              "P":"#ef3b2c",
              "TF":"#386cb0",
              "HET":"#662506",
              "PC":"#a6cee3",
              "TN":"#fb9a99",
              "L":"#dddddd"}



##
case_sei["Sample"] = case_sei.id.str.split(";").str[1].str.split("=").str[1]
ctrl_sei["Sample"] = ctrl_sei.id.str.split(";").str[1].str.split("=").str[1]

case_samples = case_sei.Sample.unique().tolist()
ctrl_samples = ctrl_sei.Sample.unique().tolist()



## For cases
# adjusted frequency
case_sei_class_freq = pd.pivot_table(data=case_sei, index='max_seqclass', columns='Sample', values='id', aggfunc='count')
case_sei_class_freq.fillna(0, inplace=True)
case_sei_class_freq = case_sei_class_freq.applymap(int)
case_sei_class_freq.insert(0, "seqclass", case_sei_class_freq.index)
case_sei_class_freq.reset_index(drop=True, inplace=True)
case_sei_class_freq.sort_values("seqclass", key=natsort_keygen(), inplace=True, ignore_index=True)
case_sei_class_freq.insert(1, "label", case_sei_class_freq.seqclass.str.split(" ").str[0])
case_sei_class_freq.insert(2, "group_label", case_sei_class_freq.label.apply(lambda x: cnames_dict[x]))


case_sei_class_adj_freq = case_sei_class_freq.copy()

for sample in case_samples:
    adj_f = float(adj_factor.loc[adj_factor.SAMPLE==sample, "AdjustFactor"])
    case_sei_class_adj_freq[sample] = case_sei_class_adj_freq[sample] * adj_f

# average scores
case_sei_class_score = pd.pivot_table(data=case_sei, index='max_seqclass', columns='Sample', values='seqclass_pred_score', aggfunc='mean')
case_sei_class_score.fillna(0, inplace=True)
case_sei_class_score = case_sei_class_score.applymap(float)
case_sei_class_score.insert(0, "seqclass", case_sei_class_score.index)
case_sei_class_score.reset_index(drop=True, inplace=True)
case_sei_class_score.sort_values("seqclass", key=natsort_keygen(), inplace=True, ignore_index=True)
case_sei_class_score.insert(1, "label", case_sei_class_score.seqclass.str.split(" ").str[0])
case_sei_class_score.insert(2, "group_label", case_sei_class_score.label.apply(lambda x: cnames_dict[x]))


case_sei_class_adj_freq.to_csv("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_proband_seq_class_adjusted_freq.tsv", sep="\t", index=False)
case_sei_class_score.to_csv("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_proband_seq_class_avg_scores.tsv", sep="\t", index=False)

## For controls
# adjusted frequency
ctrl_sei_class_freq = pd.pivot_table(data=ctrl_sei, index='max_seqclass', columns='Sample', values='id', aggfunc='count')
ctrl_sei_class_freq.fillna(0, inplace=True)
ctrl_sei_class_freq = ctrl_sei_class_freq.applymap(int)
ctrl_sei_class_freq.insert(0, "seqclass", ctrl_sei_class_freq.index)
ctrl_sei_class_freq.reset_index(drop=True, inplace=True)
ctrl_sei_class_freq.sort_values("seqclass", key=natsort_keygen(), inplace=True, ignore_index=True)
ctrl_sei_class_freq.insert(1, "label", ctrl_sei_class_freq.seqclass.str.split(" ").str[0])
ctrl_sei_class_freq.insert(2, "group_label", ctrl_sei_class_freq.label.apply(lambda x: cnames_dict[x]))


ctrl_sei_class_adj_freq = ctrl_sei_class_freq.copy()

for sample in ctrl_samples:
    adj_f = float(adj_factor.loc[adj_factor.SAMPLE==sample, "AdjustFactor"])
    ctrl_sei_class_adj_freq[sample] = ctrl_sei_class_adj_freq[sample] * adj_f

# average scores
ctrl_sei_class_score = pd.pivot_table(data=ctrl_sei, index='max_seqclass', columns='Sample', values='seqclass_pred_score', aggfunc='mean')
ctrl_sei_class_score.fillna(0, inplace=True)
ctrl_sei_class_score = ctrl_sei_class_score.applymap(float)
ctrl_sei_class_score.insert(0, "seqclass", ctrl_sei_class_score.index)
ctrl_sei_class_score.reset_index(drop=True, inplace=True)
ctrl_sei_class_score.sort_values("seqclass", key=natsort_keygen(), inplace=True, ignore_index=True)
ctrl_sei_class_score.insert(1, "label", ctrl_sei_class_score.seqclass.str.split(" ").str[0])
ctrl_sei_class_score.insert(2, "group_label", ctrl_sei_class_score.label.apply(lambda x: cnames_dict[x]))


ctrl_sei_class_adj_freq.to_csv("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_control_seq_class_adjusted_freq.tsv", sep="\t", index=False)
ctrl_sei_class_score.to_csv("/data2/sei-framework/result/Korean_ASD_WGS/hg38/data/Korean_ASD_WGS_control_seq_class_avg_scores.tsv", sep="\t", index=False)


plt.rcParams['lines.markersize'] = 10


fig, ax = plt.subplots(figsize=(14, 7))

for i in range(40):
    seq_class = case_sei_class_adj_freq.loc[i, "label"]
    group = case_sei_class_adj_freq.loc[i, "group_label"]
    case_x = case_sei_class_score.iloc[i, 3:]
    case_y = case_sei_class_adj_freq.iloc[i, 3:]

    ctrl_x = ctrl_sei_class_score.iloc[i, 3:]
    ctrl_y = ctrl_sei_class_adj_freq.iloc[i, 3:]
    
    case = ax.scatter(case_x, case_y, label=seq_class, edgecolors='red', linewidth=1, alpha=.8, c=color_dict[group])
    ctrl = ax.scatter(ctrl_x, ctrl_y, label=seq_class, edgecolors='blue', linewidth=1, alpha=.8, c=color_dict[group])


case_circle = Line2D([0], [0], marker='o', color="red", markerfacecolor="grey", linewidth=0, label="Case", markersize=10, alpha=.7)
ctrl_circle = Line2D([0], [0], marker='o', color="blue", markerfacecolor="grey", linewidth=0, label="Control", markersize=10, alpha=.7)

class_group_circle = []

for c in color_dict.keys():
    circle = Line2D([0], [0], marker='o', color="black", markerfacecolor=color_dict[c], linewidth=0, label=c, markersize=10, alpha=.7)
    class_group_circle.append(circle)



legend1 = ax.legend(handles=[case_circle, ctrl_circle], title="ASD phenotype")
legend1._legend_box.sep = 10
ax.add_artist(legend1)
legend2 = ax.legend(handles=class_group_circle, title="sequence class"+'\n'+"       group", loc='upper left')
legend2._legend_box.sep = 10
ax.add_artist(legend2)


plt.xlabel("Average of sequence class scores", fontsize=13, labelpad=20, weight="bold")
plt.ylabel("Adjusted frequency of sequence class", fontsize=13, labelpad=15, weight="bold")
plt.tight_layout()
plt.savefig("/data2/sei-framework/result/Korean_ASD_WGS/hg38/figures/Korean_ASD_WGS_874samples_seq_class_adjusted_freq_avg_scores.pdf", bbox_inches='tight')
plt.show()