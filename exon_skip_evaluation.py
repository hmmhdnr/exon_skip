# %%
import pandas as pd
import os
import math
from importlib import import_module
import argparse

from tqdm import tqdm
from scipy import stats

normalize_counts = False

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--param_name", help="a name of parameter set; i.e. param for param.py", type=str, default="params")
parser.add_argument("-s", "--summarize_method", help="specify a method to summarize counts; sum for sum of counts; mean for mean of counts", type=str, default="sum")
parser.add_argument("-a", "--p_adjust_method", help="specify a way of applying BH adjustment; all for all exons; gene for exons only in one gene", type=str, default="all")
args = parser.parse_args()

# dataset
params = import_module(args.param_name)

output_name_prefix = params.output_name_prefix
in_df = params.in_df
sample_list_dic = params.sample_list_dic
print(in_df.shape)
in_df.head()

# adding gene annotations
exon_ann_df = in_df.loc[:, ["accession", "gene_symbol", "exon_index"]]
exon_ann_df.head()
print(exon_ann_df.shape)

case_groups = ["HFD_Pio"]
control_groups = ["HFD"]

# comparison between 3 disease model mice and 3 normal ones using Data set #2 (new_n3)
def fisher_exact_(x):
    or_, pval_ = stats.fisher_exact([[x.case_skip, x.case_map], [x.cont_skip, x.cont_map]])
    return pd.Series(index = ["OddsRatio", "pValue"], data = [or_, pval_])

df_list = []

for case_group, control_group in zip(case_groups, control_groups):
    case_samples = sample_list_dic[case_group]
    control_samples = sample_list_dic[control_group]
    print(case_group, case_samples)
    case_map = in_df.loc[:, ["%s_map" % x for x in case_samples]].astype(float)
    case_skip = in_df.loc[:, ["%s_skip" % x for x in case_samples]].astype(float)
    cont_map = in_df.loc[:, ["%s_map" % x for x in control_samples]].astype(float)
    cont_skip = in_df.loc[:, ["%s_skip" % x for x in control_samples]].astype(float)
    if params.summarize_method == "sum":
        x = pd.DataFrame({
            "case_map": case_map.sum(axis = 1).astype(float),
            "case_skip": case_skip.sum(axis = 1).astype(float),
            "cont_map": cont_map.sum(axis = 1).astype(float),
            "cont_skip": cont_skip.sum(axis = 1).astype(float)
        })
    elif params.summarize_method == "mean":
        x = pd.DataFrame({
            "case_map": case_map.sum(axis = 1).astype(float).apply(lambda x: math.ceil(x / case_map.shape[1])),
            "case_skip": case_skip.sum(axis = 1).astype(float).apply(lambda x: math.ceil(x / case_skip.shape[1])),
            "cont_map": cont_map.sum(axis = 1).astype(float).apply(lambda x: math.ceil(x / cont_map.shape[1])),
            "cont_skip": cont_skip.sum(axis = 1).astype(float).apply(lambda x: math.ceil(x / cont_skip.shape[1]))
        })

    tqdm.pandas(desc = "Exon-skip in %s" % case_group)
    r = x.progress_apply(fisher_exact_, axis = 1)

    if params.p_adjust_method == "gene":
        ## BH for each gene
        gene_id = [n.split("_")[0] for n in r.index]
        r.loc[:, "gene_id"] = gene_id
        list_tmp_df = []
        for gid_ in tqdm(list(set(gene_id))):
            df_ = r.loc[r.gene_id == gid_, :].copy()
            fdr_bh = stats.false_discovery_control(df_.pValue, method = "bh")
            df_.loc[:, "fdr_bh"] = fdr_bh
            list_tmp_df.append(df_)
        r = pd.concat(list_tmp_df)
    else:
        ## BH for all
        fdr_bh = stats.false_discovery_control(r.pValue, method = "bh")
        r.loc[:, "fdr_bh"] = fdr_bh
    
    result_df_ = pd.concat([exon_ann_df, x, r], axis = 1)
    result_df_ = result_df_.rename(columns = {
        "case_map": "%s_map" % case_group,
        "case_skip": "%s_skip" % case_group,
        "OddsRatio": "%s_OddsRatio" % case_group,
        "pValue": "%s_pValue" % case_group,
        "fdr_bh": "%s_fdr_bh" % case_group
    })
    if normalize_counts:
        file_out_ = "exon_skipping_%s_in_%s_normalized.tsv" % (case_group, output_name_prefix)
    else:
        file_out_ = "exon_skipping_%s_in_%s_raw.tsv" % (case_group, output_name_prefix)

    result_df_.to_csv(file_out_, sep = "\t")
    df_list.append(result_df_)

result_merged_df = pd.concat(df_list, axis = 1)
result_merged_df = result_merged_df.loc[:, ~result_merged_df.columns.duplicated()]

if normalize_counts:
    file_out = "exon_skipping_%s_normalized.tsv" % output_name_prefix
else:
    file_out = "exon_skipping_%s_raw.tsv" % output_name_prefix

result_merged_df.to_csv(file_out, sep = "\t")
