import os
import pandas as pd

# dataset
output_name_prefix = "HFD_PQBP1_Pio_RNA-seq"

exon_skip_dir = os.path.join("result")
in_file = os.path.join(exon_skip_dir, "exon_skipping_counts.tsv")

sample_list_dic = {
    "HFD": ["HFD_1", "HFD_2", "HFD_3"],
    "HFD_Pio": ["HFD_Pio_1", "HFD_Pio_2", "HFD_Pio_3"],
    "ND": ["ND_1", "ND_2", "ND_3"]
}

# output parameters
# summarize_method:
#  'sum' for sum of counts from all samples,
#  'mean' for mean of all samples
summarize_method = 'sum'
# p_adjust_method:
#  'all' for adjusting across all exons,
#  'gene' for adjusting across exons in one gene
p_adjust_method = 'all'

# initializing input data

in_df = pd.read_table(in_file, sep = "\t", index_col="gene_exon_id")
print(in_df.shape)
in_df = in_df.loc[:, ~in_df.columns.duplicated()]
print(in_df.shape)

annotation_cols = ["accession", "gene_symbol", "exon_index"]
value_cols = ["%s_map" % x for x in sum(sample_list_dic.values(), [])]
value_cols += ["%s_skip" % x for x in sum(sample_list_dic.values(), [])]
in_df = in_df.loc[:, annotation_cols + value_cols]
