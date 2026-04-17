#%%

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from Utils import *

# ==== reading the Counts table ========

data_path = "E-GEOD-46817-raw-counts.tsv.undecorated"

counts = Load_and_filter_dataset(data_path)

#%%

"""
here i am doing a single diff. expression. I am considering
the datasets with 3 control and 4 condition. 
"""

# ==== generating the metadata for pydeseq2 ===
# the meta data are dataframe indicating what samples are condition or control
# Ctrl = control , Con = condition

# change the number of "Ctrl" and "Con" depending on your sample
meta_data = pd.DataFrame(zip(counts.index, ["Ctrl","Ctrl","Ctrl","Con","Con","Con","Con"]), columns=["Sample", "Condition"])
meta_data = meta_data.set_index(list(meta_data)[0])

# ==== extract pydeseq2 results===
# dds is the annData dataframe. used later for PCA
results, dds = Deseq2_results(counts, meta_data)

#%%

# === filtering for good expression mean value ===

results = results[results["baseMean"] > 10]
results["neg_log_p"] = -np.log(results["padj"]) # for the volcano plot

# === significant express genes. p-value adjusted < 0.05 and log2 fold change > 0.5 ====

sign_genes = results[(results["padj"] < 0.05) & (abs(results["log2FoldChange"]) > 0.5)]


#%%
import scanpy as sc

# ==== PCA plot (scanpy)====
# interesting the PCA plot hels to see if the data are similar

sc.tl.pca(dds)
sc.pl.pca(dds, color= "Condition", size=200)

#%%

import seaborn as sns

# ==== Heatmap generation ====
# we need:
# Norm counts, gene_names, sample_names


# selecting significant genes raw counts and generate the dataframe
dds_filtered = dds[:,sign_genes.index]
dds_filtered = dds_filtered[:,:100] # manual selection of elements

df_filtered = pd.DataFrame(np.log1p(dds_filtered.layers["normed_counts"]).T,
                           index = dds_filtered.var_names,
                           columns= dds_filtered.obs_names)


# === plot the heatmap (clustermap) =====
"""
clustermap groups sample and features between similar things.
the x-axis is usually the samples. it will split samples in families
depending on their similarities (in this case gene expression levels).
This is really strong to find similarties between samples.
"""

sns.clustermap(df_filtered, z_score=0)

#%%
# === Volcano plot 
Generate_volcano_plots(results)
