import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt

def Load_and_filter_dataset(path):

    counts = pd.read_csv(path, sep = "\t") # reading the tsv file

    heads = list(counts)

    counts = counts.set_index(heads[0]) # genes are now my index

    counts = counts[counts.sum(axis = 1) > 0] # fintering genes with all zeros

    return counts.T # pydseq needs the data [sample,genes]


def Deseq2_results(counts, meta_data):

    dds = DeseqDataSet(counts=counts,
                       metadata= meta_data,
                       design_factors= "Condition") # generating the Anndataset

    dds.deseq2() # start the pre-procesing of the data

    heads = list(meta_data)
    variables = meta_data[heads[0]].unique()

    stat_res = DeseqStats(dds, 
                          n_cpus = 8, 
                          contrast= ("Condition", variables[0], variables[1]))
    
    stat_res.summary()

    results = stat_res.results_df

    return results, dds


# =======================================
# Function to  generate a volcano plot
# =======================================

def Generate_volcano_plots(Final_analysis_df):

    # you have 24 subsets of your data to plot.
    # Replace this with actual subset logic if needed.
    num_rows, num_cols = 1, 1
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(6, 5))
    # axes = axes.flatten()  # To easily index subplots in a 1D loop
    ax = axes

    try:
        # Create masks
        blue_mask = (Final_analysis_df['log2FoldChange'] <= -0.5) & (Final_analysis_df['padj'] < 0.05)
        red_mask = (Final_analysis_df['log2FoldChange'] >= 0.5) & (Final_analysis_df['padj'] < 0.05)
        gray_mask = ~(blue_mask | red_mask)

        # Plot gray points
        ax.scatter(
            Final_analysis_df.loc[gray_mask, 'log2FoldChange'],
            Final_analysis_df.loc[gray_mask, 'neg_log_p'],
            color='gray', alpha=0.5, label='Not significant'
        )

        # Plot blue points
        ax.scatter(
            Final_analysis_df.loc[blue_mask, 'log2FoldChange'],
            Final_analysis_df.loc[blue_mask, 'neg_log_p'],
            color='blue', alpha=0.8, label='Significant decrease (< -5%)'
        )

        # Annotate blue points
        for idx, row in Final_analysis_df[blue_mask].iterrows():
            ax.text(row['log2FoldChange'], row['neg_log_p'], str(int(row['positions'])), 
                    fontsize=6, color='blue', ha='center', va='bottom')

        # Plot red points
        ax.scatter(
            Final_analysis_df.loc[red_mask, 'log2FoldChange'],
            Final_analysis_df.loc[red_mask, 'neg_log_p'],
            color='red', alpha=0.8, label='Significant increase (> 5%)'
        )

        # Annotate red points
        for idx, row in Final_analysis_df[red_mask].iterrows():
            ax.text(row['log2FoldChange'], row['neg_log_p'], str(int(row['positions'])), 
                    fontsize=6, color='red', ha='center', va='bottom')

    except:
        None

        # Decorations
        ax.set_xlabel('log2FoldChange', fontsize=10)
        ax.set_ylabel('log(padj)', fontsize=10)
        ax.axhline(y=0.05, color='black', linestyle='--', linewidth=1)
        ax.axvline(x=0.5, color='black', linestyle='--', linewidth=1)
        ax.axvline(x=-0.5, color='black', linestyle='--', linewidth=1)
        ax.grid(True)
        # Optional: add legend only for the first subplot

    plt.tight_layout()
    plt.show()