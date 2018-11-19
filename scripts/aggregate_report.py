#!/usr/bin/env python

import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys, os, io
import math
from scipy import stats
from matplotlib_venn import venn2
import matplotlib.dates as mdates
import datetime as dt
## Imports for clustering and heatmaps
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy import stats
from scipy import linalg

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'generate agg report'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('--samples', help='path of the file containing the samples names (one sample per line)')
    parser.add_argument('--batchdir', help='path of the directory containing the output of multiple Pandoras runs')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('-o', '--outputdir', default='aggregate_report', help='the output directory')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    vars(args)['step'] = 'report'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    return args

# -------------------------------------

def process0(args, mysamples):
    """
    Set global variables
    """

    # set global variables
    global phylo_class_levels
    phylo_class_levels = ['Accession_ver', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom']
    global source_folder
    source_folder = args.outputdir + "/../aggregate_preprocess"
    global host_reads
    host_reads = pd.DataFrame.from_csv(source_folder + "/tot_host_reads.txt", sep = "\t", header=None)
    host_reads.columns = ["Tot_reads"]
    host_reads.index.name = "Sample ID"

    global classification_table_sp
    classification_table_sp = pd.DataFrame(columns = phylo_class_levels)

    for sample in mysamples:    
        tmp_file = pd.DataFrame.from_csv(source_folder + "/" + sample + "_processed0.tsv", sep = "\t")
        tmp_file = tmp_file[phylo_class_levels]
        tmp_file = tmp_file.drop_duplicates()
        classification_table_sp = classification_table_sp.append(tmp_file)
        classification_table_sp = classification_table_sp.drop_duplicates()

    classification_table_sp.index = range(0,classification_table_sp.shape[0])

# -------------------------------------

def process1(args, mysamples):
    """
    PCA plots by classification levels
    """

    for lvl in phylo_class_levels:
        
        if lvl == "Superkingdom":
            continue

        temp_table = pd.DataFrame.from_csv(source_folder + "/Merged_table_" + lvl + "_Num_reads.tsv", sep='\t')
        temp_table = temp_table.replace(np.NaN, 0.0)

        ### Selecting **ONLY** bacterial reads ###

        temp_super_kingdoms = []

        for index_ID in temp_table.index:
            temp_loc = classification_table_sp[classification_table_sp[lvl] == index_ID.split("__")[0]].index.tolist()[0]
            temp_super_kingdoms.append(classification_table_sp.at[temp_loc, "Superkingdom"])
            
        temp_table["Superkingdom"] = temp_super_kingdoms
        temp_table = temp_table[temp_table["Superkingdom"] == "Bacteria"]
        temp_table = temp_table[mysamples]

        # Normalize with host reads
        for sample in mysamples:
            temp_table[sample] = 10**6 * (temp_table[sample]/host_reads.at[sample, "Tot_reads"])
            
        # Calculating percent variance explained by PCA axes
        [U, S, V] = linalg.svd(temp_table)
        
        lambda_1_r = S[0]/np.sum(S) *100
        lambda_2_r = S[1]/np.sum(S) *100

        # PCA projection
        pca_temp = PCA(n_components=2)
        pca_temp_proj = pca_temp.fit_transform(temp_table.T)

        with open(args.outputdir + '/' + 'PCA1vsPCA2forSamp' + lvl + 'Expression.txt', 'w') as f:
            for index in range(temp_table.shape[1]):
                f.write(temp_table.columns[index] + "\t" + str(pca_temp_proj.T[0][index]) + "\t" + str(pca_temp_proj.T[1][index]) + "\n")        

        # ax.scatter(pca_temp_proj.T[0], pca_temp_proj.T[1], c= 'b', marker = 'o', alpha = 0.7, s = 80)

        ## Point annotations
        #for index in range(temp_table.shape[1]):
        #    ax.annotate(temp_table.columns[index], (pca_temp_proj.T[0][index], pca_temp_proj.T[1][index]), fontsize = 14)

        #ax.legend(loc='best')
        #ax.set_title('PCA 1 vs PCA2 for sample '+ lvl +' expression')
        #ax.set_xlabel("PCA 1 (%2.1f %% variance)" %lambda_1_r)
        #ax.set_ylabel("PCA 2 (%2.1f %% of variance)" %lambda_2_r)

# -------------------------------------

def process2(args, mysamples):
    """
    """

    lvl = "Species"

    lvl_table = pd.DataFrame.from_csv(source_folder + "/" + "Merged_table_" + lvl + "_Num_reads.tsv", sep='\t')
    lvl_table = lvl_table.replace(np.NaN, 0.0)

    ### Selecting **ONLY** bacterial reads ###

    lvl_super_kingdoms = []

    for index_ID in lvl_table.index:
        temp_loc = classification_table_sp[classification_table_sp[lvl] == index_ID.split("__")[0]].index.tolist()[0]
        lvl_super_kingdoms.append(classification_table_sp.at[temp_loc, "Superkingdom"])

    lvl_table["Superkingdom"] = lvl_super_kingdoms
    lvl_table = lvl_table[lvl_table["Superkingdom"] == "Bacteria"]
    lvl_table = lvl_table[mysamples]

    for index in lvl_table.index:
        if ("uncultured" in index.lower()) or ("unidentified" in index.lower()):
            lvl_table = lvl_table.drop(index, axis=0)

    # Normalize with host reads
    for sample in mysamples:
        lvl_table[sample] = 10**6 * (lvl_table[sample]/host_reads.at[sample, "Tot_reads"])
        with open(args.outputdir + '/' + 'bacteria.' + sample + '.txt', 'w') as f:
            # Normalize with host reads
            f.write(str(lvl_table.sort_values(sample, ascending=False).head(10)[[sample]]))
            f.write("\n")

    lvl_table["Means"] = np.mean(lvl_table,axis=1)
    lvl_table["Stds"] = np.std(lvl_table, axis=1)

    # fix the formatting here
    with open(args.outputdir + '/lvl_table.txt', 'w') as f:
        f.write(str(lvl_table.sort_values("Stds", ascending=False)))

    species_to_consider = list(set(lvl_table.sort_values("Stds", ascending=False).index[0:10]).union(set(lvl_table.sort_values("Means", ascending=False).index[0:10])))

    for j, species in enumerate(species_to_consider):

        # Plotting normalized species reads
        # fig, ax = plt.subplots(figsize=(10,6))
        # sp_scatter = ax.scatter(range(lvl_table.shape[1] - 2), np.log10( 1 + lvl_table.loc[species].values[:-2]), color = 'r', alpha = 0.9, marker = 'o', s = 40)

        with open(args.outputdir + '/' + 'species.' + str(j) + '.txt', 'w') as f:
            f.write("Sample\t" + species + "\n")
            # Normalize with host reads
            for i in range(10):
                # print(lvl_table.columns[i])
                f.write(lvl_table.columns[i])
                f.write("\t")
                f.write(str(np.log10( 1 + lvl_table.loc[species].values[:-2])[i]))
                f.write("\n")

        # ax.set_title(species + " reads")
        # plt.xlabel("Samples")
        # plt.ylabel("Log_10 (1 + reads per 1M human reads)")
        # ax.set_xticks(range(lvl_table.shape[1] - 2))
        # ax.set_xticklabels(lvl_table.columns[:-2], rotation = '90')
        # plt.legend(loc="best")
        # plt.show()

    # print(lvl_table.loc[species_to_consider][mysamples].T)
    with open(args.outputdir + '/' + 'abundances.species.txt', 'w') as f:
        f.write(lvl_table.loc[species_to_consider][mysamples].T.to_csv(sep="\t"))

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    mysamples = []
    with open(args.samples, 'r') as f:
        # drop the one empty field
        mysamples = f.read().split('\n')[:-1]
        # print(mysamples)

    process0(args, mysamples)
    process1(args, mysamples)
    process2(args, mysamples)

    # end of step
    hp.echostep(args.step, start=0)

# -------------------------------------

if __name__ == '__main__':

    main()
