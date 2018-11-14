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

global phylo_class_levels
phylo_class_levels = ['Accession_ver', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom']

# -------------------------------------

def process1(args, mysamples):
    """
    PCA plots by classification levels
    """

    # global source_folder
    source_folder = args.outputdir + "/../aggregate_preprocess"

    host_reads = pd.DataFrame.from_csv(source_folder + "/tot_host_reads.txt", sep = "\t", header=None)
    host_reads.columns = ["Tot_reads"]
    host_reads.index.name = "Sample ID"

    # test_table = pd.DataFrame.from_csv(source_folder + "/Merged_table_Phylum_Num_reads.tsv", sep='\t')
    # 
    # test_table = test_table.replace(np.NaN, 0.0)
    # 
    # # Normalize with host reads
    # 
    # for sample in mysamples:
    #   test_table[sample] = 10**6 * (test_table[sample]/host_reads.at[sample, "Tot_reads"])
    # 
    # pca_test = PCA(n_components=2)
    # pca_test_proj = pca_test.fit_transform(test_table.T)
    # 
    # with open(args.outputdir + '/' + 'PCA1vsPCA2forSampPhylumExp.txt', 'w') as f:
    #   for index in range(test_table.shape[1]):
    #       f.write(test_table.columns[index] + "\t" + str(pca_test_proj.T[0][index]) + "\t" + str(pca_test_proj.T[1][index]) + "\n")

    classification_table_sp = pd.DataFrame(columns = phylo_class_levels)

    for sample in mysamples:    
        tmp_file = pd.DataFrame.from_csv(source_folder + "/" + sample + "_processed0.tsv", sep = "\t")
        tmp_file = tmp_file[phylo_class_levels]
        tmp_file = tmp_file.drop_duplicates()
        classification_table_sp = classification_table_sp.append(tmp_file)
        classification_table_sp = classification_table_sp.drop_duplicates()

    classification_table_sp.index = range(0,classification_table_sp.shape[0])

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

    process1(args, mysamples)

    # end of step
    hp.echostep(args.step, start=0)

# -------------------------------------

if __name__ == '__main__':

    main()
