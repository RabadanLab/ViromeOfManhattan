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

# -------------------------------------

global phylo_class_levels
phylo_class_levels = ['Accession_ver', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom']

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'preprocess'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('--taxid2names', help='location of names.dmp file mapping taxid to names')
    parser.add_argument('--taxid2nodes', help='location of nodes.dmp file')
    parser.add_argument('--samples', help='path of the file containing the samples names (one sample per line)')
    parser.add_argument('--batchdir', help='path of the directory containing the output of multiple Pandoras runs')
    parser.add_argument('--accblacklist', help='A text file containing a list of accession IDs to ignore')    
    parser.add_argument('--suffixreport', help='suffix string of the report (default: /report_ifilter/report.contig.txt)')
    parser.add_argument('--suffixstats', help='suffix string of the stats report (default: /report_ifilter/report.taxon.txt)')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('-o', '--outputdir', default='aggregate_preprocess', help='the output directory')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    vars(args)['step'] = 'preprocess'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    return args

# -------------------------------------

def get_taxid_stuff(args):
    """
    """

    print('Getting phylogeny')

    # taxid to (parent taxid, rank)
    id2parentrank = {}
    # taxid to name 
    id2name = {}

    # don't want to loop more than this many iterations through the taxonomic heirarchy
    # in case taxid not found 
    maxiterations = 50 

    # open the nodes.dmp file
    with open(args.taxid2nodes, 'r') as f:
        for line in f:
            # eliminate tabs and split on |
            l = line.replace('\t', '').split('|')
            # map taxid to (parent taxid, rank)
            id2parentrank[l[0]] = (l[1], l[2])

    # open the names.dmp file
    with open(args.taxid2names, 'r') as f:
        for line in f:
            l = line.replace('\t', '').split('|')
            if (l[3] == 'scientific name'):
                id2name[l[0]] = l[1]

    return (id2parentrank, id2name)

# -------------------------------------
  
def taxonomy(taxid, id2parentrank, id2name):

    """
    Function which returns the taxonomic tree traversal as a dictionary for a single taxid
    
    INPUT: taxid string
    OUTPUT: dictionary of string-string pairs with taxonomic levels and names
    """

    # don't want to loop more than this many iterations through the taxonomic heirarchy
    # in case taxid not found
    maxiterations = 50 

    tax_tree = {}

    myiterations = 0
    current_taxid = taxid.strip()
    while (current_taxid != '1' and myiterations < maxiterations):
        
        tax_tree[id2parentrank.get(current_taxid, ['', 'not found'])[1]] = [id2name.get(current_taxid, 'not found'), current_taxid]
        
        # get parent taxid and iterate up the tree of life
        current_taxid = id2parentrank.get(current_taxid, ['1'])[0] 
        
        myiterations += 1

    return tax_tree

# -------------------------------------

def Pandora_report_process(pandora_df, myaccblacklist, id2parentrank, id2name):

    ## Sorting input report by 'saccver':
    pandora_df = pandora_df.sort_values("saccver", ascending = True)

    ## Initializing the output dataframe
    out_df = pd.DataFrame(columns=['Accession_ver', 'Title', 'Num_reads', 'Max_contig_len', 'Contig_eval', \
                                   'TaxID', 'Species', 'Genus', 'Family', 'Order', 'Class', \
                                   'Phylum', 'Superkingdom'])
    
    last_saccver = {"Accession_ver": ""}

    ## Iterating over the rows of the report        
    for row in pandora_df.iterrows():

        row = row[1]

        # Excluding known Pandora artifacts from accession ID blacklist
        if row['saccver'] in myaccblacklist:
            continue

        # Excluding entries with unknown species or synthetic constructs (cloning vector)
        if "cloning vector" in row['stitle'].lower():
            continue  

        # If the taxid appears as list with ids separated by ";", we select the first entry for taxonomic classfication 
        row_tax = taxonomy(str(row['staxids']).split(';')[0], id2parentrank, id2name)

        ## Creating new row to append with datafields as dictionary
        temp_row = {'Accession_ver': row['saccver'], 'Title': row['stitle'], 'Num_reads': row['num_reads'],\
                   'Max_contig_len': row['qlen'], 'Contig_eval': row['evalue'], "TaxID": row['staxids']}

        # Including the phylogenetic classification terms
        for class_level in ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom']:
            if class_level.lower() in row_tax.keys():
                temp_row[class_level] = row_tax[class_level.lower()][0]
            else:
                temp_row[class_level] = np.NaN

        # Check if saccver already proceesed and update entry or introduce a new row
        if not row['saccver'] == last_saccver["Accession_ver"]:
            
            out_df = out_df.append(temp_row, ignore_index=True)

        else:

            temp_index = out_df.shape[0] - 1
            
            if last_saccver["Max_contig_len"] < temp_row["Max_contig_len"]:
                out_df.at[temp_index, "Max_contig_len"] = temp_row["Max_contig_len"]
                out_df.at[temp_index, "Contig_eval"] = temp_row["Contig_eval"]
            
            out_df.at[temp_index, "Num_reads"] += temp_row["Num_reads"]
  
        last_saccver = out_df.iloc[-1]

    return out_df

# -------------------------------------

def Pandora_report_process_class_thresh(processed_report_df, class_lvl = "Accession_ver", eval_thresh = 100, max_len_thresh = 0, num_reads_thresh = 0):  

    """
    Input: a pandora proceesed data frame, a phylogentic classification level, thresholds for e_val, max_contig and for
    number of reads; 
    Defaults at Accession_ver and 0 params thresholds
    Output: new dataframe with uniquely represented classification IDs, e_val, max_contig_len and read count
    such that the latter 3 params exceed input thresholds in the cumulative results

    note: script currently ignores NANs at the given classification level

    Assumes rows represent unique Accession_ver entries
    """

    processed_report_df = processed_report_df.sort_values(class_lvl, ascending = True)

    ## Initializing the output dataframe
    if class_lvl == 'Accession_ver': 
        out_df_2 = pd.DataFrame(columns=['Accession_ver', 'Title', 'Num_reads', 'Max_contig_len', 'Contig_eval'])
    else:
        out_df_2 = pd.DataFrame(columns=[class_lvl, 'Num_reads', 'Max_contig_len', 'Contig_eval'])
    
    if class_lvl == "Accession_ver":
        last_saccver = {"Accession_ver": "", "Title": ""}
    else:
        last_saccver = {class_lvl: ""}
        
    
    if class_lvl == "Accession_ver":
        
        ###########
        
        ## Iterating over the rows of the report        
        for row in processed_report_df.iterrows():

            row = row[1]
            
            # Excluding counts with NaN at the current classification level
            if pd.isnull(row[class_lvl]):
                continue

            ## Creating new row to append with datafields as dictionary
            temp_row = {'Accession_ver': row['Accession_ver'], 'Title': row['Title'], 'Num_reads': row['Num_reads'],\
                       'Max_contig_len': row['Max_contig_len'], 'Contig_eval': row['Contig_eval']}      


            # Check if saccver already proceesed and update entry or introduce a new row
            if not temp_row['Accession_ver'] == last_saccver["Accession_ver"]:

                out_df_2 = out_df_2.append(temp_row, ignore_index=True)

            else:

                temp_index = out_df_2.shape[0] - 1

                if last_saccver["Max_contig_len"] < temp_row["Max_contig_len"]:
                    out_df_2.at[temp_index, "Max_contig_len"] = temp_row["Max_contig_len"]
                    out_df_2.at[temp_index, "Contig_eval"] = temp_row["Contig_eval"]

                out_df_2.at[temp_index, "Num_reads"] += temp_row["Num_reads"]

            last_saccver = out_df_2.iloc[-1] 
        
        ###########
        
    else:
        ## Iterating over the rows of the report        
        for row in processed_report_df.iterrows():

            row = row[1]
            
            # Excluding counts with NaN at the current classification level
            if pd.isnull(row[class_lvl]):
                continue

            ## Creating new row to append with datafields as dictionary
            temp_row = {class_lvl: row[class_lvl], 'Num_reads': row['Num_reads'],\
                       'Max_contig_len': row['Max_contig_len'], 'Contig_eval': row['Contig_eval']}      


            # Check if saccver already proceesed and update entry or introduce a new row
            if not temp_row[class_lvl] == last_saccver[class_lvl]:

                out_df_2 = out_df_2.append(temp_row, ignore_index=True)

            else:

                temp_index = out_df_2.shape[0] - 1

                if last_saccver["Max_contig_len"] < temp_row["Max_contig_len"]:
                    out_df_2.at[temp_index, "Max_contig_len"] = temp_row["Max_contig_len"]
                    out_df_2.at[temp_index, "Contig_eval"] = temp_row["Contig_eval"]

                out_df_2.at[temp_index, "Num_reads"] += temp_row["Num_reads"]

            last_saccver = out_df_2.iloc[-1]

    out_df_2 = out_df_2[(out_df_2["Num_reads"] >= num_reads_thresh) & \
                        (out_df_2["Max_contig_len"] >= max_len_thresh) & \
                        (out_df_2["Contig_eval"] <= eval_thresh)]
        
        
    return out_df_2

# -------------------------------------

def process_batch(args, mysamples, myaccblacklist, id2parentrank, id2name):
    """
    Process Pandora reports in two stages: first by Accession_version and full phylogeny;
    then by classification level only
    """

    for sample in mysamples:
        print("Processing... sample: " + sample)
        if os.path.isfile(args.outputdir + "/" + sample + "_processed0.tsv"):
            continue
        else:
            temp_report = pd.DataFrame.from_csv(args.batchdir + '/' + sample + args.suffixreport, sep="\t")
            Pandora_report_process(temp_report, myaccblacklist, id2parentrank, id2name).to_csv(args.outputdir + '/' + sample + "_processed0.tsv", sep = "\t")

    for class_lvl in phylo_class_levels:
        print("Processing... class_lvl: " + class_lvl)
        for sample in mysamples:
            filename = sample + "_processed0.tsv"
            print(" Processing sample: " + filename)
            temp_report0 = pd.DataFrame.from_csv(args.outputdir + '/' + filename, sep="\t")
            new_file_name = sample + "_" + class_lvl + "_processed1.tsv"
            Pandora_report_process_class_thresh(temp_report0, class_lvl=class_lvl).to_csv(args.outputdir + '/' + new_file_name, sep="\t")

# -------------------------------------

def process_batch2(args, mysamples):
    """
    Looping over all the classification levels and creating merge tables for every class
    Three merge tables are created: num_reads, max_contig_length and evalue
    """

    table_stats = ['Num_reads', 'Max_contig_len', 'Contig_eval']
    destination_filname = "Merged_table_"

    for class_lvl in phylo_class_levels:

        temp_merge_reads = np.nan
        temp_merge_len = np.nan
        temp_merge_eval= np.nan
        counter = 0

        tables_dict = {'Num_reads': temp_merge_reads,
            'Max_contig_len': temp_merge_len,
            'Contig_eval': temp_merge_eval}

        correct_columns_samples = []

        for sample in mysamples:
            filename = args.outputdir + '/' + sample + '_' + class_lvl + '_processed1.tsv'
            next_table = pd.DataFrame.from_csv(filename, sep= "\t")
            correct_columns_samples.append(sample)
            if class_lvl == "Accession_ver":
                next_table["Join_title"] = next_table["Accession_ver"] + "__" + next_table["Title"]
                if counter == 0:
                    counter += 1
                    for key in tables_dict.keys():
                        tables_dict[key] = next_table[["Join_title", key]]
                else:
                    for key in tables_dict.keys():
                        tables_dict[key] = pd.merge(tables_dict[key], next_table[["Join_title", key]], on="Join_title", how = 'outer', suffixes = ["", "_" + sample])
            else:
                if counter == 0:
                    counter += 1
                    for key in tables_dict.keys():
                        tables_dict[key] = next_table[[class_lvl, key]]
                else:
                    for key in tables_dict.keys():
                        tables_dict[key] = pd.merge(tables_dict[key], next_table[[class_lvl, key]], on=class_lvl, how = 'outer', suffixes = ["", "_" + sample])

        ## Outputting the merged tables

        for key in tables_dict.keys():
            if class_lvl == "Accession_ver":
                tables_dict[key].index = tables_dict[key]["Join_title"]
                tables_dict[key] = tables_dict[key].drop("Join_title", axis = 1)
            else:
                tables_dict[key].index = tables_dict[key][class_lvl]
                tables_dict[key] = tables_dict[key].drop(class_lvl, axis = 1)

            tables_dict[key].columns = correct_columns_samples
            tables_dict[key].to_csv( args.outputdir + '/' + destination_filname + class_lvl + "_" + key + ".tsv", sep = "\t" )

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

    (id2parentrank, id2name) = get_taxid_stuff(args)

    mysamples = []
    with open(args.samples, 'r') as f:
        # drop the one empty field
        mysamples = f.read().split('\n')[:-1]
    # print(mysamples)

    myaccblacklist = []
    with open(args.accblacklist, 'r') as f:
        # drop the one empty field
        myaccblacklist = f.read().split('\n')[:-1]
    # print(myaccblacklist)

    process_batch(args, mysamples, myaccblacklist, id2parentrank, id2name)
    process_batch2(args, mysamples)

    # end of step
    hp.echostep(args.step, start=0)

# -------------------------------------

if __name__ == '__main__':

    main()
