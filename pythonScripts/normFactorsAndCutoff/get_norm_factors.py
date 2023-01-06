#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Get normalization factors for each sample of a dataset based on
the number of normal reads vs number spike-in reads.

Takes as input
(1) path to to folder with number of reads
(2) normalization_factors_out_path

Creates a csv file with the following columns:
sample_name, sample_group, time_point, primary_and_mapped_hg, primary_and_mapped_dm, normalization_factor
""" 

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import os
import time

import pandas as pd

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_number_reads_files(number_reads_folder_path):
    """Returns a list of files that have the number of reads mapped
    (on human genome and drosophila, in case of spike-ins)
    and a boolean referent to the presence of spike-ins."""
    count_reads_dm6 = []
    count_reads = []
    for f in os.listdir(number_reads_folder_path):
        if '-dm6_number_reads.csv' in f:
            count_reads_dm6.append(f)
        elif '_number_reads.csv' in f:
            count_reads.append(f)
    count_reads_dm6.sort()
    count_reads.sort()
    if len(count_reads_dm6) == 0: ## No spike-ins
        return count_reads, False
    return list(zip(count_reads, count_reads_dm6)), True

def get_normalization_factor_file(number_reads_folder_path, normalization_factors_out_path):
    """Gets list of files with the necessary info to calculate the normalization factors
    """
    list_files, spike_ins = get_number_reads_files(number_reads_folder_path)
    
    list_dfs = []
    if spike_ins:
        for pair in list_files:
            df = pd.concat([pd.read_csv(os.path.join(number_reads_folder_path, pair[0]))[['sample_name', 'sample_group', 'time_point', 'sample_number', 'primary_and_mapped']],
                            pd.read_csv(os.path.join(number_reads_folder_path, pair[1]))[['primary_and_mapped']]], axis=1)
            df.columns = ['sample_name', 'sample_group', 'time_point', 'sample_number', 'primary_and_mapped_hg', 'primary_and_mapped_dm']
            df['normalization_factor'] = df['primary_and_mapped_hg'] / df['primary_and_mapped_dm'] / ((df['primary_and_mapped_hg'] + df['primary_and_mapped_dm']) / 1000000)
            list_dfs.append(df)
    
    else: ## no spike-ins
        for file in list_files:
            df = pd.read_csv(os.path.join(number_reads_folder_path, file))[['sample_name', 'sample_group', 'time_point', 'sample_number', 'primary_and_mapped']]
            df['normalization_factor'] = 1
            ## Should explore this:
            ## if there is no spike in, the normalization factor could take into
            ## account the library size only (instead of being just 1)
            # df['normalization_factor'] = 1 / (df['primary_and_mapped'] / 1000000)
            list_dfs.append(df)
    
    final_df = pd.concat(list_dfs)
    final_df.to_csv(normalization_factors_out_path, index=False)
    return final_df


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Get normalization factors for each sample of a dataset based on
        the number of normal reads vs number spike-in reads.

        Takes as input
        (1) path to folder with number of reads
        (2) normalization_factors_out_path

        Creates a csv file with the following columns:
        sample_name, sample_group, time_point, primary_and_mapped_hg, primary_and_mapped_dm, normalization_factor
        """)

parser.add_argument('-o', '--normalizationFactorsOutFilepath',
                    nargs='?',
                    default='normalization_factors.csv',
                    help="""Path to save file with normalization factors.
                    Default: 'normalization_factors.csv'.""",
                    metavar='<normalizationFactors outFilepath>')

parser.add_argument('-d', '--numberReadsFolderPath',
                    nargs='?',
                    default='./numberReads',
                    help="""Path to folder with number of reads.
                    Default: './numberReads'.""",
                    metavar='<numberReads folderPath>')

args = parser.parse_args()

normalization_factors_out_path = args.normalizationFactorsOutFilepath
number_reads_folder_path = args.numberReadsFolderPath

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

get_normalization_factor_file(number_reads_folder_path, normalization_factors_out_path)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED getting the normalization factors file for: {}\n{} minutes and {} seconds.\n'.format(number_reads_folder_path, minutes, seconds))
