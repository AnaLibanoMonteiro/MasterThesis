#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Creates csv file with cutoff value for each sample.
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

def get_cutoff_file(norm_factors_file_path, control_cutoff, cutoff_out_file_path):
    norm_factor_df = pd.read_csv(norm_factors_file_path)

    list_cutoff_dfs = []

    for sample_group in norm_factor_df['sample_group'].unique():
        df_temp = norm_factor_df[norm_factor_df['sample_group'] == sample_group].copy()
        constant = '%temp%'

        ## Old way:
        # for sample_name in df_temp['sample_name']:
        #     if 'CTR' in sample_name:
        #         constant = df_temp[df_temp['sample_name'] == sample_name]['normalization_factor'].values[0] * control_cutoff
        # if constant == '%temp%':
        #     constant = df_temp['normalization_factor'].values[0] * control_cutoff
        #     print('There is no CTR sample. The first one will be taken as control.')

        ## New way (with sample_number):
        constant = df_temp[df_temp['sample_number'] == 1]['normalization_factor'].values[0] * control_cutoff

        df_temp['cutoff'] = df_temp.apply(lambda row: constant / row.normalization_factor, axis=1)
        list_cutoff_dfs.append(df_temp[['sample_name', 'sample_group', 'time_point', 'cutoff']])
    
    final_cutoff = pd.concat(list_cutoff_dfs)
    final_cutoff.to_csv(cutoff_out_file_path.replace('.csv', '_{}.csv'.format(control_cutoff)), index=False)
    return final_cutoff

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Get cutoff value for each sample of a dataset based on
        the normalization factors.

        Takes as input
        (1) path to to normalization factors
        (2) cutoff_out_file_path
        (3) control cutoff

        Creates a csv file with the following columns:
        sample_name, sample_group, time_point, primary_and_mapped_hg, primary_and_mapped_dm, normalization_factor
        """)

parser.add_argument('-o', '--cutoffFactorsOutFilepath',
                    nargs='?',
                    default='cutoff.csv',
                    help="""Path to save file with cutoff factors.
                    Default: 'cutoff.csv'.""",
                    metavar='<cutoffFactors outFilepath>')

parser.add_argument('-n', '--normalizationFactorsFilepath',
                    nargs='?',
                    default='normalization_factors.csv',
                    help="""Path to normalization factors file.
                    Default: 'normalization_factors.csv'.""",
                    metavar='<normalizationFactors Filepath>')

parser.add_argument('-c', '--controlCutoff',
                    nargs='?',
                    default=10,
                    help="""Cutoff value that will be applied to control
                    (and then will be adjusted for the other samples if there are spike-ins).
                    Default: 10.""",
                    metavar='<normalizationFactors Filepath>')

args = parser.parse_args()

norm_factors_file_path = args.normalizationFactorsFilepath
control_cutoff = int(args.controlCutoff)
cutoff_out_file_path = args.cutoffFactorsOutFilepath

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

get_cutoff_file(norm_factors_file_path, control_cutoff, cutoff_out_file_path)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED getting the cutoff file for: {}\n{} minutes and {} seconds.\n'.format(norm_factors_file_path, minutes, seconds))