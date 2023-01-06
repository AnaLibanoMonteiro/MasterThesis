#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input matrix (or matrices) created with deeptools computeMatrix.
Makes a plot for each matrix and then a plot with all matrices together
(in case there is more than one)
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import gzip
import json
import os
import time

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def convert_input_arg_to_list(input_arg):
    return [i for i in input_arg.split(',')]

def get_out_file_path(matrix_path, out_dir, list_sample_names=None, sample_group=None, CTR_W30=None):
    file_name = os.path.basename(matrix_path).split('.')[0]
    if sample_group:
        for sample_name in list_sample_names:
            if sample_name in file_name:
                file_name = file_name.replace(sample_name, sample_group)
    file_name += '.svg'
    if CTR_W30:
        file_name = file_name.replace(sample_group, sample_group + '_CTR_W30')
    return os.path.join(out_dir, file_name)

def read_metagenes_matrix(matrix_path):
    with gzip.open(matrix_path, 'rt') as f:
        matrix_header = json.loads(f.readline().strip('@'))
    matrix = pd.read_csv(matrix_path, sep='\t', header=None, skiprows=[0])
    return matrix_header, matrix

def get_forward_and_reverse_values(matrix_header, matrix):
    forward = matrix_header['sample_boundaries'][1]
    reverse = matrix_header['sample_boundaries'][2]
    sample_label1 = matrix_header['sample_labels'][0].replace('_woIPri_SNR', '').replace('_F', ' Forward')
    sample_label2 = matrix_header['sample_labels'][1].replace('_woIPri_SNR', '').replace('_R', ' Reverse')
    
    forward_values = matrix.iloc[:, 6:forward+6].mean(axis=0).reset_index(drop=True)
    reverse_values = matrix.iloc[:, forward+6:reverse+6].mean(axis=0).reset_index(drop=True)
    
    forward_and_reverse = pd.DataFrame({sample_label1:forward_values, sample_label2:reverse_values})
    
    return forward_and_reverse

def plot_metagenes(matrix_header, df_to_plot, out_file_path):
    wong_pallete = ['#e69f00','#d55e00','#009e73','#f0e442','#56b4e9','#0072b2','#cc79a7']
    if len(df_to_plot.columns) == 4:
        wong_pallete = ['#e69f00','#d55e00','#56b4e9','#0072b2']

    title = os.path.basename(out_file_path).split('.')[0].replace('_peaks_removed', '')

    up = matrix_header['upstream'][0]
    body = matrix_header['body'][0]
    down = matrix_header['downstream'][0]
    bin_size = matrix_header['bin size'][0]
    number_genes = matrix_header['group_boundaries'][1]
    ref_point = matrix_header['ref point'][0]
    
    with plt.style.context('seaborn-poster'):
        lp = df_to_plot.plot.line(figsize=(14,7), color=wong_pallete)
        if ref_point == None:
            lp.set_xticks([0, up/bin_size, (up+body)/bin_size, (up+body+down)/bin_size])
            lp.set_xticklabels(['- {} Kb'.format(up/1000),'TSS','TES','+ {} Kb'.format(down/1000)], fontsize=20)
        else:
            lp.set_xticks([0, up/bin_size, (up+down)/bin_size])
            lp.set_xticklabels(['- {} Kb'.format(up/1000), ref_point,'+ {} Kb'.format(down/1000)], fontsize=20)
        
        plt.yticks(fontsize=14)
        plt.title('{} Metagenes\n n={}'.format(title, number_genes))
        lp.figure.savefig(out_file_path, bbox_inches='tight')

def main_function(list_matrices_paths, out_dir, list_sample_names, sample_group):
    list_dfs_to_plot = []
    list_dfs_to_plot_CTR_and_W30 = []

    ## First plot metagenes for each matrix of this group:
    for path in list_matrices_paths:
        out_file_path = get_out_file_path(path, out_dir)
        matrix_header, matrix = read_metagenes_matrix(path)
        df_to_plot = get_forward_and_reverse_values(matrix_header, matrix)
        plot_metagenes(matrix_header, df_to_plot, out_file_path)
        
        list_dfs_to_plot.append(df_to_plot) ## Will be used later eventually
        file_name = os.path.basename(path).split('.')[0]
        if 'CTR' in file_name or 'W30' in file_name:
            list_dfs_to_plot_CTR_and_W30.append(df_to_plot) ## Will be used later eventually

    ## Then plot metagenes for all matrices of this group (if there is more than 1 matrix):
    if len(list_dfs_to_plot) > 1:
        out_file_path = get_out_file_path(path, out_dir, list_sample_names, sample_group)
        df_to_plot = pd.concat(list_dfs_to_plot, axis=1)
        plot_metagenes(matrix_header, df_to_plot, out_file_path)
        
    ## Finally, plot only CTR and W30 metagenes if they are in the list of samples:
    if len(list_dfs_to_plot_CTR_and_W30) == 2:
        out_file_path = get_out_file_path(path, out_dir, list_sample_names, sample_group, 'CTR_W30')
        df_to_plot = pd.concat(list_dfs_to_plot_CTR_and_W30, axis=1)
        plot_metagenes(matrix_header, df_to_plot, out_file_path)



###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a matrix (or matrices) from deeptools computeMatrix.
                               Makes a plot for each matrix and then a plot with all matrices together
                               (in case there is more than one).""")

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('listSampleNames',
                    help="""List with sample names from the sample group
                    (for eg: 601CTR, 601DRB, 601W30).""",
                    metavar='<sample Name>')

parser.add_argument('sampleGroup',
                    help="""Name of sample group (for eg: 601, 603, Y1P).""",
                    metavar='<sample group>')


parser.add_argument('listMatricesPaths',
                    help="""List with paths to input matrix/matrices.""",
                    metavar='<input matrix filepath>')

args = parser.parse_args()

out_dir = args.outdir
list_sample_names = convert_input_arg_to_list(args.listSampleNames)
sample_group = args.sampleGroup
list_matrices_paths = convert_input_arg_to_list(args.listMatricesPaths)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to create metagenes plots for {}\n'.format(list_matrices_paths))

main_function(list_matrices_paths, out_dir, list_sample_names, sample_group)

temp = time.time()-start
minutes = temp//60
seconds = temp - 60*minutes
print('FINISHED.\n{} minutes and {} seconds.\n'.format(minutes, seconds))
