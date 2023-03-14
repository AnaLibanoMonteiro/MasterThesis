#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input a matrix created with deeptools computeMatrix.
Removes the peaks on each bin (column) - 
goes to each column and replace the max value by Nan;
does this several times, according to the number of peaks to remove.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import gzip
import os
import time

import pandas as pd

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def str_to_bool(string):
    return string.lower() in ("true", "t", "yes", "1", "plot", "p", "replace", "r")

def get_out_file_path(replace, in_file_path):
    if replace:
        return in_file_path.replace('.mat.gz', '.mat')
    else:
        return in_file_path.replace('.mat.gz', '_peaks_removed.mat')

def replace_col_max_by_nan(matrix):
    ## Get a list with the position of the max in each column:
    idmax_per_column = matrix.drop([0,1,2,3,4,5], axis=1).abs().idxmax()
    for i in idmax_per_column.index:
        matrix.iloc[idmax_per_column[i], i] = float("nan")
    return matrix

def add_header_and_compress(input_matrix_path, out_file_path):
    matrix_header = pd.read_csv(input_matrix_path, sep='\t', header=None, nrows=1)[0].values[0]
    
    with open(out_file_path, 'r') as f: ## Read what we already had saved on the output file
        content = f.read()
    new_content_complete = matrix_header + '\n' + content ## Add header
    
    out_file_path_compressed = out_file_path.replace('.mat', '.mat.gz')
    with gzip.open(out_file_path_compressed, 'wt') as f: ## Write again but this time compressed
        f.write(new_content_complete)
    
    os.remove(out_file_path) ## Remove matrix not compressed

def remove_matrix_peaks(input_matrix_path, out_file_path, percentage_to_remove):
    matrix = pd.read_csv(input_matrix_path, sep='\t', header=None, skiprows=[0])

    number_to_remove = round(len(matrix.index) * percentage_to_remove)
    print('number_to_remove:', number_to_remove)

    temp_matrix = matrix.copy()

    for n in range(0, number_to_remove):
        temp_matrix = replace_col_max_by_nan(temp_matrix)
    
    temp_matrix.to_csv(out_file_path, sep='\t', index=False, header=False)
    add_header_and_compress(input_matrix_path, out_file_path)



###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a matrix from deeptools computeMatrix.
                               Removes the peaks of that matrix.""")

parser.add_argument('-r', '--replace',
                    nargs='?',
                    default='True',
                    help="""Boolean. True: remove original file. False: keep both files
                    (original and new one with extension '_new.mat.gz').
                    Default: True.""",
                    metavar='replace')

parser.add_argument('-p', '--percentageToRemove',
                    nargs='?',
                    default=0.005,
                    help="""Percentage to remove.
                    Default: 0.005.""",
                    metavar='percentageToRemove')

# First used percentageToRemove as 0.03 

parser.add_argument('matrixfilepath',
                    help="""Path to input matrix.""",
                    metavar='<input matrix filepath>')

args = parser.parse_args()

input_matrix_path = args.matrixfilepath
replace = str_to_bool(args.replace)
percentage_to_remove = float(args.percentageToRemove)

out_file_path = get_out_file_path(replace, input_matrix_path)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to remove peaks for {}\n'.format(input_matrix_path))

remove_matrix_peaks(input_matrix_path, out_file_path, percentage_to_remove)

temp = time.time()-start
minutes = temp//60
seconds = temp - 60*minutes
print('FINISHED.\n{} minutes and {} seconds.\n'.format(minutes, seconds))
