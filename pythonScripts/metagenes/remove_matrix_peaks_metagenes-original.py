#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input a matrix created with deeptools computeMatrix.
Removes the lines (transcripts) that have peaks.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import os
import argparse
import pandas as pd
import gzip
import time

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
        return in_file_path.replace('.mat.gz', '_peaks_removed_per_gene.mat')

def correct_matrix_header(matrix_header, new_matrix_len):
    part1 = matrix_header.split('"group_boundaries":')[0]
    part3 = matrix_header.split('"group_boundaries":')[1].split(',', 2)[-1]
    part2 = '"group_boundaries":[0,{}],'.format(new_matrix_len)
    new_header = part1 + part2 + part3
    return new_header

def add_header_and_compress(input_matrix_path, out_file_path, new_matrix_len):
    matrix_header = pd.read_csv(input_matrix_path, sep='\t', header=None, nrows=1)[0].values[0]
    print('***************\n', matrix_header, '\n***************')
    matrix_header = correct_matrix_header(matrix_header, new_matrix_len)
    
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
    
    # temp_matrix['max'] = temp_matrix.drop([0,1,2,3,4,5], axis=1).max(numeric_only=True, axis=1)
    temp_matrix['max'] = temp_matrix.drop([0,1,2,3,4,5], axis=1).abs().max(numeric_only=True, axis=1)
    temp_matrix = temp_matrix.sort_values(by='max', ascending=False)
    # temp_matrix['sum'] = temp_matrix.drop([0,1,2,3,4,5], axis=1).sum(axis=1)
    # temp_matrix['sum'] = temp_matrix.drop([0,1,2,3,4,5], axis=1).abs().sum(axis=1)
    # temp_matrix = temp_matrix.sort_values(by='sum', ascending=False)

    new_matrix = temp_matrix.copy().drop(temp_matrix.head(number_to_remove).index).drop(temp_matrix.tail(number_to_remove).index).drop('max', axis=1)
    # new_matrix = temp_matrix.copy().drop(temp_matrix.head(number_to_remove).index).drop(temp_matrix.tail(number_to_remove).index).drop('sum', axis=1)

    new_matrix.to_csv(out_file_path, sep='\t', index=False, header=False)
    add_header_and_compress(input_matrix_path, out_file_path, len(new_matrix.index))



###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a matrix from deeptools computeMatrix.
                               Removes the lines (transcripts) that have peaks.""")

parser.add_argument('-r', '--replace',
                    nargs='?',
                    default='True',
                    help="""Boolean. True: remove original file. False: keep both files
                    (original and new one with extension '_new.mat.gz').
                    Default: True.""",
                    metavar='replace')

parser.add_argument('-p', '--percentageToRemove',
                    nargs='?',
                    default=0.03,
                    help="""Percentage to remove.
                    Default: 0.03.""",
                    metavar='percentageToRemove')
# percentage_to_remove=0.03

parser.add_argument('matrixfilepath',
                    help="""Path to input matrix.""",
                    metavar='<input matrix filepath>')

args = parser.parse_args()

input_matrix_path = args.matrixfilepath
replace = str_to_bool(args.replace)
percentage_to_remove = args.percentageToRemove

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