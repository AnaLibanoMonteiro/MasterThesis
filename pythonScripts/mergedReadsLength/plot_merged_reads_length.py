#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input 1 file with 'fragment_length' and 'count'.
Creates hist plot based on the distribution of the fragments length
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import os
import time

import matplotlib.pyplot as plt
import pandas as pd

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def read_fragment_length_file(fragment_length_file_path):
    df_fragment_length = pd.read_csv(fragment_length_file_path, sep=' ',
                               header=None)
    df_fragment_length.columns = ['fragment_length', 'count']
    return df_fragment_length

def save_plot(plot, plots_out_dir, sample_name, plot_name):
    file_name = '{}_{}.png'.format(sample_name, plot_name)
    path_to_save_plot = os.path.join(plots_out_dir, file_name)
    plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

def plot_hist_fragment_length(df_fragment_length, plots_out_dir, sample_name):
    plot = df_fragment_length['fragment_length'].hist(weights=df_fragment_length[['count']],
                               bins=200, figsize=(10,5))
    plt.suptitle('Fragment lenght distribution of merged reads for {}'.format(sample_name), fontsize=15)
    plt.xlabel('Fragment Lenght', fontsize=12)
    
    plot_name = 'merged_reads_fragment_length'
    save_plot(plot, plots_out_dir, sample_name, plot_name)
    # save_plot(hp, 'merged_reads_fragment_length', plots_out_dir, sample_name)

def main_function(fragment_length_file_path, plots_out_dir, sample_name):
    df_fragment_length = read_fragment_length_file(fragment_length_file_path)
    plot_hist_fragment_length(df_fragment_length, plots_out_dir, sample_name)
    return df_fragment_length

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input 1 file with 'fragment_length' and 'count'.
        Creates hist plot based on the distribution of the fragments length.""")

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('fragmentLengthFilepath',
                    help="""Path to file with immediate splicing ratio.""",
                    metavar='<immediateSplicingRatio filepath>')

parser.add_argument('sampleName',
                    help="""Name of sample (for eg: 601CTR, 603DRB).""",
                    metavar='<sample name>')

args = parser.parse_args()

fragment_length_file_path = args.fragmentLengthFilepath
plots_out_dir = args.outdir
sample_name = args.sampleName

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

main_function(fragment_length_file_path, plots_out_dir, sample_name)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED hist plot of fragment lenght distribution for: {}\n{} minutes and {} seconds.\n'.format(sample_name, minutes, seconds))