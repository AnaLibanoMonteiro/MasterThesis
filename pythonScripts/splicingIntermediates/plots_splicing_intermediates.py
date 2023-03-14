#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input 3 files (1) cutoff file
(2) file with all events (exons last coordinate)
(3) file with splicing intermediate peaks.
Gets all the peaks from file (3) that are above the cutoff from file (1).
Creates a dataframe with total number of events, number of peaks above cutoff and percentage.
Saves that dataframe in a file and creates a barplot based on that dataframe.
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
import numpy as np
import pandas as pd
import seaborn as sns
from natsort import index_natsorted

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_cutoff_dict(cutoff_file_path):
    cutoff_df = pd.read_csv(cutoff_file_path)
    cutoff_dict = {}
    for s in cutoff_df['sample_name']:
        cutoff_dict[s] = cutoff_df[cutoff_df['sample_name'] == s]['cutoff'].values[0]
    return cutoff_dict

def get_events_above_cutoff(merged_df, cutoff_dict, sample_group):
    events_above_cutoff = merged_df.copy()
    events_above_cutoff['cutoff'] = events_above_cutoff.apply(lambda row: cutoff_dict[row.sample_name], axis=1)
    
    events_above_cutoff = events_above_cutoff.sort_values(by='sample_name',
                    key=lambda x: np.argsort(index_natsorted(events_above_cutoff['sample_name'])))

    return events_above_cutoff[events_above_cutoff['coverage_counts'] >= events_above_cutoff['cutoff']] 

def get_numb_events_and_peaks(coverage_merged_df, peaks_merged_df):
    
    numb_total_events = coverage_merged_df.groupby(['sample_name'])['name'].count().reset_index(name='numb_events')    
    numb_total_events = numb_total_events.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_total_events['sample_name'])))
    
    numb_peaks = peaks_merged_df.groupby(['sample_name'])['name'].count().reset_index(name='numb_peaks')
    numb_peaks = numb_peaks.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_peaks['sample_name'])))
    
    numb_events_and_peaks = numb_total_events.merge(numb_peaks, on='sample_name')
    numb_events_and_peaks['percentage'] = numb_events_and_peaks['numb_peaks'] / numb_events_and_peaks['numb_events'] * 100

    return numb_events_and_peaks

#####################
# Functions to Plot #
#####################

def save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff):
    file_name = '{}_{}_cutoff_{}.png'.format(sample_group, plot_name, control_cutoff)
    path_to_save_plot = os.path.join(plots_out_dir, file_name)
    plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

def add_bar_labels(list_y_positions, lits_labels):
    for i in range(len(list_y_positions)):
        plt.text(i, list_y_positions[i] + 400, '{} %'.format(round(lits_labels[i], 2)),
                 ha = 'center', fontsize=15)

def plot_numb_events_and_peaks_overlapped(numb_events_and_peaks, plots_out_dir, sample_group, control_cutoff):
    ylim_max = numb_events_and_peaks['numb_events'].max() + 4000
    with plt.style.context('seaborn-poster'):
        plot = plt.subplots(figsize=(10,7))
        plot = sns.barplot(x=numb_events_and_peaks["sample_name"], y=numb_events_and_peaks["numb_events"], color='b',
                           edgecolor="0", facecolor=(0, 0, 0, 0))
        plot = sns.barplot(x=numb_events_and_peaks["sample_name"], y=numb_events_and_peaks["numb_peaks"], color='g',
                          edgecolor="0")
        
        plt.ylim(0, ylim_max)
        plt.suptitle('Exons with splicing intermediates peaks', fontsize=25)
        add_bar_labels(numb_events_and_peaks["numb_events"].values.tolist(), numb_events_and_peaks["percentage"].values.tolist())
        plot.set(ylabel="Number Events & Number Peaks", xlabel='Time Point')
    
    plot_name = 'percentage_exons_with_peaks'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)

#################
# Main Function #
#################

def main_function(intermediates_coverage_file_path, peaks_file_path,
                  cutoff_file_path, plots_out_dir, sample_group):
    ## Get input dataframes:
    coverage_merged_df = pd.read_csv(intermediates_coverage_file_path, sep='\t')
    peaks_merged_df = pd.read_csv(peaks_file_path, sep='\t')
    ## Cutoff:
    cutoff_dict = get_cutoff_dict(cutoff_file_path)
    control_cutoff = cutoff_file_path.split('.')[0].split('_')[-1]
    print(cutoff_dict)

    ## Get events above cutoff and dataframe to plot:
    events_above_cutoff = get_events_above_cutoff(peaks_merged_df, cutoff_dict, sample_group)
    numb_events_and_peaks = get_numb_events_and_peaks(coverage_merged_df, events_above_cutoff)
      ## Save dataframe that will be ploted
    path_to_save_dataframe_to_plot = os.path.join(os.path.dirname(peaks_file_path), '{}_numb_event_and_peaks_cutoff_{}.tsv'.format(sample_group, control_cutoff))
    numb_events_and_peaks.to_csv(path_to_save_dataframe_to_plot, sep='\t', index=False)

    ## Make plot:
    plot_numb_events_and_peaks_overlapped(numb_events_and_peaks, plots_out_dir, sample_group, control_cutoff)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input 3 files (1) cutoff file, (2) file with all events (exons last coordinate),
        (3) file with splicing intermediate peaks.
        Gets all the peaks from file (3) that are above the cutoff from file (1).
        Creates a dataframe with total number of events, number of peaks above cutoff and percentage.
        Saves that dataframe in a file and creates a barplot based on that dataframe.""")

parser.add_argument('-c', '--cutoffFilepath',
                    nargs='?',
                    default='normFactorsAndCutoff/cutoff_10.csv',
                    help="""Path to file with cutoff values.
                    Default: 'normFactorsAndCutoff/cutoff_10.csv'.""",
                    metavar='<cutoff filepath>')

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('intermediatesCoverageFilepath',
                    help="""Path to file with splicing intermediates coverage.""",
                    metavar='<splicingIntermediatesCoverage filepath>')

parser.add_argument('intermediatesPeaksFilepath',
                    help="""Path to file with splicing intermediates peaks.""",
                    metavar='<intermediatesPeaks filepath>')

parser.add_argument('sampleGroup',
                    help="""Name of sample group (for eg: 601, 603, Y1P).""",
                    metavar='<sample group>')

args = parser.parse_args()

cutoff_file_path = args.cutoffFilepath
intermediates_coverage_file_path = args.intermediatesCoverageFilepath
peaks_file_path = args.intermediatesPeaksFilepath
plots_out_dir = args.outdir
sample_group = args.sampleGroup

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

main_function(intermediates_coverage_file_path, peaks_file_path, cutoff_file_path,
              plots_out_dir, sample_group)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED getting splicing intermediates peaks for: {} \n{} minutes and {} seconds.\n'.format(intermediates_coverage_file_path, minutes, seconds))
