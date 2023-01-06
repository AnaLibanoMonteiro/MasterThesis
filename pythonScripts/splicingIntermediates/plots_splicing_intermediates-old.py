#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
TODO!!!!!
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import os
import time

import numpy as np
import pandas as pd
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

#####################
# Functions to Plot #
#####################

def save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff):
    file_name = '{}_{}_cutoff_{}.png'.format(sample_group, plot_name, control_cutoff)
    path_to_save_plot = os.path.join(plots_out_dir, file_name)
    plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

# def save_plot(plot, out_file_name, plots_out_dir, sample_group, control_cutoff):
#     plot.figure.savefig('{}/{}_{}_cutoff_{}.png'.format(plots_out_dir, sample_group, out_file_name, control_cutoff),
#                         dpi=100, bbox_inches='tight')

def plot_numb_total_events(coverage_merged_df, plots_out_dir, sample_group, control_cutoff):
    numb_total_events = coverage_merged_df.groupby(['sample_name'])['name'].count().reset_index(name='count')    
    numb_total_events = numb_total_events.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_total_events['sample_name'])))
    ylim_max = numb_total_events['count'].max() + 3000
    plot = numb_total_events.plot.bar(x='sample_name', y='count', rot=0, xlabel='Time Point', ylim=(0,ylim_max),
                        title='Number of exons expressed per time point')
    plot.bar_label(plot.containers[0])
    plot_name = 'numb_exons_expressed'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, 'numb_exons_expressed', plots_out_dir, sample_group, control_cutoff)
    return numb_total_events

def plot_numb_peaks(events_above_cutoff, numb_total_events, plots_out_dir, sample_group, control_cutoff):
    numb_peaks = events_above_cutoff.groupby(['sample_name'])['name'].count().reset_index(name='count')
    numb_peaks = numb_peaks.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_peaks['sample_name'])))
    ylim_max = numb_total_events['count'].max() + 3000
    plot = numb_peaks.plot.bar(x='sample_name', y='count', rot=0, xlabel='Time Point', ylim=(0,ylim_max),
                       title='Number of exons with splicing intermediates peaks for\n{}'.format(sample_group))
    plot.bar_label(plot.containers[0])
    plot_name = 'numb_exons_with_peaks'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, 'numb_exons_with_peaks', plots_out_dir, sample_group, control_cutoff)
    return numb_peaks

def plot_numb_events_and_peaks_overlapped(coverage_merged_df, events_above_cutoff):
    
    numb_total_events = coverage_merged_df.groupby(['sample_name'])['name'].count().reset_index(name='count')    
    numb_total_events = numb_total_events.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_total_events['sample_name'])))
    
    numb_peaks = events_above_cutoff.groupby(['sample_name'])['name'].count().reset_index(name='count')
    numb_peaks = numb_peaks.sort_values(by='sample_name',
                        key=lambda x: np.argsort(index_natsorted(numb_peaks['sample_name'])))

    return

def plot_percentage(numb_total_events, numb_peaks, aka, plots_out_dir, sample_group, control_cutoff):
    percentage_df = pd.merge(numb_total_events, numb_peaks, on=['sample_name'])
    percentage_df.columns = ['sample_name', 'total_counts', 'peaks_counts']
    percentage_df['percentage'] = percentage_df['peaks_counts'] / percentage_df['total_counts'] * 100
    title = 'Percentage of expressed exons with splicing intermediates peaks for\n{}'

    if percentage_df['percentage'].max() > 40:
        ylim_max = percentage_df['percentage'].max() + 10
    else:
        ylim_max = 50

    plot = percentage_df.plot.bar(x='sample_name', y='percentage', title=title.format(aka), color=['tan'],
                                rot=0, figsize=(7,5), ylim=(0,ylim_max), xlabel='Time Point')
    plot.bar_label(plot.containers[0])
    plot_name = 'percentage_exons_with_peaks'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, 'percentage_exons_with_peaks', plots_out_dir, sample_group, control_cutoff)
    return percentage_df

#################
# Main Function #
#################

def main_function(intermediates_coverage_file_path, peaks_file_path,
                  cutoff_file_path, plots_out_dir, sample_group, aka):
    coverage_merged_df = pd.read_csv(intermediates_coverage_file_path, sep='\t')
    peaks_merged_df = pd.read_csv(peaks_file_path, sep='\t')

    cutoff_dict = get_cutoff_dict(cutoff_file_path)
    control_cutoff = cutoff_file_path.split('.')[0].split('_')[-1]
    print(cutoff_dict)

    events_above_cutoff = get_events_above_cutoff(peaks_merged_df, cutoff_dict, sample_group)
    # events_above_cutoff.to_csv(peaks_file_path.replace('.tsv', '_cutoff_{}.tsv'.format(control_cutoff)),
    #                      sep='\t', index=False)

    numb_total_events = plot_numb_total_events(coverage_merged_df, plots_out_dir, sample_group, control_cutoff)
    numb_peaks = plot_numb_peaks(events_above_cutoff, numb_total_events, plots_out_dir, sample_group, control_cutoff)
    percentage_df = plot_percentage(numb_total_events, numb_peaks, aka, plots_out_dir, sample_group, control_cutoff)

    # plot_CTR_recover_in_W30(coverage_merged_df, peaks_merged_df, aka, group_name)
    # plot_CTR_recover_in_DRB(coverage_merged_df, peaks_merged_df, aka)
    # venn_diagram = venn_diagram_CTR_vs_W30(peaks_merged_df)
    return {'numb_total_events': numb_total_events, 'numb_peaks': numb_peaks, 'percentage_df': percentage_df,
            'coverage_merged_df': coverage_merged_df, 'events_above_cutoff': events_above_cutoff}


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input 3 files (1) cutoff file, (2) file with all events (exons last coordinate),
        (3) file with splicing intermediate peaks.
        Creates a another file with the peaks in (2) that are above the cutoff from file (1).
        Creates several plots based on the splicing intermediate peaks.""")

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
              plots_out_dir, sample_group, sample_group)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED getting splicing intermediates peaks for: {} \n{} minutes and {} seconds.\n'.format(intermediates_coverage_file_path, minutes, seconds))
