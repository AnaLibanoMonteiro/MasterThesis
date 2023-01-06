#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
# Receives a file like this: (created with 'get_immediate_splicing_ratio.py' and merged after)

intron_id                  splicing_ratio      numb_spliced_reads  numb_unspliced_reads  numb_total_reads  transcript_id    intron_rank  time_point  sample_name
ENST00000382871_intron_4                       0                   0                     0                 ENST00000382871  4            CTR         601CTR
ENST00000382871_intron_5                       0                   0                     0                 ENST00000382871  5            CTR         601CTR
ENST00000382871_intron_1   0.0                 0                   1                     1                 ENST00000382871  1            CTR         601CTR
ENST00000382871_intron_13                      0                   0                     0                 ENST00000382871  13           CTR         601CTR
ENST00000382871_intron_24  1.0                 1                   0                     1                 ENST00000382871  24           CTR         601CTR
ENST00000382871_intron_12                      0                   0                     0                 ENST00000382871  12           CTR         601CTR
ENST00000382871_intron_19  0.0                 0                   2                     2                 ENST00000382871  19           CTR         601CTR
ENST00000382871_intron_14  0.3333333333333333  1                   2                     3                 ENST00000382871  14           CTR         601CTR

# Creates a file with the events above the cutoff and creates plots based on those events.
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

    return events_above_cutoff[events_above_cutoff['numb_total_reads'] >= events_above_cutoff['cutoff']] 

#####################
# Functions to Plot #
#####################

def save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff):
    file_name = '{}_{}_cutoff_{}.png'.format(sample_group, plot_name, control_cutoff)
    path_to_save_plot = os.path.join(plots_out_dir, file_name)
    plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

# def save_plot(plot, out_file_name, out_dir, sample_group, control_cutoff):
#     plot.figure.savefig('{}/{}_{}_cutoff_{}.png'.format(out_dir, sample_group, out_file_name, control_cutoff),
#                         dpi=100, bbox_inches='tight')

def box_plot_splicing_ratio(df, aka, plot_name, plots_out_dir, sample_group, control_cutoff):
    ordering = enumerate(df['sample_name'].unique())
    positions = [ind for val, ind in sorted((v, i) for (i, v) in ordering)]
    plot = df.boxplot(column=['splicing_ratio'], by='sample_name', figsize=(7,7), grid=False, fontsize=13,
                    widths=0.3, boxprops=dict(linewidth=1.5, color='b'),
                    medianprops=dict(linewidth=1.5, color='g'), positions=positions)
    plot.set_title('')
    plot.set_xlabel('Time Point', fontsize=13)
    plot.figure.suptitle('Splicing Ratio Boxplot for {}'.format(aka), fontsize=14)
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, out_file_name, plots_out_dir, sample_group, control_cutoff)

def density_plot_splicing_ratio(df, aka, plot_name, plots_out_dir, sample_group, control_cutoff):
    palette = ['#e69f00','#56b4e9','#009e73','#f0e442','#0072b2','#d55e00','#cc79a7']
    with plt.style.context('seaborn-poster'):
        plot = sns.displot(data=df, kind='kde', x='splicing_ratio',
                        hue='sample_name', aspect=1.7, bw_adjust=0.2, linewidth=3, palette=palette[:3])

        sns.move_legend(plot, loc='upper center', bbox_to_anchor=(0.45,0.8))
        plt.xlabel('Immediate Splicing Ratio')
        plt.title('Immediate Splicing Ratio Density plot for {}'.format(aka))

    # plot = sns.displot(data=df, kind='kde', x='splicing_ratio',
    #                  hue='sample_name', aspect=1.6, bw_adjust=0.2)
    # plt.xlabel('Immediate Splicing Ratio', fontsize=13)
    # plt.title('Immediate Splicing Ratio Density plot for {}'.format(aka), fontsize=14)

    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, out_file_name, plots_out_dir, sample_group, control_cutoff)


def get_CTR_W30_diff(events_above_cutoff):
    events_CTR_in_W30 = list(set(events_above_cutoff[events_above_cutoff['time_point'] == 'W30']['intron_id']) &
                           set(events_above_cutoff[events_above_cutoff['time_point'] == 'CTR']['intron_id']))

    subset_W30 = events_above_cutoff[(events_above_cutoff['time_point'] == 'W30') & 
                (events_above_cutoff['intron_id'].isin(events_CTR_in_W30))][['intron_id', 'splicing_ratio']].reset_index(drop=True)

    subset_CTR = events_above_cutoff[(events_above_cutoff['time_point'] == 'CTR') & 
                (events_above_cutoff['intron_id'].isin(events_CTR_in_W30))][['intron_id', 'splicing_ratio']].reset_index(drop=True)

    events_CTR_vs_W30 = pd.merge(subset_CTR, subset_W30, on=['intron_id'])

    events_CTR_vs_W30['diff'] = events_CTR_vs_W30['splicing_ratio_x'] - events_CTR_vs_W30['splicing_ratio_y']
    
    return events_CTR_vs_W30

def plot_CTR_W30_bx_with_lines(events_CTR_vs_W30, plots_out_dir, sample_group, control_cutoff):
    temp_CTR = events_CTR_vs_W30[['splicing_ratio_x', 'intron_id']].copy()
    temp_CTR['time_point'] = 'CTR'
    temp_CTR.columns = ['splicing_ratio', 'intron_id', 'time_point']
    temp_W30 = events_CTR_vs_W30[['splicing_ratio_y', 'intron_id']].copy()
    temp_W30['time_point'] = 'W30'
    temp_W30.columns = ['splicing_ratio', 'intron_id', 'time_point']
    CTR_W30 = pd.concat([temp_CTR, temp_W30], ignore_index=True)
    
    plt.figure(figsize=(8,11))
    plot = sns.boxplot(data=CTR_W30, y='splicing_ratio', x='time_point')
    plot = sns.swarmplot(data=CTR_W30, y='splicing_ratio', x='time_point', hue='time_point', s=1,
                       palette=['red','green'])
    plot = sns.lineplot(data=CTR_W30, y='splicing_ratio', x='time_point', units='intron_id', color='.7', estimator=None)
    plot_name = 'boxplot_with_lines'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, 'boxplot_with_lines', plots_out_dir, sample_group, control_cutoff)

def plot_numb_events(df, plots_out_dir, sample_group, control_cutoff):
    df = df.groupby(['sample_name'])['splicing_ratio'].count().reset_index(name='count')
    df = df.sort_values(by='sample_name',
            key=lambda x: np.argsort(index_natsorted(df['sample_name'])))
    plot = df.plot.bar(x='sample_name', y='count', title='Number of events > cutoff for {}'.format(sample_group))
    plot_name = 'number_events_in_boxplot'
    save_plot(plot, plots_out_dir, sample_group, plot_name, control_cutoff)
    # save_plot(plot, 'number_events_in_boxplot', plots_out_dir, sample_group, control_cutoff)

#################
# Main Function #
#################

def main_function(immediate_splicing_ratio_file_path, sample_group, aka, cutoff_file_path, plots_out_dir):
    merged_df = pd.read_csv(immediate_splicing_ratio_file_path, sep='\t')

    cutoff_dict = get_cutoff_dict(cutoff_file_path)
    control_cutoff = cutoff_file_path.split('.')[0].split('_')[-1]
    print(cutoff_dict)

    events_above_cutoff = get_events_above_cutoff(merged_df, cutoff_dict, sample_group)
    events_above_cutoff.to_csv(immediate_splicing_ratio_file_path.replace('.tsv', '_cutoff_{}.tsv'.format(control_cutoff)),
                         sep='\t', index=False)
    
    box_plot_splicing_ratio(events_above_cutoff, aka, 'boxplot_all_time_points', plots_out_dir, sample_group, control_cutoff)
    without_DRB = events_above_cutoff[events_above_cutoff['time_point'] != 'DRB']
    box_plot_splicing_ratio(without_DRB, aka, 'boxplot_splicing_ratio', plots_out_dir, sample_group, control_cutoff)
    
    density_plot_splicing_ratio(events_above_cutoff, aka, 'density_plot', plots_out_dir, sample_group, control_cutoff)

    plot_numb_events(events_above_cutoff, plots_out_dir, sample_group, control_cutoff)
    events_CTR_vs_W30 = get_CTR_W30_diff(events_above_cutoff)
    plot_CTR_W30_bx_with_lines(events_CTR_vs_W30, plots_out_dir, sample_group, control_cutoff)
    return events_above_cutoff, events_CTR_vs_W30

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input 2 files: (1) cutoff file, (2) immediate splicing ratio.
        Creates another file with the events on (2) that have total number reads > cutoff.
        Creates several plots based on the splicing ratio.""")

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

parser.add_argument('immediateSplicingRatioFilepath',
                    help="""Path to file with immediate splicing ratio.""",
                    metavar='<immediateSplicingRatio filepath>')

parser.add_argument('sampleGroup',
                    help="""Name of sample group (for eg: 601, 603, Y1P).""",
                    metavar='<sample group>')

args = parser.parse_args()

immediate_splicing_ratio_file_path = args.immediateSplicingRatioFilepath
cutoff_file_path = args.cutoffFilepath
plots_out_dir = args.outdir
sample_group = args.sampleGroup

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

main_function(immediate_splicing_ratio_file_path, sample_group, sample_group, cutoff_file_path, plots_out_dir)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED immediate splicing ratio plots for: {}\n{} minutes and {} seconds.\n'.format(immediate_splicing_ratio_file_path, minutes, seconds))



# def get_norm_factors_dict(sample_group, norm_factors_file_path):
#     df_1 = pd.read_csv(norm_factors_file_path)
#     df_2 = df_1[df_1['sample_group'] == sample_group][['sample_name', 'normalization_factor']]
#     norm_factors_dict = {}
#     for s in df_2['sample_name']:
#         norm_factors_dict[s] = df_2[df_2['sample_name'] == s]['normalization_factor'].values[0]
#     return norm_factors_dict

# def get_cutoff_based_on_control(norm_factors_dict, control_cutoff):
#     cutoff_dict = norm_factors_dict.copy()
#     constant = '%temp%'
#     for sample in norm_factors_dict:
#         if 'CTR' in sample:
#             constant = norm_factors_dict[sample] * control_cutoff
#     if constant == '%temp%':
#         constant = norm_factors_dict[list(norm_factors_dict.keys())[0]] * control_cutoff
#         print('There is no CTR sample. The first one will be taken as control.')
    
#     for sample in cutoff_dict:
#         cutoff_dict[sample] = constant / norm_factors_dict[sample]
    
#     return cutoff_dict