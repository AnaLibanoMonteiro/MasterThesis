#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
For each sample:
Creates a csv file with number reads per step counted with batch file 'count_reads.sbatch'.
Creates plots based on those values (if plot is True).
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import os
import argparse
import pandas as pd
import time

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_csv_out_path(in_file_path, sample_name):
    ## csv file will be saved on the same directory of input file:
    input_file_directory = os.path.split(number_reads_file_path)[0]
    return os.path.join(input_file_directory, sample_name + '_' + 'number_reads.csv')

def str_to_bool(string):
    return string.lower() in ("true", "t", "yes", "1", "plot", "p", "replace", "r")

def check_number(x):
    if x.isdigit():
        return int(x)
    else:
        return x

def add_columns_df(df):
    df_temp = df.copy()
    df_temp['primary_and_mapped'] = df_temp['star_bam_file-primary_and_mapped_x2']/2
    df_temp['reads_soft_clipped'] = df_temp['star_bam_file-soft_clipped_x2']/2
    df_temp['primary_and_mapped_not_sf'] = df_temp['primary_and_mapped'] - df_temp['reads_soft_clipped']
    df_temp['reads_unmapped'] = df_temp['fasta_1'] - df_temp['primary_and_mapped']
    df_temp['reads_IPri'] = df_temp['primary_and_mapped'] - df_temp['star_bam_woIPri']
    df_temp['reads_deletions_insertions'] = df_temp['star_bam_woIPri'] - df_temp['star_bam_woIPri_and_SNR']
    df_temp['soft_clipped_%'] = df_temp['reads_soft_clipped']/df_temp['primary_and_mapped']*100
    df_temp['final_%'] = df_temp['star_bam_woIPri_and_SNR']/df_temp['fasta_1']*100
    return df_temp

def read_count_read_files(file_path):
    """Returns a dataframe with values from count read files
    """
    with open(file_path, 'r') as f:
        content = f.read().strip('\n').split('\n')
    col_names = []
    for i in range(0,len(content),2):
        col_names.append(content[i])
    
    values = []
    for i in range(1,len(content),2):
        values.append([check_number(content[i])])
    
    my_dict = dict(zip(col_names, values))
    df_temp = pd.DataFrame.from_dict(my_dict)
    df = add_columns_df(df_temp)
    return df

def save_plot(plot, plots_out_dir, sample_name, plot_name):
    file_name = '{}_{}.png'.format(sample_name, plot_name)
    path_to_save_plot = os.path.join(plots_out_dir, file_name)
    plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

def plot_reads_per_step(df, sample_name, plots_out_dir):
    reads_per_step = ['fasta_1', 'fasta_2', 'primary_and_mapped', 'star_bam_woIPri', 'star_bam_woIPri_and_SNR']
    df_reads_per_step = df[reads_per_step].transpose()
    df_reads_per_step.columns = ['number_reads']
    plot = df_reads_per_step.plot.bar(legend=False, title='Numb reads per step {}'.format(sample_name), width=0.65, figsize=(8,5), rot=10)
    plot_name = 'numb_reads_per_step'
    save_plot(plot, plots_out_dir, sample_name, plot_name)
    return df_reads_per_step

def plot_lost_reads(df, sample_name, plots_out_dir):
    lost_reads = ['star_bam_woIPri_and_SNR', 'reads_unmapped', 'reads_IPri', 'reads_deletions_insertions']
    df_lost_reads = df[lost_reads].transpose()
    df_lost_reads.columns = ['number_reads']
    my_labels = ['Accepted Reads', 'Unmapped', 'Internal Primming', 'Deletions and Insertions']
    plot = df_lost_reads.plot.pie(y='number_reads', title='Accepted vs Lost reads {}'.format(sample_name),
         figsize=(6, 6), labels=my_labels, ylabel='', legend=None, autopct='%.1f%%',
         explode=[0.2, 0.1, 0.1, 0.1])
    plot_name = 'accepted_vs_lost_reads'
    save_plot(plot, plots_out_dir, sample_name, plot_name)
    return df_lost_reads

def plot_soft_clipped(df, sample_name, out_path_template):
    soft_clipped = ['primary_and_mapped_not_sf', 'reads_soft_clipped']
    df_soft_clipped = df[soft_clipped].transpose()
    df_soft_clipped.columns = ['number_reads']
    # my_labels=['Primary and mapped\nnot soft clipped', 'Soft clipped reads']
    my_labels=None
    # colors=['lightgreen', 'lightcoral']
    my_colors=['khaki', 'darkorange']
    plot = df_soft_clipped.plot.pie(y='number_reads', title='Soft clipped reads {}'.format(sample_name),
         figsize=(6, 6), explode=[0.1, 0], ylabel='', legend=None, autopct='%.1f%%', labels=my_labels,
         colors=my_colors, fontsize=16)
    plot_name = 'soft_clipped'
    save_plot(plot, plots_out_dir, sample_name, plot_name)
    return df_soft_clipped

def main_function(number_reads_file_path, plots_out_dir, plot, replace):
    df = read_count_read_files(number_reads_file_path)
    sample_name = df['sample_name'][0]
    csv_out_path = get_csv_out_path(number_reads_file_path, sample_name)
    df.to_csv(csv_out_path, index=False)
    if replace:
        os.remove(number_reads_file_path)
    if plot:
        plot_reads_per_step(df, sample_name, plots_out_dir)
        plot_lost_reads(df, sample_name, plots_out_dir)
        plot_soft_clipped(df, sample_name, plots_out_dir)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input a file with number of reads per each step for a specific sample.
        Creates a csv file with number reads per step. Creates plots based on those values (if plot is True).""")

parser.add_argument('-d','--plotsOutdir',
                    nargs='?',
                    default="./",
                    help="Directory to save plots. Default: current directory.",
                    metavar='plotsOutdir')

parser.add_argument('-p','--plot',
                    nargs='?',
                    default='True',
                    help="Boolean. True: Creates csv file and plots in the same folder of input file. False: Creates only csv file. Default: True.",
                    metavar='plot')

parser.add_argument('-r','--replace',
                    nargs='?',
                    default='True',
                    help="Boolean. True: Removes orignal file with number reads. False: Creates csv file and keep original file. Default: True.",
                    metavar='plot')

parser.add_argument('numberReadsFile',
                    help="""Path to input file. This file must have
                    sample name, folder name, number of reads for each step.""",
                    metavar='<numberReads filepath>')

args = parser.parse_args()

plot = str_to_bool(args.plot)
replace = str_to_bool(args.replace)
number_reads_file_path = args.numberReadsFile
plots_out_dir = args.plotsOutdir

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

main_function(number_reads_file_path, plots_out_dir, plot, replace)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED expressing number reads for: {} \n{} minutes and {} seconds.\n'.format(number_reads_file_path, minutes, seconds))
