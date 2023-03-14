#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Input: complete soft clipped info.
read_name | chr | mapping_start | read_size | cigar | read_seq | soft_start | soft_end | numb_soft_nucleotides | soft_nucleotide_seq | one_or_both_ends | prime_end | flag | strand | ref_nucleotide_seq | reference_>_soft
Output: several plots describing soft clipping info.
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

#####################
# Functions to Plot #
#####################

def save_plot(plot, plots_out_dir, sample_name, plot_name):
  file_name = '{}_{}.png'.format(sample_name, plot_name)
  path_to_save_plot = os.path.join(plots_out_dir, file_name)
  plot.figure.savefig(path_to_save_plot, dpi=100, bbox_inches='tight')

def save_data_frame(df, plots_out_dir, sample_name, plot_name):
  file_name = '{}_{}.csv.gz'.format(sample_name, plot_name)
  path_to_save_df = os.path.join(plots_out_dir, file_name)
  df.to_csv(path_to_save_df, index=False, compression={'method':'gzip'})

def plot_length_soft_clipped_regions(df_soft_clipped_info, plots_out_dir, sample_name):
  ## Dataframe to plot:
  numb_soft_nucleotides = df_soft_clipped_info.groupby('numb_soft_nucleotides')['numb_soft_nucleotides'].count().reset_index(name='count')
  ## Prepare variables to plot text on the plot:
  percentage_length_1 = numb_soft_nucleotides[numb_soft_nucleotides['numb_soft_nucleotides']==1]['count'][0] / numb_soft_nucleotides['count'].sum() * 100
  text_to_add = '{}% of Soft Clipped Regions have only 1 Nucleotide'.format(round(percentage_length_1,2))
  x_pos = len(numb_soft_nucleotides.index)/10
  y_pos = numb_soft_nucleotides['count'].max()/2
  ## Make plot:
  with plt.style.context('seaborn-poster'):
    plot = numb_soft_nucleotides.plot.bar(x='numb_soft_nucleotides', y='count', figsize=(13,5),
                                        xlabel='Numb soft clipped nucleotides')
    plt.xticks(fontsize=13, rotation=90)
    plot.figure.suptitle('Lenght Soft Clipped Regions', fontsize = 20)
    ## Add text:
    plt.text(x_pos, y_pos, text_to_add, fontsize = 18)

  ## Save plot:
  plot_name = 'length_soft_clipped_regions'
  save_plot(plot, plots_out_dir, sample_name, plot_name)
  save_data_frame(numb_soft_nucleotides, plots_out_dir, sample_name, plot_name)

def plot_soft_clipping_patterns(df_soft_clipped_info, plots_out_dir, sample_name, prime_end=None, top=5):
  title = 'Soft Clipping Patterns{}'
  plot_name = 'soft_clipping_patterns{}'
  if prime_end == None:
    title = title.format('')
    plot_name = plot_name.format('')
  elif int(prime_end) == 5:
    title = title.format(' (5 Prime End)')
    plot_name = plot_name.format('_5prime_end')
  elif int(prime_end) == 3:
    title = title.format(' (3 Prime End)')
    plot_name = plot_name.format('_3prime_end')
  else:
    raise NameError('Argument "prime_end" can be either "None" (default), "5" or "3".')
  
  ## Dataframe to plot:
  unique_soft_seq_and_ref = df_soft_clipped_info.groupby('reference_>_soft')['reference_>_soft'].count().reset_index(name='count').sort_values('count', ascending=False)
  ## Select only the top patterns:
  top_unique_soft_seq_and_ref = unique_soft_seq_and_ref[0:top].copy()
  ## Add space between letters:
  top_unique_soft_seq_and_ref.rename(columns={'reference_>_soft':'reference_and_soft'}, inplace=True)
  top_unique_soft_seq_and_ref['reference_and_soft'] = top_unique_soft_seq_and_ref.apply(lambda row: row.reference_and_soft.replace('>', ' > '), axis=1)

  ## Make plot:
  with plt.style.context('seaborn-poster'):
    plot = top_unique_soft_seq_and_ref.plot.bar(x='reference_and_soft', y='count', rot=0, figsize=(10,5),
                                               xlabel=('Reference Nucleotide VS Soft Clipped Nucleotide'))
    plot.figure.suptitle(title, fontsize = 20)
  
  ## Save plot:
  save_plot(plot, plots_out_dir, sample_name, plot_name)
  save_data_frame(unique_soft_seq_and_ref, plots_out_dir, sample_name, plot_name)

def plot_soft_clipping_per_chrm(df_soft_clipped_info, plots_out_dir, sample_name, top=25):
  ## Create dataframe:
  count_chr = df_soft_clipped_info.groupby('chr')['chr'].count().reset_index(name='count').sort_values('count', ascending=False)
  ## Make plot:
  with plt.style.context('seaborn-poster'):
    plot = count_chr[0:top].plot.bar(x='chr', y='count', figsize=(13,5), xlabel='Chromosome')
    plt.xticks(fontsize=13)
    plot.figure.suptitle('Soft Clipping per Chromosome', fontsize = 20)
  ## Save plot:
  plot_name = 'soft_clipping_per_chromosome'
  save_plot(plot, plots_out_dir, sample_name, plot_name)

def plot_soft_clipping_per_flag(df_soft_clipped_info, plots_out_dir, sample_name):
  ## Create dataframe:
  count_flag = df_soft_clipped_info.groupby('flag')['flag'].count().reset_index(name='count').sort_values('count', ascending=False)
  ## Make plot:
  with plt.style.context('seaborn-poster'):
    plot = count_flag.plot.bar(x='flag', y='count', figsize=(13,5), xlabel='Bam Flag')
    plt.xticks(fontsize=13)
    plot.figure.suptitle('Soft Clipping per Flag', fontsize = 20)
  ## Save plot:
  plot_name = 'soft_clipping_per_flag'
  save_plot(plot, plots_out_dir, sample_name, plot_name)
  save_data_frame(count_flag, plots_out_dir, sample_name, plot_name)

def plot_count_per_prime_end(df_5prime, df_3prime, plots_out_dir, sample_name):
  data=[['5 Prime', len(df_5prime.index)], ['3 Prime', len(df_3prime.index)]]
  count_per_prime_end_with_pairs = pd.DataFrame(data, columns=['prime_end', 'count'])
  ## Make plot
  with plt.style.context('seaborn-poster'):
    plot = count_per_prime_end_with_pairs.plot.bar(x='prime_end', y='count', rot=0, xlabel='', figsize=(8,5))
    plot.figure.suptitle('Soft Clipping per Prime End', fontsize = 20)

  ## Save plot:
  plot_name = 'soft_clipping_per_prime_end'
  save_plot(plot, plots_out_dir, sample_name, plot_name)
  save_data_frame(count_per_prime_end_with_pairs, plots_out_dir, sample_name, plot_name)

def plot_number_overlapped_reads(df_soft_clipped_info, df_5prime_only_pairs_overlapped, df_3prime_only_pairs_overlapped, plots_out_dir, sample_name):
  total_number_reads = len(set(df_soft_clipped_info['read_name']))
  overlapped_reads_5_prime = set(df_5prime_only_pairs_overlapped['read_name'])
  overlapped_reads_3_prime = set(df_3prime_only_pairs_overlapped['read_name'])
  number_overlapped = len(overlapped_reads_5_prime.union(overlapped_reads_3_prime))
  number_non_overlapped = total_number_reads - number_overlapped

  data = [['Completely Overlapped', number_overlapped], ['Others', number_non_overlapped]]
  overlapped_reads = pd.DataFrame(data, columns=['overlap', 'number'])
  ## Make plot
  with plt.style.context('seaborn-poster'):
    plot = overlapped_reads.plot.bar(x='overlap', y='number', rot=0, xlabel='', figsize=(8,5))
    plot.figure.suptitle('Completely Overlapped Reads', fontsize = 20)

  ## Save plot:
  plot_name = 'completely_overlapped_reads'
  save_plot(plot, plots_out_dir, sample_name, plot_name)
  save_data_frame(overlapped_reads, plots_out_dir, sample_name, plot_name)

def plot_match_R1_R2(df_only_pairs_overlapped, plots_out_dir, sample_name, prime_end):
  title = 'Soft Clipping Match ({})'
  plot_name = 'soft_clipping_match_{}'
  if int(prime_end) == 5:
    title = title.format('5 Prime End')
    plot_name = plot_name.format('5prime_end{}')
  elif int(prime_end) == 3:
    title = title.format('3 Prime End')
    plot_name = plot_name.format('3prime_end{}')
  else:
    raise NameError('Define prime_end argument. Either "5" or "3".')

  df_temp = df_only_pairs_overlapped.copy()
  df_temp['n_unique_reference_>_soft_per_prime_end'] = df_temp.groupby(['read_name'])['reference_>_soft'].transform('nunique')

  df_count_unique_seq = df_temp.groupby('n_unique_reference_>_soft_per_prime_end')['n_unique_reference_>_soft_per_prime_end'].count().reset_index(name='count').sort_values('count', ascending=False)
  
  ## Bar plot:
  with plt.style.context('seaborn-poster'):
    plot = df_count_unique_seq.plot.bar(x='n_unique_reference_>_soft_per_prime_end', y='count', rot=0,
                                            xlabel='Number unique soft clipped sequences', figsize=(8,5))
    plot.figure.suptitle(title, fontsize = 20)
  save_plot(plot, plots_out_dir, sample_name, plot_name.format('_bar_plot'))
  
  ## Pie plot:
  my_labels = ['R1 and R2\n matched', 'R1 and R2 \n not matched']
  plot = df_count_unique_seq.plot.pie(y='count', figsize=(6, 6), ylabel='', legend=None, labels=my_labels,
                                        colors=['lightgreen', 'lightcoral'], autopct='%.1f%%', fontsize=15)
  plot.figure.suptitle(title, fontsize = 20, y=0.88)
  save_plot(plot, plots_out_dir, sample_name, plot_name.format('_pie_plot'))

  save_data_frame(df_count_unique_seq, plots_out_dir, sample_name, plot_name.format(''))

#################
# Main Function #
#################

def soft_clipping_plots(csv_file_path, sample_name, plots_out_dir):
  ## Read file:
  df_soft_clipped_info = pd.read_csv(csv_file_path, compression={'method':'gzip'})

  ##############
  ## SUBSETS: ##
  ##############
  df_5prime = df_soft_clipped_info[(df_soft_clipped_info['prime_end'] == '5prime') &
                ((df_soft_clipped_info['flag'] == 99) | (df_soft_clipped_info['flag'] == 83))]
  
  df_3prime = df_soft_clipped_info[(df_soft_clipped_info['prime_end'] == '3prime') &
                ((df_soft_clipped_info['flag'] == 147) | (df_soft_clipped_info['flag'] == 163))]

  df_middle_5prime = df_soft_clipped_info[(df_soft_clipped_info['prime_end'] == '5prime') &
                ((df_soft_clipped_info['flag'] == 147) | (df_soft_clipped_info['flag'] == 163))]

  df_middle_3prime = df_soft_clipped_info[(df_soft_clipped_info['prime_end'] == '3prime') &
                ((df_soft_clipped_info['flag'] == 99) | (df_soft_clipped_info['flag'] == 83))]

  ## Add info from read 2 if pair overlapps:
  ### First select reads that can have a pair overlap (unique mapping start)
  unique_temp = df_soft_clipped_info[['read_name', 'mapping_start']].copy()
  unique_temp['n_unique_mapping_start'] = unique_temp.groupby(['read_name'])['mapping_start'].transform('nunique')
  set_unique_mapping_start_reads = set(unique_temp[(unique_temp['n_unique_mapping_start'] == 1)]['read_name'])
  ### Then add info for 5 prime:
  set_unique_mapping_start_5prime = set_unique_mapping_start_reads.intersection(set(df_5prime['read_name']))
  df_5prime_approved_pairs = df_middle_5prime[df_middle_5prime['read_name'].isin(set_unique_mapping_start_5prime)]
  df_5prime_with_pairs = pd.concat([df_5prime, df_5prime_approved_pairs])
  ### And for 3 prime:
  set_unique_mapping_start_3prime = set_unique_mapping_start_reads.intersection(set(df_3prime['read_name']))
  df_3prime_approved_pairs = df_middle_3prime[df_middle_3prime['read_name'].isin(set_unique_mapping_start_3prime)]
  df_3prime_with_pairs = pd.concat([df_3prime, df_3prime_approved_pairs])

  ## NOTE: THESE TWO DATASETS (df_5prime_with_pairs AND df_3prime_with_pairs)
  # HAVE ALL READ 1 INFORMATION ABOUT THE PRIME END IN QUESTION **PLUS**
  # THE INFORMATION FOR THOSE ENTRIES THAT HAVE A PAIR THAT OVERLAPS

  ## Now select ONLY the reads that have a pair:
  ### For 5 prime end:
  df_5prime_with_pairs_count = df_5prime_with_pairs.groupby(['read_name'])['read_name'].count().reset_index(name='count').sort_values('count', ascending=False)
  set_5prime_2_events = set(df_5prime_with_pairs_count[df_5prime_with_pairs_count['count']==2]['read_name'])
  df_5prime_only_pairs_overlapped = df_5prime_with_pairs[df_5prime_with_pairs['read_name'].isin(set_5prime_2_events)]
  ### And for 3 prime:
  df_3prime_with_pairs_count = df_3prime_with_pairs.groupby(['read_name'])['read_name'].count().reset_index(name='count').sort_values('count', ascending=False)
  set_3prime_2_events = set(df_3prime_with_pairs_count[df_3prime_with_pairs_count['count']==2]['read_name'])
  df_3prime_only_pairs_overlapped = df_3prime_with_pairs[df_3prime_with_pairs['read_name'].isin(set_3prime_2_events)]

  ###########
  ## PLOTS ##
  ###########

  plot_length_soft_clipped_regions(df_soft_clipped_info, plots_out_dir, sample_name)

  plot_soft_clipping_patterns(df_soft_clipped_info, plots_out_dir, sample_name)

  plot_soft_clipping_per_chrm(df_soft_clipped_info, plots_out_dir, sample_name)

  plot_soft_clipping_per_flag(df_soft_clipped_info, plots_out_dir, sample_name)

  plot_count_per_prime_end(df_5prime_with_pairs, df_3prime_with_pairs, plots_out_dir, sample_name)

  plot_number_overlapped_reads(df_soft_clipped_info, df_5prime_only_pairs_overlapped, df_3prime_only_pairs_overlapped, plots_out_dir, sample_name)

  plot_match_R1_R2(df_5prime_only_pairs_overlapped, plots_out_dir, sample_name, prime_end=5)
  plot_match_R1_R2(df_3prime_only_pairs_overlapped, plots_out_dir, sample_name, prime_end=3)

  plot_soft_clipping_patterns(df_5prime, plots_out_dir, sample_name, prime_end=5)
  plot_soft_clipping_patterns(df_3prime, plots_out_dir, sample_name, prime_end=3)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a csv file with info about soft reads.
                                                 Makes some plots based on that table.""")

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('sampleName',
                    help="""Sample Name.""",
                    metavar='<sample Name>')

parser.add_argument('csvfilepath',
                    help="""Path to input csv file.""",
                    metavar='<input csv filepath>')

args = parser.parse_args()

plots_out_dir = args.outdir
sample_name = args.sampleName
csv_file_path = args.csvfilepath

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to create soft clipped plots for {}\n'.format(csv_file_path))

soft_clipping_plots(csv_file_path, sample_name, plots_out_dir)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED this step for: {} \n{} minutes and {} seconds.\n'.format(csv_file_path, minutes, seconds))
