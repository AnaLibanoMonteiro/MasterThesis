#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
From cell 2018:

To identify spikes in the density of 3' ends of nascent transcripts we used an algorithm that finds
nucleotides where the read density is at least three standard deviations above the mean in a local region (Churchman and Weissman, 2011).
Only positions with coverage of 4 or more reads were considered.

To identify peaks at the 3' end of exons and introns indicative of splicing intermediates, the number of
reads at these positions was compared to the mean read density across the corresponding exon or intron.


To identify peaks corresponding to Pol II pause positions along exons and introns avoiding contamination
by 5' ss splicing intermediates, we removed reads aligning to the last 3 nucleotides of exons and the
first 3 nucleotides of introns.
The number of reads at each nucleotide position along the exon was then compared to the mean read density
across the entire exon.
Peaks were identified in the annotated exons and introns of expressed genes.
Exons that intersected other isoform exons were discarded and the same principle was applied for introns
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import time

import pandas as pd

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_out_file_name(intermediates_coverage_file_path):
    new_suffix = 'peaks.tsv'
    if 'coverage_filtered.tsv' in intermediates_coverage_file_path:
        return intermediates_coverage_file_path.replace('coverage_filtered.tsv', new_suffix)
    elif 'coverage.tsv' in intermediates_coverage_file_path:
        return intermediates_coverage_file_path.replace('coverage.tsv', new_suffix)
    else:
        return intermediates_coverage_file_path.replace('.tsv', new_suffix)

def get_df_local_mean_plus_3_std(df_windows_coverage):
    """Receives a data frame with the coverage of the 200 nucleotides
    surrounding each exon last coordinate (or splicing intermediates).
    Returns a data frame with the local mean plus 3 * standard deviation
    for each event.
    """
    surrounding_coverage_mean = df_windows_coverage.groupby(['name'])['coverage_counts'].mean().reset_index(name='coverage_mean')
    surrounding_coverage_std = df_windows_coverage.groupby(['name'])['coverage_counts'].std().reset_index(name='coverage_std')
    
    df_temp = pd.merge(surrounding_coverage_mean, surrounding_coverage_std, on=['name'])
    df_cutoff = pd.merge(df_temp, df_windows_coverage[['name', 'sample_name']].drop_duplicates(), on=['name'], how='left')
    df_cutoff['local_mean_plus_3_std'] = df_cutoff.apply(lambda row: row.coverage_std * 3 + row.coverage_mean, axis=1)
    
    return df_cutoff[['name', 'local_mean_plus_3_std']]

def get_splicing_intermediates_peaks(intermediates_coverage_file_path, windows_coverage_file_path,
                                     sample_group):
    df_intermediates_coverage = pd.read_csv(intermediates_coverage_file_path, sep='\t')
    df_windows_coverage = pd.read_csv(windows_coverage_file_path, sep='\t')

    df_local_mean_plus_3_std = get_df_local_mean_plus_3_std(df_windows_coverage)
    df_merged = pd.merge(df_intermediates_coverage, df_local_mean_plus_3_std, on=['name'])
    df_splicing_intermediates_peaks = df_merged[df_merged['coverage_counts'] >= df_merged['local_mean_plus_3_std']]
    out_file_name = get_out_file_name(intermediates_coverage_file_path)
    df_splicing_intermediates_peaks.to_csv(out_file_name, index=False, sep='\t')
    return df_splicing_intermediates_peaks


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input 2 files: (1) splicing intermediates coverage and
        (2) surrounding coverage, created with 'bedtools coverage'; and the sample group name.
        Returns the events on the file (1) that have coverage significantly higher than
        surrounding nucleotides in file (2).""")

parser.add_argument('intermediatesCoverageFilepath',
                    help="""Path to file with splicing intermediates coverage.""",
                    metavar='<splicingIntermediatesCoverage filepath>')

parser.add_argument('windowsCoverageFilepath',
                    help="""Path to file with surrounding nucleotides coverage.""",
                    metavar='<windowsCoverage filepath>')

parser.add_argument('sampleGroup',
                    help="""Name of sample group (for eg: 601, 603, Y1P).""",
                    metavar='<sample group>')

args = parser.parse_args()

intermediates_coverage_file_path = args.intermediatesCoverageFilepath
windows_coverage_file_path = args.windowsCoverageFilepath
sample_group = args.sampleGroup

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

print('STARTING to get splicing intermediates peaks')

get_splicing_intermediates_peaks(intermediates_coverage_file_path, windows_coverage_file_path,
                                 sample_group)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED getting splicing intermediates peaks for: {} \n{} minutes and {} seconds.\n'.format(intermediates_coverage_file_path, minutes, seconds))
