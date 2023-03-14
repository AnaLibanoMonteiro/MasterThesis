#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
COMPLETE soft clipped info (with sequence from genome reference)
Input: incomplete soft clipped info and soft clipped fasta file
(FASTA FILE WILL BE CREATED WITH bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>
FROM BED FILE CREATED WITH PREVIOUS SCRIPT)
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

def str_to_bool(string):
  return string.lower() in ("true", "t", "yes", "1", "replace")

def read_bed_fasta(fasta_file_path):
    return pd.read_csv(fasta_file_path, sep='\t', header=None, names=['chr', 'soft_start', 'soft_end', 'bed_entry', 
                                                                      'score', 'strand', 'ref_nucleotide_seq']).drop(['bed_entry', 'score'], axis=1)

def complete_soft_clipped_info_csv(csv_file_path, fasta_file_path, replace=True):
    
    ## Read files:
    df_csv = pd.read_csv(csv_file_path, compression={'method':'gzip'})
    df_fasta = read_bed_fasta(fasta_file_path)
        
    ## Merge files:
    df_merged = pd.merge(df_csv, df_fasta, on=["chr", "soft_start", "soft_end", "strand"], how='left')
    
    ## Add column with both sequences (soft and reference)
    # df_merged['soft_>_reference'] = df_merged.apply(lambda row: '{}>{}'.format(row.soft_nucleotide_seq, row.ref_nucleotide_seq), axis=1)
    df_merged['reference_>_soft'] = df_merged.apply(lambda row: '{}>{}'.format(row.ref_nucleotide_seq, row.soft_nucleotide_seq), axis=1)
    
    ## Save result:
    df_merged.to_csv(csv_file_path.replace('.csv.gz', '_complete.csv.gz'), index=False, compression={'method':'gzip'})
    if replace:
        os.remove(csv_file_path)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a csv file with info about soft reads and a fasta file.
                                                Adds two extra columns ['ref_nucleotide_seq', 'soft_>_reference'].
                                                By default saves the output on the same file.""")

parser.add_argument('-r', '--replace',
                    nargs='?',
                    default='True',
                    help="""Boolean. True: remove original file. False: keep both files (original and _complete.csv).
                    Default: True.""",
                    metavar='replace')

parser.add_argument('csvfilepath',
                    help="""Path to input csv file.""",
                    metavar='<input csv filepath>')

parser.add_argument('fastafilepath',
                    help="""Path to input fasta file.""",
                    metavar='<input fasta filepath>')

args = parser.parse_args()

csv_file_path = args.csvfilepath
fasta_file_path = args.fastafilepath
replace = str_to_bool(args.replace)


###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to complete soft clipped info for {}\n'.format(csv_file_path))

complete_soft_clipped_info_csv(csv_file_path, fasta_file_path, replace)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED this step for: {} \n{} minutes and {} seconds.\n'.format(csv_file_path, minutes, seconds))