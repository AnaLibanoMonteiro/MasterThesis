#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Call functions to create the following bed files:

1. transcripts bed file (forward and reverse)
2. exons bed file (forward and reverse)
3. exons last coordinate bed file (forward and reverse)
4. exons last coordinate window bed file (forward and reverse)

5. introns (forward and reverse)
6. introns windows (forward and reverse)
"""

"""
TODO:
- Think about the user options 
  (about bed files, forward_reverse, window_size)
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import os
import argparse
import time

###########################################
#                                         #
#            Import functions             #
#                                         #
###########################################

from transcripts_and_exons_bed_files import get_merged_df
from transcripts_and_exons_bed_files import create_transcripts_bed_file
from transcripts_and_exons_bed_files import create_exons_bed_file

from introns_and_windows_bed_files import create_introns_bed_file
from introns_and_windows_bed_files import create_windows_bed_file

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def make_sure_path_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def get_out_file_path(out_file_prefix, in_file_path, out_dir):
    suffix = '{}{}.bed' ## Brackets will be replaced with transcripts/exons/introns and F/R if it is the case
    if out_file_prefix == '%filepath%': ## If the user didn't specify a name for the out file
        out_file_prefix = os.path.basename(in_file_path).split('.')[0] ## Then get the name of the input file
    if out_file_prefix == '':
        return os.path.join(out_dir, suffix)
    else:
        return os.path.join(out_dir, out_file_prefix + '_' + suffix)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input a file with a list of expressed transcripts ids and a GTF file.
        Returns the following bed files:
        1. transcripts
        2. exons
        3. exons last coordinate
        4. exons last coordinate window
        5. introns
        6. introns windows
        """)

parser.add_argument('-o', '--outfileprefix',
                    nargs='?',
                    default='%filepath%',
                    help="""Output file prefix (no extension). A suffix "{}.bed" will be added according to the bed file.
                    Default: basename of input file.""",
                    metavar='outfileprefix')

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('gtfFilepath',
                    nargs='?',
                    default='/home/analibanomonteiro/general/gtfEnsembl_90/Homo_sapiens.GRCh38.90.gtf.csv',
                    help="""Path to GTF file. Must be converted to csv before.
                    Default:
                    '/home/analibanomonteiro/general/gtfEnsembl_90/Homo_sapiens.GRCh38.90.gtf.csv'.""",
                    metavar='<gtf filepath>')
## '/home/analibanomonteiro/general/humanGenome/gtfFile/gencode.v38.primary_assembly.annotation.gtf.checked.csv' -> First time did with this one

parser.add_argument('transcriptsFilepath',
                    help="""Path to transcripts file. This file must have
                    one column with the ids of expressed transcripts.""",
                    metavar='<transcripts filepath>')

args = parser.parse_args()


gtf_file_path = args.gtfFilepath
transcripts_file_path = args.transcriptsFilepath

out_file_prefix = args.outfileprefix
out_dir = make_sure_path_exists(args.outdir)

out_file_path = get_out_file_path(out_file_prefix, transcripts_file_path, out_dir)



###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

## Get merged df:
merged_df = get_merged_df(transcripts_file_path, gtf_file_path)
print('STARTING...\n')
############################################################

## Create transcripts bed file: (forward and reverse)
print('...transcripts bed files')
create_transcripts_bed_file(merged_df, out_file_path)

## Create exons bed file (including last coordinates): (forward and reverse)
print('\n...exons bed files')
create_exons_bed_file(merged_df, out_file_path)

############################################################

exon_file_path_F = out_file_path.format('exons', '_F')
exon_file_path_R = out_file_path.format('exons', '_R')

############################################################

## Create introns bed file (based on exons bed files): (forward and reverse)
print('\n...introns bed files')
create_introns_bed_file(exon_file_path_F, out_file_path, 'F')
create_introns_bed_file(exon_file_path_R, out_file_path, 'R')

## Create introns windows bed file (based on exons bed files): (forward and reverse)
print('\n...windows bed files')
create_windows_bed_file(exon_file_path_F, out_file_path, 'F', 10)
create_windows_bed_file(exon_file_path_R, out_file_path, 'R', 10)


end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('\nFINISHED bed files in \n{} minutes and {} seconds.\n'.format(minutes, seconds))