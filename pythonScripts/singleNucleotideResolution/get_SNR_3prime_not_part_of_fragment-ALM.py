#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Get Single Nucleotide Resolution (SNR) from a bam file.
LAST base from 147.
FIRST base from 163.
Accepts reads with soft clipping (if it is only 1 nucleotide soft clipped).
Considers the soft clipped nucleotides as part of the read.

This script was written based on https://github.com/rluis/mNET_snr.git
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import os
import time

import pysam  ## ATTENTION: pysam is always 0-based, as BAM files

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_out_file_path(out_file_name, in_file_path, out_dir):
    suffix = '_SNR3prime_sorted.bam'
    if out_file_name == '%filepath%': ## If the user didn't specify a name for the out file
        out_file_name = os.path.basename(in_file_path).split('.')[0] ## Get the name of the input file
    return os.path.join(out_dir, out_file_name + suffix)

def check_cigar(read):
    """Accept all reads as long as they don't have insertions or deletions;
       Or in case of soft clipping, they have only 1 soft clipped nucleotide
    """
    cigar = read.cigarstring
    if 'I' in cigar or 'D' in cigar:
        return False
    elif cigar.count('S') > 1:
        return False
    else:
        return True

# def get_new_start_147(read):
#     if read.cigartuples[-1][0] == 4: ## Soft clipping at the end of the read
#         return read.reference_end ## - 1 + 1 (Add 1 to the start position (will cut the previous - 1))
#     else:
#         return read.reference_end - 1 ## reference_end points to one past the last aligned residue, hence '- 1'

# def get_new_start_163(read):
#     if read.cigartuples[0][0] == 4: ## Soft clipping at the beggining of the read
#         return read.reference_start - 1
#     else:
#         return read.reference_start ## If no soft clipping, start is the same

def get_147_snr(read): ## Get last base from read (last position is the same, change start position)
    """Get last base from read. Last position is the same; change start position (read.reference_start)
    If soft clipping at the end of the read: add 1 to the start position after changing it
    """
    q = read.query_qualities
    
    ## Previous version (reads with soft clipping accepted, but soft clipped nucleotides ignored)
    end = read.query_alignment_end ## Save end of aligment in a variable
    read.query_sequence = read.query_sequence[end - 1] ## This way we are ready to accept reads with soft-clipping at pol position
    read.query_qualities = q[end - 1:end]
    read.reference_start = read.reference_end - 1
    read.flag = 99

    ## New version (reads with soft clipping accepted, soft nucleotides considered as part of the read):
    # read.query_sequence = read.query_sequence[-1]
    # read.query_qualities = q[-1:]
    # read.flag = 99
    # read.reference_start = get_new_start_147(read)
    
    read.cigarstring = '1M' ## Detail: setting cigar string changes reference_end, so reference_start has to be set first.
    read.template_length = read.template_length / abs(read.template_length) ## I think that I would change this, removing abs() - this came from original script
    return read

def get_163_snr(read): ## Get first base from read (first position is the same, it is not necessary to change last position)
    """Get first base from read. First position is the same;
    it is not necessary to change last position (when cigar string is changed to '1M', reference_end is changed automatically)
    If soft clipping at the beggining of the read: subtract 1 to the start position
    """
    q = read.query_qualities
    
    ## Previous version (reads with soft clipping accepted, but soft clipped nucleotides ignored)
    start = read.query_alignment_start
    read.query_sequence = read.query_sequence[start]
    read.query_qualities = q[start:start + 1]
    read.flag = 83

    ## New version (reads with soft clipping accepted, soft nucleotides considered as part of the read):
    # read.query_sequence = read.query_sequence[0]
    # read.query_qualities = q[0:1]
    # read.flag = 83
    # read.reference_start = get_new_start_163(read)
    
    read.cigarstring = '1M'
    read.template_length = read.template_length / abs(read.template_length)
    return read

def get_single_nucleotide_resolution(input_bam_path, snr_bam_file_path):
    input_bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    snr_bam_file = pysam.AlignmentFile(snr_bam_file_path, "wb", template=input_bam_file)
    flags = [147, 163]
    for read in input_bam_file:
        if read.flag in flags and check_cigar(read): ### NEW FILTER!!!!!
        # if read.flag in flags and (not "I" in read.cigarstring and not "D" in read.cigarstring and not "S" in read.cigarstring): ## Original filter
            if read.flag == 147:
                read = get_147_snr(read)
            elif read.flag == 163:
                read = get_163_snr(read)
            snr_bam_file.write(read)

    snr_bam_file.close()
    input_bam_file.close()

    pysam.sort('-o', snr_bam_file_path, snr_bam_file_path)
    pysam.index(snr_bam_file_path)

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description=""" Get single nucleotide resolution (with my script)
        Takes as input a bam file. Creates another bam file with the last coordinate of each read.""")

parser.add_argument('-o', '--outfilename',
                    nargs='?',
                    default='%filepath%',
                    help="""Output file prefix (no extension). A suffix "_SNR3prime_sorted.bam" will be added.
                    Default: basename of input file.""",
                    metavar='outfilename')

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('bamfilepath',
                    help="""Path to input bam file.""",
                    metavar='<input bam filepath>')

args = parser.parse_args()

input_bam_path = args.bamfilepath
out_file_name = args.outfilename
out_dir = args.outdir

out_file_path = get_out_file_path(out_file_name, input_bam_path, out_dir)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to get SNR for {}\n'.format(input_bam_path))

get_single_nucleotide_resolution(input_bam_path, out_file_path)

temp = time.time()-start
minutes = temp//60
seconds = temp - 60*minutes
print('FINISHED SNR with new filter for: {} \nWritten in: {}. \n{} minutes and {} seconds.\n'.format(input_bam_path, out_file_path, minutes, seconds))




# def check_cigar(read): ## Accept reads that have soft clipping at the opposite end of pol position
#     cigar = read.cigarstring
#     if 'I' in cigar or 'D' in cigar:
#         return False ## Reads with Insertions or Deletions will be discarted
    
#     elif 'S' in cigar:
#         if cigar.count('S') == 2:
#             return False ## Reads with soft clippings in both sides will be discarted
        
#         flag = read.flag ## If only one end has soft clipping, the read may be approved according to flag
#         if flag == 147 and cigar[-1] != 'S': ## Flag 147 cannot have soft clipping at the end
#             return True
#         else:
#             return False
#         if flag == 163 and cigar[-1] == 'S': ## Flag 163 can have soft clipping at the end
#             return True
#         else:
#             return False
    
#     else: ## If the read does not have 'I', 'D' or 'S', it will pass
#         return True
