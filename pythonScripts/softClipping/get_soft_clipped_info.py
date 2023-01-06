#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Get soft clipped info from mapped and primary reads bam file.
Creates BAM file, CSV file and BED file.
Uses pysam to parse the reads.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import os
import argparse
import pysam ## ATTENTION: pysam is always 0-based, as BAM files and BED files
import pandas as pd
import time
from Bio.Seq import Seq

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def str_to_bool(string):
    return string.lower() in ("true", "t", "yes", "1", "plot", "p", "replace", "r")

def get_out_file_path(out_file_name, in_file_path, out_dir):
    suffix = '_soft_clipped_reads{}' ## {} will be replaced with .format() according to type of file that is being created.
    if out_file_name == '%filepath%': ## If the user didn't specify a name for the out file
        out_file_name = os.path.basename(in_file_path).split('.')[0] ## Then get the name of the input file
    return os.path.join(out_dir, out_file_name + suffix)

def get_soft_clipped_bam_path(input_bam_path):
    return input_bam_path.replace('.bam', '.soft_clipped_reads.bam')

def create_bam_file_soft_reads(input_bam_path, soft_clipped_bam_path):
    print('... creating bam file with soft clipped reads. ~ {} min\n'.format((time.time() - start)//60))
    input_bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    soft_clipped_bam = pysam.AlignmentFile(soft_clipped_bam_path, "wb", template=input_bam_file)
    flags = [99, 147, 83, 163]
    for read in input_bam_file:
        if read.flag in flags and 'S' in read.cigarstring:
            soft_clipped_bam.write(read)
    soft_clipped_bam.close()
    input_bam_file.close()
    print('... bam file created. ~ {} min\n'.format((time.time() - start)//60))

def reverse_complement_cigar(cigar):
    ## TODO!!!!!
    return cigar

def get_numb_soft_nucleotides_from_cigartuples(cigartuples, left_right):
    if left_right == 'left':
        return cigartuples[0][1]
    elif left_right == 'right':
        return cigartuples[-1][1]

def get_strand_from_flag(flag):
    flag = int(flag)
    if flag == 99 or flag == 147:
        return '+'
    elif flag == 83 or flag == 163:
        return '-'
    else:
        raise NameError ('Cannot get strand from {}.'.format(flag))

def get_left_soft_clipping(read):
    numb_soft_nucleotides = get_numb_soft_nucleotides_from_cigartuples(read.cigartuples, 'left')
    read_seq = Seq(read.query_sequence)
    soft_nucleotide_seq = read_seq[:numb_soft_nucleotides]
    soft_start = read.reference_start - numb_soft_nucleotides
    soft_end = read.reference_start

    flag = read.flag
    strand = get_strand_from_flag(flag)
    if strand == '-': ## 83 or 163
        read_seq = read_seq.reverse_complement()
        soft_nucleotide_seq = soft_nucleotide_seq.reverse_complement()
        prime_end = '3prime'
        read_cigar = reverse_complement_cigar(read.cigarstring) ## TODO!!
    else:
        prime_end = '5prime'
        read_cigar = read.cigarstring
    return read_cigar, read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand

def get_right_soft_clipping(read):
    numb_soft_nucleotides = get_numb_soft_nucleotides_from_cigartuples(read.cigartuples, 'right')
    read_seq = Seq(read.query_sequence)
    soft_nucleotide_seq = read_seq[-numb_soft_nucleotides:]
    soft_start = read.reference_end
    soft_end = read.reference_end + numb_soft_nucleotides

    flag = read.flag
    strand = get_strand_from_flag(flag)
    if strand == '-':
        read_seq = read_seq.reverse_complement()
        soft_nucleotide_seq = soft_nucleotide_seq.reverse_complement()
        prime_end = '5prime'
        read_cigar = reverse_complement_cigar(read.cigarstring) ## TODO!!
    else:
        prime_end = '3prime'
        read_cigar = read.cigarstring
    return read_cigar, read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand

def create_csv_file_soft_reads(soft_clipped_bam_path, out_file_path):
    """For each read, get the info specified in the header
    """

    header = ['read_name', 'chr', 'mapping_start', 'read_size', 
              'cigar', 'read_seq', 'soft_start', 'soft_end',
              'numb_soft_nucleotides', 'soft_nucleotide_seq', 'one_or_both_ends', 'prime_end', 'flag', 'strand']
    
    get_info = 'read_cigar, read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand \
                = get_{}_soft_clipping(read)'
    
    add_info = 'soft_clipped_info.add((read.query_name, read.reference_name, read.reference_start, len(read.query_sequence), \
                                            read_cigar, read_seq, soft_start, soft_end, \
                                            numb_soft_nucleotides, soft_nucleotide_seq, one_or_both_ends, prime_end, flag, strand))'
    
    soft_clipped_bam_content = pysam.AlignmentFile(soft_clipped_bam_path, "rb")
    
    soft_clipped_info = set()
    print('... going through bam file reads. ~ {} min\n'.format((time.time() - start)//60))
    for read in soft_clipped_bam_content:
        ## ONLY IN ONE END:
        if read.cigarstring.count('S') == 1:
            one_or_both_ends = 1
            if read.cigarstring[-1] == 'S': ## Soft clippling on the right side
                exec(get_info.format('right'))
                exec(add_info)
            else: ## Soft clipping on the left side
                exec(get_info.format('left'))
                exec(add_info)
        ## BOTH ENDS:
        elif read.cigarstring.count('S') == 2:
            one_or_both_ends = 2
            exec(get_info.format('left'))
            exec(add_info)
            exec(get_info.format('right'))
            exec(add_info)
        
        else:
            print('Probably something wrong with: {}', read)
    
    print('... creating data frame and saving in csv file. ~ {} min\n'.format((time.time() - start)//60))
    soft_clipped_info_df = pd.DataFrame(soft_clipped_info, columns = header)
    soft_clipped_info_df.to_csv(out_file_path, index=False, compression={'method':'gzip'})
    return soft_clipped_info_df

def create_soft_clipped_bed_file(soft_clipped_info_df, out_file_path):
    """Creates a bed file with soft clipped regions coordinates.
       It doesn't matter from which read the soft region comes from, only the coordinates
       are necessary to get the sequence from reference genome (to see what was supposed to be there)
    """
    print('... creating bed file. ~ {} min\n'.format((time.time() - start)//60))
    soft_clipped_info_df['score'] = 0
    soft_clipped_info_df['bed_entry'] = 'soft_clipped_regions'
    
    bed_file_df = soft_clipped_info_df[['chr', 'soft_start', 'soft_end', 'bed_entry', 'score', 'strand']].drop_duplicates()
    bed_file_df = bed_file_df[bed_file_df['soft_start'] >= 0]
    bed_file_df.to_csv(out_file_path, sep='\t', header=False, index=False)


def get_soft_clipped_info(input_bam_path, out_file_path, bam_file_only):
    """Calls functions to create bam file, csv file and bed file.
       If bam_file_only is True, the last two functions will not be called.
    """

    ## First bam file:
    soft_clipped_bam_path = get_soft_clipped_bam_path(input_bam_path)
    # if not os.path.isfile(soft_clipped_bam_path): ############################## APAGAR MAIS TARDE
    #     create_bam_file_soft_reads(input_bam_path, soft_clipped_bam_path)
    create_bam_file_soft_reads(input_bam_path, soft_clipped_bam_path) ### Ficar com este
    
    if not bam_file_only:
        ## Second csv file:
        # csv_file_path = out_file_path.format('_info.csv')
        csv_file_path = out_file_path.format('_info.csv.gz')
        
        # if not os.path.isfile(csv_file_path): ############################## APAGAR MAIS TARDE
        #     soft_clipped_info = create_csv_file_soft_reads(soft_clipped_bam_path, csv_file_path)
        # else:
        #     soft_clipped_info = pd.read_csv(csv_file_path)
        soft_clipped_info = create_csv_file_soft_reads(soft_clipped_bam_path, csv_file_path) ### Ficar com este
        
        ## Third bed file:
        bed_file_path = out_file_path.format('.bed')
        create_soft_clipped_bed_file(soft_clipped_info, bed_file_path)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a bam file.
        Creates a bam file with primary and aligned reads from bam file (on the same folder as the bam file).
        Creates a csv file with soft clipped info from bam file reads (on the output folder).
        Creates a bed file with soft clipped info from csv file (on the output folder) that will be used to complete the csv file with genome reference sequence.
        If bamfileonly is set, csv and bed file won't be created, only a bam file with soft reads.""")

parser.add_argument('-o', '--outfilename',
                    nargs='?',
                    default='%filepath%',
                    help="""Output file prefix (no extension). A suffix "_soft_clipped_{}" will be added.
                    Default: basename of input file.""",
                    metavar='outfilename')

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('-b','--bamfileonly',
                    nargs='?',
                    default=False,
                    help="Boolean. False: creates all files (bam, csv, bed). True: creates only bam file with soft reads. Default: False.",
                    metavar='bamfileonly')

parser.add_argument('bamfilepath',
                    help="""Path to input bam file.""",
                    metavar='<input bam filepath>')

args = parser.parse_args()

input_bam_path = args.bamfilepath
out_file_name = args.outfilename
out_dir = args.outdir
bam_file_only = str_to_bool(args.bamfileonly)
out_file_path = get_out_file_path(out_file_name, input_bam_path, out_dir)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to get soft clipped info for {}\n'.format(input_bam_path))

get_soft_clipped_info(input_bam_path, out_file_path, bam_file_only)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED soft clipped info for: {} \nWritten in: {}. \n{} minutes and {} seconds.\n'.format(input_bam_path, out_file_path, minutes, seconds))