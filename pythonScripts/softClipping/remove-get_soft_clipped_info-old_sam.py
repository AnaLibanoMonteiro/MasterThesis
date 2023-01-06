#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Get soft clipped info from mapped and primary reads bam file.
Creates SAM file, CSV file and BED file.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import os
import argparse
import pysam
import pandas as pd
import time
from Bio.Seq import Seq

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def make_sure_path_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def get_out_file_path(out_file_name, in_file_path, out_dir):
    suffix = '_soft_clipped_reads{}'
    if out_file_name == '%filepath%': ## If the user didn't specify a name for the out file
        out_file_name = os.path.basename(in_file_path).split('.')[0] ## Then get the name of the input file
    return os.path.join(out_dir, out_file_name + suffix)



def get_sam_file_path(input_bam_path):
    return input_bam_path.replace('.bam', '.soft_clipped_reads.sam')

def create_sam_file_soft_reads(input_bam_path, sam_file_path):
    print('... creating sam file with soft clipped reads. ~ {} min\n'.format((time.time() - start)//60))
    input_bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    soft_clipped_reads = pysam.AlignmentFile(sam_file_path, "wh", template=input_bam_file)
    for read in input_bam_file:
        if not read.is_secondary and not read.is_unmapped and not read.mate_is_unmapped and 'S' in read.cigarstring:
            soft_clipped_reads.write(read)
    soft_clipped_reads.close()
    input_bam_file.close()
    print('... sam file created. ~ {} min\n'.format((time.time() - start)//60))

def read_sam_file(sam_file_path):
    file_content = []
    with open(sam_file_path, 'r') as f:
        for line in f:
            if not line.startswith("@"):
                file_content.append(line.strip("\n").split("\t"))
    return file_content


def get_numb_soft_nucleotides_from_cigar(cigar, prime_end):
    if prime_end == '5prime':
        return int(cigar.split('S')[0])
    elif prime_end == '3prime':        
        i = len(cigar)-1
        n = ''
        while cigar[i-1].isdigit():
            i -= 1
            n += cigar[i]
        return int(n[::-1])

def get_strand_from_flag(flag):
    flag = int(flag)
    if flag == 99 or flag == 147:
        return '+'
    elif flag == 83 or flag == 163:
        return '-'
    else:
        return str(flag)

def get_soft_clipped_5prime(read):
    prime_end = '5prime'
    numb_soft_nucleotides = get_numb_soft_nucleotides_from_cigar(read[5], prime_end)
    read_seq = Seq(read[9])
    soft_nucleotide_seq = read_seq[:numb_soft_nucleotides]
    soft_start = int(read[3]) - numb_soft_nucleotides - 1 # 1st position aligned - numb_soft_clipped - 1
    soft_end = int(read[3]) - 1 # 1st position aligned - 1
    flag = int(read[1])
    strand = get_strand_from_flag(flag)
    if strand == '-': ## flag 83 or 163
        read_seq = read_seq.reverse_complement()
        soft_nucleotide_seq = soft_nucleotide_seq.reverse_complement()
        prime_end = '3prime' ############## SHOULD I DO THIS ?????
    return read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand

def get_soft_clipped_3prime(read, soft_start_5prime=None): ## soft_start_5prime is None when there is only soft clipping at 3' end
    prime_end = '3prime'
    numb_soft_nucleotides = get_numb_soft_nucleotides_from_cigar(read[5], prime_end)
    read_seq = Seq(read[9])
    soft_nucleotide_seq = read_seq[-numb_soft_nucleotides:]
    if soft_start_5prime:
        soft_end = soft_start_5prime + len(read_seq) ## 1st position counting with soft clipped 5prime + read_size
        soft_start = soft_end - numb_soft_nucleotides
    else:
        soft_end = int(read[3]) + len(read_seq) - 1 # 1st position aligned + read_size - 1
        soft_start = soft_end - numb_soft_nucleotides
    flag = int(read[1])
    strand = get_strand_from_flag(flag)
    if strand == '-': ## flag 83 or 163
        read_seq = read_seq.reverse_complement()
        soft_nucleotide_seq = soft_nucleotide_seq.reverse_complement()
        prime_end = '5prime' ############## SHOULD I DO THIS ?????
    return read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand

def create_csv_file_soft_reads(sam_file_path, out_file_path):
    """
    For each read, get the info specified in the header
    """

    header = ['read_name', 'chr', 'mapping_start', 'read_size', 
              'cigar', 'read_seq', 'soft_start', 'soft_end',
              'numb_soft_nucleotides', 'soft_nucleotide_seq', 'one_or_both_ends', 'prime_end', 'flag', 'strand']
    
    get_info = 'read_seq, soft_start, soft_end, numb_soft_nucleotides, soft_nucleotide_seq, prime_end, flag, strand \
                = get_soft_clipped_{}(read{})'

    add_info = 'soft_clipped_info.add((read[0], read[2], int(read[3]) - 1, len(read[9]), \
                                            read[5], read_seq, soft_start, soft_end, \
                                            numb_soft_nucleotides, soft_nucleotide_seq, one_or_both_ends, prime_end, flag, strand))'
    
    sam_content = read_sam_file(sam_file_path)
    
    soft_clipped_info = set()
    print('... going through sam file reads. ~ {} min\n'.format((time.time() - start)//60))
    for read in sam_content:
        ## ONLY IN ONE END:
        if read[5].count('S') == 1:
            one_or_both_ends = 1
            if read[5][-1] == 'S': 
                ## 3' end
                exec(get_info.format('3prime', ''))
                exec(add_info)
            else:
                ## 5' end
                exec(get_info.format('5prime', ''))
                exec(add_info)
        ## BOTH ENDS:
        elif read[5].count('S') == 2:
            one_or_both_ends = 2
            ## 5' end
            exec(get_info.format('5prime', ''))
            exec(add_info)
            ## 3' end
            exec(get_info.format('3prime', ', soft_start'))
            exec(add_info)
        
        else:
            print('Probably something wrong with: {}', read)
    
    print('... creating data frame and saving in csv file. ~ {} min\n'.format((time.time() - start)//60))
    soft_clipped_info_df = pd.DataFrame(soft_clipped_info, columns = header)
    soft_clipped_info_df.to_csv(out_file_path, index=False)
    return soft_clipped_info_df

def create_soft_clipped_bed_file(soft_clipped_info_df, out_file_path):
    print('... creating bed files. ~ {} min\n'.format((time.time() - start)//60))
    soft_clipped_info_df['score'] = 0
    soft_clipped_info_df['bed_entry'] = 'soft_clipped_regions'
    
    bed_file_df = soft_clipped_info_df[['chr', 'soft_start', 'soft_end', 'bed_entry', 'score', 'strand']].drop_duplicates()
    bed_file_df = bed_file_df[bed_file_df['soft_start'] >= 0]
    bed_file_df.to_csv(out_file_path, sep='\t', header=False, index=False)


def get_soft_clipped_info(input_bam_path, out_file_path, sam_file_only):
    
    ## First sam file:
    sam_file_path = get_sam_file_path(input_bam_path)
    if not os.path.isfile(sam_file_path): ############################## APAGAR MAIS TARDE
        create_sam_file_soft_reads(input_bam_path, sam_file_path)
    # create_sam_file_soft_reads(input_bam_path, sam_file_path) ### Ficar com este
    
    if not sam_file_only:
        ## Second csv file:
        csv_file_path = out_file_path.format('_info.csv')
        
        # if not os.path.isfile(csv_file_path): ############################## APAGAR MAIS TARDE
        #     soft_clipped_info = create_csv_file_soft_reads(sam_file_path, csv_file_path)
        # else:
        #     soft_clipped_info = pd.read_csv(csv_file_path)
        soft_clipped_info = create_csv_file_soft_reads(sam_file_path, csv_file_path) ### Ficar com este
        
        ## Third bed file:
        bed_file_path = out_file_path.format('.bed')
        create_soft_clipped_bed_file(soft_clipped_info, bed_file_path)
    
    

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a bam file.
        Creates a sam file with primary and aligned reads from bam file (on the same folder as the bam file).
        Creates a csv file with soft clipped info from sam file reads (on the output folder).
        Creates a bed file with soft clipped info from csv file (on the output folder) that will be used to complete the csv file with genome reference sequence.
        If samfileonly is set, csv and bed file won't be created, only a sam file with soft reads.""")

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

parser.add_argument('-s','--samfileonly',
                    nargs='?',
                    default=False,
                    help="Boolean. False: creates all files (sam, csv, bed). True: creates only sam file with soft reads. Default: False.",
                    metavar='samfileonly')

parser.add_argument('bamfilepath',
                    help="""Path to input bam file.""",
                    metavar='<input bam filepath>')

args = parser.parse_args()

input_bam_path = args.bamfilepath
out_file_name = args.outfilename
out_dir = make_sure_path_exists(args.outdir)
sam_file_only = args.samfileonly

out_file_path = get_out_file_path(out_file_name, input_bam_path, out_dir)


###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to get soft clipped info for {}\n'.format(input_bam_path))

get_soft_clipped_info(input_bam_path, out_file_path, sam_file_only)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED soft clipped info for: {} \nWritten in: {}. \n{} minutes and {} seconds.\n'.format(input_bam_path, out_file_path, minutes, seconds))