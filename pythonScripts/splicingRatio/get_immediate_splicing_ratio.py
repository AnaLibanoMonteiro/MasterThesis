#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

## TODO!! FINISH DESCRIPTION!!!

"""
FORWARD > > > >
  
       Ex 1                 In 1                 Ex 2                 In 2           Ex3  
####################-------------------##########################------------##################  
               | 1 |              | 3 || 2 |                | 1 |       | 3 || 2 |  
  
REVERSE < < < <  

       Ex 3                 In 2                 Ex 2                 In 1           Ex1  
####################-------------------##########################------------##################  
               | 2 || 3 |              | 1 |                | 2 || 3 |       | 1 |

# Spliced read maps on window 1 and 2 but NOT 3
# Unspliced read has to map on window 1, 2 and 3 OR 2 and 3
"""


###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import csv
import os
import time

import pandas as pd

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_out_file_path(out_file_name, in_file_path, out_dir):
    suffix = '_splicing_ratio.tsv'
    if out_file_name == '%filepath%': ## If the user didn't specify a name for the out file
        out_file_name = os.path.basename(in_file_path).split('.')[0] ## Then get the name of the input file
    return os.path.join(out_dir, out_file_name + suffix)

def read_intersect_file(file_path, only_full_intersection = True, window_size = 10):
    """
    Example of a line from input file:
    ['chr1', '153990792', '153990802', 'ENST00000651669.1_intron_1_first_window', '0', '+', \
     'chr1', '153990763', '153990804', 'K00181:182:HFGY2BBXY:2:2214:25530:14537/2', '255', '-', '10']
    Last column refers to the number of nucleotides that overlap between the read and the window.
    If only the reads that completely overlaped are desired, than last column should be equal to window size.
    """
    with open(file_path, 'r') as my_file:
        file_content = []
        if only_full_intersection:
            [file_content.append(line) for line in csv.reader(my_file, delimiter = '\t') if int(line[-1]) == window_size]
        else:
            [file_content.append(line) for line in csv.reader(my_file, delimiter = '\t')]
    return file_content

def get_intron_name(line):
    """
    Example of entry name in bed file:
    ENST00000651669.1_intron_1_first_window OR
    ENST00000651669_intron_1_first_window
    """
    windows_names = ['_first_window', '_second_window', '_third_window']
    for i in windows_names:
        if i in line[3]:
            intron_name = line[3].replace(i, '') ## remove window name to get intron name
    return intron_name

def get_all_introns(file_content):
    introns = set()
    for line in file_content:
        introns.add(get_intron_name(line))
    return introns

def add_read_dic(line, intron, dic):
    dic[intron].append(line[9])
    # Alternative to keep number of bases overlaped:
    # dic[intron].append({ line[9]: int(line[-1]) })

def create_windows_dic(file_content, introns):
    """
    Creates one dictionary for each window, with introns as keys and list of reads mapping on that introns' window as values
    window1: {intron1: [read1, read2, read4], intron2: [read5, read7, read8]}
    window2: {intron1: [read1, read2], intron2:[read5, read8, read9]}
    window3: {intron1: [], intron2: [read5, read8]}
    """
    windows_and_reads_ids_1 = {k: [] for k in introns}
    windows_and_reads_ids_2 = {k: [] for k in introns}
    windows_and_reads_ids_3 = {k: [] for k in introns}
    for line in file_content:
        intron = get_intron_name(line)
        if 'first' in line[3]:
            add_read_dic(line, intron, windows_and_reads_ids_1)
        if 'second' in line[3]:
            add_read_dic(line, intron, windows_and_reads_ids_2)
        if 'third' in line [3]:
            add_read_dic(line, intron, windows_and_reads_ids_3)
    return windows_and_reads_ids_1, windows_and_reads_ids_2, windows_and_reads_ids_3


def get_splicing_ratio(file_path, out_file_path):
    """For each intron, get the info specified in the header
    """
    header = ['intron_id', 'splicing_ratio', 'numb_spliced_reads', 'numb_unspliced_reads', 'numb_total_reads']

    file_content = read_intersect_file(file_path)
    introns = get_all_introns(file_content)
    ## Get windows dictionaries:
    dic1, dic2, dic3 = create_windows_dic(file_content, introns)
    
    splicing_ratio_per_intron = [] ## Empty list were the results will be saved

    ## Go through all introns, and for each intron go through all reads that map
    ## at least in one of the windows of that intron
    for intron in introns:
        all_intron_reads = set(dic1[intron] + dic2[intron] + dic3[intron])
        spliced = []
        unspliced_123 = []
        unspliced_23 = []
        other = []
        ## classify the read and put it in the right list:
        for r in all_intron_reads:
            if (r in dic1[intron]) and (r in dic2[intron]) and (r in dic3[intron]):
                unspliced_123.append(r)
            elif (r not in dic1[intron]) and (r in dic2[intron]) and (r in dic3[intron]):
                unspliced_23.append(r)
            elif (r in dic1[intron]) and (r in dic2[intron]) and (r not in dic3[intron]):
                spliced.append(r)
            else:
                other.append(r)
        
        ## Calculate splicing ratio and save it in a list:
        total_spliced_unspliced = len(spliced) + len(unspliced_123) + len(unspliced_23)
        if total_spliced_unspliced == 0:
            sr = float('nan') ## Reads mapping a specific intron may not give the info necessary to calculate splicing ratio
        else:
            sr = len(spliced)/total_spliced_unspliced
        splicing_ratio_per_intron.append([intron, sr, len(spliced), len(unspliced_123 + unspliced_23), \
                                          total_spliced_unspliced])
    
    ## Convert list to dataframe and save it in output file
    splicing_ratio_per_intron_df = pd.DataFrame(splicing_ratio_per_intron, columns = header)
    splicing_ratio_per_intron_df.to_csv(out_file_path, sep='\t', index=False)



###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input a file from an intersection
        between a bed file with introns windows and a bam file with reads.
        Returns a csv file with the splicing ratio for each event (intron) considered in the bed file.""")

parser.add_argument('-o', '--outfilename',
                    nargs='?',
                    default='%filepath%',
                    help="""Output file prefix (no extension). A suffix "_splicing_ratio.csv" will be added.
                    Default: basename of input file.""",
                    metavar='outfilename')

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('filepath',
                    help="""Path to input file. This file must be a result from an intersection between a
                    bed file with introns windows and a bam file.""",
                    metavar='<input filepath>')

args = parser.parse_args()

file_path = args.filepath
out_file_name = args.outfilename
out_dir = args.outdir
# out_dir = make_sure_path_exists(args.outdir)

out_file_path = get_out_file_path(out_file_name, file_path, out_dir)


###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

get_splicing_ratio(file_path, out_file_path)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED splicing ratio for: {} \nWritten in: {}. \n{} minutes and {} seconds.\n'.format(file_path, out_file_path, minutes, seconds))
