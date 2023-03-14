#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Takes as input BAM file.
Creates another BAM file with reads that have soft clipping at 3 prime end.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse
import os
import time

import pysam  # # ATTENTION: pysam is always 0-based, as BAM files and BED files

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def get_out_file_path(out_dir, out_file_name, input_bam_path):
    if out_file_name == "%temp%":
        out_file_name = os.path.basename(input_bam_path)
    return os.path.join(out_dir, out_file_name)

def create_3prime_soft_clipped_bam(input_bam_path, output_bam_path):
    print('... creating bam file with soft clipped reads. ~ {} min\n'.format((time.time() - start)//60))
    input_bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    soft_clipped_bam = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam_file)

    for read in input_bam_file:
        if read.cigartuples[-1][0] == 4 and read.flag == 147:
            soft_clipped_bam.write(read)
        elif read.cigartuples[0][0] == 4 and read.flag == 163:
            soft_clipped_bam.write(read)

    soft_clipped_bam.close()
    input_bam_file.close()
    print('... bam file created. ~ {} min\n'.format((time.time() - start)//60))

###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""Takes as input BAM file.
    Creates another BAM file with reads that have soft clipping at 3 prime end.""")

parser.add_argument('-d','--outdir',
                    nargs='?',
                    default="./",
                    help="Output directory. Default: current directory.",
                    metavar='outdir')

parser.add_argument('-o','--outfilename',
                    nargs='?',
                    default="%temp%",
                    help="Output file name. Default: same name as input file.",
                    metavar='outfilename')

parser.add_argument('bamfilepath',
                    help="""Path to input bam file.""",
                    metavar='<input bam filepath>')

args = parser.parse_args()

input_bam_path = args.bamfilepath
out_dir = args.outdir
out_file_name = args.outfilename

output_bam_path = get_out_file_path(out_dir, out_file_name, input_bam_path)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()
print('STARTING to get bam file with 3 prime soft clipped reads for {}\n'.format(input_bam_path))

create_3prime_soft_clipped_bam(input_bam_path, output_bam_path)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED 3 prime soft clipped bam for: {}.\n{} minutes and {} seconds.\n'.format(input_bam_path, minutes, seconds))
