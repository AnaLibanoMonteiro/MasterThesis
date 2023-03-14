#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
Add suffix to sample names in
'01_sample_names.csv'
'02_metadata.csv'
files.
And change organism in '02_metadata.csv'

Use this script to create a specific folder to process Drosophila reads (spike-ins)
to get normalization factors.
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import argparse

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

def add_suffix(suffix, sample_names_file, metadata_file):
    ## First file:
    with open(sample_names_file, "r") as file_1:
        content_sample_names=file_1.read().strip('\n').split('\n')
    if suffix not in content_sample_names[0]: ## Only replace if the file does not have the suffix yet
        new_content_sample_names = content_sample_names.copy()
        for i in range(0, len(new_content_sample_names)):
            new_content_sample_names[i] += suffix
        new_content_sample_names = "\n".join(new_content_sample_names)
        ## Re-write:
        with open(sample_names_file, "w") as file_1_out:
            file_1_out.write(new_content_sample_names)

    ## Second file:
    with open(metadata_file, "r") as file_2:
        content_metadata=file_2.read().strip('\n') #.split('\n')
    if suffix not in content_metadata: ## Only replace if the file does not have the suffix yet
        new_content_metadata = content_metadata
        for i in content_sample_names:
            if i in content_metadata:
                new_content_metadata = new_content_metadata.replace(i, i + suffix)
        new_content_metadata = new_content_metadata.replace('organism,hg', 'organism,'+suffix.strip('-'))
        ## Re-write:
        with open(metadata_file, "w") as file_2_out:
            file_2_out.write(new_content_metadata)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Adds a suffix to sample names. By default adds '-dm6'.
        Changes files '01_sample_names.csv' and '02_metadata.csv'.""")

parser.add_argument('-f1','--sampleNamesFile',
                    nargs='?',
                    default="01_sample_names.csv",
                    metavar='sampleNamesFile')

parser.add_argument('-f2', '--metadataFile',
                    nargs='?',
                    default='02_metadata.csv',
                    metavar='metadataFile')

parser.add_argument('-s', '--suffix',
                    nargs='?',
                    default='-dm6',
                    help="""Suffix to be added to the sample names""",
                    metavar='suffix')

args = parser.parse_args()

sample_names_file = args.sampleNamesFile
metadata_file = args.metadataFile
suffix = args.suffix

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

add_suffix(suffix, sample_names_file, metadata_file)