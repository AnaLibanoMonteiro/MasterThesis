#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
SELECT EVENTS (EXONS LAST COORDINATE or INTRONS) ACCORDING TO TIME POINT (TRANSCRIPTION FRONT WAVE)

According to the time point, there will be different exons being transcribed.
After DRB treatment, a transcription front wave can be observed;
However, some exons downstream the front wave can be detected (due to contamination, for example)
and they should be discarded.

Possible time points: CTR, DRB, W5, W10, W15, W30.

Attention: transription wave file is only referent to replicate 1.
(In original file that Rui created, there is also info about replicate 2)
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

def checkTimePointToFilterBy(time_point_to_filter_by):
    possible = ['CTR', 'DRB', 'W5', 'W10', 'W15', 'W30']
    if time_point_to_filter_by in possible:
        return time_point
    else:
        raise NameError('Time point to filter by is not present in possible choices.')

def check_file_type(file_type):
    if file_type == '1' or file_type == 1:
        return 'intermediates'
    elif file_type == '2' or file_type == 2:
        return 'immediate'
    else:
        raise NameError('File type not inclued in available options. \
        1 for splicing intermediates, 2 for immediate splicing ratio')

def get_events_info_file_path(events_info_file_path, file_type):
    if events_info_file_path == '%eventsinfo%':
        if file_type == 'intermediates':
            return '~/general/filesFromRui/exons_info.csv'
        elif file_type == 'immediate':
            return '~/general/filesFromRui/introns_info.csv'
    else:
        return events_info_file_path

def str_to_bool(string):
    return string.lower() in ("true", "t", "yes", "1", "replace")

### splicing ratio function:
def get_transcript_id(row):
    """Get transcript_id without version from intron name from
       splicing ratio dataframe (created with script get_splicing_ratio)."""
    if '.' in row:
        return row.intron_id.split('.')[0]
    else:
        return row.intron_id.split('_')[0]

def read_input_file(input_file_path, file_type, time_point, time_point_to_filter_by, sample_name):
    if file_type == 'intermediates':
        columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'coverage_counts'] 
        input_content = pd.read_csv(input_file_path, sep='\t', header=None, names=columns)
        ## Split ids in 'name' column in individual columns:
        input_content[['transcript_id', 'exon_id', 'exon_position', 'exon_rank']] = input_content['name'].str.split('_', expand=True)
    elif file_type == 'immediate':
        input_content = pd.read_csv(input_file_path, sep='\t')
        ## Add columns with transcript_id (without version) and intron_rank to splicing ratio dataframe:
        input_content['transcript_id'] = input_content.apply(lambda row: get_transcript_id(row), axis=1)
        input_content['intron_rank'] = input_content.apply(lambda row: int(row.intron_id.split('_')[-1]), axis=1)
    input_content['time_point'] = time_point
    input_content['time_point_to_filter_by'] = time_point_to_filter_by
    input_content['sample_name'] = sample_name
    return input_content

def get_merge_on(file_type):
    if file_type == 'intermediates':
        return ['transcript_id', 'exon_id']
    if file_type == 'immediate':
        return ['transcript_id', 'intron_rank']

def get_drop_subset(file_type):
    if file_type == 'intermediates':
        return ['start', 'end', 'name']
    if file_type == 'immediate':
        return 'intron_id'

def transcription_wave_filter(input_file_path, file_type, events_info_file_path,
            transcription_wave_file_path, time_point, time_point_to_filter_by, sample_name, replace):
    ## Read files:
    input_content = read_input_file(input_file_path, file_type, time_point, time_point_to_filter_by, sample_name)
    df_events_info = pd.read_csv(events_info_file_path)
    df_transcription_wave = pd.read_csv(transcription_wave_file_path)
    df_transcription_wave.rename(columns={'time_point':'time_point_to_filter_by'}, inplace=True)

    ## Merge coverage_content with exons_info (to get distance to TSS):
    df_temp1 = pd.merge(input_content, df_events_info, on=get_merge_on(file_type))

    ## Merge coverage_content with front wave (to get info about polymerase II position):
    df_temp2 = pd.merge(df_temp1, df_transcription_wave, on=['transcript_id', 'time_point_to_filter_by'])

    ## Select events:
    df_temp3 = df_temp2.loc[df_temp2['distance_from_TSS_to_3prime'] < df_temp2['pol_position']]

    ## Select columns of interest: (will be the same as the input file + id columns and time_point) 
    df_temp4 = df_temp3[input_content.columns]

    ## Make sure there are no duplicates:
    df_final = df_temp4.drop_duplicates(subset=get_drop_subset(file_type), ignore_index=True)
    
    ## Save result:
    file_extension = os.path.splitext(input_file_path)[-1]
    new_file_extension = '_filtered' + file_extension
    df_final.to_csv(input_file_path.replace(file_extension, new_file_extension), index=False, sep='\t')
    if replace:
        os.remove(input_file_path)


###########################################
#                                         #
#             Parse Arguments             #
#                                         #
###########################################

parser = argparse.ArgumentParser(description="""
        Takes as input a file with splicing intermediates coverage created with 'bedtools coverage'
        or a file with splicing ratio per intron created with get_splicing_ratio.py.
        Returns a similiar file but only with the events expressed in that specific time point.\n
        Two auxiliary files are necessary: one with events info (to get the distance from TSS).
        Another with transcription front wave info, to discard events according to time point.
        If the time point is not specified the program will try to deduce from the file name.
        Possible time points: CTR, DRB, W5, W10, W15, W30.""")

parser.add_argument('-r', '--replace',
                    nargs='?',
                    default='True',
                    help="""Boolean. True: remove original file. False: keep both files (original and _filtered.tsv).
                    Default: True.""",
                    metavar='replace')

parser.add_argument('eventsInfoFilepath',
                    nargs='?',
                    default='%eventsinfo%',
                    help="""Path to file with events info.
                    Default: '~/general/filesFromRui/exons_info.csv' (for splicing intermediates)
                    or '~/general/filesFromRui/introns_info.csv' (for splicing ratio).""",
                    metavar='<eventsInfo filepath>')

parser.add_argument('transcriptionWaveFilepath',
                    nargs='?',
                    default='~/general/filesFromRui/transcription_wave.csv',
                    help="""Path to file with transcription front wave info.
                    Default: '~/general/filesFromRui/transcription_wave.csv'.""",
                    metavar='<transcriptionWave filepath>')

parser.add_argument('sampleName',
                    help="""Sample Name.""",
                    metavar='<sample name>')

parser.add_argument('timePoint',
                    help="""Sample Time Point.
                    It may be different from the time point to filter by""",
                    metavar='<time point>')

parser.add_argument('timePointToFilterBy',
                    help="""Time point to filter by.
                    It may be different from the sample time point.
                    Possible time points are: CTR, DRB, W5, W10, W15, W30.""",
                    metavar='<time point to filter by>')

parser.add_argument('fileType',
                    help="""Input file type. 2 options available:
                    1 for splicing intermediates;
                    2 for immediate splicing ratio.""",
                    metavar='<file type>')

parser.add_argument('inputFilepath',
                    help="""Path to input file. File with splicing intermediates coverage or with splicing ratio.""",
                    metavar='<input filepath>')

args = parser.parse_args()

replace = str_to_bool(args.replace)
file_type = check_file_type(args.fileType)
events_info_file_path = get_events_info_file_path(args.eventsInfoFilepath, file_type)
transcription_wave_file_path = args.transcriptionWaveFilepath
sample_name = args.sampleName
input_file_path = args.inputFilepath
# time_point = get_time_point(args.timePoint, input_file_path)
time_point = args.timePoint
time_point_to_filter_by = checkTimePointToFilterBy(args.timePointToFilterBy)

###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

start = time.time()

transcription_wave_filter(input_file_path, file_type, events_info_file_path,
        transcription_wave_file_path, time_point, time_point_to_filter_by, sample_name, replace)

end = time.time()
temp = end-start
minutes = temp//60
seconds = temp - 60*minutes

print('FINISHED filtering events for: {} \n{} minutes and {} seconds.\n'.format(input_file_path, minutes, seconds))