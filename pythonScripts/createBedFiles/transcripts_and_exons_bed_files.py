#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
From a file with a list of transcripts creates the following bed files:
1. transcripts bed file (forward and reverse)
2. exons bed file (forward and reverse)
3. exons last coordinate bed file (forward and reverse)
4. exons last coordinate window bed file (forward and reverse)

GTF file saved in csv format is also necessary
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import pandas as pd
# from gtfparse import read_gtf ## Python image in lobo does not have this module

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

#####################
#     Auxiliary     #
#####################

def read_transcripts_file(path):
    with open(path, 'r') as f:
        if 'ENST' in f.read()[0:5]: ## Check if the file has no header
            df = pd.read_csv(path, header=None, names=['transcript_id'])
        else:
            df = pd.read_csv(path, header=0, names=['transcript_id']) ## If the file had a header, it will be replaced
        df = check_ids(df, 'transcript_id')
        return df

def read_GTF_file(path):
    if '.csv' in path:
        df = pd.read_csv(path)
    else:
        raise NameError ('Cannot read gtf file. The gtf file has to be converted to csv previously.')
        # df = read_gtf(path) ## Cannot do like this because lobo server does not have this module
    ## Adapt dataframe before returning it:
    df['score'] = df['score'].fillna(0)
    df['exon_number'] = pd.to_numeric(df['exon_number'])
    df['start'] -= 1
    df = check_ids(df, 'transcript_id')
    df = check_ids(df, 'exon_id')
    return df

def check_ids(df, feature):
    """
    <feature> could be 'transcript_id' or 'exon_id', for example.
    Check if there is a column with the ids with no version.
    """
    if df[feature].str.contains('.', regex=False).any(): ## If the column transcript_id has the id with version
        df_temp = df.copy()
        print('Adding column "{}" with no version'.format(feature))
        df_temp['{}_with_version'.format(feature)] = df_temp[feature]
        df_temp[feature] = df_temp.apply(lambda row: str(getattr(row, feature)).split('.')[0], axis=1) # Then add col with no version
        return df_temp
    return df

def merge_transcripts_gtf(transcripts, gtf):
    return pd.merge(transcripts, gtf, on='transcript_id')

def get_merged_df(transcripts_file_path, gtf_file_path):
    transcript_df = read_transcripts_file(transcripts_file_path)
    gtf_df = read_GTF_file(gtf_file_path)
    return merge_transcripts_gtf(transcript_df, gtf_df)


def write_bed_file(df, out_file_path, feature, separated=True, together=False):
    if separated:
        df_forward = df.loc[df['strand'] == '+']
        df_reverse = df.loc[df['strand'] == '-']
        df_forward.to_csv(out_file_path.format(feature, '_F'), sep='\t', header=False, index=False)
        df_reverse.to_csv(out_file_path.format(feature, '_R'), sep='\t', header=False, index=False)
    if together:
        df.to_csv(out_file_path.format(feature, ''), sep='\t', header=False, index=False)


#####################
#     Bed Files     #
#####################

def check_bed_duplicates(df):
    """
    Check if there are any duplicated entry (based on 4th column)
    """
    colname = df.columns[3]
    duplicates_in_bed_file = df[df.duplicated(subset=colname)]
    if duplicates_in_bed_file.shape[0] > 1:
        for i in duplicates_in_bed_file[colname]:
            print('{} has a duplicated entry in bed file'.format(i))

def create_transcripts_bed_file(df, out_file_path):
    ## Select rows relative to transcripts:
    df_temp = df.loc[df['feature'] == 'transcript']
    
    ## Select columns of interest [chrom, star, end, name, score, strand]
    df_temp2 = df_temp[['seqname', 'start', 'end', 'transcript_id', 'score', 'strand']].drop_duplicates(ignore_index = True)
    
    ## Check if 4th column has unique values (Just to make sure that we don't have more than one entry for the same transcript)
    check_bed_duplicates(df_temp2)
    
    ## Write bed file:
    write_bed_file(df_temp2, out_file_path, 'transcripts')


def get_exon_rank(row):
    if row.exons_per_transcript == 1:
        return 'intronless'
    elif row.exons_per_transcript > 1 and row.exon_number == 1:
        return 'first'
    elif row.exons_per_transcript > 1 and row.exon_number == row.exons_per_transcript:
        return 'last'
    else:
        return 'middle'

def get_new_exon_id(row):
    """
    transcriptID_exonID_exonRank_exon_number
    example: ENST00000620552_ENSE00003750832_first_1
    """
    return '_'.join([row.transcript_id, row.exon_id, get_exon_rank(row), str(int(row.exon_number))])

def adjust_exon_id(df):
    """
    Input dataframe should be transcripts and GTF merged filtered by exons
    """
    df_temp = df.copy()
    # Add column with exons per transcript that will be used to get the exon rank (first, middle, last, intronless)
    df_temp['exons_per_transcript'] = df_temp.copy().groupby('transcript_id')['exon_id'].transform('count')
    # Modify the column with exon_id (after moving the original exon_id to another column)
    df_temp['exon_id_original'] = df_temp['exon_id']
    df_temp['exon_id'] = df_temp.copy().apply(lambda row: get_new_exon_id(row), axis=1)
    
    return df_temp

def get_surrounding_nucleotides(line, window_size=200):
    """If windo_size==200, it will get the 100 upstream nucleotides and 100 downstream nucleotides
    of the input nucleotide (line).
    """
    surrounding_nucleotides = []
    for i in range(1, int(window_size/2+1)):
        up_line = line.copy()
        down_line = line.copy()
        up_line[1] = int(up_line[1]) - i
        up_line[2] = int(up_line[2]) - i
        down_line[1] = int(down_line[1]) + i
        down_line[2] = int(down_line[2]) + i
        surrounding_nucleotides.append(up_line)
        surrounding_nucleotides.append(down_line)
    surrounding_nucleotides = sorted(surrounding_nucleotides)
    return surrounding_nucleotides

def create_exons_last_coordinate_window(exons_last_coordinate):
    """This function receives and returns data frames, but works with list of lists.
    Gets the content of exons last nucleotide.
    For each entry creates 100 new entries upstream and 100 new entries downstream. (The orignial entry isn't saved)
    """
    c = exons_last_coordinate.columns
    exons_last_coordinate = exons_last_coordinate.values.tolist()

    exons_windows = []
    for line in exons_last_coordinate:
        [exons_windows.append(new_line) for new_line in get_surrounding_nucleotides(line)]
    
    return pd.DataFrame(exons_windows, columns=c)

def create_exons_bed_file(df, out_file_path):
    ###########################
    # Exons normal start end: #
    ###########################
    ## Select rows relative to exons:
    df_temp = df.loc[df['feature'] == 'exon']
    
    ## Adjust exon_id (to have transcriptID_exonID_exonRank_exonNumber)
    df_temp2 = adjust_exon_id(df_temp)
    
    ## Select columns of interest [chrom, star, end, name, score, strand]
    df_temp3 = df_temp2[['seqname', 'start', 'end', 'exon_id', 'score', 'strand']].drop_duplicates(ignore_index = True)
    
    ## Check if 4th column has unique values
    check_bed_duplicates(df_temp3)
    
    ## Write exons bed file:
    write_bed_file(df_temp3, out_file_path, 'exons')
    
    def new_start(row, n):
        if row.strand == '+':
            return row.end - n
        if row.strand == '-':
            return row.start
    
    def new_end(row, n):
        if row.strand == '+':
            return row.end
        if row.strand == '-':
            return row.start + n

    ###########################
    #  Exons last coordinate: #
    ###########################
    ## Corresponds to the exons bed file but only for first and middle exons and
    # Start = End - 1 (in Forward) and End = Start + 1 (in Reverse)
    df_temp4 = df_temp3.loc[df_temp3['exon_id'].str.contains('first|middle')].copy()
    df_temp4['start'] = df_temp4.apply(lambda row: new_start(row, 1), axis=1)
    df_temp4['end'] = df_temp4.apply(lambda row: new_end(row, 1), axis=1)

    write_bed_file(df_temp4, out_file_path, 'exons_last_coordinate')

    ###############################
    #  Exons last two coordinate: #
    ###############################
    ## Corresponds to the exons bed file but only for first and middle exons and
    # Start = End - 1 (in Forward) and End = Start + 1 (in Reverse)
    df_temp5 = df_temp3.loc[df_temp3['exon_id'].str.contains('first|middle')].copy()
    df_temp5['start'] = df_temp5.apply(lambda row: new_start(row, 2), axis=1)
    df_temp5['end'] = df_temp5.apply(lambda row: new_end(row, 2), axis=1)
    
    write_bed_file(df_temp5, out_file_path, 'exons_last_two_coordinate')

    ###########################
    #      Exons windows:     #
    ###########################
    ## Bed file with an entry for each of the 200 surrounding nucleotides of the exons last coordinates
    df_temp5 = create_exons_last_coordinate_window(df_temp4)
    write_bed_file(df_temp5, out_file_path, 'exons_last_coordinate_window')

    ###############################
    #  Exons last coordinate + 1: #
    ###############################
    ## Corresponds to the exons last coordinate bed file but plus 1
    
    def plus_x(row, start_or_end, x):
        if row.strand == '+':
            return start_or_end + x
        if row.strand == '-':
            return start_or_end - x

    df_temp6 = df_temp4.copy()
    df_temp6['start'] = df_temp6.apply(lambda row: plus_x(row, row.start, 1), axis=1)
    df_temp6['end'] = df_temp6.apply(lambda row: plus_x(row, row.end, 1), axis=1)

    write_bed_file(df_temp6, out_file_path, 'exons_last_coordinate_plus_one')



###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

## Get merged df:
# merged_df = get_merged_df(transcripts_file_path, gtf_file_path)


## Create transcripst bed file:
# create_transcripts_bed_file(merged_df, out_file_path)

## Create exons bed file (including last coordinates):
# create_exons_bed_file(merged_df, out_file_path)



