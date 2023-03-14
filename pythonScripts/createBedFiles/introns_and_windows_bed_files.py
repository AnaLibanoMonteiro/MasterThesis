#!/usr/bin/env python

###########################################
#                                         #
#               Description               #
#                                         #
###########################################

"""
From an exon bed file creates the following bed files:
1. introns
2. introns windows
"""

###########################################
#                                         #
#            Import libraries             #
#                                         #
###########################################

import copy

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

####################################
#                                  #
#              Common              #
#                                  #
####################################

def read_exons_bed_file_except_intronless(exon_file_path):
    content = []
    with open(exon_file_path)as f:
        for line in f:
            if 'intronLess' not in line.strip().split()[3]:
                content.append(line.strip().split())
    return content


def write_bed_file(out_file_path, content):
    with open(out_file_path, 'w') as fileToWrite:
        for line in content:
            fileToWrite.write('\t'.join(line))
            fileToWrite.write('\n')


####################################
#                                  #
#             Introns              #
#                                  #
####################################


def create_introns_bed_file(exon_file_path, out_file_path, forward_reverse):
    """
    'forward_reverse' must be either 'F', 'R' or ''
    """
    exons_bed_file_content = read_exons_bed_file_except_intronless(exon_file_path)
    transcript = exons_bed_file_content[0][3].split('_')[0]
    sub_group = []
    new_content = []

    for i in range(0, len(exons_bed_file_content)):
        if exons_bed_file_content[i][3].split('_')[0] == transcript:
            sub_group.append(exons_bed_file_content[i])
        else:
            new_sub_group = get_introns_lines(sub_group)
            [new_content.append(line) for line in new_sub_group]
            transcript = exons_bed_file_content[i][3].split('_')[0]
            sub_group = [exons_bed_file_content[i]]
        # Correct the problem with last transcript
        if i == len(exons_bed_file_content) - 1:
            new_sub_group = get_introns_lines(sub_group)
            [new_content.append(line) for line in new_sub_group]
    
    if forward_reverse != '':
        forward_reverse = '_' + forward_reverse
    
    write_bed_file(out_file_path.format('introns', forward_reverse), new_content)


def get_coordinates(sub_group, strand):
    coordinates = []
    for i in sub_group:
        coordinates.append(i[1])
        coordinates.append(i[2])
    if strand == '+':
        return coordinates[1:-1]
    if strand == '-':
        return coordinates


def get_introns_lines(sub_group):
    introns_lines = []
    numb_lines = len(sub_group) - 1
    template = sub_group[0]
    strand = template[5]
    coordinates = get_coordinates(sub_group, strand)
    if strand == '+':
        c = 0
    if strand == '-':
        c = 3
    for i in range(0, numb_lines):
        new_line = template
        new_line[1] = coordinates[c]
        if strand == '+':
            new_line[2] = coordinates[c+1]
        if strand == '-':
            new_line[2] = coordinates[c-3]
        c += 2
        new_line[3] = template[3].split('_')[0] + '_' + 'intron_{}'.format(i+1)  # Change intron numb according exon numb
        introns_lines.append(copy.deepcopy(new_line))
    return introns_lines



####################################
#                                  #
#             Windows              #
#                                  #
####################################


def create_windows_bed_file(exon_file_path, out_file_path, forward_reverse, window_size):
    """
    'forward_reverse' must be either 'F', 'R' or ''
    """
    exons_bed_file_content = read_exons_bed_file_except_intronless(exon_file_path)
    transcript = exons_bed_file_content[0][3].split('_')[0]
    sub_group = []
    first_windows_content = []
    second_windows_content = []
    third_windows_content = []

    for i in range(0, len(exons_bed_file_content)):
        if exons_bed_file_content[i][3].split('_')[0] == transcript:
            sub_group.append(exons_bed_file_content[i])
        else:
            first, second, third = get_window_lines(sub_group, window_size)
            [first_windows_content.append(line) for line in first]
            [second_windows_content.append(line) for line in second]
            [third_windows_content.append(line) for line in third]

            transcript = exons_bed_file_content[i][3].split('_')[0]
            sub_group = [exons_bed_file_content[i]]
        # Correct the problem with last transcript
        if i == len(exons_bed_file_content) - 1:
            first, second, third = get_window_lines(sub_group, window_size)
            [first_windows_content.append(line) for line in first]
            [second_windows_content.append(line) for line in second]
            [third_windows_content.append(line) for line in third]
    
    if forward_reverse != '':
        forward_reverse = '_' + forward_reverse
    
    write_bed_file(out_file_path.format('introns_first_window', forward_reverse), first_windows_content)
    write_bed_file(out_file_path.format('introns_second_window', forward_reverse), second_windows_content)
    write_bed_file(out_file_path.format('introns_third_window', forward_reverse), third_windows_content)


def get_exon_number(line):
    return int(line[3].split('_')[-1])


def get_window_lines(sub_group, window_size):
    sub_group.sort(key=get_exon_number)  # Sort sub_group according to exons
    strand = sub_group[0][5]
    transcript = sub_group[0][3].split('_')[0]
    first_window_lines = []
    second_window_lines = []
    third_window_lines = []
    for line in sub_group:
        name_template = '{}_intron_{}_{}_window'  # .format(transcript, intron_number, window_name)
        exon_number = get_exon_number(line)

        first_window_lines.append(get_first_window_line(line, strand, window_size, transcript, exon_number, name_template))
        second_window_lines.append(get_second_window_line(line, strand, window_size, transcript, exon_number, name_template))
        third_window_lines.append(get_third_window_line(line, strand, window_size, transcript, exon_number, name_template))
    return first_window_lines[:-1], second_window_lines[1:], third_window_lines[1:]


def get_first_window_line(line, strand, window_size, transcript, exon_number, name_template, window_name='first'):
    new_line = copy.deepcopy(line)
    if strand == '+':
        new_line[1] = str(int(line[2]) - window_size)
    if strand == '-':
        new_line[2] = str(int(line[1]) + window_size)
    new_line[3] = name_template.format(transcript, exon_number, window_name)
    return new_line


def get_second_window_line(line, strand, window_size, transcript, exon_number, name_template, window_name='second'):
    new_line = copy.deepcopy(line)
    if strand == '+':
        new_line[2] = str(int(line[1]) + window_size)
    if strand == '-':
        new_line[1] = str(int(line[2]) - window_size)
    new_line[3] = name_template.format(transcript, exon_number - 1, window_name)
    return new_line


def get_third_window_line(line, strand, window_size, transcript, exon_number, name_template, window_name='third'):
    new_line = copy.deepcopy(line)
    if strand == '+':
        new_line[1] = str(int(line[1]) - window_size)
        new_line[2] = line[1]
    if strand == '-':
        new_line[1] = line[2]
        new_line[2] = str(int(line[2]) + window_size)
    new_line[3] = name_template.format(transcript, exon_number - 1, window_name)
    return new_line


###########################################
#                                         #
#                  Main                   #
#                                         #
###########################################

## Create introns bed file:
# create_introns_bed_file(exon_file_path, out_file_path, '')
# create_introns_bed_file(exon_file_path_F, out_file_path, 'F')
# create_introns_bed_file(exon_file_path_R, out_file_path, 'R')

## Create windows bed file:
# create_windows_bed_file(exon_file_path, out_file_path, '', 10)
# create_windows_bed_file(exon_file_path_F, out_file_path, 'F', 10)
# create_windows_bed_file(exon_file_path_R, out_file_path, 'R', 10)
