# get_SNR_bam_ignoreSoftClip

## Minor changes. Ready to work with samtools v1.10 or later

Python script just as [mNET_snr](https://github.com/tomasgomes/mNET_snr) modified to ignore reads with soft clipping.
Described as part of an analysis pipeline in *insert reference*
For usage, see get_SNR_bam_ignoreSoftClip.py -h

**Parameters**:

-f, One or more file paths, separated by spaces

-s, New file prefixes (no extension), separated by spaces. In the same number as the file paths.

-d, Output directory. Defaults to './'



### It operates in the following manner:

#### The input alignment file (.bam) provided together with the -f argument will be converted to .sam to ease subsequent parsing.

#### For each read in the .sam file only consider read 2 from the pair (147/163 SAM flag field) and disregard any read that contains deletions, insertions or soft clipping information in the CIGAR string.

#### Modification of SAM fields 2, 4, 6, 9, 10, 11 and strand information in order to represent the sequenced 3'OH.

##### 2
Since the original flag will be either 147 or 163 (second read in pair) this will be changed to 99 or 83 respectively. Thus, obtaining a SAM flag that points to a first read in pair with read 1 directionality.

##### 4
The 1-based position corresponds to the other end of the read for 147 flagged reads so we correct it. 83 flagged reads remain unchanged.

##### 6
The CIGAR string will be turned into "1M".

##### 9
Observed template length is divided by its absolute in order to obtain 1 and preserve its directionality.

##### 10
The nucleotide represented must be the one where the 1-based position that the 4th SAM field points to. For 147 flagged reads we extract the last nucleotide and for 83 flagged reads the first.

##### 11
The Phred-scale base quality+33 must be the one where the 1-based position that the 4th SAM field points to. For 147 flagged reads we extract the last quality value and for 83 flagged reads the first.

##### strand information
Reversed.
  
#### Convert the new single nucleotide resolution .sam file into a .bam file. Sort and index the output.

