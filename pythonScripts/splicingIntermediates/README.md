# **Python scripts for splicing intermediates analysis**

Scripts developed to identify exons with splicing intermediates using nascent RNA data.

## **get_peaks_splicing_intermediates.py**

Takes as input 2 files:  
(1) splicing intermediates coverage and  
(2) surrounding coverage (coverage of 100 nucleotides before and 100 after last exon coordinate);  
and the sample group name.
Returns the events on file (1) that have coverage significantly higher than surrounding nucleotides in file (2).


usage:  
`get_peaks_splicing_intermediates.py [-h] <splicingIntermediatesCoverage filepath> <windowsCoverage filepath> <sample group>`


### **Input files are created with `bedtools coverage`:**

**File (1)** is created with the following commands:

`bedtools coverage -a <bed_file_with_exons_last_coordinate_forward.bed> \`  
                  `-b <bam_file_forward.bam> \`  
                  `-counts \`  
                  `> <results_F.tsv>`

`bedtools coverage -a <bed_file_with_exons_last_coordinate_reverse.bed> \`  
                  `-b <bam_file_reverse.bam> \`  
                  `-counts \`  
                  `> <results_R.tsv>`

`cat <results_F.tsv> <results_R.tsv> > <results_1.tsv>`

**File (2)** is created with the following commands:

`bedtools coverage -a <bed_file_with_last_coordinate_window_forward.bed> \`  
                  `-b <bam_file_forward.bam> \`  
                  `-counts \`  
                  `> <results_F.tsv>`

`bedtools coverage -a <bed_file_with_last_coordinate_window_reverse.bed> \`  
                  `-b <bam_file_reverse.bam> \`  
                  `-counts \`  
                  `> <results_R.tsv>`

`cat <results_F.tsv> <results_R.tsv> > <results_2.tsv>`

## **plots_splicing_intermediates.py**

Takes as input 3 files:  
(1) cutoff file,  
(2) file with all events (exons last coordinate),  
(3) file with splicing intermediate peaks created by `get_peaks_splicing_intermediates.py`.

Gets all the peaks from file (3) that are above the cutoff from file (1).  
Creates a dataframe with total number of events (file (2)), number of peaks above cutoff and percentage.  
Saves that dataframe in a file and creates a barplot based on that dataframe.