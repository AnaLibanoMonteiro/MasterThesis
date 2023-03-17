# **Python scripts for splicing ratio analysis**

Scripts developed to calculate introns splicing ratio using nascent RNA data.

## **get_immediate_splicing_ratio.py**
Takes as input the result file from an intersection between a bed file with introns windows and a bam file.  
Creates a csv file with the splicing ratio for each event (intron) considered in the bed file.

Usage:  
`get_immediate_splicing_ratio.py [-h] [-o [outfilename]] [-d [outdir]] <input filepath>`

### **Input files are created with the following commands:**  

`bedtools intersect -a <introns_windows_forward.bed> \`  
                   `-b <bam_file_forward.bam> \`  
                   `-wao \`  
                   `-split \`  
                   `> <results_F.tsv>`

`bedtools intersect -a <introns_windows_reverse.bed> \`  
                   `-b <bam_file_reverse.bam> \`  
                   `-wao \`  
                   `-split \`  
                   `> <results_R.tsv>`

`cat <results_F.tsv> <results_R.tsv> > <results.tsv>`

## **plots_immediate_splicing_ratio.py**
Takes as input 2 files:  
(1) immediate splicing ratio created by `get_immediate_splicing_ratio.py` and  
(2) a cutoff file.  
Creates another file with the events on (1) that have total number reads > cutoff.  
Creates several plots based on the splicing ratio.

Usage:  
`plots_immediate_splicing_ratio.py [-h] [-c [<cutoff filepath>]] [-d [outdir]] <immediateSplicingRatio filepath> <sample group>`  
`-c` or `--cutoffFilepath` is the path for the file that has the cutoff values - minimun number of reads that an intron must have to be considered for the analysis.  
It has the following columns:  
`sample_name` | `sample_group` | `time_point` | `cutoff`