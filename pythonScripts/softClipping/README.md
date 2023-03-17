# **Python scripts for soft clipping analysys**

Scripts developed to explore unsual soft clipping patterns detected in mNET-seq data.  
Soft clipped reads have one or more nucleotide at the 5' or 3' end that do not match with the reference genome.

## **get_soft_clipped_info.py**
Takes as input a bam file.  
Creates a bam file with the soft clipped reads (on the same folder as the input).  
Creates a csv file with soft clipped reads information (saves on the output folder) based on bam file created previously.  
Creates a bed file with coordinates of soft clipped regions (saves on the output folder) that will be used to complete the csv file with genome reference sequence (to see what nucleotides were supposed to be there).

Usage:  
`get_soft_clipped_info.py [-h] [-o [outfilename]] [-d [outdir]] [-b [bamfileonly]] <input bam filepath>`

### **Output CSV file has the following columns:**

read_name | chr | mapping_start | read_size | cigar | read_seq | soft_start | soft_end | numb_soft_nucleotides | soft_nucleotide_seq | one_or_both_ends | prime_end | flag | strand

## **complete_soft_clipped_info.py**
Takes as input a csv file with soft clipped reads information and a fasta file with real sequence of soft clipeed regions.  
Adds a colum to the csv file with the real sequence of soft clipped regions.  

Usage:  
`complete_soft_clipped_info.py [-h] [-r [replace]] <input csv filepath> <input fasta filepath>`  
`-r` or `--replace` option will replace the csv file with the completed one.

### **Fastas files are created with:**  

`bedtools getfasta -fi <Human Genome FASTA file> \`  
`-bed <BED file with soft clipped regions> -s -name -tab -bedOut \`  
`> output.fasta`



## **plot_soft_clipped_info.py**
Takes as input the csv file with soft clipped reads information completed with genome reference sequence.
Creates several plots describing soft clipping info.

Usage:  
`plot_soft_clipping_info.py [-h] [-d [outdir]] <sample Name> <input csv filepath>`

## **get_soft_clipped_reads_3prime.py**
Takes as input a bam file.  
Creates a bam file with reads that have soft clipping on 3' end.  
These bam files were compared with normal bam files to check if the soft clipping nucleotides should be considered as part of the read or not (based on splicing intermediates peaks).

Usage:  
`get_soft_clipped_reads_3prime.py [-h] [-d [outdir]] [-o [outfilename]] <input bam filepath>`