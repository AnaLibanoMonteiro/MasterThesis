# **Python scripts to create several bed files based on a list of transcripts**

## **create_bed_files.py**

Imports functions from `introns_and_windows_bed_files.py` and `transcripts_and_exons_bed_files.py`.  
Takes as input a file with a list of expressed transcripts IDs and a GTF file.  
Creates the following bed files (all forward and reverse):  
(1) transcripts  
(2) exons  
(3) exons last coordinate  
(4) exons last coordinate window  
(5) introns  
(6) introns windows  

usage:  
`create_bed_files.py [-h] [-o [outfileprefix]] [-d [outdir]] \`  
`[<gtf filepath>] <transcripts filepath>`  

`gtf filepath`: Path to GTF file. Must be converted to csv before.
`transcripts filepath`: Path to transcripts file. This file must have one column with the ids of expressed transcripts.

File (1) is used for metagenes analysis - see `metagenes` folder.  
Files (3) and (4) are used for splicing intermediates analysis - see `splicingIntermediates` folder.  
Files (5) and (6) are used for splicing ratio analysis - see `splicingRatio` folder.