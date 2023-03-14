# **Python scripts for Single Nucleotide Resolution (SNR)**
Get the last nucleotide of the reads to get Pol II position (when analysing nascent RNA, the last nucleotides of the reads correspond, in principle, to the last nucleotide incorpored by Pol II).  
These scripts were based on a [previous script](https://github.com/rluis/mNET_snr), that would **discard** reads that had **soft clipping** (nucleotides on the ends of the reads that do not match with the reference genome).  
These scripts **accept** reads with **soft clipping** (if it is only 1 nucleotide soft clipped).  
`get_SNR_3prime-ALM.py` considers the soft clipped nucleotides **as part of the fragment**. (Used for the downstream analysis in this project.)  
`get_SNR_3prime_not_part_of_fragment-ALM.py` discards the soft clipped nucleotides.  

## **get_SNR_3prime-ALM.py**
Takes as input a bam file.
Creates a bam file with the last nucleotide of the reads.
LAST base from reads with flag 147.
FIRST base from reads with flag 163.

Usage:
TODO!!

## **get_SNR_3prime_not_part_of_fragment-ALM.py**
Same as previous. (The only difference is that the soft clipped nucleotides are not considered as part of the fragment)


### **Diagram of reads and flags:**

#### **If RNA molecule comes from DNA forward strand:**

    R1 - flag 99
5' ---------------> 3'  
5' ----------------------------------------------------------> 3'  
                                           5' ---------------> 3'
                                               R2 - flag 147

*Get LAST nucleotide of R2 (147)*

#### **If RNA molecule comes from DNA reverse strand:**

                                                R1 - flag 83
                                           5' ---------------> 3'
3' <---------------------------------------------------------- 5'
5' ---------------> 3'
    R2 - flag 163

*Get FIRST nucleotide of R2 (163)*