# **Python scripts developed to remove peaks and plot metagenes**

When `deeptools plotProfile` is used rigth after `deeptools computeMatrix`, there are strong peaks that prevent to see the average distribution of Pol II across expressed genes.
Therefore, those peaks are removed and the metagenes are plotted after.  
After replacing highest values by Nan, `deeptools plotProfile` no longer recognizes the matrix, so `plot_metagenes.py` was developed.

## **remove_matrix_peaks_metagenes.py**

Takes as input a matrix from `deeptools computeMatrix`.  
Removes the peaks on each bin (column):  
goes to each column and replace the max value by Nan;  
does this several times, according to the number of peaks to remove.

usage:  
`remove_matrix_peaks_metagenes.py [-h] [-r [replace]] [-p [percentageToRemove]] \`  
`<input matrix filepath>`  

`-r` or `--replace`: Boolean. True: remove original file. False: keep both files (original and new one with extension '_new.mat.gz'). Default: True.  
`-p` or `--percentageToRemove`: Percentage to remove. Default: 0.005. If intput matrix has 2000 lines, top 10 values of each column will be replaced by Nan.  


## **plot_metagenes.py**

Takes as input a matrix (or a list of matrices) from `deeptools computeMatrix` after `remove_matrix_peaks_metagenes.py`.  
Creates a metagene plot for each matrix and a plot with all matrices (in case there is more than one in the input) to be compared.

usage:  
`plot_metagenes_matrix.py [-h] [-d [outdir]] <sample Name> <sample group> <input matrix filepath>`


`samples Name`: List with sample names (for eg: 601CTR, 601DRB, 601W30).  
`samples group`: Name of sample group (for eg: 601, 603 or Y1P).  
`input matrix filepath`: List with paths to input matrix/matrices.
