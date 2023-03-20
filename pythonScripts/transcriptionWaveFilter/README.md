# **Python script to filter events according to transcription wave**

After DRB treatment, Polymerases restart transcription. According to time point after washing the drug, some downstream exons might not be transcribed yet.
POINT-seq was implemented in the same cells treated with DRB in identical conditions to define
transcription front wave (or fastest Poll II position) for a specific set of expressed genes
That information was used to filter events that were not transcribed yet in each time point.

## **transcription_wave_filter.py**

Takes as input a file with splicing intermediates coverage created with `bedtools coverage` or a file with splicing ratio per intron created with `get_splicing_ratio.py`. Returns a similiar file but only with the events expressed in that specific time point.  
Two auxiliary files are necessary:
(1) a file with events info (to get the distance from TSS).
(2) another with transcription front wave info, to discard events according to time point.  
If the time point is not specified, the program will try to deduce from the file name.  
Possible time points: CTR, DRB, W5, W10, W15, W30.

usage:  
`transcription_wave_filter.py [-h] [-r [replace]] \`  
`[<eventsInfo filepath>] [<transcriptionWave filepath>] \`  
`<sample name> <time point> <time point to filter by> <file type> \`  
`<input filepath>`


`file type`: Input file type. 2 options available: 1 for splicing intermediates; 2 for immediate splicing ratio.

`input filepath`: Path to input file. File with splicing intermediates coverage or with splicing ratio.