# MSc thesis project

The files within this repository are scripts that were produced during a MSc project (in Molecular biology) involving differential gene expression analysis. Further information can be found in the thesis report, available on request.

## Workflow summary



## Scripts

All available scripts, including the R-scripts, are non-executable. Meaning, individual commands should be run in succession.

### [all_bash_scripts.sh](./all_bash_scripts.sh)

Contains two different bash tasks:

- Salmon running on fastq files
- Gzipping quantification files for better storage

### [all_pythonscripts_thesis.py](./all_pythonscripts_thesis.py)

Contains four separate python scripts that, if ran on their own, can be executed to perform their individual tasks. All except the 'Creating metadata', takes one or more command-line argument(s) to determine input file.

- Checking mapping rates
- Creating metadata
- Heatmap genes
- Parsing GO lists

### R-scripts

Four of the R-scripts contain a variation of the same workflow, and has more or less identical code--analyzing either inter or intra-strain comparisons. The exception being (tximport_deseq2_script_multi.R)[./tximport_deseq2_script_multi.R], which features a multi-factor design executing all four previous comparisons at the same time.

Only DEGs produced using the first four scripts were included in the final report, whereas the result from the multi-factor design-script was only used to produce a heatmap.

## Author(s)

- Dorottya Ralbovszki
