#!/bin/bash

'''
Within this file, all the bash scripts that were used during the thesis
project "Transcriptome profiles associated with neuroprotection in the
Engrailed-1 mouse model for Parkinson’s disease" are available.

The scripts are listed in the order following the pipeline in which they were
used. Detailed information of the pipeline can be found in the README file.

Author(s): Dorottya Ralbovszki
Date:2020.06.25.
'''


'''
Salmon running on fastq files
This script was used to run salmon on multiple fastq files.
Date:2020.01.20.
'''

# first, preparing transcriptome indices
# the command follows the guidelines found at
# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
salmon index -t gencode.vM23.transcripts.fa.gz -i gencode_mus_index
# parsing through the fastq files found in that filepath
for fn in /Volumes/Dorottya/Seq029_RNASequencingData/Seq029_Alfredo/190312/Data/Intensities/Basecalls/*_;
do
  # printing the file names into output file to be able to keep track of loop
  echo "${fn}" >> list.txt
  # quantification step followoing the guidelines found at
  # https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode
  salmon quant -i /Volumes/Dorottya/salmon_test/gencode_mus_index -l IU -1 ${fn}/*_R1*.fastq.gz -2 ${fn}/*_R2*.fastq.gz --validateMappings -o ${fn}/quantification_gencode
  # printing the file name into strandard output to be able to keep track of loop
  echo "${fn}"
done


'''
Gzipping quantification files for better storage
This script was used to gzip quantification files for long-term storage.
Date:2019.11.26.
'''

for fn in /mnt/c/Users/Doresz/'OneDrive - Lund University'/thesis_project/salmon_data/*_;
do
  echo "${fn}"
  gzip "${fn}"/quantification2/quant.sf
done
