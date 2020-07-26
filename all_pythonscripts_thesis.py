#!/usr/bin/python3

'''
Within this file, all the python scripts that were used during the thesis
project "Transcriptome profiles associated with neuroprotection in the
Engrailed-1 mouse model for Parkinson’s disease" are available.

The scripts are listed in the order following the pipeline in which they were
used. Detailed information of the pipeline can be found in the README file.

Author(s): Dorottya Ralbovszki
Date: 2020.06.25.
'''

'''
Checking mapping rates
This script was used to compare the mapping rates of transcripts with
previous analysis after read quantification step.
Date: 2020.03.05.
'''
import sys
# creating empty dictionaries which are named after the person whose
# analysis data will populate that dictionary (Doresz is my nickname)
alfredo = dict()
doresz = dict()
# opening the two files to be compared the the output file
with open(sys.argv[1], 'r') as alf, open(sys.argv[2], 'r') as dor, open(
             sys.argv[3], 'w') as fout:
    # iterating through input file
    for line in alf:
        # splitting the line (string) by white space creating a list
        collums = line.split()
        # format of input files are the same starting with transctipt ID
        # transctipt ID becoming key in dictionary
        key = collums[0]
        # the value of mapping rate is becoming the value of
        # corresponding transcript ID wich is the key
        alfredo[key] = collums[1]
    # parsing through second input file and populating dictionary as above
    for line in dor:
        collums2 = line.split()
        key2 = collums2[0]
        doresz[key2] = collums2[1]
    # iterating through doresz dictionary both keys and values
    for key, value in doresz.items():
        # selecting transcript IDs that are part of both dictionaries
        if key in alfredo.keys():
            # comparing the values belonging to the same transcript ID
            if alfredo[key] == doresz[key]:
                # skipping the transcripts where the mapping rate is the same
                continue
            else:
                # printing the transcript ID and the mapping rate values
                # from the two dictionaries into standard output
                # to get an idea about the data
                print("Key: {}, Dval: {}, Aval: {}".format(key, value,
                                                           alfredo[key]))
    # printing keys that are only found in doresz dictionary into output file
    print(doresz.keys() - alfredo.keys(), file = fout, sep = '\n')



'''
Creating metadata
This script was used to create metadata tables so the quantification files can
be identified when imported into R environment.
Date: 2020.04.15.
'''
# opening the list containing the sample IDs and the output file
# starting with the table for SW wildtype vs C57 wildtype comparison
with open('ID_names_list.txt', 'r') as fin, open(
          'sample_metadata_WT_test.txt', 'w') as wt:
    # print headers in output
    print('Sample' + '\t' + 'Strain', file = wt)
    # parse input file (sample IDs)
    for line in fin:
        line=line.rstrip()
        # evalulate line and identify strains based on the input format
        # fourth character of the line can be used to distinguish between
        # the strains
        if line[3] == 'C':
            strain = 'C57'
            # identifying genotype and skipping the samples with En1+/-
            # (eg. mutant) genotype
            if '08_C57_2S_S10_' == line:
                # printing that line plus its strain in output
                # when the line fits the parameter
                print(line + '\t' + strain, file = wt)
            elif  '07_C57_10_S11_' == line:
                print(line + '\t' + strain, file = wt)
            elif '09_C57_43_S12_' == line:
                print(line + '\t' + strain, file = wt)
            else:
                continue
        # since there are only two strains, the ones that were not identified
        # as C57 are SW
        else:
            strain = 'SWISSOF1'
            # identifying genotype and skipping the samples with En1+/-
            # (eg. mutant) genotype
            if '02_SW_7_S8_' == line:
                print(line + '\t' + strain, file = wt)
            elif '01_SW_31_S9_' == line:
                print(line + '\t' + strain, file = wt)
            elif '03_SW_G_S7_' == line:
                print(line + '\t' + strain, file = wt)
            else:
                continue

# next, the table for SW En1+/- vs C57 En1+/- comparison
with open('ID_names_list.txt', 'r') as fin, open(
          'sample_metadata_MUT_test.txt', 'w') as mut:
    print('Sample' + '\t' + 'Strain', file = mut)
    for line in fin:
        line=line.rstrip()
        if line[3] == 'C':
            strain = 'C57'
            # identifying genotype and skipping the samples with
            # wildtype genotype
            if '08_C57_2S_S10_' == line:
                continue
            elif  '07_C57_10_S11_' == line:
                continue
            elif '09_C57_43_S12_' == line:
                continue
            else:
                # printing that line plus its strain in output
                # when the line fits the parameter
                print(line + '\t' + strain, file = mut)
        else:
            strain = 'SWISSOF1'
            if '02_SW_7_S8_' == line:
                continue
            elif '01_SW_31_S9_' == line:
                continue
            elif '03_SW_G_S7_' == line:
                continue
            else:
                print(line + '\t' + strain, file = mut)

# next, the table for SW wildtype vs SW En1+/- comparison
with open('ID_names_list.txt', 'r') as fin, open(
          'sample_metadata_SW_test.txt', 'w') as sw:
    print('Sample' + '\t' + 'Genotype', file = sw)
    for line in fin:
        line=line.rstrip()
        if line[3] == 'C':
            # skipping samples belonging to C57 strain
            continue
        else:
            strain = 'SWISSOF1'
            # identifying genotype
            if '02_SW_7_S8_' == line:
                genotype = 'wildtype'
            elif '01_SW_31_S9_' == line:
                genotype = 'wildtype'
            elif '03_SW_G_S7_' == line:
                genotype = 'wildtype'
            else:
                genotype = 'mutant'
        # since the C57 samples are skipped, only lines containing
        # SW samples will be printed
        print(line + '\t' + genotype, file = sw)

# next, the table for C57 wildtype vs C57 En1+/- comparison
with open('ID_names_list.txt', 'r') as fin, open(
          'sample_metadata_C57_test.txt', 'w') as c57:
    print('Sample' + '\t' + 'Genotype', file = c57)
    for line in fin:
        line=line.rstrip()
        if line[3] == 'C':
            strain = 'C57'
            # identifying genotype
            if '08_C57_2S_S10_' == line:
                genotype = 'wildtype'
            elif  '07_C57_10_S11_' == line:
                genotype = 'wildtype'
            elif '09_C57_43_S12_' == line:
                genotype = 'wildtype'
            else:
                genotype = 'mutant'
            # printing line and its genotype in output
            print(line + '\t' + genotype, file = c57)
        # skipping samples belonging to C57 strain
        else:
            continue

# and lastly, the table containing all samples
with open('ID_names_list.txt', 'r') as fin, open(
          'sample_metadata_all_test.txt', 'w') as all:
    print('Sample' + '\t' + 'Strain' + 'Genotype', file = all)
    for line in fin:
        line=line.rstrip()
        if line[3] == 'C':
            strain = 'C57'
            if '08_C57_2S_S10_' == line:
                genotype = 'wildtype'
            elif  '07_C57_10_S11_' == line:
                genotype = 'wildtype'
            elif '09_C57_43_S12_' == line:
                genotype = 'wildtype'
            else:
                genotype = 'mutant'
        else:
            strain = 'SWISSOF1'
            if '02_SW_7_S8_' == line:
                genotype = 'wildtype'
            elif '01_SW_31_S9_' == line:
                genotype = 'wildtype'
            elif '03_SW_G_S7_' == line:
                genotype = 'wildtype'
            else:
                genotype = 'mutant'
        # printing line with its strain and genotype in output
        print(line + '\t' + strain + genotype, file = all)



'''
Heatmap genes
This script was used to format the gene IDs so that they can be copied into R
environment for creating heatmaps. The format required is
"geneID",newline"geneID" ect.
Date:2020.05.07.
'''

import sys
# creating an empty list
original_id = []
# opening a file which contains the gene IDs in one column
with open(sys.argv[1], 'r') as fin:
    # parsing through the input file
    for line in fin:
        # stripping newline character off at the end of a line
        line = line.rstrip()
        # printing line (genen ID) into standard output in the desired format
        # the output from standard output was copied into R environment
        print('"{}",'.format(line))

'''
Parsing GO lists
This script was used to further filter the lists of Gene Ontology (GO) analysis.
Date:2020.05.22.
'''

import sys
# opening the GO list of Enrichr and the output file
with open(sys.argv[1], 'r')as fin, open(sys.argv[2], 'w') as fout:
    # parsing through input file
    for line in fin:
        # getting rid of newline characters at the end of lines
        line = line.rstrip()
        # splitting the lines by tab in order to be able to work with only
        # a certain column of input file, column variable is a list
        column = line.split("\t")
        # the list is filtered by Combined Score which is the 8th item
        # of the list so the 8th item of the list is printed into standard
        # output for checking if that is the Combined score value of the file
        print(column[7])
        # first line of the GO list is the header which is identified by
        # the first character of the 8th item of the line
        if column[7].startswith('C'):
            # the hedaer line is printed into output file
            print(line, file = fout)
        else:
            # filtering by Combined Score
            if float(column[7]) > 70:
                # lines with a greater Combined Score than 70 are printed
                # into output file
                print(line, file = fout)
            else:
                # skipping line that does not fulfill criteria
                continue
