# HLA Typing Pipeline

## HLA Pipeline:

The Python script should be run as follows:  
**`$ python HLApipeline/runner.py <options>`**
This program automates the running of HLA typing algorithms, as well as the merging of the outputs and the accuracy analysis.
(Total number of lines of code: 648)

### Options: 
**`-h` or `--help` :** displays the list of options and exits the program
**`-g` or `--input_WGS` :** specify the path to the whole genome sequence file
**`-e` or `--input_WES` :** specify the path to the whole exome sequence file
**`-r` or `--input_RNAseq` :** specify the path to both RNAseq sequence files (e.g. path/to/RNAseq1,path/to/RNAseq2)
**`-l` or `--RNAseq_min_length` :** specify RNAseq sequence length threshold for mapping
**`-G` or `--genome_dir` :** specify the path to directory containing all the genome files
**`-p` or `--population` :** specify the population origin of the sample
**`-o` or `--output_dir` :** specify the output directory for all the algorithms
**`-t` or `--num_threads` :** specify the number of threads to be used for parallel computing
**`-H` or `--dir_HLA_HD` :** specify the directory containing all the HLA-HD program files
**`-s` or `--dir_seq2HLA` :** specify the directory containing all the seq2HLA program files
**`-K` or `--dir_Kourami` :** specify the directory containing all the Kourami program files
**`-a` or `--dir_arcasHLA` :** specify the directory containing all the arcasHLA program files
**`-L` or `--dir_HLA_LA` :** specify the directory containing all the HLA-LA program files
**`-P` or `--dir_picard` :** specify the directory containing all the picard program files
**`-c` or `--correct_HLA`:** specify the CSV file containing the correct HLA typing information from PCR-Sanger sequencing

## Majority Voting:

The Python script should be run as follows:  
**`$ python majority_voting.py <path/to/results.csv>`**
This script runs the majority voting algorithm on the results of the HLA typing algorithms.
(Total number of lines of code: 55)

## Modify Sequences:

The Golang program should be run as follows:
**`$ go run goModify/modifyfasta/modifyfasta.go -w=<# of worker threads> -r=<# of reader threads> -p=<proportion of sequences to keep> -l=<length of subsequences> <path/to/input.fasta> <path/to/output.fasta`**
This program takes as input a fasta file and outputs a modified version with shorter read lengths and/or fewer reads. The program has been optimized for parallel computing.
(Total number of lines of code: 177)