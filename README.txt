## Introduction
BLAT is bioinformatics software used to produce alignments of transcripts from one species to the genome of another. In addition to large alignments that cover most of the transcript, BLAT also outputs a number of small ones where only a small portion of a transcript matches the genome, which is considered as noise. This script can be used to get rid of this noise and output clean BLAT GFF files.


## Instructions
You need a python environment with gffutils and biopython installed:

`pip install gffutils biopython`


### Usage: 
`python filter.py <gff3 filename> <fasta filename> [-h] [-l lengththreshold] [-s scorethreshold] [-o output]`

#### Flags:
`-l` specifies a length threshold, specifically the minimum allowed ratio between the length of the transcript found in the gff3 file and the length of the transcript in the FASTA file. Default: .5

`-s` specifies the minimum score threshold. Default: 0

`-o` specifies the output filename. Default: output_\<datetime>.gff3

`-h` shows the help message and exits.

### Example:

`python filter.py -o my_output.gff3 -l .4 -s 100 my_gff3_file.gff3 my_fasta_file.fasta `

This outputs a gff3 file named my_output.gff3 filtered from a gff3 file and a corresponding fasta file. The length threshold in this case is .4 (so it filters out matches that are not at least 40% of the length of the sequence in the fasta file) and the score threshold is 100 (so it filters out matches with scores that are less than 100).
