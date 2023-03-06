TITLE: motif-mark
AUTHOR: Andreas Martinez
Description: motif-mark takes in a fasta file and a list of target motifs, it determines the exons, introns, target motif sequences and all possible degenerate bases/ambiguous nucleotides and then generates a color coded PNG image depicting these features.


INSTRUCTIONS
Execute file by running "./motif-mark-oop.py -f <fasta file> -m <motif file>"
    
motif-mark will output a PNG of the sequence lengths, exons and target motifs

ENVIRONMENT
conda create -n my_pycairo pycairo
conda activate my_pycairo

INPUT
1. fasta file
2. list of target motifs

OUTPUT
1. PNG graph of equence lengths, exons and color coded target motifs along with a legend. Legend and motif colors are modular and automatically customized depending on which values are found.

MODULES USED
argparse -> soft codes input files for modularized use
re -> a common module used to include regular expressions in code
pycairo -> the module that creates the graphical representation

COLLABORATORS
Andreas Martinez, Jacob Jensen, Sam Kupp, Jason Sydes