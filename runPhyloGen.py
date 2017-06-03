#!/usr/bin/python

import sys

#This will only work if phyloGenerator.py is on your path, if not use something like os.path.append()
import phyloGenerator as pG
from Bio import SeqIO
#You will only be able to run programs like MAFFT and RAxML if you have a 'requires' folder in your current working directory
#Run the 'setupLinux.py' script that is on the pG GitHub to generate such a folder
sys.argv[1] = int(sys.argv[1])
for i in range(sys.argv[1]):
    seqs = list(SeqIO.parse(str(i+1)+"simulatedSeqs.fasta", "fasta"))
    seqs = [[x] for x in seqs]

    mafft = pG.alignSequences(seqs, method = 'mafft')

    output = pG.BEAST(mafft[0][0], runNow=False)

    print(output)
