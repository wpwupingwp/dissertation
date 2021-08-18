#!/usr/bin/python3

from Bio import SeqIO
from sys import argv
from pathlib import Path

rotated_need_reannotate = open('rotated_need_reannotate.csv', 'a')
no_rotate = open('no_rotate.csv', 'a')
same = open('same.csv', 'a')
raw_gb = Path(argv[1]).absolute()
rotated_fasta = Path(argv[1]+'-out/'+argv[1]+'.fasta')
if not rotated_fasta.exists():
    no_rotate.write(argv[1]+'\n')
else:
    raw = SeqIO.read(raw_gb, 'gb')
    rotated = SeqIO.read(rotated_fasta, 'fasta')
    if hash(raw.seq) == hash(rotated.seq):
        same.write(argv[1]+'\n')
    else:
        rotated_need_reannotate.write(argv[1]+'\n')
