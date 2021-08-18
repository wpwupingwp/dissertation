#!/usr/bin/python3

import pandas as pd
import logging
from subprocess import run
from glob import glob


FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%I:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO,
                    handlers=[logging.FileHandler('Log.txt'), ])
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger(__name__)


folder_taxon = {}
with open('list.csv', 'r') as _:
    for line in _:
        _, taxon, folder = line.strip().split(',')
        folder_taxon[folder] = taxon
handle = open('log.txt', 'w')
for i in folder_taxon:
    try:
        left = list(glob('raw/'+i+'/*1.fq'))[0]
    except IndexError:
        continue
    right = list(glob('raw/'+i+'/*2.fq'))[0]
    print(left, right)
    if ' ' in folder_taxon[i]:
        cmd = f'python3 novowrap.py -f {left} -r {right} -taxon "{folder_taxon[i]}" -reads_len 100 -mem 20'
    else:
        cmd = f'python3 novowrap.py -f {left} -r {right} -taxon {folder_taxon[i]} -reads_len 100 -mem 20'
    log.info(cmd)
    r = run(cmd, shell=True)
