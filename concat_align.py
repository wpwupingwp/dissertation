#!/usr/bin/python3

import argparse
import logging
from pathlib import Path
from subprocess import run

from Bio.Alphabet import IUPAC
from Bio.Nexus import Nexus
from Bio import SeqIO

# define logger
import coloredlogs
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)
default_level = logging.INFO
coloredlogs.install(level=default_level, fmt=FMT, datefmt=DATEFMT)
log = logging.getLogger(__name__)


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('-input', required=True, nargs='*',
                     help='input fasta files')
    arg.add_argument('-out', default='combine.nex',
                     help='output filename')
    return arg.parse_args()


def convert(tmp, files):
    py = '~/git/phd/prepare_for_tree.py'
    converted = []
    for fasta in files:
        log.info(f'Reading {fasta}')
        clean_fasta = tmp / fasta.with_suffix('.tmp')
        run(f'python3 {py} {fasta} {clean_fasta}', shell=True)

        nexus = tmp / fasta.with_suffix('.nexus')
        SeqIO.convert(clean_fasta, 'fasta', nexus, 'nexus',
                      alphabet=IUPAC.ambiguous_dna)
        converted.append(nexus)
        clean_fasta.unlink()
    return converted


def concat(files, out):
    data = [(nex.stem, Nexus.Nexus(nex)) for nex in files]
    combine = Nexus.combine(data)
    combine.write_nexus_data(filename=out)
    return out


def main():
    """
    concat fasta alignments to nexus
    Assume files were processed by prepare_for_tree
    """
    arg = parse_args()
    arg.input = [Path(i) for i in arg.input]
    arg.out = Path(arg.out)
    tmp = Path('tmp')
    tmp.mkdir()
    log.info(f'Input files: {[i.name for i in arg.input]}')
    nex_files = convert(tmp, arg.input)
    log.info('Concatenate...')
    concat(nex_files, arg.out)
    for i in nex_files:
        i.unlink()
    with open(arg.out, 'r') as _:
        # big mem
        raw = _.readlines()
    # remove invalid line
    if raw[-2].startswith('charpartition'):
        with open(arg.out, 'w') as out:
            for line in raw[:-2]:
                out.write(line)
            out.write(raw[-1])
    tmp.rmdir()
    log.info(f'Output file: {arg.out}')
    log.info('bye')


if __name__ == '__main__':
    main()
