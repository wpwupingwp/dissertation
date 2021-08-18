#!/usr/bin/python3

from pathlib import Path
import argparse
import logging

import numpy as np
from Bio import SeqIO

# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)
default_level = logging.INFO
TEMP_LOG = 'analyze_gap.log'

import coloredlogs
coloredlogs.install(level=default_level, fmt=FMT, datefmt=DATEFMT)
log_file = logging.FileHandler(TEMP_LOG)
log_file.setFormatter(formatter)
log_file.setLevel(default_level)
log = logging.getLogger(__name__)


def init_arg(arg):
    if arg.input is None or len(arg.input)==0:
        raise ValueError(f'{arg.input} does not exist.')
    else:
        arg.input = [Path(i).absolute() for i in arg.input]
    return arg


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('-input', nargs='*', help='input alignment fasta')
    # 1.25 means if one seq's len is 1.25 fold than mean len, then it's too
    # long
    arg.add_argument('-type', choices=('align', 'raw'), default='raw')
    arg.add_argument('-len_ratio', default=1.25, type=float,
                     help='length threshold')
    arg.add_argument('-out', help='output file')
    return arg.parse_args()


def fasta_to_array(aln_fasta: Path) -> tuple[np.array, np.array]:
    """
    From BarcodeFinder.evaluate
    Given fasta format alignment filename, return a numpy array for sequence:
    Faster and use smaller mem.
    Ensure all bases are capital.
    Args:
        aln_fasta(Path): aligned fasta file
    Returns:
        name(np.array): name array
        sequence(np.array): sequence array
    """
    data = []
    record = ['id', 'sequence']
    with open(aln_fasta, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('>'):
                data.append([record[0], ''.join(record[1:])])
                # remove ">" and CRLF
                name = line[1:].strip()
                record = [name, '']
            else:
                record.append(line.strip().upper())
        # add last sequence
        data.append([record[0], ''.join(record[1:])])
    # skip head['id', 'seq']
    data = data[1:]
    # check sequence length
    length_check = [len(i[1]) for i in data]
    if len(set(length_check)) != 1:
        log.error(f'Invalid alignment file {aln_fasta}')
        return None, None
    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name_array = np.array([[i[0]] for i in data], dtype=np.bytes_)
    # fromiter is faster than from list
    # S1: bytes
    sequence_array = np.array(
        [np.fromiter(i[1], dtype=np.dtype('S1')) for i in data],
        order='F')
    if name_array is None:
        log.error('Bad fasta file {}.'.format(aln_fasta))
    return name_array, sequence_array


def count_aln_len(fasta, threshold) -> tuple[np.array, np.array, float]:
    name, seq = fasta_to_array(fasta)
    # name is ndarray
    name = name[:, 0]
    name = np.array([i.decode() for i in name])
    no_gap = (seq!=b'-')
    count = np.count_nonzero(no_gap, axis=1)
    threshold = count.mean() * threshold
    too_long = (count>threshold)
    count_array = np.stack((name, count)).T
    too_long_array = np.stack((name[too_long], count[too_long])).T
    return count_array, too_long_array, threshold


def count_len(fasta, threshold) -> tuple[dict, dict, float]:
    data = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    name_len = {i: len(data[i]) for i in data}
    threshold = np.mean(list(name_len.values())) * threshold
    too_long = {i:name_len[i] for i in name_len if name_len[i]>threshold}
    return name_len, too_long, threshold


def main():
    """
    Count original length of alignment's sequences.
    """
    log.info('Start.')
    arg = parse_args()
    arg = init_arg(arg)
    log.info('Input type:'+arg.type)
    if arg.type == 'align':
        for fasta in arg.input:
            count_array, too_long_array, threshold = count_aln_len(
                fasta, arg.len_ratio)
            with open(fasta.with_suffix('.csv'), 'w') as out1:
                np.savetxt(out1, count_array, delimiter=',', fmt='%s')
            with open(fasta.with_suffix('.too_long'), 'w') as out2:
                np.savetxt(out2, too_long_array, delimiter=',', fmt='%s')
    else:
        for fasta in arg.input:
            count_dict, too_long_dict, threshold = count_len(
                fasta, arg.len_ratio)
            with open(fasta.with_suffix('.csv'), 'w') as out1:
                out1.write('Name,Len\n')
                for key, value in count_dict.items():
                    out1.write(f'{key},{value}\n')
            with open(fasta.with_suffix('.too_long'), 'w') as out2:
                out2.write('Name,Len\n')
                for key, value in too_long_dict.items():
                    out2.write(f'{key},{value}\n')
    log.info('Bye.')


if __name__ == '__main__':
    main()
