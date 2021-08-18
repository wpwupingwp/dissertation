#!/usr/bin/python3

import argparse
from pathlib import Path

from Bio import SeqIO


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    # id, input, lenght, lsc, ir, ssc, ir
    arg.add_argument('-f', help='structure table file', required=True)
    arg.add_argument('-t', dest='threshold', type=int, default=10,
                     help='overlap threshold')
    arg.add_argument('-o', default='out',
                     help='output directory')
    return arg.parse_args()


def compare_slice(a: slice, b: slice) -> tuple[str, int]:
    # b broader than b
    overlap_len = 0
    relation = 'undefined'
    if (a.start >= a.stop or b.start >= b.stop or
            b.stop-b.start <= a.stop-a.start):
        return relation, overlap_len
    # a in b
    if b.start <= a.start <= a.stop <= b.stop:
        relation = 'in'
    # a after b
    elif b.start < b.stop <= a.start < a.stop:
        relation = 'after'
    # a before b
    elif a.start < a.stop <= b.start < b.stop:
        relation = 'before'
    elif b.start < a.start < b.stop < a.stop:
        overlap_len = b.stop - a.start
        relation = 'overlap'
    elif a.start < b.start < a.stop < b.stop:
        overlap_len = a.stop - b.start
        relation = 'overlap'
    else:
        pass
    return relation, overlap_len


def get_table(csv_file) -> dict:
    # id, input, length, lsc, ir, ssc, ir
    id_loc = dict()
    with open(csv_file, 'r') as raw:
        for line in raw:
            id_, _, _, lsc, ir, ssc, _ = line.split(',')
            id_loc[id_] = (lsc, ir, ssc)
    return id_loc


def get_gb() -> dict:
    id_file = dict()
    # gb in current folder
    for gb in Path().glob('*.gb'):
        id_ = gb.stem.split('-')[0]
        id_file[id_] = gb
    return id_file

import pickle
def get_all_genes(gb_file) -> list:
    gene_features = []
    gb = SeqIO.read(gb_file, 'gb')
    for feature in gb.features:
        if feature.type != 'gene':
            continue
        if 'gene' not in feature.qualifiers:
            continue
        # feature are ordered
        gene_features.append(feature)
    return gene_features


def get_part_genes(gb_file, raw_loc) -> tuple[dict, list]:
    overlap = []
    # raw_loc: [length, length, length]
    loc_n = [int(i) for i in raw_loc]
    loc_range = [slice(0, loc_n[0]), slice(loc_n[0], loc_n[0]+loc_n[1]),
                 slice(loc_n[0]+loc_n[1], loc_n[0]+loc_n[1]+loc_n[2])]
    loc_index = 0
    loc_index_end = len(loc_range)
    name = ['LSC', 'IR', 'SSC', 'IR']
    part_gene = {i: list() for i in name}
    gene_features = get_all_genes(gb_file)
    # ignore 2nd IR
    len_genes = len(gene_features)
    gene_index = 0
    while gene_index < len_genes and loc_index < loc_index_end:
        feature = gene_features[gene_index]
        gene = feature.qualifiers['gene'][0]
        if gene == 'rps12':
            gene_index += 1
            continue
        part_slice = loc_range[loc_index]
        feature_slice = slice(feature.location.start, feature.location.end)
        cmp, overlap_len = compare_slice(feature_slice, part_slice)
        if cmp in ('in', 'before'):
            part_gene[name[loc_index]].append(gene)
        elif cmp == 'after':
            loc_index += 1
        elif cmp == 'overlap':
            part_gene[name[loc_index]].append(gene)
            overlap.append((gene, overlap_len,
                            int(feature_slice.start), int(feature_slice.stop),
                            part_slice.start, part_slice.stop))
            # print(gene, feature.location, part_slice, 'overlaped')
            part_gene[name[loc_index+1]].append(gene)
            loc_index += 1
        else:
            print(gene, feature.location, part_slice, cmp, 'bad')
            break
        gene_index += 1
    return part_gene, overlap


def get_genes(id_file, id_loc) -> dict:
    id_genes = dict()
    id_overlap = dict()
    for id_, gb_file in id_file.items():
        loc = id_loc[id_]
        part_gene, overlap = get_part_genes(gb_file, loc)
        id_genes[id_] = part_gene
        id_overlap[id_] = overlap
    return id_genes, id_overlap


def output_overlap(id_overlap, out_file) ->Path:
    with open(out_file, 'w') as out:
        out.write('ID,gene,overlap_len,gene_start,gene_end,'
                  'part_start,part_end\n')
        for id_, records in id_overlap.items():
            for record in records:
                s = '{},{},{},{},{},{}\n'.format(id_, *record)
                out.write(s)
    return out_file


def main():
    """
    Docstring.
    """
    arg = parse_args()
    arg.out = Path(arg.o)
    global threshold
    threshold = arg.threshold
    id_file = get_gb()
    id_loc = get_table(arg.f)
    id_genes, id_overlap = get_genes(id_file, id_loc)
    #print(id_genes.items())
    out_file = output_overlap(id_overlap, arg.out.with_suffix('.overlap.csv'))
    print(out_file)


if __name__ == '__main__':
    main()
