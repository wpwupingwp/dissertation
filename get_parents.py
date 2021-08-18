#!/usr/bin/python3

import argparse
import sqlite3
import os


def taxon_query(query_arg, cur):
    """
    Query taxon database using taxon_id or Name.  """
    condition = list()
    if query_arg.species_name is not None:
        condition.append(
            "Taxon_name LIKE '{0}%'".format(query_arg.species_name))
    elif query_arg.taxon_name is not None:
        condition.append(
            "Taxon_name = '{0}'".format(query_arg.taxon_name))
    if query_arg.taxon_id is not None:
        condition.append('Taxon_id = {0}'.format(query_arg.taxon_id))
    condition = ' OR '.join(condition)
    query = 'SELECT * FROM {0} WHERE {1};'.format('Taxon', condition)
    cur.execute(query)
    # % to match records like "species var."
    result = cur.fetchall()
    if len(result) == 0:
        print((query_arg.species_name or query_arg.taxon_name or
               query_arg.taxon_id), 'NOT_FOUND', sep='\t')
    return result


def print_taxon_query(cur, results):
    for result in results:
        taxon_id, taxon_name, rank, parents, sons, greatsons = result
        # remove parenthesis
        greatsons = greatsons[1:-1]
        print('{}\t{}'.format(taxon_name, id_to_text(cur, parents)))


def id_to_text(cur, taxon_id_list):
    if taxon_id_list == '':
        return 0, ''
    taxon_id_list = taxon_id_list.split(',')
    name_list = list()
    for taxon_id in taxon_id_list:
        cur.execute("""SELECT Taxon_name FROM Taxon WHERE
                    Taxon_id = {0};""".format(taxon_id))
        result = cur.fetchone()
        if result is None:
            name_list.append('')
        else:
            name_list.append('{}'.format(result[0]))
    return '\t'.join(name_list)


def parse_args():
    arg = argparse.ArgumentParser(description='Query sequence.')
    group = arg.add_mutually_exclusive_group()
    group.add_argument('-s', dest='species_name', help='species name')
    group.add_argument('-t', dest='taxon_name',
                       help='taxonomy name (higher than species)')
    group.add_argument('-tid', dest='taxon_id', help='taxonomy id')
    arg.add_argument('-f', dest='fragment_type', help='fragment type',
                     choices=('gene', 'CDS', 'spacer', 'tRNA', 'rRNA',
                              'intron', 'whole'))
    arg.add_argument('-g', dest='gene', help='gene name')
    arg.add_argument('-u', '--upstream', type=int, default=0,
                     help='upstream location')
    arg.add_argument('-d', '--downstream', type=int, default=0,
                     help='downstream location')
    arg.add_argument('-o', dest='output', default='out', help='output path')
    # arguments.print_help()
    return arg.parse_args()


def main():
    DATA = 'Data'
    SEQ = 'seq.sql'
    TAXON = 'taxon.sql'
    # global query_arg
    query_arg = parse_args()
    con_taxon = sqlite3.connect(os.path.join(DATA, TAXON))
    con_seq = sqlite3.connect(os.path.join(DATA, SEQ))
    cur_taxon = con_taxon.cursor()
    cur_seq = con_seq.cursor()

    if not os.path.exists(query_arg.output):
        os.mkdir(query_arg.output)
    result = taxon_query(query_arg, cur_taxon)
    print_taxon_query(cur_taxon, result)
    cur_taxon.close()
    con_taxon.commit()
    con_taxon.close()
    cur_seq.close()
    con_seq.commit()
    con_seq.close()


if __name__ == '__main__':
    main()
