import gffutils
from Bio import SeqIO
import csv
import sqlite3
import argparse
import os
from datetime import datetime

parser = argparse.ArgumentParser(description='Filter GFF3 files.')
parser.add_argument('gff3', help='the gff3 file to be filtered')
parser.add_argument('fasta', help='the corresponding fasta file')
parser.add_argument('-l', '--lengththreshold', type=float, 
    help='Minimum threshold for length (some % of the length of the sequences in the fasta file) (Default: .5)')
parser.add_argument('-s', '--scorethreshold', type=int, 
    help='Minimum score threshold (Default: 0)')
parser.add_argument('-o', '--output', help='Name of output file (Default: output_<datetime>.gff3)')
args = parser.parse_args()

LENGTH_THRESHOLD = args.lengththreshold if args.lengththreshold is not None else .5
if LENGTH_THRESHOLD > 1 or LENGTH_THRESHOLD < 0:
    raise ValueError('Length threshold must be a value between 0 and 1.')

SCORE_THRESHOLD = args.scorethreshold if args.scorethreshold is not None else 0
if SCORE_THRESHOLD < 0:
    raise ValueError('Score threshold must be a value greater than or equal to zero.')

gff3_file = args.gff3
fasta_file = args.fasta
output_file = args.output if args.output is not None else 'output_{}.gff3'.format(
    datetime.now().strftime('%Y%m%d_%H%M%S'))
db_name = 'gff3_{}.db'.format(datetime.now().strftime('%Y%m%d_%H%M%S'))

DB = gffutils.create_db(gff3_file, db_name, keep_order=True)

# get lengths of sequences from fasta file
lengths = {}
for record in SeqIO.parse(fasta_file, 'fasta'):
    lengths[record.id] = len(record.seq)

# filter matches
to_delete = {}
num_deleted_match = 0
num_deleted_match_part = 0
for a in DB.features_of_type('match'):
    if not a.attributes.get('m_stop', False):
        continue

    start = 0 if 'm_start' not in a.attributes else int(a.attributes['m_start'][0])
    length = int(a.attributes['m_stop'][0]) - start
    score = int(a[5])
    match = a.attributes['match_id'][0]
    gid = a.attributes['ID'][0]

    if length < LENGTH_THRESHOLD * lengths[match] or score < SCORE_THRESHOLD:
        to_delete.update({x.attributes['ID'][0]: 1 for x in [a, *list(DB.children(gid))]})
        num_deleted_match += 1
        num_deleted_match_part += len(list(DB.children(gid)))

# include gff directive at top
f = open(output_file, 'w')
f.write('##gff-version 3\n')
f.close()

# write filtered data to gff3 file
with sqlite3.connect(db_name) as connection:
    writer = csv.writer(open(output_file, 'a'), delimiter='\t')
    c = connection.cursor()
    c.execute('SELECT * FROM features')
    for row in c:
        attrs = eval(row[9])
        if row[0] in to_delete:
            continue

        new_str = []
        for k, v in attrs.items():
            if k != 'original_id' and k != 'original_name':
                new_str.append("{}={};".format(k, v[0]))
        new_str = "".join(new_str)
        writer.writerow([row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], new_str])

# print filtering statistics
original_match = DB.count_features_of_type('match')
original_match_part = DB.count_features_of_type('match_part')
filter_match = original_match - num_deleted_match
filter_match_part = original_match_part - num_deleted_match_part

table = [
    ["", "Original", "Filtered"],
    ["match", original_match, filter_match], 
    ["match_part", original_match_part, filter_match_part],  
    ["Ratio (match_part:match)", 
        "{0:.3f}".format(original_match_part/original_match),
        "{0:.3f}".format(filter_match_part/filter_match)
    ]
]

longest_cols = [(max([len(str(row[i])) for row in table]) + 3) for i in range(len(table[0]))]
row_format = "".join(["{:>" + str(longest_col) + "}" for longest_col in longest_cols])
for row in table:
    print(row_format.format(*row))

os.remove(db_name)
