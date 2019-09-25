import gffutils
from Bio import SeqIO
import csv
import sqlite3

# FIXME autogenerate db (filename)
DB = 'gff_run1.db'

db = gffutils.FeatureDB(DB, keep_order=True)
lengths = {}

# FIXME change to the parsed flag
filename = 'humpback_whale_megNov1_transcripts.final.fasta'


for record in SeqIO.parse(filename, 'fasta'):
    lengths[record.id] = len(record.seq)

# get the length of the longest match for each transcript
# if just comparing the matches to the ones in the fasta file, delete this
#for a in db.features_of_type('match'):
#    if not a.attributes.get('m_stop', False) or not a.attributes.get('m_start', False):
#        continue
#
#    length = int(a.attributes['m_stop'][0]) - int(a.attributes['m_start'][0])
#    gene = a.attributes['match_id'][0]
#
#    if not matches.get(gene, False):
#        matches[gene] = length
#    else:
#        if matches[gene] < length:
#            matches[gene] = length

# filter matches
to_delete = {}
for a in db.features_of_type('match'):
    if not a.attributes.get('m_stop', False):
        continue
    if not a.attributes.get('m_start', False):
        start = 0

    start = 0 if 'm_start' not in a.attributes else int(a.attributes['m_start'][0])
    length = int(a.attributes['m_stop'][0]) - start
    gene = a.attributes['match_id'][0]
    gid = a.attributes['ID'][0]

    if length < .5 * lengths[gene]:
        to_delete.update({x.attributes['ID'][0]: 1 for x in [a, *list(db.children(gid))]})
        #to_delete.extend([a, *list(db.children(gid))])
        #db.delete([a, *list(db.children(gid))])

#db.delete(to_delete)

fn = 'data1.csv'
with sqlite3.connect(DB) as connection:
    writer = csv.writer(open(fn, 'w'))
    c = connection.cursor()
    c.execute('SELECT * FROM features')
    rows = c.fetchall()
    for x in rows:
        writer.writerow(x)

reader = csv.reader(open(fn, 'r'))
writer = csv.writer(open('output1.csv', 'w'), delimiter='\t', lineterminator='\n')

next(reader)
for row in reader:
    new_row = eval(row[9])
    if new_row['ID'][0] in to_delete:
        continue
    new_str = []
    for k, v in new_row.items():
        if k != 'original_id' or k != 'original_name':
            new_str.append("{}={};".format(k, v[0]))
    new_str = "".join(new_str)
    writer.writerow([row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], new_str])

