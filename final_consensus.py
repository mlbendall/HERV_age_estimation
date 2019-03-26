#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
from Bio import SeqIO
import re
from collections import OrderedDict

# Load dfam seqs
ddict = {s.id.split('.')[0]:s for s in SeqIO.parse('ERV_human.dfam.gb', 'genbank')}
print('Loaded %d Dfam records' % len(ddict), file=sys.stderr)

# Load repbase seqs
rdict = {s.id.split('|')[1]:s for s in SeqIO.parse('ERV_human.repbase.fasta', 'fasta')}
print('Loaded %d RepBase ids' % len(rdict), file=sys.stderr)

# Load mapping
mapping = [l.strip('\n').split('\t') for l in open('rmsk_to_dfam.tsv', 'rU')]

conseqs = []
for row in mapping:
    if row[2]:
        assert row[2] in ddict, '%s not found in Dfam' % row[2]
        ddict[row[2]].description = '%s %s' % (ddict[row[2]].id, ddict[row[2]].description)
        ddict[row[2]].id = row[0]
        conseqs.append(ddict[row[2]])
    else:
        assert row[0] in rdict, '%s not found in RepBase' % row[0]
        rdict[row[0]].description = '%s %s' % (rdict[row[0]].id, rdict[row[0]].description)
        rdict[row[0]].id = row[0]
        conseqs.append(rdict[row[0]])

# Check sets are the same
assert set([m[0] for m in mapping]) == set([s.id for s in conseqs])

print('%d consensus sequences' % SeqIO.write(conseqs, 'ERV_human.consensus.fasta', 'fasta'), file=sys.stderr)
