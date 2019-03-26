#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

# Python 2 and 3:
from builtins import input

import sys

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import subprocess
from Bio import SeqIO
import re


rb_name = input('RepBase user name? ')
assert isinstance(rb_name, str)

rb_pw = input('RepBase password? ')
assert isinstance(rb_pw, str)

human_only = True

url = 'https://www.girinst.org/protected/repbase_extract.php'

if human_only:
    print('Taxon: Human only', file=sys.stderr)
    args = "division=Homo+sapiens&customdivision=&rank=&type=Endogenous+Retrovirus&autonomous=1&nonautonomous=1&simple=1&format=EMBL&full=true&sa=Download"
else:
    print('Taxon: ALL', file=sys.stderr)
    args = "division=&customdivision=&rank=&type=Endogenous+Retrovirus&autonomous=1&nonautonomous=1&simple=1&format=EMBL&full=true&sa=Download"

cmd = 'curl -u %s:%s "%s?%s"' % (rb_name, rb_pw, url, args)

text = subprocess.check_output(cmd, shell=True)
erec = text.split('\n//\n')

records = []
for t in erec:
    if re.match('^\s*$', t):
        continue
    t += '\n//\n'
    # Fix 1: Malformed DR line in EMBL file
    iscons = re.search('\nDR\s+\[\d+\] \(Consensus\)', t) is not None
    t = re.sub('\nDR\s+\[\d+\] \(Consensus\)', '\nDR   REPBASE; CONSENSUS', t)
    # Fix 2: End location must be greater than or equal to start location
    xt = t
    for s,e in re.findall('\nRP\s+(\d+)-(\d+)', xt):
        if int(s) > int(e):
            t = re.sub('\nRP\s+%s-%s' % (s,e), '\nRP   %s-%s' % (e,s), t)
    # Parse EMBL
    rec = SeqIO.read(StringIO(t), 'embl')
    newid = 'rn|%s|src|%s|' % (rec.name, 'con' if iscons else 'iso')
    if rec.id != '.':
        newid += 'id|%s|' % rec.id
    rec.id = newid
    records.append(rec)

if human_only:
    print('Wrote %d genbank records' % SeqIO.write(records, 'ERV_human.repbase.gb', 'genbank'), file=sys.stderr)
    print('Wrote %d fasta records' % SeqIO.write(records, 'ERV_human.repbase.fasta', 'fasta'), file=sys.stderr)
else:
    print('Wrote %d genbank records' % SeqIO.write(records, 'ERV_all.repbase.gb', 'genbank'), file=sys.stderr)
    print('Wrote %d fasta records' % SeqIO.write(records, 'ERV_all.repbase.fasta', 'fasta'), file=sys.stderr)
