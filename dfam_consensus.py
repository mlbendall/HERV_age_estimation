#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from Bio import SeqIO
import requests

""" Get all families for human """
url = "https://dfam.org/api/families"
params = {
    "format": "summary",
    "clade": "9606",
    "clade_relatives": "both",
}
response = requests.get(url, params=params)
results = response.json()["results"]

records = []
for r in results:
    if r['repeat_type_name'] == 'LTR':
        nurl = url + '/' + r['accession'] + '/sequence'
        response2 = requests.get(nurl, params={'format':'embl'})
        rec = SeqIO.read(StringIO(response2.text.encode('ascii','ignore')), 'embl')
        records.append(rec)

SeqIO.write(records, 'ERV_human.dfam.gb', 'genbank')
SeqIO.write(records, 'ERV_human.dfam.fasta', 'fasta')
