#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
from builtins import *
import subprocess

from Bio import SeqIO
from Bio import AlignIO

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


def align_mafft(seqs, infile='tmp_unaligned.fasta', quiet=True):
    SeqIO.write(seqs, infile, 'fasta')
    if quiet:
        output = subprocess.check_output('mafft %s' % infile, shell=True, 
                                         stderr=subprocess.PIPE)
    else:
        output = subprocess.check_output('mafft %s' % infile, shell=True)    
    align = AlignIO.read(StringIO(output), 'fasta')
    return align


def select_ltrs(clust):
    # Get the longest exon. Everything to the left is considered 5', right is 3'
    exons = [m for m in clust.members if m.feature == 'exon']
    lexon = sorted(exons, key=lambda x:x.length())[-1]
    mid = clust.members.index(lexon)
    
    opt5 = [m for m in clust.members[:mid] if m.attr.get('geneRegion', '') == 'ltr']
    opt3 = [m for m in clust.members[mid:] if m.attr.get('geneRegion', '') == 'ltr']
    
    # Handle internal and one side
    if len(opt5) == 0 and len(opt3) == 0:
        # print('internal\t%s' % clust.attr['locus'])
        return 'internal', None, None
    elif len(opt5) == 0:
        # Use the rightmost
        return 'oneside_3', None, opt3[-1]
    elif len(opt3) == 0:
        # Use the leftmost
        return 'oneside_5', opt5[0], None
    
    # Handle one LTR for each side
    if len(opt5) == 1 and len(opt3) == 1:
        ltr5 = opt5[0]
        ltr3 = opt3[0]
        if ltr5.attr['repName'] == ltr3.attr['repName']:
            return 'prototype', ltr5, ltr3
        else:
            return 'prototype_mis', ltr5, ltr3
    else:
        # Handle multiple LTR
        allpairs = [(l5,l3) for l5 in opt5 for l3 in opt3]
        allpairs.sort(key=lambda t: t[0].score + t[1].score)
        allpairs.sort(key=lambda t: t[0].attr['repName'] == t[1].attr['repName'])
        ltr5 = allpairs[-1][0]
        ltr3 = allpairs[-1][1]
        if ltr5.attr['repName'] == ltr3.attr['repName']:
            return 'prototype', ltr5, ltr3
        else:
            return 'prototype_mis', ltr5, ltr3
    
    assert False
    return None, None, None
