#! /usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
# from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
from builtins import *

import sys
# import os
import argparse

from collections import OrderedDict, defaultdict, Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet, IUPAC

from utils import select_ltrs
from telebuilder.utils import gtfutils
import subprocess
from dist import (pdistance, JCdistance, TNdistance, 
                  K2Pdistance, Tamuradistance, Coffindistance)
from utils import align_mafft

DISTS = ['p.dist', 'jc.dist', 'tn.dist', 'k2p.dist', 't3p.dist', 'cof.dist',]

def pairwise_distance(seqs):
    assert len(seqs) == 2
    align = align_mafft(seqs)
    s1, s2 = (str(align[0].seq).upper(), str(align[1].seq).upper())
    return {
        'p.dist': pdistance(s1,s2),
        'jc.dist': JCdistance(s1,s2),
        'tn.dist': TNdistance(s1,s2),
        'k2p.dist': K2Pdistance(s1,s2),
        't3p.dist': Tamuradistance(s1,s2),
        'cof.dist': Coffindistance(s1,s2),
    }

def pairwise_distance_R(seqs):
    assert len(seqs) == 2
    align = align_mafft(seqs)
    AlignIO.write(align, 'tmp_aligned.fasta', 'fasta')
    out = subprocess.check_output('Rscript calcdist.R tmp_aligned.fasta', shell=True, stderr=subprocess.PIPE)
    dists = out.strip('\n').split('\t')
    s1, s2 = (str(align[0].seq).upper(), str(align[1].seq).upper())
    return {
        'p.dist': dists[0],
        'jc.dist': dists[1],
        'tn.dist': dists[2],
        'k2p.dist': dists[3],
        't3p.dist': dists[4],
        'cof.dist': Coffindistance(s1,s2),
    }


def age_estimates(gtffile, seqfile, consfile, outfile, pydist):
    paired_ltrs = OrderedDict()      # Paired LTRs from prototype
    ltr_subtypes = defaultdict(list) # Single LTRs by subfamily
    subfam = {}                      # {locus -> subfamily}
    category = {}                    # {locus -> category}

    seqs = OrderedDict((s.id, s) for s in SeqIO.parse(seqfile, 'fasta'))
    cons = OrderedDict((s.id, s) for s in SeqIO.parse(consfile, 'fasta'))

    # Iterate over loci and identify best LTRs
    header = ['locus', 'category', 'subfamily',] + DISTS
    distances = []
    clusters = gtfutils.read_gtf_clusters(gtffile)
    for clust in clusters:
        loc = clust.attr['locus']
        cat, ltr5, ltr3 = select_ltrs(clust)
        if ltr5 is None and ltr3 is None:        
            distances.append([loc, cat, 'NA',] + [float('inf') for _ in DISTS])
            continue
        repnames = set()
        if ltr5 is not None:
            rec5 = SeqRecord(
                Seq(str(seqs[loc].seq[(ltr5.start-1):(ltr5.end-1)]).upper(), SingleLetterAlphabet),
                id = loc + '[5]',
                description = ltr5.attr['repName']
            )
            ltr_subtypes[ltr5.attr['repName']].append(rec5)
            repnames.add(ltr5.attr['repName'])
        if ltr3 is not None:
            rec3 = SeqRecord(
                Seq(str(seqs[loc].seq[(ltr3.start-1):(ltr3.end-1)]).upper(), SingleLetterAlphabet),
                id = loc + '[3]',
                description = ltr3.attr['repName']
            )
            ltr_subtypes[ltr3.attr['repName']].append(rec3)
            repnames.add(ltr3.attr['repName'])
        
        if ltr5 is not None and ltr3 is not None:
            paired_ltrs[loc] = (rec5, rec3)        
            best_ltr = rec5 if len(rec5) > len(rec3) else rec3
            # print('%d %d -> %d' % (len(rec5), len(rec3), len(best_ltr)))
        elif ltr5 is not None:
            assert ltr3 is None
            best_ltr = rec5
        elif ltr3 is not None:
            assert ltr5 is None
            best_ltr = rec3
        else:
            assert False
        
        subfam = best_ltr.description
        assert subfam in cons
        curcons = cons[subfam]
        if pydist:
            d = pairwise_distance([curcons, best_ltr])
        else:
            d = pairwise_distance_R([curcons, best_ltr])
        # distances.append((loc, cat, subfam, d['k2p.dist']))
        distances.append([loc, cat, subfam,] + [d[_] for _ in DISTS])
        category[loc] = cat
    
    print('\t'.join(header), file=outfile)
    for row in distances:
        print('\t'.join(map(str, row)), file=outfile)

def console():
    parser = argparse.ArgumentParser(
        description='''Age estimation'''
    )
    parser.add_argument('--pydist', action='store_true',
                        help='Use distances calculated with python functions. Default is to use the "ape" package in R'
                        )
    parser.add_argument('gtf', help='Extracted GTF file')    
    parser.add_argument('seqs', help='Extracted sequence file')    
    parser.add_argument('cons', help='Consensus sequences for each LTR')
    parser.add_argument('outfile',
                        nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output file. Default is stdout.')
    args = parser.parse_args()
    age_estimates(args.gtf, args.seqs, args.cons, args.outfile, args.pydist)


if __name__ == '__main__':
    console()
    

# age_estimates('HML2/HML2_extracted.gtf', 'HML2/HML2.fna', 'ERV_human.consensus.fasta', sys.stdout)





