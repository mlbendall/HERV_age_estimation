# -*- coding: utf-8 -*-

def zero(v):
    return abs(v) if v == -0.0 else v
 
def trim_end_gaps(s1,s2):
    """
    Trim to first position where both sequences have non-gap
    """
    for spos in range(len(s1)):
        if s1[spos] != '-' and s2[spos] != '-':
            break
    for epos in range(len(s1), 0, -1):
        if s1[epos-1] != '-' and s2[epos-1] != '-':
            break
    return s1[spos:epos], s2[spos:epos]


def estimate_nucleotide_frequencies(seq):
    seq = seq.replace('-','').upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]

def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
    #for (x,y) in zip(seq1,seq2):
    for (x,y) in pairs:
        if x != y:
            p += 1
    #length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    return zero(float(p) / length)

def JCdistance(seq1, seq2):
    """ 
    distance = -b log(1 - p / b)
    where:
    b = 3/4
    and p = p-distance, i.e. uncorrected distance between seq1 and seq2
    """
    from math import log
    b = 0.75
    p = pdistance(seq1,seq2)
    try: d = -b * log( 1-p/b )
    except ValueError: 
        print "Tried to take log of a negative number"
        return None
    return zero(d)

def TNdistance(seq1, seq2):
    """ 
    Tajima-Nei distance = -b log(1 - p / b)
    where:
    b = 0.5 * [ 1 - Sum i from A to T(Gi^2+p^2/h) ]
    h = Sum i from A to G( Sum j from C to T (Xij^2/2*Gi*Gj))
    p = p-distance, i.e. uncorrected distance between seq1 and seq2
    Xij = frequency of pair (i,j) in seq1 and seq2, with gaps removed
    Gi = frequency of base i over seq1 and seq2 """
    from math import log

    ns = ['A','C','G','T']
    G = estimate_nucleotide_frequencies(seq1 + seq2)
    p = pdistance(seq1,seq2)
    pairs = []
    h = 0

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
       
    #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns)-1):
        for j in range(i+1,len(ns)):
            if i != j: paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
            Xij_sq = (float(paircount)/len(pairs))**2
            GiGj = G[i]*G[j]
            h += 0.5*Xij_sq/GiGj  #h value used to calculate b
    try:
        b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
        return zero(-b * log(1 - p/b))
    except ZeroDivisionError:
        assert h == 0.0
        return 0.0
    except ValueError:
        assert (1 - p/b) < 0.0
        print "Tried to take log of a negative number"
        return None
#     
#     try:
#         d = -b * log(1 - p/b)
#     except ValueError:
#         print "Tried to take log of a negative number"
#         return None
#     return d

def K2Pdistance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError: 
        print "Tried to take log of a negative number"
        return None
    return zero(d)

def Tamuradistance(seq1,seq2):
    """
    Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
    where:
    P = transition frequency
    Q = transversion frequency
    C = GC1 + GC2 - 2 * GC1 * GC2
    GC1 = GC-content of sequence 1
    GC2 = GC-coontent of sequence 2
    """
    from math import log
    pairs = []
    
    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    gc1 = sum(estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(estimate_nucleotide_frequencies(seq2)[1:3])
    c = gc1 + gc2 - 2 * gc1 * gc2

    try: d = -c * log( 1 - p/c - q) - 0.5 * ( 1 - c ) * log ( 1 - 2*q )
    except ValueError: 
        print "Tried to take log of a negative number"
        return None
    return zero(d)


def Coffindistance(seq1,seq2):
    """
    Indels longer than 2 bases are treated as a single substitution
    """
    s1, s2 = trim_end_gaps(seq1.upper(), seq2.upper())
    
    match = mismatch = opengap = 0
    ingap = False
    for c1,c2 in zip(s1, s2):
        if c1 == '-' or c2 == '-':
            if not ingap:
                opengap += 1
            ingap = True
        else:
            ingap = False
            if c1 == c2:
                match += 1
            else:
                mismatch += 1
    tl = match + mismatch + opengap   
    return zero(float(mismatch + opengap) / tl)

import re
def Coffindistance0(seq1,seq2):
    S1,S2 = seq1.upper(), seq2.upper()
    gspan1 = [m.span() for m in re.finditer('-+', S1)]
    g1 = len([t for t in gspan1 if t[0] != 0 and t[1] != len(S1)])
    gspan2 = [m.span() for m in re.finditer('-+', S2)]
    g2 = len([t for t in gspan2 if t[0] != 0 and t[1] != len(S2)])
    #collect ungapped pairs
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)

    diff = sum(c1!=c2 for c1,c2 in pairs)
    d = float(g1 + g2 + diff) / (g1 + g2 + len(pairs))
    return d            
