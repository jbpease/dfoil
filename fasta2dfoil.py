#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
http://www.github.com/jbpease/dfoil

fasta2dfoil - fasta site counting script
@author: James B. Pease

If you use this software please cite:
Pease JB and MW Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny"
Systematic Biology. Online.
http://www.dx.doi.org/10.1093/sysbio/syv023

v.2015-06-11: Initial Version

This file is part of DFOIL.

DFOIL is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DFOIL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DFOIL.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
from io import open
import sys
import argparse
from itertools import groupby

def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name)
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
        

def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="""
    This script takes one or more FASTA files containing
    5 or 4 taxa and counts site patterns for use in DFOIL/Dstat analysis.
    To combine multiple FASTA files, the names of sequences must be
    identical in all files to be combined.     
    """)
    parser.add_argument("fastafile", nargs='*',
                        help="""one or more input fasta files""")
    parser.add_argument("--out", "-o", required=True,
                        help="""output count file""")
    parser.add_argument("--names", "-n", nargs='*', required=True,
                       help="""Order of the 5 (or 4) taxa, names must be 
                               consistent in all input files,
                               outgroup should be last""")
    args = parser.parse_args(args=arguments)
  
    position = 0
    NTAXA = len(args.names)
    if NTAXA is 4:
        headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
                   'BAAA', 'BABA', 'BBAA', 'BBBA']  
    elif NTAXA is 5:
        headers =  ['AAAAA', 'AAABA', 'AABAA', 'AABBA',
                    'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
                    'BAAAA', 'BAABA', 'BABAA', 'BABBA',
                    'BBAAA', 'BBABA', 'BBBAA', 'BBBBA']
    else:
        raise RuntimeError("Invalid number of taxa, use 5 or 4")
    with open(args.out, 'wb') as outfile:
        outfile.write("#chrom\tposition\t{}\n".format('\t'.join(headers)))
    for infilename in args.fastafile:
        site_count = {}
        seqs = dict(list(fasta_iter(infilename)))
        for i in range(len(seqs.values()[0])):
            site = [str(seqs[name][i]) for name in args.names]
            if len(set(site)) > 2:
                continue
            if set(site) - set('ATGC'):
                continue
            site_code = ''.join([x == site[-1] and 'A' or 'B'
                                 for x in site])
            site_count[site_code] = site_count.get(site_code, 0) + 1
                
        with open(args.out, 'ab') as outfile:
            outfile.write("FASTA\t{}\t{}".format(position, 
                          '\t'.join([str(site_count.get(x, 0)) 
                                     for x in headers])))
        position += 1
    return ''
    
if __name__ == "__main__":
    main()