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
v.2015-11-21: Fixes for Python2to3

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
import sys
import argparse
from itertools import groupby


def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def main(arguments=None):
    if arguments is None:
        arguments = sys.argv[1:]
    parser = argparse.ArgumentParser(description="""
    This script takes one or more FASTA files containing
    5 or 4 taxa and counts site patterns for use in DFOIL/Dstat analysis.
    To combine multiple FASTA files, each file should be sequences
    from one locus (i.e., one entry in the final table) and
    the names of sequences must be identical in all files.
    """)
    parser.add_argument("fastafile", nargs='*',
                        help="""one or more input fasta
                                files for each locus""")
    parser.add_argument("--out", "-o", required=True,
                        help="""output count file, one entry per fasta""")
    parser.add_argument("--names", "-n", nargs='*', required=True,
                        help="""Order of the 5 (or 4) taxa,
                                names must be
                                consistent in all input files,
                                outgroup should be last""")
    parser.add_argument("--version", action="store_true",
                        help="""display version info""")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("fasta2dfoil: v.2015-11-21")
        sys.exit()
    position = 0
    NTAXA = len(args.names)
    if NTAXA is 4:
        headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
                   'BAAA', 'BABA', 'BBAA', 'BBBA']
    elif NTAXA is 5:
        headers = ['AAAAA', 'AAABA', 'AABAA', 'AABBA',
                   'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
                   'BAAAA', 'BAABA', 'BABAA', 'BABBA',
                   'BBAAA', 'BBABA', 'BBBAA', 'BBBBA']
    else:
        raise RuntimeError("Invalid number of taxa, use 5 or 4")
    with open(args.out, 'w') as outfile:
        outfile.write("#chrom\tposition\t{}\n".format('\t'.join(headers)))
    for infilename in args.fastafile:
        site_count = {}
        seqs = {}
        for header, seq in fasta_iter(infilename):
            seqs[header] = seq
        if list(sorted(seqs.keys())) != list(sorted(args.names)):
            raise RuntimeError(
                "Error: Labels from {} ({}) do not match --names ({})".format(
                    infilename, seqs.keys(), args.names))
        for label in seqs:
            if len(seqs[label]) != len(seqs[args.names[0]]):
                raise RuntimeError(
                    "Error: Sequences in {} are of unequal length ({})".format(
                        infilename, ",".join(["{}={}".format(
                            x, len(seqs[x])) for x in seqs])))
        for i in range(len(list(seqs.values())[0])):
            site = [str(seqs[name][i]) for name in args.names]
            if len(set(site)) > 2:
                continue
            if set(site) - set('ATGC'):
                continue
            site_code = ''.join([x == site[-1] and 'A' or 'B'
                                 for x in site])
            site_count[site_code] = site_count.get(site_code, 0) + 1

        with open(args.out, 'a') as outfile:
            outfile.write("{}\t{}\t{}\n".format(infilename, position,
                          '\t'.join([str(site_count.get(x, 0))
                                     for x in headers])))
        position += 1
    return ''

if __name__ == "__main__":
    main()
