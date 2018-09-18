#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
James B. Pease
http://www.github.com/jbpease/dfoil

fasta2dfoil -
This script takes one or more FASTA files containing
5 or 4 taxa and counts site patterns for use in DFOIL/Dstat analysis.
To combine multiple FASTA files, each file should be sequences
from one locus (i.e., one entry in the final table) and
the names of sequences must be identical in all files.

Example Usage:

python3 fasta2dfoil.py INPUT.fasta --out OUTPUT.fasta \
        --names TAXA1,TAXA2,TAXA3,TAXA4

"""

from __future__ import print_function
import sys
import argparse
from itertools import groupby

_LICENSE = """
If you use this software please cite:
Pease JB and MW Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny"
Systematic Biology. 64 (4): 651-662.
http://www.dx.doi.org/10.1093/sysbio/syv023

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


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="fasta2dfoil.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("fastafile", nargs='*',
                        help="""one or more input fasta
                                files for each locus""")
    parser.add_argument("--out", "-o", required=True,
                        help="""output count file, one entry per fasta""")
    parser.add_argument("--names", "-n", nargs=1, required=True,
                        help="""Order of the 5 (or 4) taxa separated by commas.
                                Names must be  consistent in all input files,
                                outgroup should be last""")
    parser.add_argument("--version", action="version", version="2017-11-07",
                        help="display version information and quit")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    position = 0
    taxa = args.names[0].split(',')
    ntaxa = len(taxa)
    if ntaxa is 4:
        headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
                   'BAAA', 'BABA', 'BBAA', 'BBBA']
    elif ntaxa is 5:
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
        if list(sorted(seqs.keys())) != list(sorted(taxa)):
            raise RuntimeError(
                "Error: Labels from {} ({}) do not match --names ({})".format(
                    infilename, seqs.keys(), taxa))
        for label in seqs:
            if len(seqs[label]) != len(seqs[taxa[0]]):
                raise RuntimeError(
                    "Error: Sequences in {} are of unequal length ({})".format(
                        infilename, ",".join(["{}={}".format(
                            x, len(seqs[x])) for x in seqs])))
        for i in range(len(list(seqs.values())[0])):
            site = [str(seqs[name][i]).upper() for name in taxa]
            if len(set(site)) > 2:
                continue
            if set(site) - set('ATGC'):
                continue
            site_code = ''.join(['A' if x == site[-1].upper() else 'B'
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
