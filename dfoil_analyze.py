#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
dfoil_analyze: Given a dfoil output file, gives summary statistics to stdout
James B. Pease
http://www.github.com/jbpease/dfoil
"""

from __future__ import print_function
import sys
import argparse
from numpy import mean, percentile, var, std


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


def printlist(entry, delim="\t", ndigits=3):
    """Print a tab-separated list of items
        Arguments:
            entry: list of elements in line
            delim: field delimeter (default=\t)
    """
    decstring = "{0:." + str(ndigits) + "f}"
    newlist = []
    for elem in entry:
        try:
            if int(elem) == elem:
                newlist.append(int(elem))
            else:
                raise ValueError
        except ValueError as errstr:
            try:
                newlist.append(decstring.format(elem))
            except ValueError as errstr:
                newlist.append(str(elem))
    print(delim.join(["{}".format(x) for x in newlist]))


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="dfoil_analyze.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("infile", help="dfoil output file")
    parser.add_argument("--ndigits", type=int, default=3,
                        help="number of decimal places")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    dfo = [[], [], []]
    dil = [[], [], []]
    dfi = [[], [], []]
    dol = [[], [], []]
    headers = []
    introg_call = {}
    with open(args.infile) as infile:
        for line in infile:
            if line[0] == '#':
                headers = line[1:].rstrip().split("\t")
            else:
                entry = line.split()
                introg_call[entry[31]] = introg_call.get(entry[31], 0) + 1
                for name, dstat in (('DFO', dfo), ('DIL', dil),
                                    ('DFI', dfi), ('DOL', dol)):
                    chisq = float(entry[headers.index(name + '_chisq')])
                    dval = float(entry[headers.index(name + '_stat')])
                    pval = float(entry[headers.index(name + '_Pvalue')])
                    if chisq != 0:
                        dstat[0].append(dval)
                        dstat[1].append(chisq)
                        dstat[2].append(pval)
    statobj = (('DFO', dfo), ('DIL', dil), ('DFI', dfi), ('DOL', dol))
    statfields = ['D', 'chisq', 'Pvalue']
    print("\n# DFOIL component summary statistics:\n")
    for name, dstat in statobj:
        print("\t".join([
            "stat", "min", "mean", "max",
            "5%ile", "25%ile", "50%ile", "75%ile", "95%ile",
            "var", "sd  ", 'ab<0.01', '>0  ',
            '<0.05', 'cr0.5', 'cr0.05', 'cr0.01',
            ]))
        for i in range(len(dstat)):
            total = len(dstat[i])
            entry = (
                [(name if i == 0 else statfields[i]),
                 min(dstat[i]), mean(dstat[i]), max(dstat[i])] +
                [percentile(dstat[i], x) for x in (5, 25, 60, 75, 95)] +
                [dstat[i] and var(dstat[i]) or 0,
                 dstat[i] and std(dstat[i]) or 0,
                 float(sum([int(abs(x) <= 0.01) for x in dstat[i]])) / total,
                 float(sum([int(x > 0) for x in dstat[i]])) / total,
                 float(sum([int(x <= 0.05) for x in dstat[i]])) / total,
                 float(sum([int(x >= 0.46) for x in dstat[i]])) / total,
                 float(sum([int(x >= 3.84) for x in dstat[i]])) / total,
                 float(sum([int(x >= 6.64) for x in dstat[i]])) / total])
            printlist(entry, ndigits=args.ndigits)
    print("\n# Introgression Calls:\n")
    for (key, value) in iter(introg_call.items()):
        print(key, value)
    return ''


if __name__ == "__main__":
    main()
