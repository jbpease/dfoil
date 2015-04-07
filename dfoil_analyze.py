#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
http://www.github.com/jbpease/dfoil

dfoil_analyzed - summary statistics of dfoil output
@author: James B. Pease

Version: 2015-02-07 - Re-release on GitHub

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
import sys, argparse
from numpy import mean, percentile, var, std

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
        except ValueError:
            try:
                newlist.append(decstring.format(elem))
            except ValueError:
                newlist.append(str(elem))
    print(delim.join(["{}".format(x) for x in newlist]))


def main(arguments=sys.argv[1:]):
    """Main dfoil_analyze method"""
    parser = argparse.ArgumentParser(description="""
    Given a dfoil output file, gives summary statistics to stdout
    """)
    parser.add_argument("infile", help="dfoil output file")
    parser.add_argument("--ndigits", type=int, default=3,
                        help="number of decimal places")
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
        for i in xrange(len(dstat)):
            total = len(dstat[i])
            entry = (
                [(i == 0 and name or statfields[i]),
                 min(dstat[i]), mean(dstat[i]), max(dstat[i])]
                + [percentile(dstat[i], x) for x in (5, 25, 60, 75, 95)]
                + [dstat[i] and var(dstat[i]) or 0,
                   dstat[i] and std(dstat[i]) or 0,
                   float(sum([int(abs(x) <= 0.01) for x in dstat[i]])) / total,
                   float(sum([int(x > 0) for x in dstat[i]])) / total,
                   float(sum([int(x <= 0.05) for x in dstat[i]])) / total,
                   float(sum([int(x >= 0.46) for x in dstat[i]])) / total,
                   float(sum([int(x >= 3.84) for x in dstat[i]])) / total,
                   float(sum([int(x >= 6.64) for x in dstat[i]])) / total])
            printlist(entry, ndigits=args.ndigits)
    print ("\n# Introgression Calls:\n")
    for key, value in introg_call.iteritems():
        print(key, value)
    return ''


if __name__ == "__main__":
    main()
