#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
http://www.github.com/jbpease/dfoil

pre-dfoil - site count checking testt
@author: James B. Pease

If you use this software please cite:
Pease JB and MW Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny"
Systematic Biology. 64 (4): 651-662.
http://www.dx.doi.org/10.1093/sysbio/syv023

@version 2015-11-13 - Change to input specification, allow characters for zeros

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

from __future__ import print_function, unicode_literals
import sys
import argparse
from warnings import warn
from scipy.stats import chi2


SIGNCODES = {'dfoil': {'+++0': 2, '--0+': 3, '++-0': 4, '--0-': 5,
                       '+0++': 6, '-0++': 7, '0+--': 8, '0---': 9,
                       '++00': 10, '--00': 11},
             'partitioned': {'-00': 2, '0-0': 3, '+00': 4, '0+0': 5,
                             '-0+': 6, '0-+': 7, '+0-': 8, '0+-': 9},
             'dstat': {'+': 2, '-': 3}}
STATNAMES = {'dfoil': ['DFO', 'DIL', 'DFI', 'DOL'],
             'partitioned': ['D1', 'D2', 'D12'],
             'dstat': ['D']}

DIVNAMES = {'dfoil': ['T12', 'T34', 'T1234'],
            'dstat': ['T12', 'T123']}
DIVNAMES['partitioned'] = DIVNAMES['dfoil']

INTROGPATTERNS = {'dfoil': ['na', 'none', '13', '14', '23', '24',
                            '31', '41', '32', '42', '123', '124'],
                  'dstat': ['na', 'none', '23', '13']}
INTROGPATTERNS['partitioned'] = INTROGPATTERNS['dfoil']

INTROGLABELS = {'dfoil': ['N/A', 'None',
                          '1$\\Rightarrow$3', '1$\\Rightarrow$4',
                          '2$\\Rightarrow$3', '2$\\Rightarrow$4',
                          '3$\\Rightarrow$1', '4$\\Rightarrow$1',
                          '3$\\Rightarrow$2', '4$\\Rightarrow$2',
                          '1+2$\\Leftrightarrow$3', '1+2$\\Leftrightarrow$4'],
                'dstat': ['N/A', 'None',
                          '2$\\Leftrightarrow$3', '1$\\Leftrightarrow$3']}
INTROGLABELS['partitioned'] = INTROGPATTERNS['dfoil']


SITECODES = dict([(x, "{}{}".format('0' * (7 - len(bin(x))),
                  str(bin(x))[2:].replace('0', 'A').replace('1', 'B')))
                  for x in range(0, 32, 2)])


class DataWindow(object):
    """Basic Handler of each data entry"""
    def __init__(self, counts=None, meta=None, stats=None):
        self.counts = counts or {}
        self.meta = meta or {}
        self.stats = stats or {}

    def dcalc(self, mincount=0):
        """Calculate D-statistics
            Arguments:
                mincount: minimum total count to calculate P-values
        """
        (beta0, beta1, beta2) = self.meta['beta']
        if self.meta['mode'] == 'dfoil':
            self.stats['DFO'] = dcrunch(
                (self.counts[2] * beta0 + self.counts[10] * beta1 +
                 self.counts[20] * beta1 + self.counts[28] * beta2),
                (self.counts[4] * beta0 + self.counts[12] * beta1 +
                 self.counts[18] * beta1 + self.counts[26] * beta2),
                mincount=mincount)
            self.stats['DIL'] = dcrunch(
                (self.counts[2] * beta0 + self.counts[12] * beta1 +
                 self.counts[18] * beta1 + self.counts[28] * beta2),
                (self.counts[4] * beta0 + self.counts[10] * beta1 +
                 self.counts[20] * beta1 + self.counts[26] * beta2),
                mincount=mincount)
            self.stats['DFI'] = dcrunch(
                (self.counts[8] * beta0 + self.counts[10] * beta1 +
                 self.counts[20] * beta1 + self.counts[22] * beta2),
                (self.counts[16] * beta0 + self.counts[12] * beta1 +
                 self.counts[18] * beta1 + self.counts[14] * beta2),
                mincount=mincount)
            self.stats['DOL'] = dcrunch(
                (self.counts[8] * beta0 + self.counts[12] * beta1 +
                 self.counts[18] * beta1 + self.counts[22] * beta2),
                (self.counts[16] * beta0 + self.counts[10] * beta1 +
                 self.counts[20] * beta1 + self.counts[14] * beta2),
                mincount=mincount)
            self.stats['Dtotal'] = (
                (sum([self.counts[x] for x in (2, 4, 8, 16)]) * beta0) +
                (sum([self.counts[x] for x in (10, 12, 18, 20)]) * beta1) +
                (sum([self.counts[x] for x in (14, 22, 26, 28)]) * beta2))
            self.stats['Tvalues'] = self.calculate_5taxon_tvalues()
        elif self.meta['mode'] == "partitioned":
            self.stats['D1'] = dcrunch(self.counts[12], self.counts[20],
                                       mincount=mincount)
            self.stats['D2'] = dcrunch(self.counts[10], self.counts[18],
                                       mincount=mincount)
            self.stats['D12'] = dcrunch(self.counts[14], self.counts[22],
                                        mincount=mincount)
            self.calculate_5taxon_tvalues()
            self.stats['Dtotal'] = sum([self.counts[x] for x in
                                        [10, 12, 14, 18, 20, 22]])
        elif self.meta['mode'] == 'dstat':
            self.stats['D'] = dcrunch(
                self.counts[6] * beta1 + self.counts[8] * beta0,
                self.counts[10] * beta1 + self.counts[4] * beta0,
                mincount=mincount)
            self.stats['Dtotal'] = (
                (sum([self.counts[x] for x in (4, 8)]) * beta0) +
                (sum([self.counts[x] for x in (6, 10)]) * beta1))
            self.calculate_4taxon_tvalues(self.counts)
        return ''

    def calculate_5taxon_tvalues(self, counts=None):
        """Calculate approximate divergence times for a five-taxon tree
            Arguments:
                counts: dict of site counts (int keys)
        """
        counts = counts or self.counts
        total = float(sum(counts.values()))
        if not total:
            return {'T34': 0.0, 'T12': 0.0, 'T1234': 0.0}
        self.stats['T34'] = float(counts.get(2, 0) +
                                  counts.get(4, 0)) / (2 * total)
        self.stats['T12'] = float(counts.get(8, 0) +
                                  counts.get(16, 0)) / (2 * total)
        self.stats['T1234'] = 0.5 * (((float(counts.get(24, 0)) / total) +
                                      self.stats['T12']) +
                                     ((float(counts.get(6, 0)) / total) +
                                      self.stats['T34']))
        return ''

    def calculate_4taxon_tvalues(self, counts=None):
        """Calculate approximate divergence times for a four-taxon tree
            Arguments:
                counts: dict of site counts (int keys)
        """
        counts = counts or self.counts
        total = float(sum(counts.values()))
        if not total:
            return {'T12': 0.0, 'T123': 0.0}
        self.stats['T12'] = (float(counts.get(8, 0) + counts.get(4, 0)) /
                             (2.0 * total))
        self.stats['T123'] = 0.5 * (((float(counts.get(12, 0)) / total) +
                                     self.stats['T12']) +
                                    float(counts.get(2, 0)) / total)
        return ''

    def calc_signature(self, pvalue_cutoffs=None):
        """Determine D/DFOIL signature
            Arguments:
                pvalue_cutoffs = list of 1 or 2 P-value cutoffss
        """
        mode = self.meta['mode']
        if mode == "dfoil":
            pvalues = [pvalue_cutoffs[0], pvalue_cutoffs[0],
                       pvalue_cutoffs[1], pvalue_cutoffs[1]]
        elif mode == "partitioned":
            pvalues = [pvalue_cutoffs[0], pvalue_cutoffs[0],
                       pvalue_cutoffs[1]]
        elif mode in ["dstat"]:
            pvalues = [pvalue_cutoffs[0]]
        dfoil_signature = []
        for j, statname in enumerate(STATNAMES[mode]):
            if self.stats[statname]['Pvalue'] == 1:
                self.stats['signature'] = 0
                return ''
            if self.stats[statname]['Pvalue'] <= pvalues[j]:
                if self.stats[statname]['D'] > 0:
                    dfoil_signature.append('+')
                else:
                    dfoil_signature.append('-')
            else:
                dfoil_signature.append('0')
        self.stats['signature'] = SIGNCODES[mode].get(''.join(dfoil_signature),
                                                      1)
        return ''


def dcrunch(left_term, right_term, mincount=0):
    """Calculate D-statistic
        Arguments:
            left_term: left term of D
            right_term: right term of D
            mincount: minimum total to calculate P-values, otherwise P=1
    """
    result = {}
    result['left'] = left_term
    result['right'] = right_term
    result['Dtotal'] = left_term + right_term
    if not left_term + right_term:
        result['Pvalue'] = 1.0
        result['chisq'] = 0
        result['D'] = 0
    elif left_term + right_term < mincount:
        result['chisq'] = 0
        result['Pvalue'] = 1.0
        result['D'] = (float(left_term - right_term) /
                       (left_term + right_term))
    else:
        (val, pval) = chi2_test(left_term, right_term)
        result['chisq'] = val
        result['Pvalue'] = pval
        result['D'] = (float(left_term - right_term) /
                       (left_term + right_term))

    return result


def chi2_test(val0, val1):
    """Calculate Pearson Chi-Squared for the special case of
       two values that are expected to be equal
       Arguments:
           val0: first value
           val1: second value
    """
    try:
        chisq = float((val0 - val1)**2) / float(val0 + val1)
        if not chisq:
            return (0, 1)
        pval = 1.0 - chi2.cdf(chisq, 1)
        return (chisq, pval)
    except ZeroDivisionError as errstr:
        return (0, 1)


def make_header(mode):
    """Create Column Headers for Various Modes
        Arguments:
            mode: dfoil statistical mode
    """

    return ("{}\n".format('\t'.join(
        ['#chrom', 'coord', 'total', 'dtotal'] +
        DIVNAMES[mode] +
        ['{}_{}'.format(x, y)
         for x in STATNAMES[mode]
         for y in ('left', 'right', 'total',
                   'stat', 'chisq', 'Pvalue')] +
        ['introgression'] +
        ['introg{}'.format(x) for x in INTROGPATTERNS[mode]]))
        ).encode('utf-8')


def main(arguments=sys.argv[1:]):
    """Main pre-dfoil method"""
    parser = argparse.ArgumentParser(description=("""
    Check your site counts before running DFOIL to ensure you've
    specified the right taxon order and your counts indicate a tree that
    is correct and does not violate assumptions of the DFOIL model
    USAGE: pre-dfoil.py INPUTFILE1 INPUTFILE2 ... """))
    parser.add_argument('--infile', help="input tab-separated counts file",
                        nargs='*', required=True)
    parser.add_argument('--out', help="outputs tab-separated DFOIL stats",
                        nargs='*')
    parser.add_argument('--mincount', type=int, default=10,
                        help="minium number of D denominator sites per window")
    parser.add_argument("--mintotal", type=int, default=50,
                        help="minimum total number of sites in a region")
    parser.add_argument('--pvalue', type=float, default=[0.01, 0.01],
                        nargs='*',
                        help="""minimum P-value cutoff for regions,
                                can specify one P-value for all four tests
                                or two separate ones for DFO/DIL and DFI/DOL
                                (or D1/D2 and D12 for 'partitioned')""")
    parser.add_argument('--mode', default="dfoil",
                        choices=["dfoil", "dfoilalt", "partitioned"],
                        help="""dfoil = DFOIL,
                                dfoilalt = DFOIL without single-B patterns,
                                partitioned = Partitioned D-statistics""")
    parser.add_argument("--beta1", type=float,
                        help="""beta1 coefficient for single-B patterns,
                                 defaults: DFOIL/Dstatalt=1.0,
                                 DFOILalt,Dstat=0,Dpart=N.A.""")
    parser.add_argument("--beta2", type=float,
                        help="""beta2 coefficient for double-B patterns,
                                defaults: Dpart=N.A., others=1.0""")
    parser.add_argument("--beta3", type=float,
                        help="""beta3 coefficient for triple-B patterns
                                defaults: DFOIL/DFOILalt=1.0,
                                Dstat/Dpart=N.A.""")
    parser.add_argument("--zerochar", default=[".", "NA"], nargs='*',
                        help="""list of strings used in place of zeros
                                in the input file default is [".", "NA"]""")
    parser.add_argument("--version", action="store_true",
                        help="display version information and quit")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("2015-11-23")
        sys.exit()
    # ===== INITIALIZE =====
    if len(args.pvalue) == 1:
        args.pvalue = [args.pvalue[0], args.pvalue[0]]
    # Set beta parameters for presets
    if not args.beta1:
        args.beta1 = args.mode in ['dfoil', 'dstatalt'] and 1. or 0.
    if not args.beta2:
        args.beta2 = 1.
    if not args.beta3:
        args.beta3 = 1.
    if args.mode == 'dstatalt':
        args.mode = 'dstat'
    elif args.mode == 'dfoilalt':
        args.mode = 'dfoil'
    # ===== PARSE COUNT FILE  =====
    for ifile, infilename in enumerate(args.infile):
        window_data = []
        with open(infilename) as infile:
            for line in infile:
                if line[0] == '#':
                    continue
                try:
                    arr = line.rstrip().split()
                    window = DataWindow(meta=dict(
                        chrom=arr[0], position=int(arr[1]),
                        mode=args.mode,
                        beta=(args.beta1, args.beta2, args.beta3)))
                    if args.mode in ["dfoil", "partitioned"]:
                        window.counts = dict([
                            (j - 2) * 2,
                            arr[j] not in args.zerochar and int(arr[j]) or 0]
                                             for j in range(2, 18))
                    elif args.mode == 'dstat':
                        window.counts = dict([
                            (j - 2) * 2,
                            arr[j] not in args.zerochar and int(arr[j]) or 0]
                                             for j in range(2, 9))
                    if sum(window.counts.values()) < args.mintotal:
                        continue
                    window.meta['total'] = sum(window.counts.values())
        #            window.dcalc(mincount=args.mincount)
        #            window.calc_signature(pvalue_cutoffs=args.pvalue)
                    window_data.append(window)
                except:
                    warn(
                        "line invalid, skipping...\n{}".format(line))
                    continue
            # ===== ANALYZE WINDOWS AND DETERMINE RUNS ====-
        # ===== WRITE TO OUTPUT =====
#        with open(args.out[ifile], 'wb') as outfile:
#           outfile.write(make_header(args.mode))
        sum_data = {}
        for window in window_data:
            for code in window.counts:
                sum_data[code] = sum_data.get(
                    code, 0) + window.counts[code]
        # Check1
        checkok = True
        print("="*79)
        print("""\
Checking that concordant patterns are more common that discordant
(Note that this is normal when introgression is extreme, but generally
indicates the taxa are out of order):""")
        print("-"*79)
        for concode in (2, 4, 8, 16, 6, 24):
            for discode in (10, 12, 14, 18, 20, 22, 26, 28):
                if sum_data[concode] < sum_data[discode]:
                    print("""\
Total count of discordant pattern {}={} is higher than \
concordant pattern {}={}""".format(
                        SITECODES[discode], sum_data[discode],
                        SITECODES[concode], sum_data[concode]))
                    checkok = False
        if checkok:
            print("Pass")
        # Check2
        checkok = True
        print("="*79)
        print("""\
Checking that divergences are correctly ordered
(P1 and P2 should diverge AFTER P3 and P4)""")
        print("-"*79)
        for abcode in (8, 16):
            for cdcode in (2, 4):
                if sum_data[abcode] > sum_data[cdcode]:
                    print("""\
Total count of A/B terminal substitutions {}={} is higher than \
C/D terminal substitutions {}={}""".format(
                        SITECODES[abcode], sum_data[abcode],
                        SITECODES[cdcode], sum_data[cdcode]))
                    checkok = False
        if checkok:
            print("Pass")
        print("="*79)
        # Check3
        checkok = True
        print("="*79)
        print("Checking that terminal branch pairs are"
              "proportionate approximately")
        print("-"*79)
        abratio = float(sum_data[16]) / sum_data[8]
        print("BAAAA/ABAAA ratio = {}".format(abratio))
        if 0.8 < abratio > 1.25:
            checkok = False
            print("Warning: P1/P2 ratio is somewhat high")
        cdratio = float(sum_data[4]) / sum_data[2]
        print("AABAA/AAABA ratio = {}".format(cdratio))
        if 0.8 < cdratio > 1.25:
            checkok = False
            print("Warning: P3/P4 ratio is somewhat high")
        if checkok:
            print("Pass")
        print("="*79)

    return ''
if __name__ == "__main__":
    main()
