#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
dfoil - Calculate DFOIL and D-statistics stats from one or more count files.
James B. Pease
http://www.github.com/jbpease/dfoil

USAGE: dfoil.py INPUTFILE1 ... --out OUTPUTFILE1 ...
"""

from __future__ import print_function, unicode_literals
import sys
import argparse
from warnings import warn
# from numpy import mean
from scipy.stats import chi2
from Bio import Phylo
# import matplotlib
from itertools import combinations
# from precheck import pre_check
from fasta2dfoil import fasta_iter

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


PLOTFORMATS = ("eps", "jpeg", "jpg", "pdf", "pgf", "png",
               "ps", "raw", "rgba", "svg", "svgz", "tif", "tiff")


class DfoilTree(object):

    def __init__(self, quintets=None, treepath=None):
        self.quintets = quintets if quintets is not None else {}
        self.tree = Phylo.read(treepath, 'newick')
        self.parents = {}
        self._process_tree(self.tree)

    def _process_tree(self, tree):
        # nodes = [x for x in tree.get_terminals()]
        for node in tree.get_nonterminals():
            if node.count_terminals() < 5:
                continue
            if node[0].count_terminals() < 4 and node[1].count_terminals() < 4:
                continue
            # clades0 = [x.name for x in node[0].get_terminals()]
            # clades1 = [x.name for x in node[1].get_terminals()]
            for (i, j) in [(0, 1), (1, 0)]:
                if (node[i].count_terminals() >= 4 and
                    node[i][0].count_terminals() > 1 and
                        node[i][1].count_terminals() > 1):
                    for outgroup in node[j].get_terminals():
                        for first_taxa_pair in combinations(
                                node[i][0].get_terminals(), 2):
                            for second_taxa_pair in combinations(
                                    node[i][1].get_terminals(), 2):
                                if (node.distance(first_taxa_pair[0]) +
                                    node.distance(first_taxa_pair[1])) > (
                                        node.distance(second_taxa_pair[0]) +
                                        node.distance(second_taxa_pair[1])):
                                    quintet = tuple([
                                        second_taxa_pair[0],
                                        second_taxa_pair[1],
                                        first_taxa_pair[0],
                                        first_taxa_pair[1],
                                        outgroup])
                                else:
                                    quintet = tuple([
                                        first_taxa_pair[0],
                                        first_taxa_pair[1],
                                        second_taxa_pair[0],
                                        second_taxa_pair[1],
                                        outgroup])
                                self.quintets[quintet] = {}
            return ""


class DataWindow(object):
    """Basic Handler of each data entry"""
    def __init__(self, counts=None, meta=None, stats=None):
        self.counts = counts or {}
        self.meta = meta or {}
        self.stats = stats or {}

    def _getcount(self, bits):
        try:
            return sum([self.counts.get(x, 0) for x in bits])
        except TypeError as te:
            return self.counts.get(bits, 0)

    def dcalc(self, mincount=0):
        """Calculate D-statistics
            Arguments:
                mincount: minimuim total count to calculate P-values
        """
        (beta0, beta1, beta2) = self.meta['beta']
        if self.meta['mode'] == 'dfoil':
            self.stats['DFO'] = dcrunch(
                (self._getcount(2) * beta0 +
                 self._getcount((10, 20)) * beta1 +
                 self._getcount(28) * beta2),
                (self._getcount(4) * beta0 +
                 self._getcount((12, 18)) * beta1 +
                 self._getcount(26) * beta2),
                mincount=mincount)
            self.stats['DIL'] = dcrunch(
                (self._getcount(2) * beta0 +
                 self._getcount((12, 18)) * beta1 +
                 self._getcount(28) * beta2),
                (self._getcount(4) * beta0 +
                 self._getcount((10, 20)) * beta1 +
                 self._getcount(26) * beta2),
                mincount=mincount)
            self.stats['DFI'] = dcrunch(
                (self._getcount(8) * beta0 +
                 self._getcount((10, 20)) * beta1 +
                 self._getcount(22) * beta2),
                (self._getcount(16) * beta0 +
                 self._getcount((12, 18)) * beta1 +
                 self._getcount(14) * beta2),
                mincount=mincount)
            self.stats['DOL'] = dcrunch(
                (self._getcount(8) * beta0 +
                 self._getcount((18, 22)) * beta1 +
                 self._getcount(22) * beta2),
                (self._getcount(16) * beta0 +
                 self._getcount((10, 20)) * beta1 +
                 self._getcount(14) * beta2),
                mincount=mincount)
            self.stats['Dtotal'] = (
                (self._getcount((2, 4, 8, 16)) * beta0) +
                (self._getcount((10, 12, 18, 20)) * beta1) +
                (self._getcount((14, 22, 26, 28)) * beta2))
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
    except ZeroDivisionError as exc:
        return (0, 1)


def count_patterns(sequences):
    patterns = {}
    for i in range(len(sequences)):
        alleles = [x[i] for x in sequences]
        if len(set(alleles)) != 2:
            continue
        if any(x not in 'ATGC' for x in alleles):
            continue
        pattern = sum([2 ** (4 - i) for i, x in enumerate(alleles[:4]) if x !=
                       alleles[-1]])
        patterns[pattern] = patterns.get(pattern, 0) + 1
    return patterns


def make_header(mode):
    """Create Column Headers for Various Modes
        Arguments:
            mode: dfoil statistical mode
    """

    return ("{}\n".format('\t'.join(
        ['#P1', 'P2', 'P3', 'P4', 'PO', 'total', 'dtotal'] +
        DIVNAMES[mode] +
        ['{}_{}'.format(x, y)
         for x in STATNAMES[mode]
         for y in ('left', 'right', 'total',
                   'stat', 'chisq', 'Pvalue')] +
        ['introgression']))
        # ['introg{}'.format(x) for x in INTROGPATTERNS[mode]])
           ).encode('utf-8')


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="dfoil.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument('--fasta', help="alignment in FASTA format")
    parser.add_argument('--tree', help="tree in newick format")
    parser.add_argument('--out', help="outputs tab-separated DFOIL stats",
                        required=True)
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
    parser.add_argument('--runlength', type=int, default=0,
                        help="""if two introgressing windows are separated
                                by this many windows of non-introgression
                                color in the intervening windows to
                                create a more continuous visual appearance""")
    parser.add_argument('--mode', default="dfoil",
                        choices=["dfoil", "dfoilalt", "partitioned",
                                 "dstat", "dstatalt"],
                        help="""dfoil = DFOIL,
                                dfoilalt = DFOIL without single-B patterns,
                                partitioned = Partitioned D-statistics,
                                dstat = Four-Taxon D-statistic,
                                dstatalt = Four-Taxon D-statistic
                                with single-B patterns""")
    parser.add_argument("--beta1", type=float,
                        help="""beta1 coefficient for single-B patterns,
                                 defaults: DFOIL/Dstatalt=1.0,
                                 DFOILalt,Dstat=0,Dpart=N.A.""")
    parser.add_argument("--beta2", type=float, default=1.0,
                        help="""beta2 coefficient for double-B patterns,
                                defaults: Dpart=N.A., others=1.0""")
    parser.add_argument("--beta3", type=float, default=1.0,
                        help="""beta3 coefficient for triple-B patterns
                                defaults: DFOIL/DFOILalt=1.0,
                                Dstat/Dpart=N.A.""")
    parser.add_argument("--zerochar", default=[".", "NA"], nargs='*',
                        help="""list of strings used in place of zeros
                                in the input file default is [".", "NA"]""")
    parser.add_argument("--plot", nargs='*',
                        help="""write plot to file path(s) given""")
    parser.add_argument("--plot_labels", nargs='*', help="taxon labels")
    parser.add_argument("--plot_color", default="color",
                        choices=["color", "colordark",
                                 "bw", "bwdark"],
                        help="""choose color mode""")
    parser.add_argument("--plot_noanc", action="store_true",
                        help="""do not plot background for
                                ancestral introgression""")
    parser.add_argument("--plot_lineweight", type=float, default=1.,
                        help="line weight for dplots (default=1pt)")
    parser.add_argument("--plot_yscale", type=float, default=1.0,
                        help="Y-axis min-max value, default is 1")
    parser.add_argument("--plot_smooth", type=int,
                        help="average D-stats over this number of points")
    parser.add_argument("--plot_background", type=float, default=0.3,
                        help="""0-1.0 background intensity
                                (0=none, default=0.3)""")
    parser.add_argument("--plot_totals", action="store_true",
                        help="add a background plot of total site counts")
    parser.add_argument("--plot_hidekey", action="store_true",
                        help="hide plot key")
    parser.add_argument("--plot_hideaxes", action="store_true",
                        help="hide axes labels")
    parser.add_argument("--plot_height", type=float, default=8.,
                        help="height of plot (in cm)")
    parser.add_argument("--plot_width", type=float, default=24.,
                        help="width of plot (in cm)")
    parser.add_argument("--pre-check-only", action="store_true",
                        help=("Only run the data pre-check "
                              "(formely pre-dfoil.py)"))
    parser.add_argument("--skip-pre-check", action="store_true",
                        help=("Skip running the data pre-check "
                              "(formely pre-dfoil)"))
    parser.add_argument("--version", action="version", version="2017-011-25",
                        help="display version information and quit")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    # ===== INITIALIZE =====
#    if set(args.infile) & set(args.out):
#        raise NameError("input and output file have same path")
#    if len(args.pvalue) == 1:
#        args.pvalue = [args.pvalue[0], args.pvalue[0]]
    # ===== Process Tree =====
    tree = DfoilTree(treepath=args.tree)
    sequences = dict((hdr, seq) for hdr, seq in fasta_iter(args.fasta))

    # Set beta parameters for presets
    if args.beta1 is None:
        args.beta1 = 1.0 if args.mode in ['dfoil', 'dstatalt'] else 0.
    if args.mode == 'dstatalt':
        args.mode = 'dstat'
    elif args.mode == 'dfoilalt':
        args.mode = 'dfoil'
    # ===== MAIN DFOIL CALC =========
    window_data = []
    for quintet in tree.quintets:
        try:
            window = DataWindow(meta=dict(
                chrom='0', position=0,
                quintet=[x.name for x in quintet],
                mode=args.mode,
                beta=(args.beta1, args.beta2, args.beta3)))
            window.counts = count_patterns(
                [sequences[x.name] for x in quintet])
            # print(quintet, window.counts)
            if sum(window.counts.values()) < args.mintotal:
                continue
            window.meta['total'] = sum(window.counts.values())
            window.dcalc(mincount=args.mincount)
            window.calc_signature(pvalue_cutoffs=args.pvalue)
            window_data.append(window)
        except Exception as exc:
            warn(
                "quintet invalid, skipping...\n{}".format(quintet))
            continue
# ===== WRITE TO OUTPUT =====
    with open(args.out, 'wb') as outfile:
        outfile.write(make_header(args.mode))
        for window in window_data:
            entry = window.meta['quintet'][:]
            entry.append(str(window.meta['total']))
            entry.append(str(window.stats['Dtotal']))
            entry.extend([str(window.stats[x])
                          for x in DIVNAMES[args.mode]])
            for dname in STATNAMES[args.mode]:
                entry.extend([str(window.stats[dname][y])
                              for y in ('left', 'right', 'Dtotal', 'D',
                                        'chisq', 'Pvalue')])
            entry.append(
                INTROGPATTERNS[args.mode][window.stats['signature']])
            outfile.write(('\t'.join(entry) + '\n').encode('utf-8'))
    return ''


if __name__ == "__main__":
    main()
