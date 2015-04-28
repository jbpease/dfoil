#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
http://www.github.com/jbpease/dfoil

dfoil - main introgression testing script
@author: James B. Pease

If you use this software please cite:
Pease JB and MW Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny"
Systematic Biology. Online.
http://www.dx.doi.org/10.1093/sysbio/syv023

v.2015-02-07 - Re-release on GitHub
v.2015-04-28 - Upgrades and Python3 compatibility fixes

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
import argparse, sys
from numpy import mean
from scipy.stats import chi2
import matplotlib


SIGNCODES = {'dfoil': {'+++0': 2, '--0+': 3, '++-0': 4, '--0-': 5, '+0++': 6,
                       '-0++': 7, '0+--': 8, '0---': 9, '++00': 10, '--00': 11},
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

class DataWindow(object):
    """Basic Handler of each data entry"""
    def __init__(self, counts=None, meta=None, stats=None):
        self.counts = counts or {}
        self.meta = meta or {}
        self.stats = stats or {}

    def dcalc(self, mincount=0):
        """Calculate D-statistics
            Arguments:
                mincount: minimuim total count to calculate P-values
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
                 self.counts[20] * beta1 +  self.counts[14] * beta2),
                mincount=mincount)
            self.stats['Dtotal'] = (
                (sum([self.counts[x] for x in (2, 4, 8, 16)]) * beta0)
                + (sum([self.counts[x] for x in (10, 12, 18, 20)]) * beta1)
                + (sum([self.counts[x] for x in (14, 22, 26, 28)]) * beta2))
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
                (sum([self.counts[x] for x in (4, 8)]) * beta0)
                + (sum([self.counts[x] for x in (6, 10)]) * beta1))
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
        self.stats['T34'] = float(counts.get(2, 0)
                                  + counts.get(4, 0)) / (2 * total)
        self.stats['T12'] = float(counts.get(8, 0)
                                  + counts.get(16, 0)) / (2 * total)
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
        self.stats['T12'] = (float(counts.get(8, 0) + counts.get(4, 0))
                             / (2.0 * total))
        self.stats['T123'] = 0.5 * (((float(counts.get(12, 0)) / total)
                                     + self.stats['T12'])
                                    + float(counts.get(2, 0)) /total)
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
        result['D'] = (float(left_term - right_term)
                       / (left_term + right_term))
    else:
        (val, pval) = chi2_test(left_term, right_term)
        result['chisq'] = val
        result['Pvalue'] = pval
        result['D'] = (float(left_term - right_term)
                       / (left_term + right_term))

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
    return "{}\n".format('\t'.join(
        ['#chrom', 'coord', 'total', 'dtotal']
        + DIVNAMES[mode]
        + ['{}_{}'.format(x, y)
           for x in STATNAMES[mode]
           for y in ('left', 'right', 'total',
                     'stat', 'chisq', 'Pvalue')]
        + ['introgression']
        + ['introg{}'.format(x) for x in INTROGPATTERNS[mode]]))


def plot_colors(colormode="color", linealpha=1, bgalpha=0.7):
    """Generate Plotting Line and Background Colors for Various Modes"""
    if colormode in ["color", "colornoanc"]:
        dstat_colors = [[(r, g, b, linealpha), '-'] for (r, g, b) in
                        [(0.11, 0.62, 0.47), (0.85, 0.37, 0),
                         (0.46, 0.44, 0.70), (0.91, 0.16, 0.54)]]
        bin_colors = [(0.60, 0.96, 0.85, bgalpha), (0.99, 0.55, 0.38, bgalpha),
                      (0.5, 0.5, 0.90, bgalpha), (1.0, 0.65, 0.86, bgalpha),
                      (0.20, 0.56, 0.45, bgalpha), (0.69, 0.25, 0.08, bgalpha),
                      (0.35, 0.35, 0.80, bgalpha), (0.71, 0.35, 0.56, bgalpha),
                      (0.2, 0.2, 0.2, bgalpha), (0.7, 0.7, 0.7, bgalpha)]
        bgcolors = ['w', 'k', bgalpha]
    elif colormode in ["colordark", "colornoancdark"]:
        bgalpha = bgalpha != 0 and 1 - bgalpha or 0
        dstat_colors = [[(r, g, b, linealpha), '-'] for (r, g, b) in
                        [(0.11, 0.62, 0.47), (0.85, 0.37, 0.00),
                         (0.46, 0.44, 0.70), (0.91, 0.16, 0.54)]]
        bgcolors = ['k', 'w', bgalpha]
        bin_colors = [(0.40, 0.76, 0.65, bgalpha), (0.99, 0.55, 0.38, bgalpha),
                      (0.55, 0.63, 0.80, bgalpha), (0.91, 0.54, 0.76, bgalpha),
                      (0.65, 0.85, 0.33, bgalpha), (1.00, 0.85, 0.18, bgalpha),
                      (0.90, 0.77, 0.58, bgalpha), (0.70, 0.70, 0.70, bgalpha),
                      (0.87, 0.80, 0.47, bgalpha), (0.79, 0.70, 0.84, bgalpha)]

    elif colormode == "bw":
        dstat_colors = [[(0, 0, 0, 1), "-"], [(0, 0, 0, 1), "--"],
                        [(0.6, 0.6, 0.6, 1), "-"], [(0.6, 0.6, 0.6, 1), "--"]]
        bgcolors = ['w', 'k', bgalpha]
        bin_colors = [(0, 0, 0, x) for x in [0.1]*10]
    elif colormode == "bwdark":
        dstat_colors = [[(1, 1, 1, 1), "-"], [(1, 1, 1, 1), "--"],
                        [(0.6, 0.6, 0.6, 1), "-"], [(0.6, 0.6, 0.6, 1), "--"]]
        bgcolors = ['k', 'w', bgalpha]
        bin_colors = [(0.9, 0.9, 0.9, x) for x in [0.1]*10]
    if colormode in ["colornoanc", "colornoancdark"]:
        bin_colors[8] = (0, 0, 0, 0)
        bin_colors[9] = (0, 0, 0, 0)
    return (dstat_colors, bin_colors, bgcolors)

def plot_dfoil(path, params, window_data, bool_data):
    """Plot DFOIL stats
        Arguments:
            path: file path for output or '' for interactive
            params: dictionary conversion of main params
            window_data: data from windows
            bool_data: introgression pres/absence binary data
    """
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    ## Set up labels
    dstat_names = STATNAMES[params['mode']]
    introgression_labels = INTROGLABELS[params['mode']]
    if params['plot_labels']:
        for i, elem in enumerate(introgression_labels):
            for (j, label) in enumerate(params['plot_labels']):
                introgression_labels[i] = (elem.replace(str(j+1), label))
    for i, elem in enumerate(introgression_labels):
        if '+' in elem:
            introgression_labels[i] = elem.replace('$\\Rightarrow$',
                                                   '$\\Leftrightarrow$')
    ## Establish Colors
    (dstat_colors, bin_colors,
     bgcolors) = plot_colors(colormode=params['plot_color'],
                             bgalpha=params['plot_background'])
    ## Calcuate plot values
    xdstat = [int(window.meta['position']) for window in window_data]
    dplots = [[window.stats[dstat]['D'] for window in window_data]
              for dstat in dstat_names]
    ## Begin PLot
    if params['plot_smooth']:
        dplots = [[mean(dplot[x:x + params['plot_smooth']])
                   for x in range(0, len(dplot), params['plot_smooth'])]
                  for dplot in dplots]
        xdstat = xdstat[::params['plot_smooth']]
    xbin = [int(window.meta['position']) for window in window_data]
    if params['plot_color'] == 'bw':
        if params['mode'] == 'dstat':
            binplot = [int('1' in bool_data[x][2:]) * 10000 - 5000
                       for x in range(len(bool_data))]
    else:
        binplots = [[int(window[x])*10000 - 5000 for window in bool_data]
                    for x in range(2, len(bool_data[0]))]
    if params['plot_smooth']:
        totalplot = [mean([x[3] for x in
                           window_data[y:y + params['plot_smooth']]])
                     for y in range(0, len(window_data),
                                     params['plot_smooth'])]
    else:
        totalplot = [int(window.meta['total'])
                     for window in window_data]
    fig, host = plt.subplots(figsize=(params['plot_width'],
                                      params['plot_height']))
    if params['plot_totals']:
        par1 = host.twinx()
    for i, dplot in enumerate(dplots):
        dashlen = dstat_colors[i][1] == '--' and (2, 2) or ''
        if params['mode'] == 'dstat' and params['plot_color'] == 'bw':
            dlabels = params['plot_labels'] or ['$P+1$', '$P_2$', '$P_3$']
            host.plot(xdstat, dplot, color=dstat_colors[i][0],
                      linewidth=params['plot_lineweight'],
                      linestyle=dstat_colors[i][1], dashes=dashlen,
                      label=("$D(+)$=" + dlabels[2] + "$\\Leftrightarrow$"
                             + dlabels[1] +  ";   $D(-)$="
                             + dlabels[2] + "$\\Leftrightarrow$" + dlabels[0]),
                      drawstyle="steps-pre")
        else:
            host.plot(xdstat, dplot, color=dstat_colors[i][0],
                      linewidth=params['plot_lineweight'],
                      linestyle=dstat_colors[i][1], dashes=dashlen,
                      label=("$" + dstat_names[i][0] + "_{"
                             + dstat_names[i][1:] + "}$"),
                      drawstyle="steps-pre")
    if params['plot_background']:
        if params['plot_color'] == 'bw':
            host.fill_between(xbin, binplot, -1, edgecolor='none',
                              facecolor=(0.7, 0.7, 0.7,
                                         params['plot_background']),
                              label=("Significant (P<{})").format(
                                  params['pvalue'][0]))
        else:
            for i, binplot in enumerate(binplots):
                host.fill_between(xbin, binplot, -1, edgecolor='none',
                                  facecolor=bin_colors[i],
                                  label=(introgression_labels[i+2]))
    host.set_ylim(-1 * params['plot_yscale'], params['plot_yscale'])
    if params['plot_totals']:
        par1.fill_between(xdstat, y1=totalplot, y2=0, facecolor=(0, 0, 0, 0.5),
                          edgecolor="none")
        par1.set_ylim(0, max(totalplot) * 5)
        par1.set_ylabel("Total Count of Sites")
    host.set_xlim(xmin=min(xbin), xmax=max(xbin))
    host.set_axis_bgcolor(bgcolors[0])
    host.tick_params(direction="in", length=5, width=0.5)
    if params['plot_hideaxes']:
        host.set_xticklabels([], visible=False)
        host.set_yticklabels([], visible=False)
        fig.tight_layout()
    else:
        host.set_xlabel("Position")
        host.set_ylabel("$D$")
    fig.patch.set_alpha(0.0)
    if not params['plot_hidekey']:
        leg = host.legend(loc=9, ncol=7)
        leg.get_frame()
        if params['plot_background']:
            for label in leg.get_texts():
                label.set_fontsize(10)
                label.set_color(bgcolors[1])
            for i, patch in enumerate(leg.get_patches()):
                patch.set_alpha(bin_colors[i][3])
    host.axhline(y=0, color=bgcolors[1], linewidth=0.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        host.spines[axis].set_linewidth(0.5)
    if not path:
        plt.show()
    else:
        plt.savefig(path)
    return ''

def fill_windows(window_data, run_length):
    """Fill background color between introgressing windows that
       are only separated by a certain number of non-introgressing windows
    """
    i = 0
    n_windows = len(window_data)
    while i < n_windows:
        if window_data[i].stats['signature'] > 1:
            j = i + 1
            while j < n_windows:
                if j == n_windows - 1:
                    i = j
                    break
                if j - i == run_length:
                    i = j - 1
                    break
                if window_data[j].stats['signature'] > 1:
                    if (window_data[j].stats['signature'] ==
                            window_data[i].stats['signature']):
                        for k in range(i + 1, j):
                            window_data[k].stats['signature'] = (
                                window_data[i].stats['signature'] + 0)
                    else:
                        i = j - 1
                        break
                j += 1
        i += 1
    return window_data

def main(arguments=sys.argv[1:]):
    """Main dfoil method"""
    parser = argparse.ArgumentParser(description=("""
    Calculate DFOIL and D-statistics stats from one or more count files.
    USAGE: dfoiler.py INPUTFILE1 ... --out OUTPUTFILE1 ..."""))
    parser.add_argument('infile', help="input tab-separated counts file",
                        nargs='*')
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
    parser.add_argument("--beta2", type=float,
                        help="""beta2 coefficient for double-B patterns,
                                defaults: Dpart=N.A., others=1.0""")
    parser.add_argument("--beta3", type=float,
                        help="""beta3 coefficient for triple-B patterns
                                defaults: DFOIL/DFOILalt=1.0,
                                Dstat/Dpart=N.A.""")
    parser.add_argument("--plot", choices=["show", "write", "none"],
                        default="none",
                        help="""show=interactive,
                                write=write to file (--plot_path)
                                or none=no plot""")
    parser.add_argument("--plot_path",
                        help="output filepath for plot",
                        nargs='*')
    parser.add_argument('--plot_labels', nargs='*', help="taxon labels")
    parser.add_argument('--plot_color', default="color",
                        choices=["color", "colordark",
                                 "bw", "bwdark",
                                 "colornoanc", "colornoancdark"],
                        help="""choose color mode""")
    parser.add_argument("--plot_lineweight", type=float, default=1.,
                        help="line weight for dplots (default=1pt)")
    parser.add_argument("--plot_yscale", type=float, default=1.0,
                        help="Y-axis min-max value, default is 1")
    parser.add_argument("--plot_smooth", type=int,
                        help="average D-stats over this number of points")
    parser.add_argument("--plot_background", type=float, default=0.3,
                        help="0-1.0 background intensity, default = 0.3")
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
    parser.add_argument("--version", action="store_true",
                        help="display version information and quit")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("DFOIL v. 2015-02-07")
        sys.exit()
#### ===== INITIALIZE =====
    if not args.out:
        raise RuntimeError("--out path not specified")
    if args.plot_path:
        if set(args.infile) & set(args.plot_path):
            raise NameError
    if set(args.infile) & set(args.out):
        raise NameError("input and output file have same path")
    if args.plot == 'write' and not args.plot_path:
        raise NameError(
            "write mode specified for plot but --plot_path not specified")
    if len(args.pvalue) == 1:
        args.pvalue = [args.pvalue[0], args.pvalue[0]]
    ## Set beta parameters for presets
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
#### ===== PARSE COUNT FILE  =====
    for ifile, infilename in enumerate(args.infile):
        window_data = []
        with open(infilename) as infile:
            for line in infile:
                if line[0] == '#':
                    continue
                arr = line.rstrip().split()
                window = DataWindow(meta=dict(
                    chrom=arr[0], position=int(arr[1]),
                    mode=args.mode,
                    beta=(args.beta1, args.beta2, args.beta3)))
                if args.mode in ["dfoil", "partitioned"]:
                    window.counts = dict([(j - 2) * 2, int(arr[j])]
                                         for j in range(2, 18))
                elif args.mode == 'dstat':
                    window.counts = dict([(j - 2) * 2, int(arr[j])]
                                         for j in range(2, 9))
                if sum(window.counts.values()) < args.mintotal:
                        continue
                window.meta['total'] = sum(window.counts.values())
                window.dcalc(mincount=args.mincount)
                window.calc_signature(pvalue_cutoffs=args.pvalue)
                window_data.append(window)
######## ===== ANALYZE WINDOWS AND DETERMINE RUNS ====-
        if args.runlength:
            window_data = fill_windows(window_data, args.runlength)
######## ===== WRITE TO OUTPUT =====
        with open(args.out[ifile], 'w') as outfile:
            outfile.write(make_header(args.mode))
            bool_data = []
            for window in window_data:
                bool_flags = [
                    '0' for x in range(len(INTROGPATTERNS[args.mode]))]
                bool_flags[window.stats['signature']] = '1'
                bool_data.append(bool_flags)
                entry = [str(window.meta[k]) for k in [
                    'chrom', 'position', 'total']]
                entry.append(str(window.stats['Dtotal']))
                entry.extend([str(window.stats[x])
                              for x in DIVNAMES[args.mode]])
                for dname in STATNAMES[args.mode]:
                    entry.extend([str(window.stats[dname][y])
                                  for y in ('left', 'right', 'Dtotal', 'D',
                                            'chisq', 'Pvalue')])
                entry.append(
                    INTROGPATTERNS[args.mode][window.stats['signature']])
                entry.extend(bool_flags)
                outfile.write('\t'.join(entry) + '\n')
######## ==== PLOT GRAPHS =====
        if args.plot != "none":
            plot_dfoil(args.plot == 'write' and args.plot_path[ifile] or '',
                       vars(args), window_data, bool_data)
    return ''
if __name__ == "__main__":
    main()
