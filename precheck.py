#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DFOIL: Directional introgression testing a five-taxon phylogeny
pre-dfoil - site count checking testt
James B. Pease
http://www.github.com/jbpease/dfoil
"""

from __future__ import print_function, unicode_literals

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


SITECODES = dict([(x, "{}{}".format('A' * (7 - len(bin(x))),
                  str(bin(x))[2:].replace('0', 'A').replace('1', 'B')))
                  for x in range(0, 32, 2)])


def pre_check(window_data, mode='dfoil', verbose=True):
    sum_data = dict([(x, 0) for x in range(0, 32, 2)])
    for window in window_data:
        for code in window.counts:
            sum_data[code] = sum_data.get(
                code, 0) + window.counts[code]
    # Check 1
    check_concordant(sum_data, mode, verbose=verbose)
    # Check2
    divergence_order(sum_data, mode, verbose=verbose)
    # Check3
    check_terminal(sum_data, mode, verbose=verbose)
    return ''


def check_concordant(sum_data, mode, verbose=True):
    checkok = True
    if verbose is True:
        print("="*79)
        print("""Checking that concordant patterns are more common that discordant
(Note that this is normal when introgression is extreme, but could also
indicate the taxa are out of order):""")
        print("-"*79)
    if mode == 'dfoil':
        concordant_patterns = (2, 4, 8, 16, 6, 24)
        discordant_patterns = (10, 12, 14, 18, 20, 22, 26, 28)
    elif mode == 'dstat':
        concordant_patterns = (2, 4, 8, 12)
        discordant_patterns = (6, 10)
    for concode in concordant_patterns:
        for discode in discordant_patterns:
            if sum_data.get(concode, 0) < sum_data.get(discode, 0):
                print(("WARNING: Total count of "
                       "discordant pattern {}={} is higher than "
                       "concordant pattern {}={}").format(
                    SITECODES.get(discode, 0), sum_data.get(discode, 0),
                    SITECODES.get(concode, 0), sum_data.get(concode, 0)))
                checkok = False
    if checkok and verbose is True:
        print("Pass")
    print("="*79)
    return ''


def divergence_order(sum_data, mode, verbose=True):
    if mode != 'dfoil':
        return ''
    checkok = True
    if verbose is True:
        print("="*79)
        print("Checking that divergences are correctly ordered "
              "(P1 and P2 should diverge AFTER P3 and P4)")
        print("-"*79)
    for abcode in (8, 16):
        for cdcode in (2, 4):
            if sum_data.get(abcode, 0) > sum_data.get(cdcode, 0):
                print(("WARNING: Total count of P1/P2 terminal substitutions "
                       "{}={} is higher than P3/P4 terminal substitutions "
                       "{}={}"
                       ).format(SITECODES.get(abcode, 0),
                                sum_data.get(abcode, 0),
                                SITECODES.get(cdcode, 0),
                                sum_data.get(cdcode, 0)))
                checkok = False
    if checkok:
        print("Pass")
    print("="*79)
    return ''


def check_terminal(sum_data, mode, verbose=True):
    checkok = True
    if verbose is True:
        print("="*79)
        print("Checking that terminal branch pairs are "
              "proportionate approximately")
        print("-"*79)
    if mode == 'dfoil':
        abratio = (float(sum_data[16]) / sum_data[8]
                   if sum_data[8] > 0 else "inf")
        print("BAAAA/ABAAA ratio = {} ({}/{})".format(abratio, sum_data[16],
                                                      sum_data[8]))
        if abratio == "inf" or (0.8 < abratio > 1.25):
            checkok = False
            print("WARNING: P1/P2 ratio deviates more than 25% from 1.0")
        cdratio = (float(sum_data[4]) / sum_data[2]
                   if sum_data[2] > 0 else "inf")
        print("AABAA/AAABA ratio = {} ({}/{})".format(cdratio, sum_data[4],
                                                      sum_data[2]))
        if cdratio == "inf" or (0.8 < cdratio > 1.25):
            checkok = False
            print("WARNING: P3/P4 ratio deviates more than 25% from 1.0")
    elif mode == 'dstat':
        cdratio = (float(sum_data[8]) / sum_data[4]
                   if sum_data[4] > 0 else "inf")
        print("BAAA/ABAA ratio = {} ({}/{})".format(cdratio, sum_data[8],
                                                    sum_data[4]))
        if cdratio == "inf" or (0.8 < cdratio > 1.25):
            checkok = False
            print("WARNING: P1/P2 ratio deviates more than 25% from 1.0")
        if checkok:
            print("Pass")
    print("="*79)


if __name__ == "__main__":
    print("pre-dfoil can no longer be run on its own.")
    print("Please use dfoil.py with the --pre-check-only "
          "flag to just run a pre-check")
