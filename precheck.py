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
    sum_data = {}
    for window in window_data:
        for code in window.counts:
            sum_data[code] = sum_data.get(
                code, 0) + window.counts[code]
    # Check 1
    checkok = True
    if verbose is True:
        print("="*79)
        print("""Checking that concordant patterns are more common that discordant
(Note that this is normal when introgression is extreme, but could also
indicate the taxa are out of order):""")
        print("-"*79)
    for concode in (2, 4, 8, 16, 6, 24):
        for discode in (10, 12, 14, 18, 20, 22, 26, 28):
            if sum_data[concode] < sum_data[discode]:
                print(("WARNING: Total count of "
                       "discordant pattern {}={} is higher than "
                       "concordant pattern {}={}").format(
                    SITECODES[discode], sum_data[discode],
                    SITECODES[concode], sum_data[concode]))
                checkok = False
    if checkok:
        print("Pass")
    # Check2
    checkok = True
    if verbose is True:
        print("="*79)
        print("Checking that divergences are correctly ordered "
              "(P1 and P2 should diverge AFTER P3 and P4)")
        print("-"*79)
    for abcode in (8, 16):
        for cdcode in (2, 4):
            if sum_data[abcode] > sum_data[cdcode]:
                print(("WARNING: Total count of A/B terminal substitutions "
                       "{}={} is higher than C/D terminal substitutions {}={}"
                       ).format(SITECODES[abcode], sum_data[abcode],
                                SITECODES[cdcode], sum_data[cdcode]))
                checkok = False
    if checkok:
        print("Pass")
    print("="*79)
    # Check3
    checkok = True
    if verbose is True:
        print("="*79)
        print("Checking that terminal branch pairs are "
              "proportionate approximately")
        print("-"*79)
    abratio = float(sum_data[16]) / sum_data[8] if sum_data[8] > 0 else "inf"
    print("BAAAA/ABAAA ratio = {} ({}/{})".format(abratio, sum_data[16],
                                                  sum_data[8]))
    if abratio == "inf" or (0.8 < abratio > 1.25):
        checkok = False
        print("WARNING: P1/P2 ratio deviates more than 25% from 1.0")
    cdratio = float(sum_data[4]) / sum_data[2] if sum_data[2] > 0 else "inf"
    print("AABAA/AAABA ratio = {} ({}/{})".format(cdratio, sum_data[4],
                                                  sum_data[2]))
    if cdratio == "inf" or (0.8 < cdratio > 1.25):
        checkok = False
        print("WARNING: P3/P4 ratio deviates more than 25% from 1.0")
    if checkok:
        print("Pass")
    print("="*79)
    return ''

if __name__ == "__main__":
    print("pre-dfoil can no longer be run on its own.")
    print("Please use dfoil.py with the --pre-check-only "
          "flag to just run a pre-check")
