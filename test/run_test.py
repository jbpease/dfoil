#!/usr/bin/env python
"""Run dfoil test dataset"""

import sys
import subprocess

# Run fasta2dfoil
args = ("../fasta2dfoil.py", sys.argv[1],
        "-o", sys.argv[1] + ".counts",
        "--names P1,P2,P3,P4,PO")
print("Running ", " ".join(args))
proc = subprocess.Popen(' '.join(args), stdout=sys.stdout, stderr=sys.stderr,
                        shell=True)
proc.communicate()
del proc

# Run dfoil

args = ("../dfoil.py",
        "--infile", sys.argv[1] + ".counts",
        "--out", sys.argv[1] + ".dfoil",
        "--mode dfoil",
        "--plot", sys.argv[1] + ".dfoil.pdf")
print("Running ", " ".join(args))
proc = subprocess.Popen(' '.join(args), stdout=sys.stdout, stderr=sys.stderr,
                        shell=True)
proc.communicate()
del proc
