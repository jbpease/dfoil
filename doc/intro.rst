###############
Getting Started
###############

What is dfoil?
==============
*D*\ :sub:`FOIL`\ is a method for testing introgression in a five-taxon symmetric phylogeny.  


How do I cite this sofware?
===========================

If you use this program, please cite:

``James B Pease, Matthew W. Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny" 
Systematic Biology. 64 (4): 651-662.
http://www.dx.doi.org/10.1093/sysbio/syv023
doi: 10.1093/sysbio/syv023``

Please also include the link <https://www.github.com/jbpease/dfoil> in your publication.

Requirements
============
* Python 2.7.x or 3.x
* Scipy: http://www.scipy.org/
* Numpy: http://www.numpy.org/
* matplotlib 1.5.3+: http://www.matplotlib.org/

Optional
--------
To run simulations within dfoil_sim.py you will also need:
* ms: http://home.uchicago.edu/rhudson1/source/mksamples.htm

Installation
============
No installation is necessary, just download dfoil and run scripts through Python.
``git clone https://www.github.com/jbpease/dfoil``
Other required software should be installed according to their individual instructions. 


Preparing your data
===================

Generating an AB-pattern count file from a FASTA file.
------------------------------------------------------
Use the script fasta2dfoil (see options below) to generate the count file from a FASTA.

Using a prepared count file:
----------------------------

One or more input files can be specified. These files can have any number of header lines at the beginning (including none), but **all header lines must start with '#'**

Data fields should be tab/space separated files with two starting fields:

* ``CHROMOSOME`` (this can be dummy values) 
* ``POSITION`` (also can be dummy values)

followed by the pattern counts in this order:

``AAAAA AAABA AABAA AABBA ABAAA ABABA ABBAA ABBBA, BAAAA BAABA BABAA BABBA BBAAA BBABA BBBAA BBBBA``

for the four-taxon test the patterns are
``AAAA AABA ABAA ABBA BAAA BABA BBAA BBBA``

.. hint::  these patterns are in 'binary' order 0000, 0010, 0100...

.. important:: 
        The order of taxa must be P1 P2 P3 P4 O, such that:
        * "O" is the outgroup
        * P1 and P2 are a monophyletic pair of taxa
        * P3 and P4 are a monophyletic pair of taxa
        * P3 and P4 divergence occurs before (in forward time) the divergence of P1 and P2
        (The choice of P1/P2 and P3/P4 within the pairings is arbitrary)


Running dfoil.py
================

Basic usage
------------

```python dfoil.py --infile INPUTFILE1 [INPUTFILE2 ...] --out OUTPUTFILE1 [OUTPUTFILE2 ...]```

Pre-check
---------
The data will undergo a precheck (use ``--skip-pre-check`` to turn off, or ``--pre-check-only`` to only run the pre-check.
This will check for common issues in count data that might affect the result or violate the assumptions.

Common issues include:
* An accelerated rate of substitutions (i.e. an excess of B's) on a specific branch relative to its sister taxon
* Mis-labeling of P1/P2 and P3/P4 (remember P3/P4 divergence should come first, in forward time)

Modes
-----
* dfoil = standard DFOIL for five-taxa (this is default)
* dfoilalt = dfoil without single-B patterns (ABAAA, etc.)
* partitioned = Partitioned D-statistics (Eaton & Ree 2013)
* dstat = four-taxon D-statistic (Green et al. 2010)
* dstatalt = four-taxon D-statistic with inverse patterns added (use with caution)

Advanced Weight Parameters
--------------------------
The parameters `--beta1`, `--beta2`, and `--beta3` are weighting factors (0 <= b <= 1), for single-B, double-B, and triple-B patterns, respectively. Ordinarily you will not need to set these.

By default these are set as:

* dfoil = 1,1,1 
* dfoilalt = 0,1,1
* dstat: 0,1,N/A (no triple-B)
* dstatalt: 1,1,N/A (no triple-B)
* partitioned: N/A (does not use weighting parameters)

Output Format
=============
One or more output files are specified (equal to number of inputs).
The outputs will have fields:

* ``CHROMOSOME``
* ``POSITION``

then for each D-statistic:

* ``Dxx_left`` (left term value)
* ``Dxx_right`` (right term value)
* ``Dxx_stat`` (D-statistic value)
* ``Dxx_chisq`` (Chi_Squared value)
* ``Dxx_pval`` (Chi_Squared P-value)

(where "xx" will be replaced with FO, IL, FI, OL)

-----------------------------------------------------------
### Plotting Options:

Releases
========

2017-06-14
----------
Major upgrade to Sphinx documentation. Integrated the pre-check (formerly pre-dfoil.py) into the main script. Minor fixes to syntax.

2017-01-29
----------
Fixes for visual graph outputs, `--plot_labels` has been fixed, replaced `colornoanc` and `colornoancdark` with `--plot_noanc` option, changes to color pallette for code comprehension. Replaced `--plot-path` with just `--plot`, and `--plot show` is now deprecated.

2015-11-23 
----------
More fixes for Python3 compability, added 'pre-dfoil.py' that checks for issues in count files before running dfoil.py

2015-04-17 
----------
Minor updates and citation information, Publication Release Version

2014-04-28 
----------
Fixes for Python3 compatibility

2014-02-07 
----------
Re-release version on GitHub


License
=======

This file is part of dfoil.

dfoil is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as publihed by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

dfoil is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).
