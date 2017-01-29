# dfoil
**Software for detection of introgression in a five-taxon symmetric phylogeny** 

##Author
**James Pease** - [website](http://jbpease.github.io) - [@jamesbpease](http://www.twitter.com/jamesbpease/)

## Citation Information
If you use this program, please cite:
```
James B Pease, Matthew W. Hahn. 2015.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny" 
Systematic Biology. 64 (4): 651-662.
http://www.dx.doi.org/10.1093/sysbio/syv023
doi: 10.1093/sysbio/syv023
```
## Version History
* v. 2014-02-07 Re-release version on GitHub
* v. 2015-04-17 Minor updates and citation information, Publication Release Version
* v. 2014-04-28 Fixes for Python3 compatibility
* v. 2015-11-23 More fixes for Python3 compability, added 'pre-dfoil.py' that checks for issues in count files before running dfoil.py
* v. 2017-01-29 Fixes for visual graph outputs, `--plot_labels` has been fixed, replaced `colornoanc` and `colornoancdark` with `--plot_noanc` option, changes to color pallette for code comprehension. Replaced `--plot-path` with just `--plot`, and `--plot show` is now deprecated.

=======

## Requirements:
* Python 2.7+ or 3.x
* [Scipy](http://www.scipy.org/)
* [Numpy](http://www.numpy.org/)
* [matplotlib](http://www.matplotlib.org/)

To run simulations within dfoil_sim.py you will also need:
* [ms](http://home.uchicago.edu/rhudson1/source/mksamples.htm)

## Installation
No installation is necessary, just download dfoil and run scripts through Python.
Dependencies should be installed according to their individual instructions. 

## Usage
`python dfoil.py --infile INPUTFILE1 [INPUTFILE2 ...] --out OUTPUTFILE1 [OUTPUTFILE2 ...]`

### Input Format
One or more input files can be specified.
These files can have any number of header lines at the beginning (including none), but **all header lines must start with '#'**

Data fields should be tab/space separated files with two starting fields:

* `CHROMOSOME` (this can be dummy values) 
* `POSITION` (also can be dummy values)

followed by the pattern counts in this order:

`AAAAA AAABA AABAA AABBA ABAAA ABABA ABBAA ABBBA, 
BAAAA BAABA BABAA BABBA BBAAA BBABA BBBAA BBBBA`

for the four-taxon test the patterns are
`AAAA AABA ABAA ABBA BAAA BABA BBAA BBBA`

*Note: these patterns are in 'binary' order 0000, 0010, 0100...*

**IMPORTANT**
The order of taxa must be P1 P2 P3 P4 O, such that:
*"O" is the outgroup
*P1 and P2 are a monophyletic pair of taxa
*P3 and P4 are a monophyletic pair of taxa
*P3 and P4 divergence >= P1 and P2 divergence
(The choice of P1/P2 and P3/P4 within the pairings is arbitrary)

### Output Format
One or more output files are specified (equal to number of inputs).
The outputs will have fields:

* `CHROMOSOME`
* `POSITION`

then for each D-statistic:

* Dxx_left (left term value)
* Dxx_right (right term value)
* Dxx_stat (D-statistic value)
* Dxx_chisq (Chi_Squared value)
* Dxx_pval (Chi_Squared P-value)

-----------------------------------------------------------
### Main Options
`--mincount INT`: mininum number of total pattern counts in the D-statistic denominator (default=10)

`--mintotal INT`: minimum number of total sites in the given region (i.e. per line, specified in column 3 of input)

`--pvalue FLOAT [FLOAT]`: minimum P-value(s) for significance 
                   (use two values for separate P for DFO/DIL and DFI/DOL, or D1/D2 and D12)

`--mode {dfoil (default), dfoilalt, partitioned, dstat, dstatalt}`

* dfoil = standard DFOIL for five-taxa
* dfoilalt = dfoil without single-B patterns (ABAAA, etc.)
* partitioned = Partitioned D-statistics (Eaton & Ree 2013)
* dstat = four-taxon D-statistic (Green et al. 2010)
* dstatalt = four-taxon D-statistic with inverse patterns added (untested, use with caution)

`--beta1, --beta2, --beta3`: weighting factors (0 <= b <= 1), for single-B, double-B, and triple-B patterns (respectively) by default these are set as:

* dfoil = 1,1,1; 
* dfoilalt = 0,1,1;
* dstat: 0,1,N/A (no triple-B);
* dstatalt: 1,1,N/A (no triple-B);
* partitioned: N/A (does not use weighting parameters);

### Plotting Options:

`--plot PLOTFILE1 [PLOTFILE2 ...]`: Turns on plotting, outputs to one or more space-separated plot file paths (path extensions used to determine plot file format)

`--plot_labels LABEL1 LABEL2 LABEL3 LABEL4`: Space-separated Labels for the figure key in the order P1 P2 P3 P4 (use only 3 for `--mode dstat` or `dstatalt` )

`--plot_color {color,dark,bw,bwdark}`:

* color: standard full color plot
* colordark: standard plot with dark background and light text
* bw: grayscale plot 
* bwdark: greyscale plot with dark background and light text

`--plot_noanc`: boolean to remove ancestral introgressions background coloring and removes from key

`--plot_yscale FLOAT`: max value of y-axis (default = 1), y-min is complementary value (default = -1)

`--plot_smooth INT`: number of windows to average DFOIL statistics over for the plot (default = none)

`--plot_totals`: boolean to add a grey background plot of the total site pattern counts for each window

`--plot_hidekey`: boolean to hide the key in the plot (useful for making your own figures for publication, when fonts need to conform)

`--plot_hideaxes`: boolean to hide the axis labels and numbers

`--plot_lineweight FLOAT`: weight of lines in plot (in pt, default=1.0)

`--plot_height FLOAT`: height of plot (in cm)

`--plot_width FLOAT`: width of plot (in cm)

`--plot_background FLOAT (0-1)`: alpha level of the introgression background shading (0 = transparent, 1=opaque, default=0.3)

## License
This file is part of dfoil.

dfoil is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as publihed by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

dfoil is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).
