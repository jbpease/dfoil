# dfoil
**Software for detection of introgression in a five-taxon symmetric phylogeny** 

##Author
**James Pease** - [website](http://pages.iu.edu/~jbpease/) - [@jamesbpease](http://www.twitter.com/jamesbpease/)

## Citation Information
If you use this program, please cite:

```
James B Pease, Matthew W. Hahn.
"Detection and Polarization of Introgression in a Five-taxon Phylogeny" In revision.
http://biorxiv.org/content/early/2014/05/01/004689.short
doi: 10.1101/004689 
```
## Version History
* v. 2014-02-07 Re-release version on GitHub

## Requirements:
* Python 2.6+ (2.6+ and 3.x compatible)
* [Scipy](http://www.scipy.org/)
* [Numpy](http://www.numpy.org/)
* [matplotlib](http://www.matplotlib.org/)

## Installation
No installation is necessary, just download dfoil and run scripts through Python.
Dependencies should be installed according to their individual instructions. 

## Usage
`python dfoil.py INPUTFILE1 [, INPUTFILE2...] --out OUTPUTFILE1 [, OUTPUTFILE2]`

### Input Format
One or more input files can be specified.
Input files must be tab/space separated files with two starting fields:

* `CHROMOSOME` (this can be dummy values) 
* `POSITION` (also can be dummy values)

followed by the pattern counts in this order:

`AAAAA AAABA AABAA AABBA ABAAA ABABA ABBAA ABBBA, 
BAAAA BAABA BABAA BABBA BBAAA BBABA BBBAA BBBBA`

for the four-taxon test the patterns are
`AAAA AABA ABAA ABBA BAAA BABA BBAA BBBA`

*Note: these patterns are in 'binary' order 0000, 0010, 0100...*

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

` --plot_mode {none (default), write, show}`: show = show DFOIL/Dstat graph interactively, write = write to file(s) as specified 
	   
### Additional Plotting Options:

`--plot_path PLOTFILE1 [, PLOTFILE2, ...]`: One or more plot file paths (PNG format, must use "--plot_mode write")

`--plot_labels LABEL1 [, LABEL2, ...]`: Labels for the figure key in the order P1, P2, P3, P4

`--plot_color {color, dark, colornoanc, darknoanc, bw}`:

* color: standard full color plot
* colordark: standard plot with dark background and light text
* colornoanc/colornoancdark: color or dark plot without ancestral introgression shown
* bw: greyscale plot 
* bwdark: black background greyscale plot

`--plot_yscale FLOAT`: max value of y-axis (default = 1), y-min is complementary value (default = -1)

`--plot_smooth INT`: number of windows to average DFOIL statistics over for the plot (default = none)

`--plot_totals`: add a grey background plot of the total site pattern counts for each window

`--plot_hidekey`: hide the key in the plot (useful for making your own figures for publication, when fonts need to conform)

`--plot_hideaxes`: hide the axis labels and numbers

`--plot_lineweight`: weight of lines in plot (in pt, default=1.0)

`--plot_height FLOAT`: height of plot (in cm)

`--plot_width FLOAT`: width of plot (in cm)

`--plot_background FLOAT (0-1)`: alpha level of the introgression background shading (0 = transparent, 1=opaque, default=0.3)

## License
This file is part of dfoil.

dfoil is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as publihed by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

dfoil is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).