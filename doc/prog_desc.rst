Program Parameter Descriptions
##############################

.. dfoil:

dfoil
=====

Description
-----------

DFOIL: Directional introgression testing a five-taxon phylogeny
dfoil - Calculate DFOIL and D-statistics stats from one or more count files.
James B. Pease
http://www.github.com/jbpease/dfoil

USAGE: dfoil.py INPUTFILE1 ... --out OUTPUTFILE1 ...


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--infile`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** input tab-separated counts file

**Type:** None; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** outputs tab-separated DFOIL stats

**Type:** None; **Default:** None



``--beta1``
^^^^^^^^^^^

**Description:** beta1 coefficient for single-B patterns,
                                 defaults: DFOIL/Dstatalt=1.0,
                                 DFOILalt,Dstat=0,Dpart=N.A.

**Type:** float; **Default:** None



``--beta2``
^^^^^^^^^^^

**Description:** beta2 coefficient for double-B patterns,
                                defaults: Dpart=N.A., others=1.0

**Type:** float; **Default:** 1.0



``--beta3``
^^^^^^^^^^^

**Description:** beta3 coefficient for triple-B patterns
                                defaults: DFOIL/DFOILalt=1.0,
                                Dstat/Dpart=N.A.

**Type:** float; **Default:** 1.0



``--mincount``
^^^^^^^^^^^^^^

**Description:** minium number of D denominator sites per window

**Type:** integer; **Default:** 10



``--mintotal``
^^^^^^^^^^^^^^

**Description:** minimum total number of sites in a region

**Type:** integer; **Default:** 50



``--mode``
^^^^^^^^^^

**Description:** dfoil = DFOIL,
                                dfoilalt = DFOIL without single-B patterns,
                                partitioned = Partitioned D-statistics,
                                dstat = Four-Taxon D-statistic,
                                dstatalt = Four-Taxon D-statistic
                                with single-B patterns

**Type:** None; **Default:** dfoil

**Choices:** ['dfoil', 'dfoilalt', 'partitioned', 'dstat', 'dstatalt']


``--plot``
^^^^^^^^^^

**Description:** write plot to file path(s) given

**Type:** None; **Default:** None



``--plot_background``
^^^^^^^^^^^^^^^^^^^^^

**Description:** 0-1.0 background intensity
                                (0=none, default=0.3)

**Type:** float; **Default:** 0.3



``--plot_color``
^^^^^^^^^^^^^^^^

**Description:** choose color mode

**Type:** None; **Default:** color

**Choices:** ['color', 'colordark', 'bw', 'bwdark']


``--plot_height``
^^^^^^^^^^^^^^^^^

**Description:** height of plot (in cm)

**Type:** float; **Default:** 8.0



``--plot_hideaxes``
^^^^^^^^^^^^^^^^^^^

**Description:** hide axes labels

**Type:** boolean flag



``--plot_hidekey``
^^^^^^^^^^^^^^^^^^

**Description:** hide plot key

**Type:** boolean flag



``--plot_labels``
^^^^^^^^^^^^^^^^^

**Description:** taxon labels

**Type:** None; **Default:** None



``--plot_lineweight``
^^^^^^^^^^^^^^^^^^^^^

**Description:** line weight for dplots (default=1pt)

**Type:** float; **Default:** 1.0



``--plot_noanc``
^^^^^^^^^^^^^^^^

**Description:** do not plot background for
                                ancestral introgression

**Type:** boolean flag



``--plot_smooth``
^^^^^^^^^^^^^^^^^

**Description:** average D-stats over this number of points

**Type:** integer; **Default:** None



``--plot_totals``
^^^^^^^^^^^^^^^^^

**Description:** add a background plot of total site counts

**Type:** boolean flag



``--plot_width``
^^^^^^^^^^^^^^^^

**Description:** width of plot (in cm)

**Type:** float; **Default:** 24.0



``--plot_yscale``
^^^^^^^^^^^^^^^^^

**Description:** Y-axis min-max value, default is 1

**Type:** float; **Default:** 1.0



``--pre-check-only``
^^^^^^^^^^^^^^^^^^^^

**Description:** Only run the data pre-check (formely pre-dfoil.py)

**Type:** boolean flag



``--pvalue``
^^^^^^^^^^^^

**Description:** minimum P-value cutoff for regions,
                                can specify one P-value for all four tests
                                or two separate ones for DFO/DIL and DFI/DOL
                                (or D1/D2 and D12 for 'partitioned')

**Type:** float; **Default:** [0.01, 0.01]



``--runlength``
^^^^^^^^^^^^^^^

**Description:** if two introgressing windows are separated
                                by this many windows of non-introgression
                                color in the intervening windows to
                                create a more continuous visual appearance

**Type:** integer; **Default:** 0



``--skip-pre-check``
^^^^^^^^^^^^^^^^^^^^

**Description:** Skip running the data pre-check (formely pre-dfoil)

**Type:** boolean flag



``--zerochar``
^^^^^^^^^^^^^^

**Description:** list of strings used in place of zeros
                                in the input file default is [".", "NA"]

**Type:** None; **Default:** ['.', 'NA']


.. dfoil_analyze:

dfoil_analyze
=============

Description
-----------

DFOIL: Directional introgression testing a five-taxon phylogeny
dfoil_analyze: Given a dfoil output file, gives summary statistics to stdout
James B. Pease
http://www.github.com/jbpease/dfoil


Parameters
----------

infile
^^^^^^

**Description:** dfoil output file

**Type:** None; **Default:** None



``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--ndigits``
^^^^^^^^^^^^^

**Description:** number of decimal places

**Type:** integer; **Default:** 3


.. dfoil_sim:

dfoil_sim
=========

Description
-----------

DFOIL: Directional introgression testing a five-taxon phylogeny
dfoil_sim - simulation of sequences for testing dfoil
James B. Pease
http://www.github.com/jbpease/dfoil


Parameters
----------

outputfile
^^^^^^^^^^

**Description:**  output site count filename

**Type:** file path; **Default:** None



``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--coaltimes``
^^^^^^^^^^^^^^^

**Description:** coalescent times in 4Ne units

**Type:** float; **Default:** (3, 2, 1, 1)



``--mdest``
^^^^^^^^^^^

**Description:** 1-based index of migration recipient population

**Type:** integer; **Default:** None



``--mrate``
^^^^^^^^^^^

**Description:** per individual per generation migration
                                rate (default=5e-4)

**Type:** float; **Default:** 0.0005



``--msfile``
^^^^^^^^^^^^

**Description:** use pre-computed ms output
                                             file instead of running ms.

**Type:** None; **Default:** None



``--msource``
^^^^^^^^^^^^^

**Description:** 1-based index of migration source population

**Type:** integer; **Default:** None



``--mspath``
^^^^^^^^^^^^

**Description:** path to ms executable

**Type:** None; **Default:** ms



``--mtimes``
^^^^^^^^^^^^

**Description:** time bounds for the migration period

**Type:** float; **Default:** None



``--mu``
^^^^^^^^

**Description:** per site per generation mutation rate
                                (default=7e-9)

**Type:** float; **Default:** 7e-09



``--nconverge``
^^^^^^^^^^^^^^^

**Description:** number of convergent sites per window

**Type:** integer; **Default:** 0



``--nloci``
^^^^^^^^^^^

**Description:** number of windows to simulate

**Type:** integer; **Default:** 100



``--popsize``
^^^^^^^^^^^^^

**Description:** Ne, effective population size (default=1e6)

**Type:** integer; **Default:** 1000000.0



``--quiet``
^^^^^^^^^^^

**Description:** suppress screen output

**Type:** boolean flag



``--recomb``
^^^^^^^^^^^^

**Description:** per site per generation recombination rate
                                (default=0)

**Type:** float; **Default:** 0.0



``--rho``
^^^^^^^^^

**Description:** specific rho = 4*Ne*mu instead of using
                                --recomb

**Type:** float; **Default:** None



``--window``
^^^^^^^^^^^^

**Description:** length (bp) of windows

**Type:** integer; **Default:** 100000


.. fasta2dfoil:

fasta2dfoil
===========

Description
-----------

DFOIL: Directional introgression testing a five-taxon phylogeny
James B. Pease
http://www.github.com/jbpease/dfoil

fasta2dfoil -
This script takes one or more FASTA files containing
5 or 4 taxa and counts site patterns for use in DFOIL/Dstat analysis.
To combine multiple FASTA files, each file should be sequences
from one locus (i.e., one entry in the final table) and
the names of sequences must be identical in all files.


Parameters
----------

fastafile
^^^^^^^^^

**Description:** one or more input fasta
                                files for each locus

**Type:** None; **Default:** None



``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--names/-n`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Order of the 5 (or 4) taxa,
                                names must be
                                consistent in all input files,
                                outgroup should be last

**Type:** None; **Default:** None



``--out/-o`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output count file, one entry per fasta

**Type:** None; **Default:** None


