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


