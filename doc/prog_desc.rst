Program Parameter Descriptions
##############################

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


