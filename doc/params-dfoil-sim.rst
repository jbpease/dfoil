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


