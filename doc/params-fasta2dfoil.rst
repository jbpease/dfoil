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

Example Usage:

python3 fasta2dfoil.py INPUT.fasta --out OUTPUT.fasta         --names TAXA1,TAXA2,TAXA3,TAXA4



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

**Description:** Order of the 5 (or 4) taxa separated by commas.
                                Names must be  consistent in all input files,
                                outgroup should be last

**Type:** None; **Default:** None



``--out/-o`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output count file, one entry per fasta

**Type:** None; **Default:** None


