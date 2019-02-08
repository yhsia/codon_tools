=====
Usage
=====

.. argparse::
   :module: codon_harmony.codon_harmony
   :func: get_parser
   :prog: codon_harmony

Executing Codon Harmony as a script
-----------------------------------

    python codon_harmony/codon_harmony.py --input misc/INPUT_LIST.fasta --output out.fasta

To get started, create a conda environment from the ``environment.yml`` file::

    conda env create -f environment.yml

contents of ``misc/INPUT_LIST.fasta``:

.. code-block:: text

  >SEQ_1
  ACDEFGHIKLMNPQRSTVWY
  >SEQ_2
  ACDEFGHIKLMNPQRSTVWY

Using Codon Harmony in a project
--------------------------------

.. code-block:: python

  import codon_harmony
  codon_harmony.runner()

The ``runner`` function will handle parsing all command line arguments.
