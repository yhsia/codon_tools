=====
Usage
=====

To execute Codon Harmony as a script::

    python codon_harmony/codon_harmony.py --input misc/INPUT_LIST.fasta --output out.fasta

To get started, create a conda environment from the ``environment.yml`` file::

    conda env create -f environment.yml

contents of ``misc/INPUT_LIST.fasta``:

.. code-block:: text

  >SEQ_1
  ACDEFGHIKLMNPQRSTVWY
  >SEQ_2
  ACDEFGHIKLMNPQRSTVWY

To use Codon Harmony in a project::

    import codon_harmony
