===========
Codon Tools
===========


.. image:: https://img.shields.io/pypi/v/codon_tools.svg
        :target: https://pypi.python.org/pypi/codon_tools

.. image:: https://img.shields.io/travis/weitzner/codon_tools.svg
        :target: https://travis-ci.org/weitzner/codon_tools

.. image:: https://readthedocs.org/projects/codon-tools/badge/?version=latest
        :target: https://codon-tools.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/weitzner/codon_tools/shield.svg
     :target: https://pyup.io/repos/github/weitzner/codon_tools/
     :alt: Updates



Amino acid reverse translation and DNA optimization tool based on species-specific codon-use distributions.
Species-specifc data can be found on the [Codon Usage Database](http://www.kazusa.or.jp) using the [NCBI Taxonomy database](http://www.ncbi.nlm.nih.gov/taxonomy) id (e.g. 413997) or the organism's Latin name (e.g. _Escherichia coli_ B). Mapping species names to Taxonomy IDs can be done [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi).


* Free software: MIT license
* Documentation: https://codon-tools.readthedocs.io.


Use
---

```sh
$python codon_tools.py --input misc/INPUT_LIST.fasta --output out.fasta
```

To get started, create a conda environment from the `environment.yml` file.

```sh
conda env create -f environment.yml
```

contents of `misc/INPUT_LIST.fasta`:

```
>SEQ_1
ACDEFGHIKLMNPQRSTVWY
>SEQ_2
ACDEFGHIKLMNPQRSTVWY
```

Features
--------

1. Reverse translates input amino acid sequence to DNA.
2. Calculates the host's per-AA codon usage profile – codons used less than a specified threshold (defaults to 10%) are dropped.
3. Compares the reverse-translated DNA sequence to the host profile, determines which codons are overused/underused.
4. Stochastically mutates codons according to host profile.
5. Processes DNA to remove unwanted features:
    * high GC content within a sliding window and across the entire sequence
    * unwanted restriction sites
    * alternate start positions (GA-rich regions 18 bp upstream of ATG/GTG/TTG)
    * 3-consecutive identical codons and 9-mer repeat chunks
    * areas with more than 4 (variable) consecutive identical bps ("local homopolymers")

The process is repeated from step 3 for a specified number of cycles (defaults to 1000) OR until the per-AA codon profile of current DNA and host profile matches (within tolerance).

To do
-----

- [x] remove RNA structure from sequence
    * ~[CONTRAfold](http://contra.stanford.edu/contrafold/)~ overkill for now
    * ~[nupack](http://nupack.org)~ overkill for now
    * restrict structure to hairpins, detected by looking for 10-mers with reverse complements (including wobble bases) in the sequence
- [ ] remove predicted splice sites
- [x] store "best" sequence by codon adaptation index relative to host and return that to the user

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
