#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_tools.util` module."""


import unittest

import codon_harmony.util as chu


class TestCodon_tools_data(unittest.TestCase):
    """Tests for `codon_tools.data` module."""

    def setUp(self):
        """Set up test fixtures, if any."""
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        self.test_aa = Seq("TESTTESTTEST", IUPAC.protein)
        self.test_dna = Seq(
            "ACCGAGTCTACCACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna
        )

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_003_back_translate(self):
        """Test `codon_harmony.util.seq` -- unoptimized reverse translation"""
        assert self.test_aa.back_translate() == self.test_dna

    def test_004_codon_count(self):
        known_counts = {"ACC": 6, "GAG": 3, "TCT": 3}
        cdn_count = chu.codon_use.count_codons(self.test_dna)
        for codon, count in cdn_count.items():
            if codon in known_counts:
                assert known_counts[codon] == count
            else:
                assert not count
