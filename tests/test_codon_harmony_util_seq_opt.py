#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony.util.seq` module."""


import unittest
from codon_harmony.util import seq_opt


class TestCodon_harmony_util_seq_opt(unittest.TestCase):
    """Tests for `codon_harmony.util.seq` module."""

    def setUp(self):
        """Set up test fixtures, if any."""
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        self.test_aa = Seq("TESTTESTTEST", IUPAC.protein)
        self.test_dna = Seq(
            "ACCGAGTCTACCACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna
        )

        from codon_harmony.util import codon_use

        self.codon_use_table, self.host_profile, _ = codon_use.host_codon_usage(
            9606, 0.0
        )

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_mutate_codon(self):
        """Test `codon_harmony.util.seq` -- unoptimized reverse translation"""
        # Phe -- two codons
        codon_in = "TTT"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out == "TTC"

        # Leu -- six codns
        codon_in = "TTA"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out in ["TTG", "CTT", "CTC", "CTA", "CTG"]

        # Met -- one codon
        codon_in = "ATG"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out == codon_in
