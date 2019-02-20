#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony` package."""


import unittest

from codon_harmony import codon_harmony


class TestCodon_tools(unittest.TestCase):
    """Tests for `codon_tools` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_args(self):
        """Test args parsing."""
        args_to_parse = ["--input", "dummy.fasta", "--output", "out.fasta"]
        parser = codon_harmony.get_parser()
        parsed_args = parser.parse_args(args_to_parse)
        print()
        print(parsed_args)
        assert parsed_args.input == "dummy.fasta"
        assert parsed_args.output == "out.fasta"
        assert parsed_args.cycles == 100
        assert parsed_args.host == "413997"
        assert parsed_args.host_threshold == 0.1
        assert parsed_args.inner_cycles == 10
        assert parsed_args.local_homopolymer_threshold == 4
        assert parsed_args.max_relax == 0.1
        assert parsed_args.restriction_enzymes == [
            "NdeI",
            "XhoI",
            "HpaI",
            "PstI",
            "EcoRV",
            "NcoI",
            "BamHI",
        ]
        assert parsed_args.verbose == 0
