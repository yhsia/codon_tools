#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony` package."""


import unittest

from codon_harmony import codon_harmony


class TestCodon_tools(unittest.TestCase):
    """Tests for `codon_tools` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.args_to_parse = [
            "--input",
            "misc/INPUT_LIST.fasta",
            "--output",
            "out.fasta",
        ]

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_args(self):
        """Test args parsing."""
        parser = codon_harmony.get_parser()
        parsed_args = parser.parse_args(self.args_to_parse)
        assert parsed_args.input == "misc/INPUT_LIST.fasta"
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

    def test_main(self):
        """Test `codon_harmony` main function"""
        # all this does is make sure that nothing crashes --
        # more of an integration test than a unit test with
        # each unit being separately tested.

        # test three sequences:
        #   one that can be optimized
        #   one that cannot be optimized (given max_relax = 0.1)
        #   one that triggers a warning about GC content
        codon_harmony.main(self.args_to_parse)
        assert True
