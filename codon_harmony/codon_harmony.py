#! /usr/bin/env python
import argparse
import random
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from .util import codon_use, logging, log_levels, seq_opt
from .data import GC_content, RibosomeBindingSites, RestrictionEnzymes


def get_parser():
    parser = argparse.ArgumentParser(
        description="Reverse translate your amino acid sequence harmoniously with "
        + "a host's condon usage.",
        epilog="v0.9.4 (contact yhsia@uw.edu or bweitzner@lyellbio.com if you "
        + "encounter errors)",
    )
    parser.add_argument(
        "--input", type=str, required=True, help="input file with sequence"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="output file to write DNA sequence(s)"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="413997",
        help="host table code: http://www.kazusa.or.jp/codon/, default is "
        + '"Escherichia coli B"',
    )
    parser.add_argument(
        "--host_threshold",
        type=float,
        default="0.10",
        help="lowest codon fraction per AA in the host that is allowed",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="verbose output level (0=only result, 1=standard output, "
        + "2=extra output 3=debugging)",
    )
    parser.add_argument(
        "--local_homopolymer_threshold",
        type=int,
        default="4",
        help="number of consecutive NT repeats allowed",
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=100,
        help="number of independent codon samples to run. 0 means 1 pass",
    )
    parser.add_argument(
        "--inner_cycles",
        type=int,
        default=10,
        help="number of times to iteratively optimize each independent codon"
        + " sample. 0 means 1 pass",
    )
    parser.add_argument(
        "--max_relax",
        type=float,
        default="0.1",
        help="maximum percent deviation from host profile",
    )
    parser.add_argument(
        "--restriction_enzymes",
        nargs="*",
        type=str,
        default="NdeI XhoI HpaI PstI EcoRV NcoI BamHI".split(),
        help="list of restriction enzyme sites to remove "
        + "(e.g. --restriction_enzymes NdeI XhoI HpaI). ",
    )

    return parser


def main(argv=None):
    """Read in a fasta-formatted file containing amino acid sequences and
    reverse translate each of them in accordance with a specified host's
    codon usage frequency. The DNA sequence is then processed to remove
    unwanted features.
    """
    args = get_parser().parse_args(argv)
    logging.basicConfig(level=log_levels[args.verbose])
    logger = logging.getLogger(__name__)

    random.seed()
    logger.info("Beginning codon use optimization")

    # stores the final sequences
    out_seqs = []

    # generate host profile
    codon_use_table, host_profile, codon_relative_adativeness = codon_use.host_codon_usage(
        args.host, args.host_threshold
    )

    # initialize the restriction sites of interest
    rest_enz = RestrictionEnzymes(args.restriction_enzymes)

    # process through all supplied sequences
    for seq_no, record in enumerate(SeqIO.parse(args.input, "fasta", IUPAC.protein)):
        logger.info(
            "Processing sequence number {}:\n{}".format(
                seq_no + 1, record.format("fasta")
            )
        )

        dna = record.seq.back_translate()
        logger.detail(
            "Initial DNA sequence:\n{}".format(
                SeqRecord(dna, id=record.id).format("fasta")
            )
        )

        # intialize bookkeeping variables
        best_cai, best_dna = 0.0, ""

        args.cycles = 1 if args.cycles == 0 else args.cycles
        args.inner_cycles = 1 if args.inner_cycles == 0 else args.inner_cycles

        # run `args.cycles` independent trials
        for sample_no in range(args.cycles):
            logger.info("Current sample no: {}/{}".format(sample_no + 1, args.cycles))

            # relax harmonization requirement
            relax = 1 + (args.max_relax * ((sample_no + 1) / (args.cycles)))
            dna = seq_opt.resample_codons_and_enforce_host_profile(
                dna, codon_use_table, host_profile, relax
            )

            # go through a few cycles with the same starting sequence to
            # allow iterative improvements to the same sample of codons
            for _ in range(args.inner_cycles):
                # identify and remove undesirable features
                for gc_content in GC_content:
                    # check various GC content requirements
                    dna = seq_opt.gc_scan(dna, codon_use_table, gc_content)

                dna = seq_opt.remove_start_sites(
                    dna, codon_use_table, RibosomeBindingSites
                )
                dna = seq_opt.remove_repeating_sequences(dna, codon_use_table, 9)
                dna = seq_opt.remove_local_homopolymers(
                    dna,
                    codon_use_table,
                    n_codons=2,
                    homopolymer_threshold=args.local_homopolymer_threshold,
                )
                dna = seq_opt.remove_hairpins(dna, codon_use_table, stem_length=10)
                if len(rest_enz):
                    dna = seq_opt.remove_restriction_sites(
                        dna, codon_use_table, rest_enz
                    )

            # measure the deviation from the host profile post-cleanup
            # only move forward if we haven't deviated too much from host
            _, difference = seq_opt.compare_profiles(
                codon_use.count_codons(dna), host_profile, relax
            )
            if difference >= args.max_relax:
                continue

            # if the codon adaptation index is better than what we've
            # seen so far, store this sequence
            cai = codon_relative_adativeness.cai_for_gene(str(dna))
            if cai > best_cai:
                best_cai = cai
                best_dna = SeqRecord(
                    dna, id=record.id, name=record.name, description=record.description
                )

        logger.info(
            "Completed {} independent codon samples with optimization!".format(
                args.cycles
            )
        )

        if isinstance(best_dna, str):
            logger.warning(
                "Unable to create suitable DNA sequence for input sequence {}.\n{}".format(
                    seq_no + 1, record.format("fasta")
                )
            )
            continue

        logger.output("Optimized gene metrics and sequence")
        # check GC content
        gc_frac = GC(best_dna.seq) / 100
        logger.output("Final overall GC content is {:.0%}".format(gc_frac))
        if gc_frac < 0.3 or gc_frac > 0.65:
            logger.warning(
                "The sequence's GC content ({:.2f}) is beyond normal ranges (0.3 > GC < 0.65)!".format(
                    gc_frac
                )
            )

        # measure the final deviation from the host profile
        _, difference = seq_opt.compare_profiles(
            codon_use.count_codons(best_dna.seq), host_profile, relax
        )

        logger.output(
            "Final codon-use difference between host and current sequence: {:.2f}".format(
                difference
            )
        )

        out_seqs.append(best_dna.format("fasta"))
        logger.output(
            "The designed gene's CAI is: {:.2f}\n{}".format(best_cai, out_seqs[-1])
        )

    # write sequences to file
    with open(args.output, "w") as f:
        f.write("".join(out_seqs))
