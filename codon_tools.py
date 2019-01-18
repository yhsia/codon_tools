#! /usr/bin/env python
import argparse
import random
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from codon_tools.util import codon_use, logging, log_levels, seq_opt
from codon_tools.data import GC_content, RibosomeBindingSites, RestrictionEnzymes


def get_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Reverse translate your amino acid sequence harmoniously with "
        + "a host's condon usage.",
        epilog="v0.9 (contact yhsia@uw.edu or bweitzner@lyellbio.com if you "
        + "encounter errors)",
    )
    parser.add_argument(
        "--input", type=str, required=True, help="input file with sequence"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="413997",
        help='host table code: http://www.kazusa.or.jp/codon/, default is "Escherichia coli B"',
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
        help="verbose output level (0=only result, 1=standard output, 2=extra output 3=debugging)",
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
        default=1000,
        help="max number of cycles to run optimization, 0=unlimited",
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

    args = parser.parse_args()

    return args


def main(argv):
    """Read in a fasta-formatted file containing amino acid sequences and
    reverse translate each of them in accordance with a specified host's
    codon usage frequency. The DNA sequence is then processed to remove
    unwanted features.

    Args:
        argv (list[str]): Command line arguments passed into the function.
    """
    args = get_args(argv)
    global logger
    logging.basicConfig(level=log_levels[args.verbose])
    logger = logging.getLogger(__name__)

    random.seed()
    logger.info("===== SCRIPT START =====")

    # generate host profile
    codon_use_table, host_profile, codoon_relative_adativeness = codon_use.host_codon_usage(
        args.host, args.host_threshold
    )

    # initialize the restriction sites of interest
    rest_enz = RestrictionEnzymes(args.restriction_enzymes)

    # process through all supplied sequences
    for seq_no, record in enumerate(SeqIO.parse(args.input, "fasta", IUPAC.protein)):
        logger.info(
            "===== PROCESSING SEQUENCE {0} ===== {1}".format(seq_no + 1, record.id)
        )

        # check input seq style
        logger.info("INPUT IS AA SEQUENCE")
        dna = record.seq.back_translate()

        logger.detail("===== DUMPING SEQUENCE =====")
        logger.detail(str(dna))

        # intialize bookkeeping variables
        difference, current_cycle, relax, best_cai, best_dna = 1.0, 1, 1.0, 0.0, ""

        # keep running while there are cycles AND difference between current and host is less than the % relax allowed
        while ((current_cycle <= args.cycles) or (args.cycles == 0)) and (
            difference >= (relax - 1)
        ):
            logger.info(
                "~~~~~~~~~~ Current cycle: {0}/{1} ~~~~~~~~~~".format(
                    current_cycle, args.cycles
                )
            )

            # relax harmonization requirement
            relax = 1 + (args.max_relax * ((current_cycle - 1) / args.cycles))
            dna = seq_opt.resample_codons_and_enforce_host_profile(
                dna, codon_use_table, host_profile, relax
            )

            # identify and remove undesirable features
            for gc_content in GC_content:
                # check various GC content requirements
                dna = seq_opt.gc_scan(dna, gc_content, codon_use_table)

            dna = seq_opt.remove_start_sites(dna, RibosomeBindingSites, codon_use_table)
            dna = seq_opt.remove_repeating_sequences(dna, 9, codon_use_table)
            dna = seq_opt.remove_local_homopolymers(
                dna,
                codon_use_table,
                n_codons=2,
                homopolymer_threshold=args.local_homopolymer_threshold,
            )
            if len(rest_enz):
                dna = seq_opt.remove_restriction_sites(dna, rest_enz, codon_use_table)

            # if the codon adaptation index is better than what we've
            # seen so far, store this sequence
            cai = codoon_relative_adativeness.cai_for_gene(str(dna))
            if cai > best_cai:
                best_dna = SeqRecord(
                    dna, id=record.id, name=record.name, description=record.description
                )

            # tick cycle
            current_cycle += 1

        # hit the max number of cycles?
        if current_cycle > args.cycles:
            logger.info("You hit the max number of cycles: {0}".format(args.cycles))

        # check GC content
        logger.info("===== GC CONTENT =====")
        gc_percent = round(GC(best_dna.seq) / 100, 3)
        if gc_percent < 0.3 or gc_percent > 0.65:
            logger.warning("Overall GC content is {0}!".format(gc_percent))

        # measure the final deviation from the host profile
        _, difference = seq_opt.compare_profiles(
            codon_use.count_codons(best_dna.seq), host_profile, relax
        )

        logger.output(
            "Final codon-use difference between host and current sequence "
            + "(0.00 is ideal): {0}".format(round(difference, 2))
        )

        logger.output("===== NAME AND SEQUENCE =====")
        print(best_dna.format("fasta"))


if __name__ == "__main__":
    main(sys.argv)
