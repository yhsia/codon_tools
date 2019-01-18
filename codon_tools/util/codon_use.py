from Bio.SeqUtils import CodonUsage

from . import Seq, logging
from ..data.codon_use_by_host import codon_tables, latin_name_to_tax_id

logger = logging.getLogger(__name__)


def count_codons(dna_sequence):
    """Count the number of times each codon appears in a DNA sequence.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.

    Returns:
        dict({str: int}): A dictionary with codons as keys and the
            corresponding number of occurences as values.
    """
    logger.info("===== COUNTING CODONS =====")

    codons_dict = CodonUsage.CodonsDict.copy()
    for codon_start in range(0, len(dna_sequence), 3):
        codon_idx = slice(codon_start, codon_start + 3)
        codons_dict[str(dna_sequence[codon_idx])] += 1

    return codons_dict


def calc_profile(codons_count):
    """Calculate the frequency of usage of each synonymous codon from an
    input dictionary of counts.

    Args:
        codons_count (dict({str: int})): A dictionary with codons as keys
            and the corresponding number of occurences as values.

    Returns:
        dict({str: int}): A dictionary with codons as keys and the
            corresponding frequency of occurences as values.
    """
    logger.info("===== CALCULATING PROFILE =====")
    codons_freq = CodonUsage.CodonsDict.copy()

    # loop through all amino acids
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # add up number of times each AA appears
        tot_usage = sum([codons_count[codon] for codon in synonymous_codons])

        # calculate the distribution of synonymous codon use
        if tot_usage == 0:
            continue
        codons_freq.update(
            {codon: codons_count[codon] / tot_usage for codon in synonymous_codons}
        )

    return codons_freq


def calc_codon_relative_adaptiveness(codons_count):
    """Calculate the relative adaptiveness of each synonymous codon from an
    input dictionary of counts.

    Note:
        The claculation and some nomenclature is taken from Sharp and
        Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).

    Args:
        codons_count (dict({str: int})): A dictionary with codons as keys
            and the corresponding number of occurences as values.

    Returns:
        Bio.SeqUtils.CodonUsage.CodonAdaptationIndex: A CodonAdaptationIndex
            instance configured to calculate CAI for a target gene.
    """
    logger.info("===== CALCULATING RELATIVE ADAPTIVENESS =====")
    codons_rel_adapt = CodonUsage.CodonsDict.copy()

    # loop through all amino acids
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # get the number of occurrences of the most frequently used codon
        # `X_max` is the notation from Sharp & Li, Nucleic Acids Res. 1987
        X_max = max([codons_count[codon] for codon in synonymous_codons])

        # calculate the relative adaptiveness of each synonymous codon
        if X_max == 0:
            continue
        codons_rel_adapt.update(
            {codon: codons_count[codon] / X_max for codon in synonymous_codons}
        )

    codon_adaptation_index = CodonUsage.CodonAdaptationIndex()
    codon_adaptation_index.set_cai_index(codons_rel_adapt)

    return codon_adaptation_index


def _load_host_table(host):
    logger.info("===== PROCESSING HOST TABLE: {0} =====".format(host))
    table = CodonUsage.CodonsDict.copy()

    if host not in codon_tables and host in latin_name_to_tax_id:
        host = latin_name_to_tax_id[host]

    try:
        raw_table = codon_tables[host]
    except KeyError:
        raise ValueError(
            '"{}" is not a valid host id. '.format(host)
            + "Currently supported hosts (Latin and NCBI taxonomy IDs) are {}".format(
                host, ", ".join(latin_name_to_tax_id.keys() + codon_tables.keys())
            )
        )
    else:
        for codon, usage in raw_table.items():
            # codons are stored with RNA alphabet
            table[str(Seq(codon).back_transcribe())] = usage["count"]

    return calc_profile(table)


def process_host_table(host, threshold):
    """Load the codon usage table for the desired host, filter codons with
    a lower occurence than the threshold, and renormalize the frequency of
    usage of each synonymous codon.

    Args:
        host (str): Latin name or NCBI taxonomy ID of the host organism.
        threshold (float): Lowest fraction of codon usage to keep.

    Returns:
        dict({str: int}): A dictionary with codons as keys and the
            corresponding frequency of occurences as values.
    """
    table = _load_host_table(host)

    logger.info("HOST THRESHOLD SET TO: {0}".format(threshold))
    logger.detail("pre-threshold host table:")

    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        logger.detail("== {0} ==".format(AA))
        for codon in synonymous_codons:
            logger.detail("{0}: {1}".format(codon, table[codon]))
            if table[codon] < threshold:
                table[codon] = 0

    # recalculate profile after dropping rare codons
    table = calc_profile(table)

    if logger.isEnabledFor(logging.DETAIL):
        logger.detail("post-threshold host table:")
        for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
            logger.detail("== {0} ==".format(AA))
            for codon in synonymous_codons:
                logger.detail("{0}: {1}".format(codon, table[codon]))

    return table


def host_codon_usage(host, threshold=0.10):
    """Load and process the per amino acid codon usage for the desired host in
    accordance with the supplied threshold and configure a CodonAdaptationIndex
    instance to calculate CAI for a target gene.

    Note:
        The relative adaptiveness used in the CodonAdaptationIndex is based
        on the filtered codon use frequencies, not the raw counts.

    Args:
        host (str): Latin name or NCBI taxonomy ID of the host organism.
        threshold (float, optional): Lowest fraction of codon usage to keep.
            Defaults to 0.10.

    Returns:
        dict({str : list[list, list]}): A dictionary with each amino acid
            three-letter code as keys, and a list of two lists as values.
            The first list is the synonymous codons that encode the amino
            acid, the second is the frequency with which each synonymous
            codon is used.
        dict({str: int}): A dictionary with codons as keys and the
            corresponding frequency of occurences as values.
        Bio.SeqUtils.CodonUsage.CodonAdaptationIndex: A `CodonAdaptationIndex`
            instance configured to calculate CAI for a target gene.
    """
    host_profile = process_host_table(host, threshold)
    cra = calc_codon_relative_adaptiveness(host_profile)

    codon_use_by_aa = {}
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        codon_use_by_aa[AA] = [
            synonymous_codons,
            [host_profile[codon] for codon in synonymous_codons],
        ]

    return codon_use_by_aa, host_profile, cra
