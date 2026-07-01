"""Interleave forward and reverse reads to improve SPOA performance."""

from pathlib import Path
import sys

import pysam
import spoa

from .util import get_named_logger, wf_parser  # noqa: ABS101

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def sw_score(s1, s2):
    """Compute a simple Smith-Waterman-like score using Python.

    Uses a basic match/mismatch scoring without gaps for speed.
    Returns the number of matching positions when sequences are aligned
    by comparing same-length substrings.
    """
    # Use the shorter sequence
    l = min(len(s1), len(s2))
    s1 = s1[:l]
    s2 = s2[:l]
    return sum(a == b for a, b in zip(s1, s2))


def get_seqs_from_fastx_and_check_lengths(file, logger, max_len=None):
    """Get sequences from FASTx file while making sure they are not too long.

    SPOA was not intended for very long reads and will take an exorbitant amount
    of memory.

    :param file: FASTx file
    :param logger: logger for error message
    :param max_len: maximum allowed read length, defaults to None
    :yield: sequences from FASTx file
    """
    with pysam.FastxFile(file, "r") as f:
        for entry in f:
            seq = entry.sequence
            if max_len is not None and len(seq) > max_len:
                logger.error(
                    f"Tried to run SPOA with reads longer than {max_len}. "
                    "This is what assemblers are for. Aborting..."
                )
                sys.exit(0)
            yield seq


def interleave_lists(l1, l2):
    """Uniformly interleave two lists.

    :param l1: first list
    :param l2: second list
    :return: list containing items of ``l1`` and ``l2`` uniformly interleaved
    """
    interleaved = []
    rate_1 = 1 / len(l1)
    rate_2 = 1 / len(l2)

    c_1, c_2 = 0, 0

    itr_1 = iter(l1)
    itr_2 = iter(l2)

    while True:
        try:
            if c_1 > c_2:
                interleaved.append(next(itr_2))
                c_2 += rate_2
            else:
                interleaved.append(next(itr_1))
                c_1 += rate_1
        except StopIteration:
            break

    return interleaved


def main(args):
    """Run the entry point."""
    logger = get_named_logger("runSPOA")

    logger.info("Get read orientations...")
    seqs = get_seqs_from_fastx_and_check_lengths(
        args.fastq, logger, args.max_allowed_read_length
    )
    try:
        first_seq = next(seqs)
    except StopIteration:
        logger.error(f"Input file '{args.fastq}' appears to be empty. Aborting...")
        sys.exit(0)

    fwd = []
    rev = []
    for seq in seqs:
        rc = reverse_complement(seq)
        fwd_score = sw_score(seq, first_seq)
        rev_score = sw_score(rc, first_seq)
        if fwd_score >= rev_score:
            fwd.append(seq)
        else:
            rev.append(rc)

    if fwd and rev:
        interleaved_reads = interleave_lists(fwd, rev)
        logger.info("Finished interleaving reads.")
    else:
        logger.info("Only got reads from one strand; not interleaving...")
        interleaved_reads = fwd or rev

    min_cov = None
    if args.relative_min_coverage is not None:
        min_cov = int(round(len(interleaved_reads) * args.relative_min_coverage))
        logger.info(f"SPOA min coverage: {min_cov}.")

    cons, _ = spoa.poa(interleaved_reads, genmsa=False, min_coverage=min_cov)
    cons, _ = spoa.poa([cons, *interleaved_reads], genmsa=False, min_coverage=min_cov)

    with open(args.output, "w") as outfile:
        outfile.write(f">consensus\n{cons}\n")

    logger.info("Wrote consensus to output file.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("run_spoa")
    parser.add_argument("fastq", type=Path, help="Path to input FASTQ file")
    parser.add_argument(
        "-o", required=True, dest="output", help="Output FASTA file name"
    )
    parser.add_argument(
        "--relative-min-coverage",
        type=float,
        help="Minimum relative coverage of POA graph",
    )
    parser.add_argument(
        "--max-allowed-read-length",
        type=int,
        help=(
            "SPOA can't deal with very long reads; "
            "don't run it when encountering reads longer than this"
        ),
    )
    return parser
