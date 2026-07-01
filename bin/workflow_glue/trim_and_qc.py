"""Select a contig if multiple were created, do some QC, trim low-coverage ends."""
from pathlib import Path
import sys

import pandas as pd
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("TrimAndQC")

    logger.info("Read input files.")

    flagstat = pd.read_csv(args.flagstat, sep="\t", index_col=0)

    depths = pd.read_csv(
        args.depth, sep="\t", header=None, names=["ref", "start", "end", "depth"]
    )

    with pysam.FastxFile(args.consensus) as f:
        fastx = {entry.name: entry for entry in f}

    seq_ids = [x for x in flagstat.index if x != "*"]
    if set(fastx.keys()) != set(seq_ids):
        raise ValueError("Sequence IDs in consensus file and flagstat file don't agree.")

    # if all reads in the BAM are unmapped mosdepth outputs an empty per-base depth file
    if depths.empty:
        mean_depths = pd.Series(0, index=seq_ids)
    else:
        mean_depths = depths.groupby("ref").apply(
            lambda df: df.eval("(end - start) * depth").sum() / df["end"].iloc[-1]
        )

    # perform QC checks on each contig
    qc_stats = pd.DataFrame(
        columns=[
            "method",
            "length",
            "mean_depth",
            "primary",
            "secondary",
            "supplementary",
            "primary_ratio",
            "status",
            "fail_reason",
        ]
    )
    qc_stats.index.name = "contig"
    for contig, mean_depth in mean_depths.items():
        status = "passed"
        contig_len = len(fastx[contig].sequence)
        fail_reasons = []
        if mean_depth < args.minimum_depth:
            status = "failed"
            fail_reasons.append("low depth")
        primary_ratio = flagstat.loc[contig, "primary"] / flagstat.loc[contig, "total"]
        if primary_ratio < args.primary_threshold:
            status = "failed"
            fail_reasons.append("low primary ratio")

        qc_stats.loc[contig] = [
            args.asm_method,
            contig_len,
            round(mean_depth, 3),
            *flagstat.loc[contig, ["primary", "secondary", "supplementary"]],
            round(primary_ratio, 3),
            status,
            ", ".join(fail_reasons),
        ]

    qc_passed = qc_stats.query("status == 'passed'")
    failed = qc_passed.empty

    outdir = Path(args.outdir_fail if failed else args.outdir_pass)
    outdir.mkdir(parents=True, exist_ok=True)

    qc_stats.to_csv(outdir / args.qc_summary_tsv, sep="\t")

    if failed:
        logger.warning(
            "Consensus generation failed for this sample. No consensus candidate "
            "passed the quality checks."
        )
        sys.exit()

    # select the contig with the largest mean depth
    selected = qc_passed["mean_depth"].sort_values(ascending=False).index[0]
    logger.info(f"Selected contig '{selected}'.")

    depths.query("ref == @selected", inplace=True)

    # determine trim limits based on depth threshold
    above_threshold = depths.query(
        "depth > depth.max() * @args.relative_depth_trim_threshold"
    )
    if above_threshold.empty:
        trim_start = 0
        trim_end = depths.iloc[-1]["end"]
    else:
        first_and_last_idx = above_threshold.index[[0, -1]]
        depths = depths.loc[slice(*first_and_last_idx)]
        trim_start = int(depths.iloc[0]["start"])
        trim_end = int(depths.iloc[-1]["end"])

    # extract and write the selected consensus sequence
    logger.info("Extract and write selected consensus sequence.")
    with pysam.FastxFile(args.consensus) as f:
        for entry in f:
            if entry.name == selected:
                break
        else:
            raise ValueError(
                f"Selected contig '{selected}' not found in '{args.consensus}'"
            )

    entry.name = args.alias
    entry.sequence = entry.sequence[trim_start:trim_end]
    # handle FASTA (no quality) or FASTQ (has quality)
    if entry.quality is not None:
        entry.quality = entry.quality[trim_start:trim_end]

    out_name = Path(args.consensus).name
    with open(outdir / out_name, "w") as f:
        f.write(str(entry) + "\n")

    logger.info("Done")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("trim_and_qc")
    parser.add_argument("--alias", type=str, required=True, help="Sample alias.")
    parser.add_argument(
        "--asm-method", type=str, required=True,
        help="Assembly method (added to QC summary TSV)."
    )
    parser.add_argument(
        "--depth", type=Path, required=True,
        help="Path to mosdepth per-base BED file."
    )
    parser.add_argument(
        "--flagstat", type=Path, required=True,
        help="Path to bamstats flagstat TSV file."
    )
    parser.add_argument(
        "--consensus", type=Path, required=True,
        help="Path to consensus FASTA or FASTQ file."
    )
    parser.add_argument(
        "--outdir-pass", type=Path, required=True,
        help="Output directory if QC passes."
    )
    parser.add_argument(
        "--outdir-fail", type=Path, required=True,
        help="Output directory if QC fails."
    )
    parser.add_argument(
        "--relative-depth-trim-threshold", type=float, required=True,
        help="Trim ends where relative coverage drops below this fraction."
    )
    parser.add_argument(
        "--minimum-depth", type=int, required=True,
        help="Minimum mean depth required to trust a consensus sequence."
    )
    parser.add_argument(
        "--primary-threshold", type=float, required=True,
        help="Minimum fraction of primary alignments required."
    )
    parser.add_argument(
        "--qc-summary-tsv", type=str, required=True,
        help="Name of TSV file for QC summary."
    )
    return parser
