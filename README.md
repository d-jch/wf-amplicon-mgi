# wf-amplicon-mgi

**MGI long-read haploid amplicon de novo consensus workflow**

A Nextflow DSL2 workflow for generating de novo consensus sequences from MGI
long-read amplicon data. The workflow assembles a consensus per sample using
miniasm + racon, falling back to SPOA for difficult samples, and applies
depth-based QC and trimming to produce a final high-confidence FASTA.

---

## Overview

```
Input FASTQ
  │
  ├─ Length / quality filter (seqkit)
  ├─ Subsample reads (by length or randomly)
  │
  ├─ minimap2 all-vs-all overlap → miniasm draft
  │     ├─ [success] → racon polish
  │     └─ [failed / short contig] → SPOA fallback
  │
  ├─ Re-align reads to candidate consensus (minimap2)
  ├─ mosdepth per-base depth
  │
  ├─ QC: filter by mean depth & primary alignment fraction
  ├─ Select contig with highest mean depth
  ├─ Trim low-coverage ends
  │
  └─ Final consensus FASTA + aligned BAM
```

---

## Requirements

The following tools must be available in the execution environment:

| Tool | Purpose |
|------|---------|
| `minimap2` | Read overlap and consensus alignment |
| `miniasm` | All-vs-all assembly |
| `racon` | One-round polishing |
| `spoa` (+ Python bindings) | Fallback POA consensus |
| `samtools` | BAM/FASTA indexing and manipulation |
| `mosdepth` | Per-base depth calculation |
| `seqkit` | Read filtering and statistics |
| `bamstats` | Alignment statistics |
| `bgzip` | FASTQ compression |
| Python ≥ 3.8 with `pysam`, `pandas` | Python helper scripts |

When using Docker/Singularity, a suitable container is specified in
`nextflow.config`. Adjust the container image as needed.

---

## Quick start

```bash
# Clone the repository
git clone https://github.com/d-jch/wf-amplicon-mgi.git
cd wf-amplicon-mgi

# Single FASTQ file
nextflow run main.nf \
    --fastq /path/to/sample.fastq.gz \
    --out_dir results \
    -profile standard

# Directory of FASTQ files (treated as one sample)
nextflow run main.nf \
    --fastq /path/to/fastq_dir/ \
    --sample MySample \
    --out_dir results

# Barcode sub-directory layout (one sample per barcode)
nextflow run main.nf \
    --fastq /path/to/run_dir/ \
    --out_dir results
```

---

## Parameters

### I/O

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq` | *(required)* | FASTQ file, flat directory, or barcode directory |
| `--out_dir` | `output` | Output directory |
| `--sample` | *(auto)* | Override sample name (single-sample runs) |

### Read filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_read_length` | `200` | Discard reads shorter than this (bp) |
| `--max_read_length` | *(off)* | Discard reads longer than this (bp) |
| `--min_read_qual` | *(off)* | Minimum mean base quality |
| `--min_n_reads` | `20` | Skip samples with fewer reads after filtering |

### Read subsampling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads_downsampling_size` | `1500` | Target read count for assembly |
| `--drop_frac_longest_reads` | `0.0` | Drop this fraction of longest reads (e.g. `0.05`) |
| `--take_longest_remaining_reads` | `true` | Take the longest reads from the remainder |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--overlap_preset` | `ava-pb` | minimap2 preset for all-vs-all overlap |
| `--minimap2_preset` | `map-pb` | minimap2 preset for read → consensus alignment |
| `--force_spoa_length_threshold` | `2000` | Fall back to SPOA if longest miniasm contig < this |

### SPOA (fallback)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--spoa_minimum_relative_coverage` | `0.15` | Minimum relative POA graph coverage |
| `--spoa_max_allowed_read_length` | `5000` | Abort SPOA if any read exceeds this length |

### QC thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--minimum_mean_depth` | `30` | Minimum mean depth to accept a consensus |
| `--primary_alignments_threshold` | `0.7` | Minimum primary alignment fraction |
| `--number_depth_windows` | `100` | Windows for depth profile in QC output |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `4` | CPU threads per process |

---

## Output

Results are written to `--out_dir` (default: `output/`). The directory structure is:

```
output/
├── <sample>/
│   ├── consensus/
│   │   └── consensus.fasta       # final trimmed consensus sequence
│   ├── alignments/
│   │   ├── <sample>.aligned.sorted.bam
│   │   └── <sample>.aligned.sorted.bam.bai
│   └── qc/
│       ├── per-window-depth.tsv.gz
│       └── qc-summary.tsv
├── versions.txt
├── params.json
└── execution/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

`qc-summary.tsv` contains per-contig QC metrics including mean depth, primary
alignment fraction, and pass/fail status.

---

## minimap2 presets

Choose `--overlap_preset` and `--minimap2_preset` based on your MGI long-read
platform:

| Platform | Overlap preset | Alignment preset |
|----------|---------------|-----------------|
| MGI CycloneSEQ / PacBio CLR-like | `ava-pb` | `map-pb` |
| MGI with ONT-like error profile | `ava-ont` | `map-ont` |
| PacBio HiFi | `ava-hifi` | `map-hifi` |

---

## How it works

### 1. Read ingress and filtering

Reads are collected from the input path, filtered by length and quality using
`seqkit`, and optionally subsampled.

### 2. miniasm assembly

`minimap2` is run in all-vs-all overlap mode, piped into `miniasm` to produce
a GFA assembly. If the longest contig is shorter than
`--force_spoa_length_threshold`, the sample is routed to the SPOA fallback.

### 3. racon polishing

For successful miniasm assemblies, one round of `racon` polishing is performed
using the original reads.

### 4. SPOA fallback

For samples where miniasm fails, `spoa` is used to compute a partial order
alignment (POA) consensus. Reads are interleaved forward/reverse before being
passed to SPOA.

### 5. QC and trimming

Reads are re-aligned to the candidate consensus. `mosdepth` computes per-base
depth, and the consensus is trimmed at the ends where depth drops below the
threshold. Contigs failing mean-depth or primary-alignment-fraction thresholds
are discarded.

---

## License

MIT

---

## Acknowledgements

Inspired by [epi2me-labs/wf-amplicon](https://github.com/epi2me-labs/wf-amplicon).
This is an independent project focused on MGI long-read data and removes all
ONT-specific dependencies (Medaka, Porechop, basecaller model handling).
