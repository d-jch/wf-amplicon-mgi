#!/usr/bin/env nextflow
// wf-amplicon-mgi: MGI long-read haploid amplicon de novo consensus workflow

nextflow.enable.dsl = 2

include { catFastqDir } from "./lib/ingress"
include { getParams } from "./lib/common"
include {
    pipeline as deNovoPipeline_asm;
    pipeline as deNovoPipeline_spoa;
} from "./modules/local/de-novo"


process getVersions {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    output: path "versions.txt"
    script:
    """
    python --version | tr ' ' ',' | sed 's/P/p/' > versions.txt
    seqkit version | tr ' ' ',' >> versions.txt
    samtools --version | head -n1 | tr ' ' ',' >> versions.txt
    printf "minimap2,%s\\n" \$(minimap2 --version) >> versions.txt
    mosdepth --version | tr ' ' ',' >> versions.txt
    printf "miniasm,%s\\n" \$(miniasm -V) >> versions.txt
    printf "racon,%s\\n" \$(racon --version) >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import spoa; print(f'spoa,{spoa.__version__}')" >> versions.txt 2>/dev/null || true
    """
}


process downsampleReads {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    memory "2 GB"
    input:
        tuple val(meta), path("reads.fastq.gz")
        val n_reads
    output: tuple val(meta), path("downsampled.fastq.gz")
    script:
    int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    seqkit sample -2 reads.fastq.gz -n $n_reads \
    | bgzip -@ $bgzip_threads > downsampled.fastq.gz
    """
}


process filterReads {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    memory "2 GB"
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta), path("filtered.fastq.gz"), env(N_READS)
    script:
    String len_args = ""
    if (params.min_read_length) { len_args += " -m ${params.min_read_length}" }
    if (params.max_read_length) { len_args += " -M ${params.max_read_length}" }
    String qual_arg = params.min_read_qual ? " -Q ${params.min_read_qual}" : ""
    """
    seqkit seq ${len_args}${qual_arg} reads.fastq.gz | bgzip > filtered.fastq.gz
    N_READS=\$(seqkit stats -T filtered.fastq.gz | tail -1 | cut -f4)
    """
}


/*
Subset reads by length: optionally drop the longest fraction and then
keep the longest of the remaining reads. This helps avoid concatemers
and anomalously long reads.
*/
process subsetReads {
    label "wfamplicon"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta), path("reads.fastq.gz")
        val drop_longest_frac
        val take_longest
        val n_reads
    output: tuple val(meta), path("subset.fastq.gz")
    script:
    // Sort by length descending; optionally drop longest fraction; take top-N
    int bgzip_threads = Math.min(params.threads, 2)
    """
    # Sort reads by length (longest first)
    seqkit sort -lr reads.fastq.gz > sorted_reads.fastq

    total=\$(seqkit stats -T sorted_reads.fastq | tail -1 | cut -f4)

    # Calculate how many reads to drop from the front (longest reads)
    drop_n=0
    if awk "BEGIN { exit !($drop_longest_frac > 0) }"; then
        drop_n=\$(python3 -c "print(int(\$total * $drop_longest_frac))")
    fi

    keep_start=\$(( drop_n + 1 ))

    # Extract remaining reads (skip dropped ones) then select n_reads
    if [ "${take_longest}" = "true" ]; then
        # Take longest of remaining
        seqkit range -r \${keep_start}:-1 sorted_reads.fastq \
            | seqkit head -n $n_reads \
            | bgzip -@ $bgzip_threads > subset.fastq.gz
    else
        # Random subsample of remaining
        seqkit range -r \${keep_start}:-1 sorted_reads.fastq \
            | seqkit sample -n $n_reads -2 \
            | bgzip -@ $bgzip_threads > subset.fastq.gz
    fi
    """
}


process concatQCSummaries {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("input/f*.tsv")
    output: tuple val(meta), path("qc-summary.tsv")
    script:
    """
    head -1 \$(ls input/f*.tsv | head -1) > qc-summary.tsv
    tail -n +2 -q input/f*.tsv >> qc-summary.tsv
    """
}


process publish {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}


// Main workflow
workflow {

    // Validate required params
    if (!params.fastq) {
        error "Please provide --fastq pointing to a FASTQ file or directory."
    }

    def workflow_params = getParams()
    def software_versions = getVersions()

    // Ingest reads
    // Detect input layout: single file, flat directory, or barcode subdirectories
    def input_path = file(params.fastq, checkIfExists: true)

    if (input_path.isFile()) {
        // Single FASTQ file
        def alias = params.sample ?: input_path.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, "")
        def meta  = [alias: alias, id: alias]
        ch_reads = Channel.of([meta, input_path])
    } else {
        // Directory input
        def children = input_path.listFiles() ?: []
        def barcode_dirs = children.findAll {
            it.isDirectory() && (it.name =~ /^barcode\d+$/)
        }

        if (barcode_dirs) {
            // Barcode sub-directory layout: one sample per barcode
            ch_reads = Channel.fromPath("${params.fastq}/*", type: 'dir')
            | filter { it.name =~ /^barcode\d+$/ }
            | map { dir -> [[alias: dir.name, id: dir.name], dir] }
            | catFastqDir
        } else {
            // Flat directory: all FASTQ files → one sample
            def alias = params.sample ?: input_path.name
            def meta  = [alias: alias, id: alias]
            ch_reads = Channel.of([meta, input_path]) | catFastqDir
        }
    }

    // Filter reads by length and quality
    ch_reads = filterReads(ch_reads)

    // Warn if no reads survive filtering
    ch_reads
    | map { meta, reads, n_reads -> n_reads as int }
    | collect
    | map { counts ->
        if (counts.sum() == 0) {
            log.warn "No reads survived pre-processing filters. " +
                "Consider relaxing --min_read_length, --max_read_length, --min_read_qual."
        }
    }

    // Drop samples with fewer than min_n_reads
    ch_reads = ch_reads
    | filter { meta, reads, n_reads -> (n_reads as int) >= params.min_n_reads }
    | map { meta, reads, n_reads -> [meta, reads] }

    // Subset reads (by length) or downsample randomly
    if (params.drop_frac_longest_reads || params.take_longest_remaining_reads) {
        ch_reads = subsetReads(
            ch_reads,
            params.drop_frac_longest_reads ?: 0,
            params.take_longest_remaining_reads,
            params.reads_downsampling_size
        )
    } else if (params.reads_downsampling_size) {
        ch_reads = downsampleReads(ch_reads, params.reads_downsampling_size)
    }

    // Run de novo pipeline: try miniasm -> racon first
    deNovoPipeline_asm(ch_reads, "miniasm")

    // Re-run with spoa for samples that failed miniasm
    deNovoPipeline_spoa(
        deNovoPipeline_asm.out.metas_failed
        | combine(ch_reads, by: 0),
        "spoa",
    )

    // Merge results from both pipelines
    ch_de_novo_results = deNovoPipeline_asm.out.passed
    | mix(deNovoPipeline_spoa.out.passed)
    | multiMap { meta, cons, bam, bai, bamstats, flagstat, depth ->
        consensus: [meta, cons]
        mapped:    [meta, bam, bai]
        depth:     [meta, depth]
    }

    // Concatenate QC summary TSVs (miniasm + spoa may both have run)
    ch_qc_summaries = concatQCSummaries(
        deNovoPipeline_asm.out.qc_summaries
        | mix(deNovoPipeline_spoa.out.qc_summaries)
        | groupTuple
    )

    // Collect files to publish
    ch_to_publish = Channel.empty()
    | mix(
        ch_de_novo_results.consensus
        | map { meta, cons -> [cons, "${meta.alias}/consensus"] },
        ch_de_novo_results.mapped
        | map { meta, bam, bai -> [[bam, bai], "${meta.alias}/alignments"] }
        | transpose,
        ch_de_novo_results.depth
        | map { meta, depth -> [depth, "${meta.alias}/qc"] },
        ch_qc_summaries
        | map { meta, tsv -> [tsv, "${meta.alias}/qc"] },
        software_versions | map { [it, null] },
        workflow_params | map { [it, null] },
    )

    ch_to_publish | toList | flatMap | publish
}
