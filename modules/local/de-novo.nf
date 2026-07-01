include {
    alignReads as alignDraft;
    alignReads as alignPolished;
    bamstats as bamstatsDraft;
    bamstats as bamstatsPolished;
    mosdepthPerBase;
    mosdepthWindows;
} from "./common"


/*
Get draft consensus with SPOA.
The Python script interleaves reads so that forward and reverse reads are
uniformly distributed before being passed to SPOA.
*/
process spoa {
    label "wfamplicon"
    cpus 1
    memory "8 GB"
    input: tuple val(meta), path("reads.fastq.gz")
    output: tuple val(meta), path("reads.fastq.gz"), path("asm.fasta"), optional: true
    script:
    String min_cov_args = ""
    if (params.spoa_minimum_relative_coverage) {
        min_cov_args = "--relative-min-coverage $params.spoa_minimum_relative_coverage"
    }
    """
    echo "${meta.alias}"  # makes some debugging easier

    workflow-glue run_spoa reads.fastq.gz \
        $min_cov_args \
        --max-allowed-read-length $params.spoa_max_allowed_read_length \
        -o asm.fasta
    """
}

/*
Assemble a draft consensus with miniasm.
Uses minimap2 with params.overlap_preset for all-vs-all read overlap.
If none of the assembled contigs is longer than params.force_spoa_length_threshold,
STATUS is set to 'failed' and the sample will fall back to SPOA.
*/
process miniasm {
    label "wfamplicon"
    cpus params.threads
    memory "${[8, 15, 31][task.attempt - 1]} GB"
    errorStrategy { task.exitStatus in [137, 140] ? "retry" : "terminate" }
    maxRetries 2
    input: tuple val(meta), path("reads.fastq.gz")
    output: tuple val(meta), path("reads.fastq.gz"), path("asm.fasta"), env(STATUS)
    script:
    int mapping_threads = Math.max(1, task.cpus - 1)
    """
    STATUS=failed
    echo "${meta.alias}"  # makes some debugging easier

    (
        set +eo pipefail

        minimap2 -L -x ${params.overlap_preset} -t $mapping_threads \
            --cap-kalloc 100m --cap-sw-mem 50m \
            reads.fastq.gz reads.fastq.gz \
        | miniasm -s 100 -e 3 -f reads.fastq.gz - \
        | awk '/^S/{print ">"\$2"\\n"\$3}' > asm.fasta

        exit_codes=("\${PIPESTATUS[@]}")
        echo "pipe exit codes: \${exit_codes[*]}"
        memory_fail=0
        other_fail=0
        for exit_code in "\${exit_codes[@]}"; do
            case \$exit_code in
                0) : ;;
                137 | 140) memory_fail=\$exit_code ;;
                *) other_fail=\$exit_code ;;
            esac
        done
        [[ \$memory_fail == 0 ]] || exit \$memory_fail
        [[ \$other_fail == 0 ]] || exit \$other_fail
    )

    if [[ -s asm.fasta ]]; then
        samtools faidx asm.fasta
        longest=\$(sort -k2rn asm.fasta.fai | head -n1 | cut -f2)
        if [[ \$longest -gt $params.force_spoa_length_threshold ]]; then
            STATUS=passed
        fi
    fi
    """
}

/*
Polish draft consensus with one round of racon.
Uses minimap2 with params.overlap_preset for draft alignment.
`--no-trimming` is added because racon sometimes trims too aggressively;
we trim downstream based on depth instead.
*/
process racon {
    label "wfamplicon"
    cpus params.threads
    memory "8 GB"
    input: tuple val(meta), path("reads.fastq.gz"), path("draft.fasta")
    output: tuple val(meta), path("reads.fastq.gz"), path("polished.fasta")
    script:
    int mapping_threads = Math.max(1, task.cpus - 1)
    """
    echo "${meta.alias}"  # makes some debugging easier

    minimap2 -L -x ${params.overlap_preset} -t $mapping_threads \
        --cap-kalloc 100m --cap-sw-mem 50m \
        draft.fasta reads.fastq.gz \
    | bgzip > pre-racon.paf.gz

    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus -q -1 --no-trimming \
        reads.fastq.gz pre-racon.paf.gz draft.fasta \
        > polished.fasta
    """
}

/*
Use mosdepth to calculate per-base depths from draft alignment.
This is used only for trimming the consensus; see mosdepthWindows for the
windowed version used in QC reports.
*/
process trimAndQC {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta),
            path("consensus.fasta"),
            path("flagstat.tsv"),
            path("depth.per-base.bed.gz")
        val asm_method
    output:
        tuple val(meta), path("passed/consensus.fasta"), path("passed/qc-summary.tsv"),
            optional: true, emit: passed
        tuple val(meta), path("failed/qc-summary.tsv"), optional: true, emit: failed
    script:
    """
    echo "${meta.alias}"  # makes some debugging easier

    workflow-glue trim_and_qc \
        --alias ${meta.alias} \
        --asm-method $asm_method \
        --depth depth.per-base.bed.gz \
        --flagstat flagstat.tsv \
        --consensus consensus.fasta \
        --outdir-pass passed \
        --outdir-fail failed \
        --minimum-depth $params.minimum_mean_depth \
        --primary-threshold $params.primary_alignments_threshold \
        --relative-depth-trim-threshold $params.spoa_minimum_relative_coverage \
        --qc-summary-tsv qc-summary.tsv

    if [[ ( -f passed/qc-summary.tsv ) && ( -f failed/qc-summary.tsv ) ]]; then
        echo "Found 'qc-summary.tsv' in both 'passed' and 'failed' directories." >&2
        exit 1
    elif [[ ( ! -f passed/qc-summary.tsv ) && ( ! -f failed/qc-summary.tsv ) ]]; then
        echo "Found 'qc-summary.tsv' in neither 'passed' nor 'failed' directories." >&2
        exit 1
    fi
    """
}


// De novo consensus workflow
// Runs miniasm -> racon (primary) or spoa (fallback), then QC + trimming.
workflow pipeline {
    take:
        // expected shape: [meta, reads]
        ch_reads
        method
    main:
        if (method !in ["miniasm", "spoa"]) {
            error "Invalid de novo method '$method' (must be 'miniasm' or 'spoa')."
        }

        if (method == "miniasm") {
            ch_branched = miniasm(ch_reads)
            | branch { meta, reads, asm, status ->
                passed: status == "passed"
                failed: status == "failed"
                err: error "Post-assembly status is neither 'passed' nor 'failed'."
            }

            ch_draft = ch_branched.passed
            | map { meta, reads, asm, status -> [meta, reads, asm] }
            | racon
            | map { meta, reads, polished -> [meta, reads, polished] }

            ch_failed = ch_branched.failed
        } else {
            ch_draft = spoa(ch_reads)
            ch_failed = Channel.empty()
        }

        // re-align reads against the draft/polished consensus
        alignDraft(ch_draft)
        bamstatsDraft(alignDraft.out)
        mosdepthPerBase(alignDraft.out)

        // QC and trim the consensus
        trimAndQC(
            ch_draft
            | map { meta, reads, asm -> [meta, asm] }
            | join(
                bamstatsDraft.out | map { meta, bamstats, flagstat -> [meta, flagstat] }
            )
            | join(mosdepthPerBase.out),
            method
        )

        // re-align reads against the final selected and trimmed consensus
        ch_reads
        | join(trimAndQC.out.passed, remainder: true)
        | filter { null !in it }
        | map { meta, reads, cons, qc_summary -> [meta, reads, cons] }
        | alignPolished

        bamstatsPolished(alignPolished.out)
        mosdepthWindows(alignPolished.out, params.number_depth_windows)

    emit:
        passed = trimAndQC.out.passed
        | map { meta, cons, qc_summary -> [meta, cons] }
        | join(alignPolished.out)
        | join(bamstatsPolished.out)
        | join(mosdepthWindows.out)

        metas_failed = ch_failed
        | map { meta, reads, asm, status -> meta }
        | mix(trimAndQC.out.failed | map { meta, qc_summary -> meta })

        qc_summaries = trimAndQC.out.passed
        | map { meta, cons, qc_summary -> [meta, qc_summary] }
        | mix(trimAndQC.out.failed)
}
