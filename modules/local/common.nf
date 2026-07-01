/*
Align reads to a reference (or draft consensus) with minimap2.
The mapping preset is controlled by params.minimap2_preset (default: map-pb).
*/
process alignReads {
    label "wfamplicon"
    cpus params.threads
    memory "4 GB"
    input: tuple val(meta), path("reads.fastq.gz"), path("reference.fasta")
    output: tuple val(meta), path("*.bam"), path("*.bai")
    script:
    """
    minimap2 -t $task.cpus -ax ${params.minimap2_preset} \
        --cap-kalloc 100m --cap-sw-mem 50m \
        reference.fasta reads.fastq.gz \
        -R '@RG\\tID:${meta.alias}\\tSM:${meta.alias}' \
    | samtools sort -@ $task.cpus -o "${meta.alias}.aligned.sorted.bam" -

    samtools index "${meta.alias}.aligned.sorted.bam"
    """
}

/*
Calculate alignment statistics (bamstats) for a BAM file.
Outputs a per-read TSV and a flagstat TSV.
*/
process bamstats {
    label "wfamplicon"
    cpus Math.min(params.threads, 2)
    memory "2 GB"
    input: tuple val(meta), path("input.bam"), path("input.bam.bai")
    output: tuple val(meta), path("bamstats.tsv"), path("bamstats-flagstat.tsv")
    script:
    """
    bamstats -u input.bam -s ${meta.alias} -f bamstats-flagstat.tsv -t $task.cpus \
    > bamstats.tsv
    """
}

/*
Run mosdepth to calculate per-base depths and windowed summary depths.
Used for trimming the consensus and for QC reporting.
*/
process mosdepthPerBase {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    memory "4 GB"
    input: tuple val(meta), path("input.bam"), path("input.bam.bai")
    output: tuple val(meta), path("depth.per-base.bed.gz")
    script:
    int mosdepth_extra_threads = task.cpus - 1
    """
    mosdepth -t $mosdepth_extra_threads depth input.bam
    """
}

/*
Run mosdepth over a fixed number of windows for a per-sample depth report.
Emits a TSV with header.
*/
process mosdepthWindows {
    label "wfamplicon"
    cpus Math.min(params.threads, 3)
    memory "4 GB"
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai")
        val n_windows
    output: tuple val(meta), path("per-window-depth.tsv.gz")
    script:
    int mosdepth_extra_threads = task.cpus - 1
    """
    samtools idxstats input.bam | grep -v '^\\*' > idxstats

    if [[ \$(wc -l < idxstats) -ne 1 ]]; then
        echo "Unexpected number of references in input BAM." >&2
        exit 1
    fi

    REF_LENGTH=\$(cut -f2 idxstats)

    window_length=1
    if [ "\$REF_LENGTH" -gt "${n_windows}" ]; then
        window_length=\$(expr \$REF_LENGTH / ${n_windows})
    fi

    mosdepth -t $mosdepth_extra_threads -b \$window_length -n depth input.bam

    cat <(printf "ref\\tstart\\tend\\tdepth\\n" | gzip) depth.regions.bed.gz \
    > per-window-depth.tsv.gz

    rm depth.regions.bed.gz
    """
}
