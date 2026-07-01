/**
 * Concatenate all FASTQ files found in a directory into a single compressed FASTQ.
 * Handles .fastq, .fq, .fastq.gz, and .fq.gz files.
 */
process catFastqDir {
    label "wfamplicon"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("input_dir")
    output:
        tuple val(meta), path("${meta.alias}.fastq.gz")
    script:
    """
    (
        find -L input_dir -maxdepth 1 \\( -name '*.fastq.gz' -o -name '*.fq.gz' \\) \\
            -print0 | sort -z | xargs -0 -r zcat
        find -L input_dir -maxdepth 1 \\( -name '*.fastq' -o -name '*.fq' \\) \\
            -print0 | sort -z | xargs -0 -r cat
    ) | bgzip > ${meta.alias}.fastq.gz
    """
}
