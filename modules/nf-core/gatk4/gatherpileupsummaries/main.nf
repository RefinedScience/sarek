process GATK4_GATHERPILEUPSUMMARIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(pileup)
    path dict

    output:
    tuple val(meta), path("*.pileups.table"), emit: table
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = pileup.collect { pileup_ -> "--I ${pileup_}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GatherPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GatherPileupSummaries \\
        ${input_list} \\
        --O ${prefix}.pileups.table.unsorted \\
        --sequence-dictionary ${dict} \\
        --tmp-dir . \\
        ${args}

    # GatherPileupSummaries concatenates its --I inputs in the order given, with no re-sort.
    # When the inputs come from scattered intervals that were merged out of genomic order
    # (e.g. grouped in task-completion order rather than by position), the merged table can
    # contain backward jumps in position. CalculateContamination's ContaminationSegmenter
    # assumes genome-sorted input and builds segments off consecutive rows, so a backward
    # jump crashes it with "Invalid interval". Re-sort the merged table's data rows here by
    # (contig rank, position), using the .dict's authoritative @SQ order (not lexical
    # chr1,chr10,chr2 sorting), so CalculateContamination always gets a sorted table
    # regardless of how the merge assembled it. The two header lines (#<METADATA>... and the
    # column header) are preserved as-is.
    awk 'BEGIN { FS=OFS="\\t" } /^@SQ/ { n = split(\$0, a, "\\t"); for (i = 1; i <= n; i++) { if (a[i] ~ /^SN:/) { print substr(a[i], 4), ++r } } }' ${dict} > contig_order.tsv

    {
        head -n 2 ${prefix}.pileups.table.unsorted
        tail -n +3 ${prefix}.pileups.table.unsorted \\
            | awk -v OFS='\\t' 'NR == FNR { rank[\$1] = \$2; next } { r = (\$1 in rank) ? rank[\$1] : 999999; print r, \$2, \$0 }' contig_order.tsv - \\
            | sort -k1,1n -k2,2n -s \\
            | cut -f3-
    } > ${prefix}.pileups.table

    rm -f ${prefix}.pileups.table.unsorted contig_order.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileups.table
    """
}
