// Local module: Snippy on an assembly instead of reads (snippy --ctgs). Snippy cuts
// the contigs into simulated paired reads and maps them to the reference, so the
// outputs match the nf-core snippy/run module and feed snippy-core unchanged.
// Used when the samplesheet gives assemblies (--input_type assembly).
process SNIPPY_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2':
        'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(contigs)
    path reference

    output:
    tuple val(meta), path("${prefix}/${prefix}.vcf")       , emit: vcf
    tuple val(meta), path("${prefix}/${prefix}.aligned.fa"), emit: aligned_fa
    tuple val(meta), path("${prefix}/${prefix}.tab")       , emit: tab
    tuple val(meta), path("${prefix}/${prefix}.txt")       , emit: txt
    tuple val(meta), path("${prefix}/${prefix}.log")       , emit: log
    tuple val("${task.process}"), val('snippy'), eval("snippy --version 2>&1 | sed 's/snippy //'"), emit: versions_snippy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    snippy \\
        ${args} \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        --outdir ${prefix} \\
        --reference ${reference} \\
        --prefix ${prefix} \\
        --ctgs ${contigs}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.vcf ${prefix}/${prefix}.aligned.fa ${prefix}/${prefix}.tab ${prefix}/${prefix}.txt ${prefix}/${prefix}.log
    """
}
