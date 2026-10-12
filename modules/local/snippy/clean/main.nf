// Local module: snippy-clean_full_aln, which replaces non-ACGT characters
// with N so Gubbins can read snippy-core's full alignment. The output is named
// after the phylogeny group (meta.id) because the nf-core gubbins module has
// no meta map and names its outputs after the input file.
process SNIPPY_CLEAN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2':
        'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(full_aln)

    output:
    tuple val(meta), path("${prefix}.aln"), emit: aln
    tuple val("${task.process}"), val('snippy'), eval("snippy --version 2>&1 | sed 's/snippy //'"), emit: versions_snippy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    snippy-clean_full_aln $full_aln > ${prefix}.aln
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln
    """
}
