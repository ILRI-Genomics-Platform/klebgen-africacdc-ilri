// Local module: E. coli Clermont phylogroup assignment with EzClermont
// (no nf-core module exists for EzClermont or ClermonTyping).
process EZCLERMONT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // No container pinned: with -profile wave, a container is built from environment.yml

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.ezclermont.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.ezclermont.log"), emit: log
    tuple val("${task.process}"), val('ezclermont'), eval("python -c 'import importlib.metadata as m; print(m.version(\"ezclermont\"))'"), emit: versions_ezclermont, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # EzClermont can exit non-zero when no phylogroup is assigned; keep the
    # result as long as it wrote one, so an 'unknown' call does not fail the run.
    ezclermont $args $fasta > ${prefix}.ezclermont.tsv 2> ${prefix}.ezclermont.log || [ -s ${prefix}.ezclermont.tsv ]
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ezclermont.tsv ${prefix}.ezclermont.log
    """
}
