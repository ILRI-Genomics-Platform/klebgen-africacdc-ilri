// Local module: short-read polishing with Polypolish, including alignment.
// Polypolish needs each read file aligned separately with 'bwa mem -a'
// (all alignments kept) and then 'polypolish filter' on the pair; the
// nf-core polypolish/polish module covers only the final step.
process POLYPOLISH_SHORTREAD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // No container pinned: with -profile wave, a container is built from environment.yml

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${prefix}.polypolish.fasta"), emit: fasta
    tuple val("${task.process}"), val('polypolish'), eval("polypolish polish --version |& sed '1!d;s/Polypolish-polish //'"), emit: versions_polypolish, topic: versions
    tuple val("${task.process}"), val('bwa')       , eval('bwa 2>&1 | sed -n "s/^Version: //p"')                , emit: versions_bwa       , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    if [[ $assembly == *.gz ]]; then gzip -cd $assembly > draft.fasta; else cp $assembly draft.fasta; fi

    bwa index draft.fasta
    bwa mem -t $task.cpus -a draft.fasta ${reads[0]} > alignments_1.sam
    bwa mem -t $task.cpus -a draft.fasta ${reads[1]} > alignments_2.sam

    polypolish filter \\
        --in1 alignments_1.sam --in2 alignments_2.sam \\
        --out1 filtered_1.sam  --out2 filtered_2.sam

    polypolish polish $args draft.fasta filtered_1.sam filtered_2.sam > ${prefix}.polypolish.fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.polypolish.fasta
    """
}
