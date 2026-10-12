// Local module: Autocycler consensus long-read assembly, end to end.
// nf-core has modules for Autocycler's individual sub-commands (subsample,
// compress, cluster, trim, resolve, combine) but not for the step that builds
// the input assemblies with several assemblers, so the whole recipe is one
// process here. Trainees can split it into the nf-core sub-modules as an exercise.
process AUTOCYCLER_ASSEMBLE {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // No container pinned: with -profile wave, a container is built from environment.yml

    input:
    tuple val(meta), path(longreads)

    output:
    tuple val(meta), path("${prefix}.autocycler.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.autocycler.gfa")  , emit: gfa
    tuple val("${task.process}"), val('autocycler'), eval("autocycler --version 2>&1 | sed 's/^[^ ]* //'"), emit: versions_autocycler, topic: versions
    tuple val("${task.process}"), val('flye')      , eval("flye --version 2>&1")                      , emit: versions_flye      , topic: versions
    tuple val("${task.process}"), val('raven')     , eval("raven --version 2>&1")                     , emit: versions_raven     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix          = task.ext.prefix ?: "${meta.id}"
    def assemblers  = task.ext.assemblers ?: 'flye raven'   // add miniasm etc. via ext.assemblers
    """
    genome_size=\$(autocycler helper genome_size --reads $longreads --threads $task.cpus)

    autocycler subsample \\
        --reads $longreads \\
        --out_dir subsampled_reads \\
        --genome_size \$genome_size

    mkdir -p assemblies
    for assembler in $assemblers; do
        for read_set in subsampled_reads/sample_*.fastq; do
            i=\$(basename \$read_set .fastq | sed 's/sample_//')
            autocycler helper \$assembler \\
                --reads \$read_set \\
                --out_prefix assemblies/\${assembler}_\${i} \\
                --threads $task.cpus \\
                --genome_size \$genome_size
        done
    done

    autocycler compress -i assemblies -a autocycler_out
    autocycler cluster -a autocycler_out
    for cluster in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c \$cluster
        autocycler resolve -c \$cluster
    done
    autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

    cp autocycler_out/consensus_assembly.fasta ${prefix}.autocycler.fasta
    cp autocycler_out/consensus_assembly.gfa   ${prefix}.autocycler.gfa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.autocycler.fasta ${prefix}.autocycler.gfa
    """
}
