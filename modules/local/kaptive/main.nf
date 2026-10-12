// Local module: standalone Kaptive 3 on K. pneumoniae species complex assemblies.
// Kleborate already runs Kaptive and reports the K and O loci in its table; this
// step adds Kaptive's full report (every gene found in and outside the locus,
// with identity and coverage) and the locus sequences, for capsule (K) and
// LPS (O) typing in detail. It uses the Kleborate environment, so the Kaptive
// version and databases are the same as in Kleborate's results.
process KAPTIVE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // No container pinned: with -profile wave, a container is built from environment.yml;
    // with the course's shared images, Kleborate's image is used (conf/singularity_images.config)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.kaptive_K.tsv"), emit: k_locus
    tuple val(meta), path("${prefix}.kaptive_O.tsv"), emit: o_locus
    tuple val(meta), path("${prefix}_{K,O}_locus/*"), emit: loci, optional: true
    tuple val("${task.process}"), val('kaptive'), eval("kaptive --version 2>&1 | sed 's/^[^0-9]*//'"), emit: versions_kaptive, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Kaptive 3.3 renamed the locus-FASTA option from --fasta to --loci
    if kaptive assembly --help 2>&1 | grep -q -- '--loci'; then fasta_opt=--loci; else fasta_opt=--fasta; fi

    # kpsc_k: capsule (K) loci; kpsc_o: LPS (O) loci, both for the K. pneumoniae species complex
    for locus in K O; do
        db=kpsc_\$(echo \$locus | tr 'KO' 'ko')
        mkdir -p ${prefix}_\${locus}_locus
        kaptive assembly \$db $fasta \\
            $args \\
            -o ${prefix}.kaptive_\${locus}.tsv \\
            \$fasta_opt ${prefix}_\${locus}_locus
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.kaptive_K.tsv ${prefix}.kaptive_O.tsv
    """
}
