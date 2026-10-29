// Local module: SKA2 (split k-mer analysis) for one phylogeny group, used with --phylo_method ska2.
//
//   ska build    : split k-mers of every isolate in the group (trimmed reads or contigs; both may be mixed)
//   ska map      : place the split k-mers on the group's reference -> a reference-length alignment.
//                  The reference sequence is added as 'Reference', so the alignment has the same layout
//                  as snippy-core's <group>.full.aln and goes through the same steps
//                  (snippy-clean_full_aln -> Gubbins -> IQ-TREE -> snp-dists).
//   ska align    : the reference-free alternative, used when no reference is staged (--skip_gubbins with
//                  --phylo_method ska2). Writes the variant split-k-mer sites only: no 'Reference' row and
//                  no genome coordinates, so Gubbins cannot run on it and the tree comes straight from
//                  this alignment. It keeps the invariant sites, so IQ-TREE needs no ASC correction.
//   ska distance : reference-free pairwise SNP distances between the isolates (all shared split k-mers,
//                  not only those on the reference)
//
// Split k-mers find SNPs only when the 31 bp around them match: SNPs closer than ~15 bp to another
// variant, and indels, are missed. That suits closely related isolates (one ST), as here.
process SKA2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // No container pinned: built from environment.yml (bin/prepare_shared_containers.sh), or -profile wave

    input:
    tuple val(meta), val(ids), val(nfiles), path(seqs, stageAs: 'in???/*')   // nfiles: 2 for reads (R1, R2), 1 for contigs
    path reference                                                           // [] for a reference-free run (ska align)

    output:
    tuple val(meta), path("${prefix}.full.aln")         , emit: full_aln
    tuple val(meta), path("${prefix}.ska_distances.tsv"), emit: distances
    tuple val(meta), path("${prefix}.skf")              , emit: skf
    tuple val(meta), path("${prefix}.ska_inputs.tsv")   , emit: inputs
    tuple val("${task.process}"), val('ska2'), eval("ska --version 2>&1 | sed 's/^ska //'"), emit: versions_ska2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''      // ska build
    def args2 = task.ext.args2 ?: ''      // ska align (reference-free) or ska map
    def args3 = task.ext.args3 ?: ''      // ska distance
    prefix    = task.ext.prefix ?: "${meta.id}"
    def files = seqs instanceof List ? seqs : [seqs]
    // One line per isolate: name <tab> file [<tab> file2]
    def lines = []
    def k = 0
    ids.eachWithIndex { id, i ->
        def n = nfiles[i] as int
        lines << ([id] + files[k..<(k + n)].collect { f -> f.toString() }).join('\t')
        k += n
    }
    // No reference staged: align the split k-mers against each other instead of onto a reference
    def refFree = !reference || (reference instanceof List && reference.isEmpty())
    def alignCmd = refFree ?
        """
    # Reference-free: no 'Reference' row and no genome coordinates. ska align's own defaults are
    # --min-freq 0.9 and --filter no-const; ext.args2 (conf/modules.config) changes them.
    ska align ${args2} --threads ${task.cpus} -o ${prefix}.full.aln ${prefix}.skf
        """ :
        """
    ska map ${args2} ${reference} ${prefix}.skf --threads ${task.cpus} -o ${prefix}.ska_map.aln

    # Reference first (contigs joined in file order, as ska map does), then the isolates
    {
        printf '>Reference\\n'
        grep -v '^>' ${reference} | tr -d '\\r\\n ' | tr 'acgtn' 'ACGTN'
        printf '\\n'
        cat ${prefix}.ska_map.aln
    } > ${prefix}.full.aln

    # Every sequence must have the reference's length
    awk '/^>/ { if (seq != "") print length(seq); seq = ""; next } { seq = seq \$0 } END { print length(seq) }' ${prefix}.full.aln \\
        | sort -u | awk 'END { if (NR != 1) { print "ERROR: sequences in ${prefix}.full.aln differ in length" > "/dev/stderr"; exit 1 } }'
        """

    """
    printf '%s\\n' ${lines.collect { line -> "'${line}'" }.join(' ')} > ${prefix}.ska_inputs.tsv

    ska build \\
        ${args} \\
        --threads ${task.cpus} \\
        -f ${prefix}.ska_inputs.tsv \\
        -o ${prefix}
    ${alignCmd}

    ska distance ${args3} --threads ${task.cpus} -o ${prefix}.ska_distances.tsv ${prefix}.skf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.skf ${prefix}.ska_inputs.tsv
    printf '>Reference\\nACGTACGTAC\\n' > ${prefix}.full.aln
    for id in ${ids.join(' ')}; do printf '>%s\\nACGTACGTAC\\n' \$id >> ${prefix}.full.aln; done
    printf 'Sample1\\tSample2\\tDistance\\tMismatches\\n' > ${prefix}.ska_distances.tsv
    """
}
