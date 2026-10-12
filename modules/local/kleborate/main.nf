// Local module: Kleborate v3.
// The nf-core 'kleborate' module is pinned to v2.1.0, which predates the v3
// rewrite (minimap2, modular presets, E. coli MLST) used in this course.
process KLEBORATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // No container pinned: with -profile wave, a container is built from environment.yml

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.kleborate.txt"), emit: txt
    tuple val(meta), path("${prefix}.kleborate_hamronization.txt"), emit: hamronization, optional: true
    tuple val("${task.process}"), val('kleborate'), eval("kleborate --version 2>&1 | sed 's/^.*Kleborate v//; s/ .*//'"), emit: versions_kleborate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    // kpsc = Klebsiella pneumoniae species complex (species, MLST, AMR, virulence, K/O via Kaptive)
    // escherichia = E. coli MLST (Achtman and Pasteur schemes)
    def preset = task.ext.preset ?: (meta.species == 'ecoli' ? 'escherichia' : 'kpsc')
    """
    kleborate \\
        $args \\
        -a $fasta \\
        -p $preset \\
        -o kleborate_out

    # Kleborate names its table after the species it detected, and for the K. pneumoniae species complex it
    # also writes the same AMR results in hAMRonization format. Take the main table by name.
    main_out=""
    for f in klebsiella_pneumo_complex_output.txt klebsiella_oxytoca_complex_output.txt escherichia_output.txt; do
        if [ -s "kleborate_out/\$f" ]; then main_out="kleborate_out/\$f"; break; fi
    done
    if [ -z "\$main_out" ]; then
        echo "ERROR: Kleborate wrote no results table for ${fasta} (preset ${preset}): the assembly's species" >&2
        echo "       does not match the preset. Check the species column of the samplesheet and the assembly QC." >&2
        ls kleborate_out >&2
        exit 1
    fi
    cp "\$main_out" ${prefix}.kleborate.txt

    if [ -s kleborate_out/klebsiella_pneumo_complex_hAMRonization_output.txt ]; then
        cp kleborate_out/klebsiella_pneumo_complex_hAMRonization_output.txt ${prefix}.kleborate_hamronization.txt
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // One-row table: species from the samplesheet, ST from phylo_group (if it looks like an ST)
    def stub_species = meta.species == 'ecoli' ? 'Escherichia coli' : 'Klebsiella pneumoniae'
    def stub_st      = (meta.phylo_group ?: '') ==~ /ST\d+/ ? meta.phylo_group : 'NA'
    """
    printf 'strain\\tenterobacterales__species__species\\tklebsiella_pneumo_complex__mlst__ST\\n%s\\t%s\\t%s\\n' \\
        '${meta.id}' '${stub_species}' '${stub_st}' > ${prefix}.kleborate.txt
    touch ${prefix}.kleborate_hamronization.txt
    """
}
