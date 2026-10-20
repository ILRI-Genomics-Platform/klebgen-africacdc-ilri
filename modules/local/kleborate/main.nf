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
    run_kleborate() {   # \$1 = preset, \$2 = output folder
        kleborate $args -a $fasta -p "\$1" -o "\$2"
    }
    main_table() {      # the species table Kleborate wrote in folder \$1, if any
        for f in klebsiella_pneumo_complex_output.txt klebsiella_oxytoca_complex_output.txt escherichia_output.txt; do
            if [ -s "\$1/\$f" ]; then echo "\$1/\$f"; return; fi
        done
    }

    run_kleborate ${preset} kleborate_out
    main_out=\$(main_table kleborate_out)
    note=""

    # Kleborate writes no table when the genome is not the species of the preset (e.g. the lab said
    # K. pneumoniae, the genome is E. coli). Try the other presets, so the isolate is still typed.
    if [ -z "\$main_out" ]; then
        for other in kpsc kosc escherichia; do
            [ "\$other" = "${preset}" ] && continue
            run_kleborate "\$other" "kleborate_\$other" || true
            main_out=\$(main_table "kleborate_\$other")
            if [ -n "\$main_out" ]; then
                note="species does not match samplesheet (${meta.species}); typed with Kleborate preset \$other"
                break
            fi
        done
    fi

    if [ -n "\$main_out" ]; then
        if [ -n "\$note" ]; then
            awk -F'\\t' -v OFS='\\t' -v n="\$note" 'NR==1 {print \$0, "pipeline_note"; next} {print \$0, n}' \
                "\$main_out" > ${prefix}.kleborate.txt
            echo "WARNING: ${meta.id}: \$note" >&2
        else
            cp "\$main_out" ${prefix}.kleborate.txt
        fi
    else
        # Not K. pneumoniae complex, K. oxytoca complex or E. coli: record that and let the run continue.
        # Without an ST the isolate is left out of the phylogeny (with a warning); AMR and plasmid results
        # are still produced.
        note="not KpSC, KoSC or E. coli according to Kleborate: check the species (Kraken2, QUAST, lab ID) and assembly QC"
        printf 'strain\\tenterobacterales__species__species\\tpipeline_note\\n%s\\t%s\\t%s\\n' \
            '${meta.id}' 'not KpSC, KoSC or E. coli' "\$note" > ${prefix}.kleborate.txt
        echo "WARNING: ${meta.id}: \$note" >&2
    fi

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
