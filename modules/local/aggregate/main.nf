// Local module: one table per unit of observation, from the results of this run
// (bin/aggregate_results.py; see README.md, "Aggregating the results").
//
// Every source is STAGED into in/<source>/ and passed with --dir, so this step reads its
// inputs from the work directory like any other process. It never reads params.outdir, and
// so cannot race with publishDir, and -resume works normally.
//
// Two sources are renamed by the caller before they get here, because their work-directory
// names do not identify the isolate or group: snp-dists writes <group>.tsv and MOB-suite
// writes mobtyper_results.txt for every isolate.
//
// No container and no image in conf/singularity_images.config: the script uses only the
// Python standard library, so python3 on the node is enough.
process AGGREGATE_RESULTS {
    label 'process_single'

    conda "conda-forge::python=3.12"

    input:
    path fastp        , stageAs: 'in/fastp/*'
    path kraken2      , stageAs: 'in/kraken2/*'
    path bracken      , stageAs: 'in/bracken/*'
    path quast        , stageAs: 'in/quast/*'
    path checkm2      , stageAs: 'in/checkm2/*'
    path bakta        , stageAs: 'in/bakta/*'
    path prokka       , stageAs: 'in/prokka/*'
    path kleborate    , stageAs: 'in/kleborate/*'
    path kaptive      , stageAs: 'in/kaptive/*'
    path ectyper      , stageAs: 'in/ectyper/*'
    path ezclermont   , stageAs: 'in/ezclermont/*'
    path amrfinderplus, stageAs: 'in/amrfinderplus/*'
    path mobsuite     , stageAs: 'in/mobsuite/*'
    path phylogeny    , stageAs: 'in/phylogeny/*'
    path samplesheet
    path metadata

    output:
    path "out/master_table.tsv"     , emit: master
    path "out/*_all.tsv"            , emit: tool_tables, optional: true
    path "out/*_long.tsv"           , emit: long_tables, optional: true
    path "out/aggregate_summary.txt", emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | sed 's/^Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // One --dir per source that actually has files; an empty input stages no directory
    def sources = [fastp: fastp, kraken2: kraken2, bracken: bracken, quast: quast,
                   checkm2: checkm2, bakta: bakta, prokka: prokka, kleborate: kleborate,
                   kaptive: kaptive, ectyper: ectyper, ezclermont: ezclermont,
                   amrfinderplus: amrfinderplus, mobsuite: mobsuite, phylogeny: phylogeny]
    def dir_args = sources
        .findAll { _name, files -> files }
        .collect { name, _files -> "--dir ${name}=in/${name}" }
        .join(' \\\n        ')
    """
    aggregate_results.py \\
        ${dir_args} \\
        ${samplesheet ? "--samplesheet ${samplesheet}" : ''} \\
        ${metadata ? "--metadata ${metadata}" : ''} \\
        --outdir out \\
        ${args}
    """

    stub:
    """
    mkdir -p out
    printf 'sample\\n' > out/master_table.tsv
    printf 'stub run: no tables aggregated\\n' > out/aggregate_summary.txt
    """
}
