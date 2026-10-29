//
// NEONATAL_AMR: bacterial genomic surveillance of neonatal K. pneumoniae
// species complex and E. coli isolates (Illumina + optional ONT reads).
//
// Course module map:
//   READ_QC     -> Module 1   (QC, contamination)
//   ASSEMBLY    -> Module 2   (short-read, long-read and hybrid assembly)
//   PROKKA      -> optional genome annotation (--run_prokka)
//   BAKTA       -> optional genome annotation (--run_bakta --bakta_db <db>)
//   TYPING_AMR  -> Modules 3, 4, 5, 7 (typing, AMR/virulence, plasmids, K/O)
//   PHYLOGENY   -> Module 6   (recombination-aware trees, SNP distances)
//   AGGREGATE   -> Module 8   (one table per unit of observation; --skip_aggregate)
//   MULTIQC     -> Module 8   (summary report; integration in R is an exercise)
//

include { READ_QC    } from '../subworkflows/local/read_qc/main'
include { ASSEMBLY   } from '../subworkflows/local/assembly/main'
include { TYPING_AMR } from '../subworkflows/local/typing_amr/main'
include { PHYLOGENY  } from '../subworkflows/local/phylogeny/main'
include { PROKKA     } from '../modules/nf-core/prokka/main'
include { BAKTA_BAKTA } from '../modules/nf-core/bakta/bakta/main'
include { AGGREGATE_RESULTS } from '../modules/local/aggregate/main'
include { MULTIQC    } from '../modules/nf-core/multiqc/main'

workflow NEONATAL_AMR {

    take:
    ch_samplesheet   // reads: [ meta, [ fastq_1, fastq_2 ], long_fastq | [] ]; assembly: [ meta, [ fasta ], [] ]

    main:
    // Samples given as assemblies (--input_type assembly) skip read QC, Kraken2 and assembly
    ch_input = ch_samplesheet.branch { meta, _files, _long ->
        assembly: meta.input_type == 'assembly'
        reads:    true
    }
    ch_given_assemblies = ch_input.assembly.map { meta, files, _long -> [meta, files[0]] }

    READ_QC(ch_input.reads, params.kraken2_db)

    ASSEMBLY(READ_QC.out.short_reads, READ_QC.out.long_reads, ch_given_assemblies)

    // Optional genome annotation (genes, proteins, GenBank files) for every assembly
    ch_prokka_mqc = channel.empty()
    if (params.run_prokka) {
        PROKKA(ASSEMBLY.out.assembly, [], [])
        ch_prokka_mqc = PROKKA.out.txt.map { _meta, txt -> txt }
    }
    ch_bakta_mqc = channel.empty()
    if (params.run_bakta) {
        BAKTA_BAKTA(ASSEMBLY.out.assembly, file(params.bakta_db, checkIfExists: true), [], [], [], [])
        ch_bakta_mqc = BAKTA_BAKTA.out.txt.map { _meta, txt -> txt }
    }

    TYPING_AMR(ASSEMBLY.out.assembly, params.amrfinderplus_db, params.skip_kaptive)

    ch_snp_distances = channel.empty()
    ch_ska_distances = channel.empty()
    if (!params.skip_phylogeny) {
        PHYLOGENY(
            // what Snippy (or SKA2) compares: trimmed reads, or the given contigs
            READ_QC.out.short_reads.mix(ch_given_assemblies),
            ASSEMBLY.out.assembly,
            TYPING_AMR.out.kleborate,
            params.phylo_grouping,
            params.phylo_reference_sheet,
            params.phylo_reference_dir,
            params.phylo_references,
            params.phylo_min_group_size,
            params.phylo_method,
            params.skip_gubbins
        )
        ch_snp_distances = PHYLOGENY.out.snp_distances
        ch_ska_distances = PHYLOGENY.out.ska_distances
    }

    //
    // One table per unit of observation (bin/aggregate_results.py). Every source is staged
    // into the task's work directory, so nothing here reads params.outdir.
    //
    if (!params.skip_aggregate) {
        // Two sources need renaming first: snp-dists writes <group>.tsv and MOB-suite writes
        // mobtyper_results.txt for every isolate, so the file name alone cannot identify them.
        ch_snpdists_named = ch_snp_distances
            .collectFile { meta, tsv -> ["${meta.id}.snp_distances.tsv", tsv.text] }
        ch_mobtyper_named = TYPING_AMR.out.plasmids
            .collectFile { meta, txt -> ["${meta.id}.mobtyper.txt", txt.text] }

        def gather = { ch -> ch.map { _meta, f -> f }.collect().ifEmpty([]) }

        AGGREGATE_RESULTS(
            gather(READ_QC.out.fastp_json),
            gather(READ_QC.out.kraken2),
            gather(READ_QC.out.bracken),
            gather(ASSEMBLY.out.quast),
            ASSEMBLY.out.checkm2.collect().ifEmpty([]),
            ch_bakta_mqc.collect().ifEmpty([]),
            ch_prokka_mqc.collect().ifEmpty([]),
            gather(TYPING_AMR.out.kleborate),
            gather(TYPING_AMR.out.kaptive_k.mix(TYPING_AMR.out.kaptive_o)),
            gather(TYPING_AMR.out.serotype),
            gather(TYPING_AMR.out.phylogroup),
            gather(TYPING_AMR.out.amrfinder),
            ch_mobtyper_named.collect().ifEmpty([]),
            ch_snpdists_named.mix(ch_ska_distances.map { _meta, f -> f }).collect().ifEmpty([]),
            file(params.input, checkIfExists: true),
            params.aggregate_metadata ? file(params.aggregate_metadata, checkIfExists: true) : []
        )
    }

    //
    // Software versions: every module publishes [process, tool, version] to the
    // 'versions' topic; collect them into one YAML file (also shown in MultiQC).
    //
    ch_versions_yml = channel.topic('versions')
        .map { process, tool, version -> [process.tokenize(':').last(), tool, version.toString().trim()] }
        .unique()
        .collect(flat: false)
        .map { rows -> versionsToYaml(rows) }
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'neonatal_amr_software_mqc_versions.yml',
            newLine: true
        )

    //
    // MultiQC summary report
    //
    ch_multiqc_files = READ_QC.out.multiqc
        .mix(ASSEMBLY.out.multiqc)
        .mix(ch_prokka_mqc)
        .mix(ch_bakta_mqc)
        .mix(ch_versions_yml)
        .collect()
    ch_multiqc_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []
    MULTIQC(
        ch_multiqc_files.map { files -> [[id: 'neonatal_amr'], files, ch_multiqc_config, [], [], []] }
    )

    emit:
    assembly  = ASSEMBLY.out.assembly
    kleborate = TYPING_AMR.out.kleborate
    amrfinder = TYPING_AMR.out.amrfinder
    multiqc   = MULTIQC.out.report
}

//
// Turn [process, tool, version] rows into MultiQC "software versions" YAML
//
def versionsToYaml(List rows) {
    def byProcess = [:].withDefault { [:] }
    rows.each { row -> byProcess[row[0]][row[1]] = row[2] }
    return byProcess
        .sort()
        .collect { process, tools ->
            "${process}:\n" + tools.sort().collect { tool, version -> "    ${tool}: \"${version}\"" }.join('\n')
        }
        .join('\n')
}
