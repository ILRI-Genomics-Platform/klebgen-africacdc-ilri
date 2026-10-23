//
// ASSEMBLY (Module 2): short-read assembly for every isolate; for isolates
// with ONT reads, a long-read-first hybrid assembly:
//   Autocycler consensus -> Medaka (long-read polish)
//   -> Polypolish -> Pypolca (short-read polish)
// Each isolate's best assembly (hybrid if available) goes downstream.
// Assemblies given in the samplesheet (--input_type assembly) join here: they
// get QUAST statistics like the others and go downstream unchanged.
//

include { SHOVILL              } from '../../../modules/nf-core/shovill/main'
include { MEDAKA               } from '../../../modules/nf-core/medaka/main'
include { PYPOLCA_RUN          } from '../../../modules/nf-core/pypolca/run/main'
include { QUAST                } from '../../../modules/nf-core/quast/main'
include { CHECKM2_PREDICT      } from '../../../modules/nf-core/checkm2/predict/main'
include { AUTOCYCLER_ASSEMBLE  } from '../../../modules/local/autocycler/assemble/main'
include { POLYPOLISH_SHORTREAD } from '../../../modules/local/polypolish/shortread/main'

workflow ASSEMBLY {

    take:
    ch_short_reads   // [ meta, [ R1, R2 ] ]
    ch_long_reads    // [ meta, long_fastq ]
    ch_given         // [ meta, fasta ] assemblies from the samplesheet (meta.assembly_type = 'provided')

    main:
    // Short-read assembly (all isolates)
    SHOVILL(ch_short_reads)

    // Hybrid assembly (isolates with ONT reads)
    AUTOCYCLER_ASSEMBLE(ch_long_reads)
    MEDAKA(ch_long_reads.join(AUTOCYCLER_ASSEMBLE.out.fasta))
    POLYPOLISH_SHORTREAD(MEDAKA.out.assembly.join(ch_short_reads))
    PYPOLCA_RUN(POLYPOLISH_SHORTREAD.out.fasta.join(ch_short_reads))

    // Best assembly per isolate: hybrid if it exists, otherwise short-read
    ch_best = SHOVILL.out.contigs
        .join(PYPOLCA_RUN.out.polished, remainder: true)
        .map { meta, short_asm, hybrid_asm ->
            [meta + [assembly_type: hybrid_asm ? 'hybrid' : 'short'], hybrid_asm ?: short_asm]
        }
        .mix(ch_given)

    // Assembly statistics per isolate (no reference genome)
    QUAST(ch_best, [[:], []], [[:], []])

    // Completeness and contamination (CheckM2), when a database is given. Given contigs
    // (--input_type assembly) are included, so a country run gets a contamination check even though
    // Kraken2 (read-based) is skipped.
    //
    // One run per isolate, as for QUAST: the nf-core shovill module emits every assembly as 'contigs.fa',
    // so one run over all of them would collide on the staged file name. CheckM2 also names each row
    // after its input file, which would be 'contigs' for every short-read assembly, so the sample ID is
    // put back into the first column when the per-isolate reports are merged into one
    // 02_assembly/checkm2/quality_report.tsv (the file name the course pages use).
    ch_checkm2 = channel.empty()
    if (params.checkm2_db) {
        ch_checkm2_db = channel.value([[:], file(params.checkm2_db, checkIfExists: true)])
        CHECKM2_PREDICT(ch_best, ch_checkm2_db)
        ch_checkm2 = CHECKM2_PREDICT.out.checkm2_tsv
            .map { meta, tsv -> checkm2Rows(meta, tsv) }
            // Drop isolates with no report (a stub run, or a failed CheckM2): an empty entry sorts first
            // and would leave collectFile taking its missing header, so the merged table would lose it
            .filter { text -> text as Boolean }
            .collectFile(name: 'quality_report.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/02_assembly/checkm2")
        // Name the assemblies that look contaminated; the run carries on, the trainer decides
        ch_checkm2.subscribe { tsv -> warnOnContamination(tsv, params.checkm2_max_contamination) }
    }

    emit:
    assembly = ch_best                                          // [ meta, fasta ]
    graph    = AUTOCYCLER_ASSEMBLE.out.gfa                      // [ meta, gfa ] for Bandage
    checkm2  = ch_checkm2                                       // quality_report.tsv (empty without --checkm2_db)
    multiqc  = QUAST.out.results.map { _meta, dir -> dir }
}

//
// One isolate's CheckM2 report as text, with the sample ID in the first column: CheckM2 names each row
// after the input file ('contigs' for every short-read assembly), which is not the sample. Returns the
// header plus that isolate's rows, for collectFile to merge; '' when there is nothing (e.g. a stub run).
//
def checkm2Rows(meta, tsv) {
    if (!tsv.exists() || tsv.size() == 0) {
        return ''
    }
    def lines = tsv.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return ''
    }
    def rows = lines.drop(1).collect { line ->
        def values = line.split('\t', -1) as List
        values[0] = meta.id
        values.join('\t')
    }
    return ([lines[0]] + rows).join('\n') + '\n'
}

//
// Warn about assemblies whose CheckM2 contamination is above a threshold (%), and about low completeness.
// Stub runs write an empty placeholder, so there is nothing to read there.
//
def warnOnContamination(tsv, maxContamination) {
    if (workflow.stubRun || !tsv.exists() || tsv.size() == 0) {
        return
    }
    def lines = tsv.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return
    }
    def header = lines[0].split('\t', -1) as List
    def nameCol = header.indexOf('Name')
    def contCol = header.indexOf('Contamination')
    def compCol = header.indexOf('Completeness')
    if (nameCol < 0 || contCol < 0) {
        log.warn "CheckM2: unexpected columns in ${tsv.name} (${header.join(', ')}); no contamination check made"
        return
    }
    def flagged = []
    lines.drop(1).each { line ->
        def values = line.split('\t', -1)
        def contamination = values.size() > contCol ? (values[contCol] as Double) : null
        def completeness = (compCol >= 0 && values.size() > compCol) ? (values[compCol] as Double) : null
        if (contamination != null && contamination > maxContamination) {
            flagged << "${values[nameCol]} (${contamination}% contamination" +
                (completeness != null ? ", ${completeness}% complete)" : ")")
        }
    }
    if (flagged) {
        log.warn "CheckM2: ${flagged.size()} assembly/assemblies above ${maxContamination}% contamination: " +
            "${flagged.join('; ')}. Typing and AMR results for these are unreliable, and a contaminated " +
            "assembly adds missing data to the phylogeny: check them before reporting, and consider leaving " +
            "them out of the samplesheet."
    }
}
