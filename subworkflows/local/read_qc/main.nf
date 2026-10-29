//
// READ_QC (Module 1): trim Illumina reads, filter ONT reads, report quality,
// and check species composition / contamination with Kraken2 + Bracken.
//

include { FASTP           } from '../../../modules/nf-core/fastp/main'
include { FILTLONG        } from '../../../modules/nf-core/filtlong/main'
include { NANOPLOT        } from '../../../modules/nf-core/nanoplot/main'
include { KRAKEN2_KRAKEN2 } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN } from '../../../modules/nf-core/bracken/bracken/main'

workflow READ_QC {

    take:
    ch_samples     // [ meta, [ fastq_1, fastq_2 ], long_fastq | [] ]
    kraken2_db     // path to Kraken2 database directory, or null to skip

    main:
    ch_multiqc = channel.empty()

    // Illumina: adapter and quality trimming
    FASTP(
        ch_samples.map { meta, sr, _lr -> [meta, sr, []] },
        false,   // discard_trimmed_pass
        false,   // save_trimmed_fail
        false    // save_merged
    )
    ch_multiqc = ch_multiqc.mix(FASTP.out.json.map { _meta, json -> json })

    // ONT: length/quality filtering, then a QC report on the filtered reads
    FILTLONG(
        ch_samples
            .filter { meta, _sr, _lr -> meta.has_long }
            .map { meta, _sr, lr -> [meta, [], lr] }
    )
    NANOPLOT(FILTLONG.out.reads)
    ch_multiqc = ch_multiqc.mix(NANOPLOT.out.txt.map { _meta, txt -> txt })

    // Species composition (catches mixed cultures and mislabelled isolates)
    ch_kraken2_report = channel.empty()
    ch_bracken_report = channel.empty()
    if (kraken2_db) {
        def db = file(kraken2_db, checkIfExists: true)
        KRAKEN2_KRAKEN2(FASTP.out.reads, db, false, false)
        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, db)
        ch_kraken2_report = KRAKEN2_KRAKEN2.out.report
        ch_bracken_report = BRACKEN_BRACKEN.out.txt
        ch_multiqc = ch_multiqc.mix(KRAKEN2_KRAKEN2.out.report.map { _meta, report -> report })
    } else {
        // warn only when there are reads to check (not for an assembly samplesheet)
        ch_samples.first().subscribe { _sample -> log.warn "No --kraken2_db given: skipping the Kraken2/Bracken contamination check" }
    }

    emit:
    short_reads = FASTP.out.reads       // [ meta, [ R1, R2 ] ]
    long_reads  = FILTLONG.out.reads    // [ meta, long_fastq ]
    multiqc     = ch_multiqc
    fastp_json  = FASTP.out.json        // [ meta, json ]      for bin/aggregate_results.py
    kraken2     = ch_kraken2_report     // [ meta, report ]    empty without --kraken2_db
    bracken     = ch_bracken_report     // [ meta, report ]    empty without --kraken2_db
}
