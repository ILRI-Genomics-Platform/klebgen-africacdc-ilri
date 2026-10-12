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

    emit:
    assembly = ch_best                                          // [ meta, fasta ]
    graph    = AUTOCYCLER_ASSEMBLE.out.gfa                      // [ meta, gfa ] for Bandage
    multiqc  = QUAST.out.results.map { _meta, dir -> dir }
}
