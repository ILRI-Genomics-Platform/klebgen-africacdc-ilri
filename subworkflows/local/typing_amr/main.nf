//
// TYPING_AMR (Modules 3, 4, 5, 7): species, MLST, AMR, virulence, K/O
// antigens (Kleborate + Kaptive), AMRFinderPlus, plasmid reconstruction
// (MOB-suite), and E. coli serotype / phylogroup.
//

include { AMRFINDERPLUS_UPDATE } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN    } from '../../../modules/nf-core/amrfinderplus/run/main'
include { MOBSUITE_RECON       } from '../../../modules/nf-core/mobsuite/recon/main'
include { ECTYPER              } from '../../../modules/nf-core/ectyper/main'
include { KLEBORATE            } from '../../../modules/local/kleborate/main'
include { EZCLERMONT           } from '../../../modules/local/ezclermont/main'
include { KAPTIVE              } from '../../../modules/local/kaptive/main'
include { parseKleborate       } from '../phylogeny/main'

workflow TYPING_AMR {

    take:
    ch_assembly        // [ meta, fasta ]
    amrfinderplus_db   // path to an AMRFinderPlus database (.tar.gz or directory), or null to download
    skip_kaptive       // true: no standalone Kaptive report (Kleborate's K/O columns remain)

    main:
    // Species, MLST, AMR, virulence and K/O loci (KpSC); MLST (E. coli)
    KLEBORATE(ch_assembly)

    // Full Kaptive report and locus sequences, after Kleborate, for isolates that
    // Kleborate places in the K. pneumoniae species complex (whatever the samplesheet says)
    ch_kaptive_k = channel.empty()
    ch_kaptive_o = channel.empty()
    if (!skip_kaptive) {
        ch_kpsc = ch_assembly
            .join(KLEBORATE.out.txt)
            .filter { _meta, _fasta, txt -> isKpsc(parseKleborate(txt).species) }
            .map { meta, fasta, _txt -> [meta, fasta] }
        KAPTIVE(ch_kpsc)
        ch_kaptive_k = KAPTIVE.out.k_locus
        ch_kaptive_o = KAPTIVE.out.o_locus
    }

    // AMR genes and point mutations (dated database recorded in versions)
    if (amrfinderplus_db) {
        ch_amr_db = channel.value(file(amrfinderplus_db, checkIfExists: true))
    } else {
        AMRFINDERPLUS_UPDATE()
        ch_amr_db = AMRFINDERPLUS_UPDATE.out.db
    }
    AMRFINDERPLUS_RUN(
        ch_assembly.map { meta, fasta -> [meta + [organism: amrfinderOrganism(meta.species)], fasta] },
        ch_amr_db
    )

    // Plasmid reconstruction and typing (most informative on hybrid assemblies)
    MOBSUITE_RECON(ch_assembly)

    // E. coli only: O:H serotype and Clermont phylogroup
    ch_ecoli = ch_assembly.filter { meta, _fasta -> meta.species == 'ecoli' }
    ECTYPER(ch_ecoli)
    EZCLERMONT(ch_ecoli)

    emit:
    kleborate     = KLEBORATE.out.txt
    kaptive_k     = ch_kaptive_k
    kaptive_o     = ch_kaptive_o
    amrfinder     = AMRFINDERPLUS_RUN.out.report
    plasmids      = MOBSUITE_RECON.out.mobtyper_results
    contig_report = MOBSUITE_RECON.out.contig_report
    serotype      = ECTYPER.out.tsv
    phylogroup    = EZCLERMONT.out.tsv
}

//
// AMRFinderPlus organism names; these enable species-specific point mutations
//
def amrfinderOrganism(String species) {
    return [kpneumoniae: 'Klebsiella_pneumoniae', ecoli: 'Escherichia'][species]
}

//
// K. pneumoniae species complex members, as named by Kleborate (Kaptive's kpsc databases cover these)
//
def isKpsc(species) {
    return species && species ==~ /^Klebsiella (pneumoniae|quasipneumoniae|variicola|quasivariicola|africana)\b.*/
}
