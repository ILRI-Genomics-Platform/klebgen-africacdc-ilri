//
// PHYLOGENY (Module 6): one recombination-aware phylogeny and SNP-distance
// matrix per group of related isolates, for transmission analysis.
//
//   reads --snippy--> core alignment --clean--> Gubbins --> IQ-TREE
//                                                     \--> snp-dists
//
// Method (--phylo_method):
//   snippy (default) map each isolate's reads (or contigs: snippy --ctgs) to the reference; snippy-core
//   ska2             split k-mers of the whole group (SKA2): ska build, then ska map onto the same reference.
//                    Its alignment has the same layout as snippy-core's, so the steps after it are the same.
//                    Also writes reference-free SNP distances (<group>.ska_distances.tsv).
//
// Groups and references (--phylo_grouping):
//   kleborate   (default) group = species + ST reported by Kleborate, e.g. Kpneumoniae_ST39.
//               Reference from the sheet --phylo_reference_sheet (conf/phylo_references.csv):
//               an ST-specific row if there is one, otherwise the species row.
//   samplesheet group = the samplesheet's phylo_group column. Reference from
//               --phylo_references (CSV group,fasta) if given, otherwise the
//               group's best own assembly (hybrid first, then fewest contigs).
//
// Gubbins (RAxML) and IQ-TREE bootstrap need at least 4 sequences: the reference
// plus 3 isolates. A group of 2 isolates (--phylo_min_group_size 2) gets SNP
// distances from the core alignment only, without recombination filtering or a tree.
//

include { SNIPPY_RUN   } from '../../../modules/nf-core/snippy/run/main'
include { SNIPPY_CORE  } from '../../../modules/nf-core/snippy/core/main'
include { GUBBINS      } from '../../../modules/nf-core/gubbins/main'
include { IQTREE       } from '../../../modules/nf-core/iqtree/main'
include { SNPDISTS     } from '../../../modules/nf-core/snpdists/main'
include { SNIPPY_CLEAN } from '../../../modules/local/snippy/clean/main'
include { SNIPPY_CONTIGS } from '../../../modules/local/snippy/contigs/main'
include { SKA2          } from '../../../modules/local/ska2/main'

workflow PHYLOGENY {

    take:
    ch_short_reads      // [ meta, [ R1, R2 ] ], or [ meta, contigs ] for assemblies (meta.input_type == 'assembly')
    ch_assembly         // [ meta, fasta ]
    ch_kleborate        // [ meta, kleborate.txt ]
    grouping            // 'kleborate' or 'samplesheet'
    reference_sheet     // CSV species,st,name,accession (kleborate mode)
    reference_dir       // folder with <accession>.fasta (kleborate mode)
    phylo_references    // CSV group,fasta or null (samplesheet mode)
    min_group_size      // groups smaller than this are skipped
    method              // 'snippy' or 'ska2'

    main:
    if (grouping == 'kleborate') {
        // ---- group by Kleborate species + ST; reference by species (or ST) from the sheet
        def refRows = readReferenceSheet(reference_sheet)

        ch_grouped = ch_short_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_kleborate.map { meta, txt -> [meta.id, parseKleborate(txt)] })
            .flatMap { id, meta, reads, typing ->
                def st = typing.st
                if (!typing.species || !st || st in ['NA', '-', '']) {
                    log.warn "Phylogeny: ${id} has no species/ST from Kleborate (species '${typing.species}', ST '${st}'); left out"
                    return []
                }
                def ref = chooseReference(refRows, typing.species, st)
                if (!ref) {
                    log.warn "Phylogeny: ${id} is '${typing.species}', which has no reference in ${reference_sheet}; left out"
                    return []
                }
                def group = "${speciesLabel(typing.species)}_${st}".replaceAll(/[^A-Za-z0-9_-]/, '')
                def fasta = file("${reference_dir}/${ref.accession}.fasta", checkIfExists: true)
                return [[meta + [phylo_group: group, kleborate_species: typing.species, kleborate_st: st], reads, ref, fasta]]
            }

        ch_reads_grouped = ch_grouped.map { meta, reads, _ref, _fasta -> [meta, reads] }
        ch_refs = ch_grouped
            .map { meta, _reads, ref, fasta -> [meta.phylo_group, ref.name, ref.accession, fasta] }
            .unique { item -> item[0] }
            .map { group, name, acc, fasta ->
                log.info "Phylogeny group ${group}: reference = ${name} (${acc})"
                [group, fasta]
            }
    } else {
        // ---- samplesheet phylo_group column
        ch_reads_grouped = ch_short_reads.filter { meta, _reads -> meta.phylo_group }

        if (phylo_references) {
            ch_refs = channel
                .fromPath(phylo_references, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row -> [row.group, file(row.fasta, checkIfExists: true)] }
        } else {
            ch_refs = ch_assembly
                .filter { meta, _fasta -> meta.phylo_group }
                .map { meta, fasta ->
                    // rank: hybrid assemblies first, then fewest contigs
                    def rank = [meta.assembly_type == 'hybrid' ? 0 : 1, fasta.countFasta(), meta.id]
                    [meta.phylo_group, rank, fasta]
                }
                .groupTuple()
                .map { group, ranks, fastas ->
                    // order by assembly type, then contig count, then sample ID
                    def best = (0..<ranks.size()).toList().sort { a, b ->
                        (ranks[a][0] <=> ranks[b][0]) ?: (ranks[a][1] <=> ranks[b][1]) ?: (ranks[a][2] <=> ranks[b][2])
                    }.first()
                    log.info "Phylogeny group ${group}: reference = ${ranks[best][2]} (${fastas[best].name})"
                    [group, fastas[best]]
                }
        }
    }

    // ---- keep groups with enough isolates
    ch_group_sizes = ch_reads_grouped
        .map { meta, _reads -> [meta.phylo_group, meta.id] }
        .groupTuple()
        .map { group, ids -> [group, ids.size()] }
    ch_kept_groups = ch_group_sizes
        .filter { group, n ->
            if (n < min_group_size) {
                log.warn "Phylogeny group ${group} has ${n} isolate(s) (< ${min_group_size}); skipped"
            }
            n >= min_group_size
        }
        .map { group, _n -> group }

    ch_full_aln      = channel.empty()
    ch_core_aln      = channel.empty()
    ch_ska_distances = channel.empty()
    if (method == 'ska2') {
        // ---- SKA2: one run per group, on all its isolates (reads and/or contigs) and its reference
        ch_ska_in = ch_reads_grouped
            .map { meta, files -> [meta.phylo_group, meta, files] }
            .combine(ch_kept_groups, by: 0)
            .map { group, meta, files ->
                def list = files instanceof List ? files : [files]
                [group, meta.id, list.size(), list]
            }
            .groupTuple()
            .combine(ch_refs, by: 0)
            .multiMap { group, ids, nfiles, files, ref ->
                ska: [[id: group], ids, nfiles, files.flatten()]
                ref: ref
            }
        SKA2(ch_ska_in.ska, ch_ska_in.ref)
        ch_full_aln      = SKA2.out.full_aln
        ch_ska_distances = SKA2.out.distances
    } else {
        // ---- map each isolate to its group reference: reads with Snippy, assemblies with snippy --ctgs
        ch_to_map = ch_reads_grouped
            .map { meta, reads -> [meta.phylo_group, meta, reads] }
            .combine(ch_kept_groups, by: 0)
            .combine(ch_refs, by: 0)
            .branch { _group, meta, _reads, _ref ->
                contigs: meta.input_type == 'assembly'
                reads:   true
            }
        ch_snippy_reads = ch_to_map.reads.multiMap { _group, meta, reads, ref ->
            reads: [meta, reads]
            ref:   ref
        }
        SNIPPY_RUN(ch_snippy_reads.reads, ch_snippy_reads.ref)

        ch_snippy_contigs = ch_to_map.contigs.multiMap { _group, meta, contigs, ref ->
            contigs: [meta, contigs]
            ref:     ref
        }
        SNIPPY_CONTIGS(ch_snippy_contigs.contigs, ch_snippy_contigs.ref)

        // ---- one core alignment per group (reads and contigs can share a group)
        ch_core_in = SNIPPY_RUN.out.vcf.mix(SNIPPY_CONTIGS.out.vcf)
            .join(SNIPPY_RUN.out.aligned_fa.mix(SNIPPY_CONTIGS.out.aligned_fa))
            .map { meta, vcf, aligned -> [meta.phylo_group, vcf, aligned] }
            .groupTuple()
            .combine(ch_refs, by: 0)
            .multiMap { group, vcfs, aligned, ref ->
                core: [[id: group], vcfs, aligned]
                ref:  ref
            }
        SNIPPY_CORE(ch_core_in.core, ch_core_in.ref)
        ch_full_aln      = SNIPPY_CORE.out.full_aln
        ch_core_aln      = SNIPPY_CORE.out.aln
    }
    SNIPPY_CLEAN(ch_full_aln)

    // ---- groups big enough for a tree: alignment = reference + isolates >= 4 sequences
    ch_clean = SNIPPY_CLEAN.out.aln
        .map { meta, aln -> [meta.id, meta, aln] }
        .combine(ch_group_sizes, by: 0)
        .branch { group, meta, aln, n ->
            tree: n >= 3
                return [meta, aln]
            pair: true
                log.warn "Phylogeny group ${group} has ${n} isolates: too few for Gubbins and IQ-TREE (they need 3 or more). " +
                    "Reporting SNP distances from the core alignment only, without recombination filtering."
                return [meta, aln]
        }

    // ---- remove recombination, then build the tree on the filtered SNP sites.
    // nf-core gubbins has no meta map; outputs are named after the input file
    // (<group>.aln), so the group ID is recovered from the file name.
    GUBBINS(ch_clean.tree.map { _meta, aln -> aln })
    ch_snps = GUBBINS.out.fasta.map { fasta -> [[id: fasta.name.tokenize('.')[0]], fasta] }

    IQTREE(ch_snps.map { meta, fasta -> [meta, fasta, []] }, [], [], [], [], [], [], [], [], [], [], [], [])
    SNPDISTS(ch_snps.mix(ch_clean.pair))

    emit:
    core_aln      = ch_core_aln                    // [ meta(group), aln ] (snippy only)
    gubbins_tree  = GUBBINS.out.tree
    tree          = IQTREE.out.phylogeny           // [ meta(group), treefile ]
    snp_distances = SNPDISTS.out.tsv               // [ meta(group), tsv ]
    ska_distances = ch_ska_distances               // [ meta(group), tsv ] (ska2 only)
}

//
// Reference sheet: CSV species,st,name,accession; lines starting with '#' are comments
//
def readReferenceSheet(sheet) {
    def lines = file(sheet, checkIfExists: true).readLines()
        .collect { line -> line.trim() }
        .findAll { line -> line && !line.startsWith('#') }
    if (!lines) {
        error("Reference sheet ${sheet} has no rows")
    }
    def header = lines[0].split(',', -1).collect { col -> col.trim() }
    def required = ['species', 'st', 'name', 'accession']
    def missing = required.findAll { col -> !(col in header) }
    if (missing) {
        error("Reference sheet ${sheet}: missing column(s) ${missing.join(', ')}; expected ${required.join(',')}")
    }
    return lines.drop(1).collect { line ->
        def values = line.split(',', -1).collect { value -> value.trim() }
        def row = [header, values].transpose().collectEntries { kv -> [kv[0], kv[1]] }
        if (!row.species || !row.accession) {
            error("Reference sheet ${sheet}: every row needs species and accession: '${line}'")
        }
        row
    }
}

//
// First ST-specific row for this species, otherwise the first species-wide row
//
def chooseReference(List rows, String species, String st) {
    def bySpecies = rows.findAll { row -> species.toLowerCase().startsWith(row.species.toLowerCase()) }
    return bySpecies.find { row -> row.st && row.st == st } ?: bySpecies.find { row -> !row.st }
}

//
// 'Klebsiella pneumoniae subsp. pneumoniae' -> 'Kpneumoniae'
//
def speciesLabel(String species) {
    def words = species.tokenize(' ')
    return words.size() >= 2 ? "${words[0][0].toUpperCase()}${words[1]}" : words[0]
}

//
// Species and ST from a Kleborate v3 table (full or --trim_headers column names)
//
def parseKleborate(txt) {
    def lines = txt.readLines().findAll { line -> line.trim() }
    if (lines.size() < 2) {
        return [species: null, st: null]
    }
    def header = lines[0].split('\t', -1) as List
    def values = lines[1].split('\t', -1) as List
    def row = [header, values].transpose().collectEntries { kv -> [kv[0], kv[1]] }
    def speciesCol = header.find { col -> col == 'species' || col.endsWith('__species__species') }
    def stCol = header.find { col -> col == 'ST' || col ==~ /.*__mlst(_achtman)?__ST/ }
    return [
        species: speciesCol ? row[speciesCol]?.trim() : null,
        st:      stCol ? row[stCol]?.trim() : null,
    ]
}
