//
// PHYLOGENY (Module 6): one recombination-aware phylogeny and SNP-distance
// matrix per group of related isolates, for transmission analysis.
//
//   Kleborate species (+ ST) --> one group per species (default) or per ST, each with a reference genome
//
//   reads (or contigs) of the group's isolates
//      |-- snippy: snippy per isolate (--ctgs for contigs) --> snippy-core  --> <group>.full.aln
//      `-- ska2:   ska build (whole group) --> ska map onto the reference  --> <group>.full.aln
//                    `--> ska distance: reference-free SNP distances           (<group>_ska2 with 'both')
//   each full alignment --clean--> Gubbins (3+ isolates) --> IQ-TREE
//                                          `--> snp-dists
//
// Method (--phylo_method):
//   snippy (default) map each isolate's reads (or contigs: snippy --ctgs) to the reference; snippy-core
//   ska2             split k-mers of the whole group (SKA2): ska build, then ska map onto the same reference.
//                    Its alignment has the same layout as snippy-core's, so the steps after it are the same.
//                    Also writes reference-free SNP distances (<group>.ska_distances.tsv).
//   both             both of the above for every group, to compare them: Snippy's results are named after the
//                    group (Kpneumoniae_ST39), SKA2's get '_ska2' (Kpneumoniae_ST39_ska2).
//                    bin/compare_snp_methods.py puts the two sets of SNP distances side by side.
//
// Groups and references (--phylo_grouping):
//   species     (default) group = species reported by Kleborate, whatever the ST, e.g. Kpneumoniae, Ecoli:
//               every isolate of a species is mapped to the same reference (the species row of
//               --phylo_reference_sheet; ST rows are not used) and all of them get one SNP-distance matrix and
//               one tree. Isolates without an ST from Kleborate are left out, as in kleborate mode.
//   kleborate   group = species + ST reported by Kleborate, e.g. Kpneumoniae_ST39.
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
    grouping            // 'species', 'kleborate' or 'samplesheet'
    reference_sheet     // CSV species,st,name,accession (species and kleborate modes)
    reference_dir       // folder with <accession>.fasta (species and kleborate modes)
    phylo_references    // CSV group,fasta or null (samplesheet mode)
    min_group_size      // groups smaller than this are skipped
    method              // 'snippy', 'ska2' or 'both'
    skip_gubbins        // true: no recombination filtering; with 'ska2' alone, SKA2 runs
                        //       'ska align' (reference-free) instead of 'ska map'

    main:
    if (grouping in ['species', 'kleborate']) {
        // ---- group by Kleborate species (species mode) or species + ST (kleborate mode);
        //      reference by species (or, in kleborate mode, ST) from the sheet
        def refRows = readReferenceSheet(reference_sheet)
        def bySpecies = grouping == 'species'

        ch_grouped = ch_short_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_kleborate.map { meta, txt -> [meta.id, parseKleborate(txt)] })
            .flatMap { id, meta, reads, typing ->
                def st = typing.st
                def noSt = !st || st in ['NA', '-', '']
                if (!typing.species || noSt) {
                    log.warn "Phylogeny: ${id} has no species/ST from Kleborate (species '${typing.species}', ST '${st}'); left out"
                    return []
                }
                def ref = chooseReference(refRows, typing.species, bySpecies ? null : st)
                if (!ref) {
                    log.warn "Phylogeny: ${id} is '${typing.species}', which has no reference in ${reference_sheet}; left out"
                    return []
                }
                def label = bySpecies ? speciesLabel(typing.species) : "${speciesLabel(typing.species)}_${st}"
                def group = label.replaceAll(/[^A-Za-z0-9_-]/, '')
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
    if (method in ['ska2', 'both']) {
        // ---- SKA2: one run per group, on all its isolates (reads and/or contigs) and its reference.
        // Reference-free only when SKA2 is the sole method: in 'both' mode SKA2 keeps 'ska map' so that its
        // alignment stays on the same coordinates as Snippy's (bin/compare_snp_methods.py).
        def skaRefFree = skip_gubbins && method == 'ska2'
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
                ska: [[id: method == 'both' ? "${group}_ska2" : group, group: group, method: 'ska2',
                       aln_type: skaRefFree ? 'align' : 'map'], ids, nfiles, files.flatten()]
                ref: skaRefFree ? [] : ref
            }
        SKA2(ch_ska_in.ska, ch_ska_in.ref)
        ch_full_aln      = ch_full_aln.mix(SKA2.out.full_aln)
        ch_ska_distances = SKA2.out.distances
    }
    if (method in ['snippy', 'both']) {
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
                core: [[id: group, group: group, method: 'snippy'], vcfs, aligned]
                ref:  ref
            }
        SNIPPY_CORE(ch_core_in.core, ch_core_in.ref)
        ch_full_aln      = ch_full_aln.mix(SNIPPY_CORE.out.full_aln)
        ch_core_aln      = SNIPPY_CORE.out.aln
    }
    // A reference-free 'ska align' alignment holds variant sites only, with no 'Reference' row and no
    // genome coordinates, so snippy-clean_full_aln and Gubbins do not apply to it.
    ch_aln = ch_full_aln.branch { meta, _aln ->
        variant: meta.aln_type == 'align'
        whole:   true
    }
    SNIPPY_CLEAN(ch_aln.whole)

    // ---- groups big enough for a tree: alignment = reference + isolates >= 4 sequences.
    // Gubbins drops every isolate with more than 25% missing data (N, '-', '?': its default --filter-percentage)
    // and stops if fewer than 3 sequences are left. Typical cause: isolates mapped to a reference of another
    // species (about half the genome is then missing), or low-coverage assemblies. Those isolates are named
    // here; if fewer than 3 usable isolates remain, the group gets SNP distances only instead of failing the run.
    ch_clean = SNIPPY_CLEAN.out.aln
        .map { meta, aln -> [meta.group, meta, aln] }
        .combine(ch_group_sizes, by: 0)
        .map { _group, meta, aln, n ->
            def usable = n
            if (!workflow.stubRun) {
                def gappy = alignmentMissing(aln).drop(1).findAll { rec -> rec.missing > gubbinsMaxMissing() }
                if (gappy) {
                    def named = gappy.collect { rec -> "${rec.name} (${String.format('%.1f', rec.missing)}%)" }.join(', ')
                    if (skip_gubbins) {
                        log.warn "Phylogeny group ${meta.id}: ${gappy.size()} of ${n} isolates have more than " +
                            "${gubbinsMaxMissing()}% missing data: ${named}. With --skip_gubbins they are kept in the tree " +
                            "and SNP distances, so values involving them are unreliable. Check that the reference is the " +
                            "right species (conf/phylo_references.csv) and the isolates' assembly quality."
                    } else {
                        usable = n - gappy.size()
                        log.warn "Phylogeny group ${meta.id}: ${gappy.size()} of ${n} isolates have more than " +
                            "${gubbinsMaxMissing()}% missing data, so Gubbins leaves them out of the tree and SNP distances: " +
                            named +
                            ". Check that the reference is the right species (conf/phylo_references.csv) and the isolates' assembly quality."
                    }
                }
            }
            [meta, aln, n, usable]
        }
        .branch { meta, aln, n, usable ->
            tree: usable >= 3
                return [meta, aln]
            pair: true
                if (usable < n) {
                    log.warn "Phylogeny group ${meta.id}: only ${usable} isolate(s) with enough data for Gubbins and IQ-TREE " +
                        "(they need 3 or more). Reporting SNP distances from the core alignment only (all ${n} isolates, " +
                        "without recombination filtering); distances involving the isolates above are unreliable."
                } else {
                    log.warn "Phylogeny group ${meta.id} has ${n} isolates: too few for Gubbins and IQ-TREE (they need 3 or more). " +
                        "Reporting SNP distances from the core alignment only, without recombination filtering."
                }
                return [meta, aln]
        }

    // ---- remove recombination, then build the tree on the filtered SNP sites.
    // nf-core gubbins has no meta map; outputs are named after the input file
    // (<group>.aln), so the group ID is recovered from the file name.
    // Reference-free alignments cannot go through Gubbins, so they go straight to the tree. They keep
    // their invariant sites (ska align does not drop them), so IQ-TREE needs no ASC correction.
    ch_variant = ch_aln.variant
        .map { meta, aln -> [meta.group, meta, aln] }
        .combine(ch_group_sizes, by: 0)
        .branch { _group, meta, aln, n ->
            tree: n >= 3
                return [meta, aln]
            pair: true
                log.warn "Phylogeny group ${meta.id} has ${n} isolate(s): too few for IQ-TREE (3 or more needed). " +
                    "Reporting SNP distances only."
                return [meta, aln]
        }

    // With --skip_gubbins, Gubbins is given an empty channel and never runs; the tree is then built on the
    // cleaned alignment itself, which is NOT recombination-corrected.
    GUBBINS(skip_gubbins ? channel.empty() : ch_clean.tree.map { _meta, aln -> aln })
    ch_snps = (skip_gubbins
            ? ch_clean.tree
            : GUBBINS.out.fasta.map { fasta -> [[id: fasta.name.tokenize('.')[0]], fasta] })
        .mix(ch_variant.tree)

    // A group whose isolates differ by no SNPs outside recombination has an empty filtered alignment, or one
    // whose few columns hold only N, '-' or '?': IQ-TREE stops with "Unknown sequence type". Such a group gets no tree (with a warning), but still its SNP
    // distances (all 0). Stub runs write empty placeholder files, so the check is skipped there.
    ch_tree_input = ch_snps
        .branch { meta, fasta ->
            tree: workflow.stubRun || alignmentVariableSites(fasta) > 0
                return [meta, fasta]
            no_sites: true
                log.warn "Phylogeny group ${meta.id}: no variable A/C/G/T sites left after Gubbins (the isolates do not " +
                    "differ outside recombination and missing data), so no tree is built. Check the group's core " +
                    "alignment, the isolates' coverage, and Gubbins' recombination predictions."
                return [meta, fasta]
        }
    IQTREE(ch_tree_input.tree.map { meta, fasta -> [meta, fasta, []] }, [], [], [], [], [], [], [], [], [], [], [], [])
    SNPDISTS(ch_snps.mix(ch_clean.pair).mix(ch_variant.pair))

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
// First ST-specific row for this species (when an ST is given), otherwise the first species-wide row
//
def chooseReference(List rows, String species, String st) {
    def bySpecies = rows.findAll { row -> species.toLowerCase().startsWith(row.species.toLowerCase()) }
    return (st ? bySpecies.find { row -> row.st && row.st == st } : null) ?: bySpecies.find { row -> !row.st }
}

//
// 'Klebsiella pneumoniae subsp. pneumoniae' -> 'Kpneumoniae'
//
def speciesLabel(String species) {
    def words = species.tokenize(' ')
    return words.size() >= 2 ? "${words[0][0].toUpperCase()}${words[1]}" : words[0]
}

//
// Gubbins' default --filter-percentage: sequences with more missing data (%) are left out
//
def gubbinsMaxMissing() { return 25.0 }

//
// Percentage of missing data (N, '-', '?', X) per sequence of an alignment, in file order (reference first):
// [[name: , missing: ], ...]. Read line by line, so long single-line sequences are fine.
//
def alignmentMissing(fasta) {
    def recs = []
    def name = null
    long total = 0
    long gaps = 0
    fasta.eachLine { line ->
        if (line.startsWith('>')) {
            if (name != null) { recs << [name: name, missing: total ? 100.0 * gaps / total : 100.0] }
            name = line.substring(1).trim().tokenize(' ')[0]
            total = 0
            gaps = 0
        } else if (name != null) {
            def s = line.trim()
            total += s.length()
            gaps += s.length() - s.replaceAll(/[NnXx?-]/, '').length()
        }
    }
    if (name != null) { recs << [name: name, missing: total ? 100.0 * gaps / total : 100.0] }
    return recs
}

//
// Number of alignment columns where at least two sequences carry different bases (A, C, G or T).
// Columns of only N, '-' or '?' do not count: IQ-TREE cannot detect the sequence type from them.
//
def alignmentVariableSites(fasta) {
    def seqs = []
    def current = null
    fasta.eachLine { line ->
        if (line.startsWith('>')) {
            if (current != null) { seqs << current.toString().toUpperCase() }
            current = new StringBuilder()
        } else if (current != null) {
            current.append(line.trim())
        }
    }
    if (current != null) { seqs << current.toString().toUpperCase() }
    if (seqs.size() < 2) { return 0 }
    def length = seqs*.length().min()
    def variable = 0
    for (int i = 0; i < length; i++) {
        def bases = [] as Set
        for (seq in seqs) {
            def c = seq.charAt(i)
            if (c == ('A' as char) || c == ('C' as char) || c == ('G' as char) || c == ('T' as char)) { bases << c }
        }
        if (bases.size() >= 2) { variable++ }
    }
    return variable
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
