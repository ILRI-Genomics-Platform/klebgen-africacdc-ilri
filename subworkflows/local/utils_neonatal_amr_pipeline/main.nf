//
// PIPELINE_INITIALISATION: read and validate the samplesheet.
//
// Two kinds of samplesheet (CSV, header required); --input_type auto (default) tells them apart by the header:
//
//   reads:     sample,species,fastq_1,fastq_2,long_fastq,phylo_group
//   assembly:  sample,species,fasta,phylo_group
//
//   - species     : kpneumoniae | ecoli (also accepted: 'K. pneumoniae', 'Klebsiella pneumoniae',
//                   'E. coli', 'Escherichia coli', ...); chooses the Kleborate preset. Kleborate then
//                   reports the species from the genome itself.
//   - fastq_1/2   : paired-end Illumina reads
//   - long_fastq  : Oxford Nanopore reads (optional; triggers hybrid assembly)
//   - fasta       : an assembly (contigs), e.g. from another assembler; read QC, Kraken2 and assembly
//                   are skipped for these samples, and the phylogeny maps the contigs (snippy --ctgs)
//   - phylo_group : optional label (e.g. ST39); with --phylo_grouping samplesheet,
//                   isolates sharing a label are analysed together in the phylogeny
//                   step (the default groups by Kleborate species + ST instead)
//
// Validation is written out in Groovy rather than with the nf-schema plugin
// so that trainees can read every rule and the pipeline runs offline.
//

include { readReferenceSheet } from '../phylogeny/main'

workflow PIPELINE_INITIALISATION {

    take:
    input    // path to samplesheet CSV

    main:
    if (!input) {
        error("Please provide a samplesheet with --input (see assets/samplesheet.csv)")
    }
    def sheet = file(input, checkIfExists: true)
    def inputType = detectInputType(sheet, params.input_type)
    log.info "Input: ${inputType == 'assembly' ? 'assemblies (contigs) - read QC, Kraken2 and assembly are skipped' : 'sequencing reads'} (${sheet.name})"
    if (!(params.phylo_min_group_size instanceof Integer) || params.phylo_min_group_size < 2) {
        error("--phylo_min_group_size must be a whole number of 2 or more (got '${params.phylo_min_group_size}'); " +
            "groups of 3 or more also get a recombination-filtered tree")
    }

    // Phylogeny references: every reference genome in the sheet must already be downloaded
    if (!params.skip_phylogeny) {
        if (!(params.phylo_grouping in ['kleborate', 'samplesheet'])) {
            error("--phylo_grouping must be 'kleborate' or 'samplesheet' (got '${params.phylo_grouping}')")
        }
        if (params.phylo_grouping == 'kleborate') {
            def missingRefs = readReferenceSheet(params.phylo_reference_sheet)
                .collect { row -> "${params.phylo_reference_dir}/${row.accession}.fasta" }
                .unique()
                .findAll { path -> !file(path).exists() }
            if (missingRefs) {
                error("Missing phylogeny reference genome(s):\n  ${missingRefs.join('\n  ')}\n" +
                    "Download them once with:\n" +
                    "  bin/fetch_references.py --dir ${params.phylo_reference_dir}\n" +
                    "(the references come from ${params.phylo_reference_sheet}), or run with " +
                    "--phylo_grouping samplesheet, or --skip_phylogeny.")
            }
        }
    }

    ch_samplesheet = channel
        .fromPath(sheet)
        .splitCsv(header: true, strip: true)
        .map { row -> inputType == 'assembly' ? validateAssemblyRow(row, sheet.parent) : validateRow(row, sheet.parent) }

    // Duplicate sample IDs would silently overwrite each other's results
    ch_samplesheet
        .map { meta, _sr, _lr -> meta.id }
        .collect()
        .map { ids ->
            def dups = ids.countBy { it }.findAll { _id, n -> n > 1 }.keySet()
            if (dups) {
                error("Samplesheet: duplicate sample IDs: ${dups.join(', ')}")
            }
        }

    emit:
    samplesheet = ch_samplesheet   // reads: [ meta, [ fastq_1, fastq_2 ], long_fastq | [] ]; assembly: [ meta, [ fasta ], [] ]
}

//
// Reads or assemblies? --input_type reads|assembly, or auto: from the samplesheet header
//
def detectInputType(Path sheet, String requested) {
    def header = sheet.readLines().find { line -> line.trim() }?.split(',')?.collect { col -> col.trim() } ?: []
    def hasReads = 'fastq_1' in header
    def hasFasta = 'fasta' in header
    def type = requested ?: 'auto'
    if (!(type in ['auto', 'reads', 'assembly'])) {
        error("--input_type must be auto, reads or assembly (got '${type}')")
    }
    if (type == 'auto') {
        if (hasReads && hasFasta) {
            error("Samplesheet ${sheet.name} has both 'fastq_1' and 'fasta' columns: use one kind per run, or set --input_type")
        }
        if (!hasReads && !hasFasta) {
            error("Samplesheet ${sheet.name}: no 'fastq_1' (reads) or 'fasta' (assemblies) column. Header: ${header.join(',')}")
        }
        type = hasFasta ? 'assembly' : 'reads'
    }
    def needed = type == 'assembly' ? ['sample', 'species', 'fasta'] : ['sample', 'species', 'fastq_1', 'fastq_2']
    def missing = needed.findAll { col -> !(col in header) }
    if (missing) {
        error("Samplesheet ${sheet.name} (--input_type ${type}): missing column(s) ${missing.join(', ')}")
    }
    return type
}

//
// Species the pipeline knows how to type: 'kpneumoniae' (K. pneumoniae species complex) or 'ecoli'.
// Common spellings from lab metadata are accepted: 'K. pneumoniae', 'Klebsiella pneumoniae', 'E. coli', ...
//
def normaliseSpecies(String sample, String value) {
    def s = (value ?: '').toLowerCase().replaceAll(/[^a-z]/, '')
    if (s in ['kpneumoniae', 'klebsiellapneumoniae', 'kpn', 'kp'] || s.contains('pneumoniae') || s.contains('variicola')) {
        return 'kpneumoniae'
    }
    if (s in ['ecoli', 'escherichiacoli', 'eco'] ) {
        return 'ecoli'
    }
    error("Samplesheet: sample '${sample}' has species '${value}'; expected K. pneumoniae (kpneumoniae) or E. coli (ecoli)")
}

def checkSampleAndGroup(Map row) {
    if (!(row.sample ==~ /[A-Za-z0-9_.-]+/)) {
        error("Samplesheet: sample '${row.sample}' may only contain letters, digits, '_', '.' and '-'")
    }
    def group = row.phylo_group ?: null
    if (group && !(group ==~ /[A-Za-z0-9_-]+/)) {
        error("Samplesheet: phylo_group '${group}' (sample ${row.sample}) may only contain letters, digits, '_' and '-'")
    }
    return group
}

//
// Validate one assembly row: [ meta, [ fasta ], [] ]
//
def validateAssemblyRow(Map row, Path sheetDir) {
    def missing = ['sample', 'species', 'fasta'].findAll { col -> !row[col] }
    if (missing) {
        error("Samplesheet row ${row}: missing value(s) for ${missing.join(', ')}")
    }
    def group = checkSampleAndGroup(row)
    def fasta = resolvePath(row.fasta, sheetDir)
    if (fasta.name.endsWith('.gz')) {
        error("Samplesheet: sample '${row.sample}': ${fasta.name} is compressed; not every tool reads gzipped contigs, " +
            "so decompress it first (gunzip -k ${fasta.name}) and list the .fasta file")
    }
    if (!(fasta.name ==~ /.*\.(fa|fasta|fna|fas)$/)) {
        error("Samplesheet: sample '${row.sample}': '${fasta.name}' does not look like a FASTA file (.fa, .fasta, .fna)")
    }
    def meta = [
        id            : row.sample,
        species       : normaliseSpecies(row.sample, row.species),
        single_end    : false,
        phylo_group   : group,
        has_long      : false,
        input_type    : 'assembly',
        assembly_type : 'provided',
    ]
    return [meta, [fasta], []]
}

//
// Validate one reads row and build the meta map
//
def validateRow(Map row, Path sheetDir) {
    def required = ['sample', 'species', 'fastq_1', 'fastq_2']
    def missing = required.findAll { col -> !row[col] }
    if (missing) {
        error("Samplesheet row ${row}: missing value(s) for ${missing.join(', ')}")
    }
    def group = checkSampleAndGroup(row)
    def meta = [
        id          : row.sample,
        species     : normaliseSpecies(row.sample, row.species),
        single_end  : false,
        phylo_group : group,
        has_long    : row.long_fastq ? true : false,
        input_type  : 'reads',
    ]
    def short_reads = [resolvePath(row.fastq_1, sheetDir), resolvePath(row.fastq_2, sheetDir)]
    def long_reads  = row.long_fastq ? resolvePath(row.long_fastq, sheetDir) : []
    return [meta, short_reads, long_reads]
}

//
// Relative paths in the samplesheet are read relative to the samplesheet's folder
//
def resolvePath(String path, Path sheetDir) {
    def isRemoteOrAbsolute = path ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || path.startsWith('/')
    return file(isRemoteOrAbsolute ? path : "${sheetDir}/${path}", checkIfExists: true)
}
