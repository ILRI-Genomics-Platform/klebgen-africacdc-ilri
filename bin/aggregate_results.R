#!/usr/bin/env Rscript
#
# aggregate_results.R - collect neonatal-amr-nf results into one table per unit of observation.
#
# The R counterpart of bin/aggregate_results.py: same sources, same options, same output file names,
# so the two can be run on one folder and their master tables diffed against each other. Use whichever
# fits your workflow; course/day08.md (section 8.1) builds the same table by hand, one join at a time.
#
#   bin/run_R.sh aggregate_results.R --results results --long --tool-tables
#   bin/run_R.sh aggregate_results.R --results results --metadata training_data/metadata_compact.tsv
#   bin/run_R.sh aggregate_results.R --dir kleborate=/path/a --dir amrfinderplus=/path/b --outdir combined
#
# Options (defaults in brackets):
#   --results DIR      pipeline output folder (the run's --outdir). Optional if every source
#                      you want is named with --dir                                       [none]
#   --dir SOURCE=PATH  read one source from this folder instead of under --results; repeatable
#   --outdir DIR       where the tables are written      [<results>/09_aggregate; required without --results]
#   --metadata FILE    TSV/CSV of extra per-isolate columns (ward, year, AST ...)              [none]
#   --id-col NAME      sample ID column in --metadata                                        [sample]
#   --samplesheet CSV  defines the expected isolates; extras and missing ones are reported      [none]
#   --sections LIST    subset of qc,assembly,typing,amr,plasmids,annotation,phylogeny,all        [all]
#   --snp-threshold N  for the n_within_<N>_snps column                                          [20]
#   --amr-flag NAME=REGEX  add or replace an AMR flag column; NAME= with no regex removes one
#   --kleborate-extra REGEX  extra Kleborate columns to keep, matched on the column name
#   --prefix NAME      base name of the wide table                                     [master_table]
#   --na STR           text written for missing values                                           [NA]
#   --tool-tables      also write one lossless table per tool (kleborate_all.tsv, ...)
#   --beside-source    with --tool-tables, write each next to its source folder
#   --long             also write amr_genes_long.tsv, plasmids_long.tsv, snp_distances_long.tsv
#   --summaries        also write summary_overall.tsv, summary_by_<col>.tsv and the per-gene,
#                      per-K-locus, per-O-type and per-plasmid-cluster counts
#   --group-by LIST    with --summaries, comma-separated columns to summarise by; a column that
#                      is not in the table is skipped with a warning       [ST,kleborate_species]
#   --strict           exit non-zero if a samplesheet isolate produced nothing
#   --list-sources     print the source names and where they are looked for, then exit
#   --describe         print what can and cannot be combined, then exit

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "r_setup.R"))
use_packages(cran = "tidyverse")
suppressPackageStartupMessages(library(tidyverse))

DESCRIBE <- "
What can and cannot be combined
===============================

The question is always: what is one ROW? Two tables can be joined only when they share a unit of
observation (or when one is a lookup for a key the other holds).

1. ONE ROW PER ISOLATE - joined directly into master_table.tsv
   fastp, Kraken2/Bracken, QUAST, CheckM2, Kleborate, Kaptive, ECTyper/EzClermont.
   One file (or row) per isolate, so the join cannot change the row count.

2. MANY ROWS PER ISOLATE - can only enter the wide table as a SUMMARY you choose
   AMRFinderPlus (gene hits), MOB-suite (plasmids), Bakta/Prokka (features).
   Summarising loses information BY DEFINITION: 'n_plasmids = 3' and 'amr_genes = a;b;c' cannot say
   WHICH plasmid carries WHICH gene. Use --long for that.

3. CANNOT BE COMBINED - a different unit of observation entirely
   SNP distances   the unit is a PAIR of isolates (n(n-1)/2 of them per group)
   Trees           the unit is a GROUP, and the content is Newick text
   Gubbins         the unit is a BRANCH or a recombinant BLOCK
   Plasmid clusters the unit is a PLASMID, which occurs in many isolates (many-to-many)
   Contig reports  the unit is a CONTIG

4. LEGITIMATE DERIVED COLUMNS - pair or group data reduced to per-isolate facts, named as such
   phylo_group, phylo_group_n, nearest_sample, nearest_snps, n_within_<threshold>_snps.
   'The SNP distance' is not an isolate-level fact; 'the distance to this isolate's nearest
   neighbour' is.

5. PER-TOOL TABLES (--tool-tables) - one tool, every row, its own columns
   Concatenation, not a join: nothing is summarised and the tool's columns are kept verbatim with a
   'sample' column added. kleborate_all.tsv is one row per isolate, amrfinderplus_all.tsv one row per
   gene hit, mobtyper_all.tsv one row per plasmid - same tool family, different units, hence separate
   files.

A rule of thumb: if adding a column would change the number of rows, it does not belong in the wide
table. Put it in a long or per-tool table and join when you need it.
"

# Where each source lives under --results; override any with --dir NAME=PATH
SOURCE_DIRS <- c(fastp = "01_qc/fastp", kraken2 = "01_qc/kraken2", bracken = "01_qc/bracken",
                 quast = "02_assembly/quast", checkm2 = "02_assembly/checkm2",
                 bakta = "02_assembly/bakta", prokka = "02_assembly/prokka",
                 kleborate = "03_typing/kleborate", kaptive = "03_typing/kaptive",
                 ectyper = "03_typing/ectyper", ezclermont = "03_typing/ezclermont",
                 amrfinderplus = "04_amr/amrfinderplus", mobsuite = "05_plasmids/mobsuite",
                 phylogeny = "06_phylogeny")
SECTIONS <- list(qc = c("fastp", "kraken2", "bracken"), assembly = c("quast", "checkm2"),
                 typing = c("kleborate", "kaptive", "ectyper", "ezclermont"),
                 amr = "amrfinderplus", plasmids = "mobsuite",
                 annotation = c("bakta", "prokka"), phylogeny = "phylogeny")
# tool-table key -> c(file name, source folder it is written beside with --beside-source)
TOOL_TABLES <- list(fastp = c("fastp_all.tsv", "fastp"),
                    kraken2 = c("kraken2_species_all.tsv", "kraken2"),
                    bracken = c("bracken_all.tsv", "bracken"),
                    quast = c("quast_all.tsv", "quast"),
                    checkm2 = c("checkm2_all.tsv", "checkm2"),
                    bakta = c("bakta_all.tsv", "bakta"),
                    prokka = c("prokka_all.tsv", "prokka"),
                    kleborate = c("kleborate_all.tsv", "kleborate"),
                    kaptive_K = c("kaptive_K_all.tsv", "kaptive"),
                    kaptive_O = c("kaptive_O_all.tsv", "kaptive"),
                    ectyper = c("ectyper_all.tsv", "ectyper"),
                    ezclermont = c("ezclermont_all.tsv", "ezclermont"),
                    amrfinderplus = c("amrfinderplus_all.tsv", "amrfinderplus"),
                    mobsuite = c("mobtyper_all.tsv", "mobsuite"))
DEFAULT_AMR_FLAGS <- c(esbl_or_ampc = "CEPHALOSPORIN", carbapenemase = "CARBAPENEM",
                       gentamicin_r = "GENTAMICIN", colistin_r = "COLISTIN")
LEADING <- c("sample", "phylo_group", "phylo_group_n", "kleborate_species", "ST", "kleborate_K_locus",
             "kleborate_O_type", "kaptive_K_locus", "kaptive_K_type", "kaptive_K_confidence",
             "kaptive_O_locus", "kaptive_O_type", "kaptive_O_confidence", "resistance_score",
             "virulence_score", "amr_n_genes", "n_plasmids", "completeness", "contamination",
             "contigs", "assembly_length", "n50")

# ------------------------------------------------------------------------------------------------ options
args  <- commandArgs(trailingOnly = TRUE)
opt   <- list(results = NA, outdir = NA, metadata = NA, id_col = "sample", samplesheet = NA,
              sections = "all", snp_threshold = "20", prefix = "master_table", na = "NA",
              group_by = "ST,kleborate_species")
flags <- list(tool_tables = FALSE, beside_source = FALSE, long = FALSE, summaries = FALSE,
              strict = FALSE, list_sources = FALSE, describe = FALSE)
multi <- list(dir = character(), amr_flag = character(), kleborate_extra = character())

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (any(args %in% c("-h", "--help"))) {
    lines <- readLines(self)
    cat(sub("^# ?", "", lines[2:(which(lines == "")[1] - 1)]), sep = "\n")
    quit(status = 0)
}
for (f in names(flags)) {                       # boolean flags take no value: remove them first
    hit <- args == paste0("--", gsub("_", "-", f))
    if (any(hit)) { flags[[f]] <- TRUE; args <- args[!hit] }
}
if (flags$describe)     { cat(trimws(DESCRIBE), "\n"); quit(status = 0) }
if (flags$list_sources) {
    section_of <- unlist(lapply(names(SECTIONS), function(s) setNames(rep(s, length(SECTIONS[[s]])), SECTIONS[[s]])))
    cat(sprintf("%-16s %-24s %-12s %s\n", "source", "default path under --results", "--sections", "tool table"))
    for (nm in names(SOURCE_DIRS)) {
        tt <- vapply(TOOL_TABLES, function(x) if (x[2] == nm) x[1] else NA_character_, character(1))
        tt <- paste(tt[!is.na(tt)], collapse = ", ")
        cat(sprintf("%-16s %-24s %-12s %s\n", nm, SOURCE_DIRS[[nm]], section_of[[nm]], if (nzchar(tt)) tt else "-"))
    }
    cat("\nOverride any of them with --dir SOURCE=PATH (repeatable).\n")
    quit(status = 0)
}
if (length(args) %% 2 != 0) stop("Options come in pairs: --name value (see --help)")
for (k in seq_len(length(args) %/% 2)) {
    key <- gsub("-", "_", sub("^--", "", args[2 * k - 1]))
    val <- args[2 * k]
    if (key %in% names(multi))      multi[[key]] <- c(multi[[key]], val)
    else if (key %in% names(opt))   opt[[key]]   <- val
    else stop("Unknown option: ", args[2 * k - 1], " (see --help)")
}
if (is.na(opt$results) && !length(multi$dir)) {
    stop("give --results, or --dir SOURCE=PATH for at least one source (--list-sources shows the names)")
}
if (!is.na(opt$results) && !dir.exists(opt$results)) stop("--results is not a directory: ", opt$results)
if (is.na(opt$results) && is.na(opt$outdir)) stop("--outdir is required when --results is not given")
if (flags$beside_source && !flags$tool_tables) stop("--beside-source only applies with --tool-tables")
if (is.na(opt$outdir)) opt$outdir <- file.path(opt$results, "09_aggregate")
snp_threshold <- suppressWarnings(as.integer(opt$snp_threshold))
if (is.na(snp_threshold)) stop("--snp-threshold must be a whole number")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

problems <- character()
note <- function(...) { msg <- paste0(...); problems <<- c(problems, msg); warning(msg, call. = FALSE) }
say  <- function(...) cat(..., "\n", sep = "")

# --dir overrides
overrides <- list()
for (item in multi$dir) {
    if (!str_detect(item, "=")) { note("ignoring --dir '", item, "' (expected NAME=PATH)"); next }
    nm <- str_trim(str_split_fixed(item, "=", 2)[1]); pth <- str_trim(str_split_fixed(item, "=", 2)[2])
    if (!nm %in% names(SOURCE_DIRS)) { note("ignoring --dir '", nm, "': unknown source"); next }
    if (!dir.exists(pth)) { note("--dir ", nm, "=", pth, " is not a directory; skipping that source"); next }
    overrides[[nm]] <- pth
}
src_dir <- function(nm) {
    if (!is.null(overrides[[nm]])) return(overrides[[nm]])
    if (!is.na(opt$results)) {
        p <- file.path(opt$results, SOURCE_DIRS[[nm]])
        if (dir.exists(p)) return(p)
    }
    NA_character_
}

sections <- str_trim(str_split(opt$sections, ",")[[1]])
if ("all" %in% sections) sections <- names(SECTIONS)
unknown <- setdiff(sections, names(SECTIONS))
if (length(unknown)) stop("unknown --sections: ", paste(unknown, collapse = ","))
for (s in names(SECTIONS)) if (any(SECTIONS[[s]] %in% names(overrides))) sections <- union(sections, s)
want <- function(nm) any(vapply(sections, function(s) nm %in% SECTIONS[[s]], logical(1)))

amr_spec <- DEFAULT_AMR_FLAGS
for (item in multi$amr_flag) {
    if (!str_detect(item, "=")) { note("ignoring --amr-flag '", item, "' (expected NAME=REGEX)"); next }
    nm <- str_trim(str_split_fixed(item, "=", 2)[1]); pat <- str_trim(str_split_fixed(item, "=", 2)[2])
    if (!nzchar(pat)) amr_spec <- amr_spec[names(amr_spec) != nm] else amr_spec[nm] <- pat
}

# ------------------------------------------------------------------------------------------------ helpers
clean_id <- function(x) {
    x <- basename(as.character(x))
    x <- sub("\\.(fasta|fa|fna|gfa)(\\.gz)?$", "", x)
    sub("(_contigs|\\.contigs|\\.hybrid|\\.scaffolds|\\.shovill)$", "", x)
}
rd <- function(p, ...) suppressWarnings(read_tsv(p, show_col_types = FALSE, progress = FALSE, ...))
need <- function(d, cols, label) {
    if (is.null(d) || !nrow(d)) return(NULL)
    miss <- setdiff(cols, names(d))
    if (length(miss)) { note(label, ": no column(s) ", paste(miss, collapse = ", "), " - skipped"); return(NULL) }
    d
}
col_or_na <- function(d, pattern) {
    hit <- grep(pattern, names(d), value = TRUE)
    if (length(hit)) as.character(d[[hit[1]]]) else NA_character_
}
read_each <- function(dir, glob, suffix) {        # one file per isolate; ID from the file name
    if (is.na(dir)) return(NULL)
    files <- Sys.glob(file.path(dir, glob))
    if (!length(files)) return(NULL)
    map_dfr(files, function(f) {
        d <- rd(f)
        if (!nrow(d)) return(NULL)
        if (!"sample" %in% names(d)) d$sample <- sub(paste0(suffix, "$"), "", basename(f))
        d
    })
}

wide <- list()      # per-isolate tibbles, joined at the end
raw  <- list()      # tool tables
rows <- list()      # report lines
add_wide <- function(d) if (!is.null(d) && nrow(d)) wide[[length(wide) + 1]] <<- distinct(d, sample, .keep_all = TRUE)
add_raw  <- function(key, d) if (!is.null(d) && nrow(d)) raw[[key]] <<- bind_rows(raw[[key]], d)
report   <- function(nm, nfiles, samples, note_txt = "", where = NA) {
    rows[[length(rows) + 1]] <<- tibble(source = nm, files = nfiles, samples = samples,
                                        note = note_txt, directory = if (is.na(where)) "-" else where)
}

# ------------------------------------------------------------------------------------------------ readers
if (want("fastp")) {
    d <- src_dir("fastp"); files <- if (is.na(d)) character() else Sys.glob(file.path(d, "*.json"))
    if (length(files) && requireNamespace("jsonlite", quietly = TRUE)) {
        tab <- map_dfr(files, function(f) {
            j <- try(jsonlite::fromJSON(f), silent = TRUE)
            if (inherits(j, "try-error")) return(NULL)
            tibble(sample = sub("\\.fastp\\.json$|\\.json$", "", basename(f)),
                   reads_before = j$summary$before_filtering$total_reads %||% NA,
                   reads_after  = j$summary$after_filtering$total_reads  %||% NA,
                   q30_rate_after = j$summary$after_filtering$q30_rate   %||% NA,
                   duplication_rate = j$duplication$rate                 %||% NA)
        })
        add_wide(tab); add_raw("fastp", tab)
        report("fastp", length(files), n_distinct(tab$sample), where = d)
    } else if (length(files)) {
        note("fastp: jsonlite is not installed - skipped")
    } else report("fastp", 0, 0, "directory not found")
}
if (want("kraken2")) {
    d <- src_dir("kraken2"); files <- if (is.na(d)) character() else Sys.glob(file.path(d, "*report.txt"))
    tab <- map_dfr(files, function(f) {
        k <- rd(f, col_names = c("percent", "clade", "taxon", "rank", "taxid", "name"))
        if (!nrow(k)) return(NULL)
        k$name <- str_trim(k$name)
        s <- sub("\\.kraken2\\.report\\.txt$|\\.report\\.txt$|\\.txt$", "", basename(f))
        sp <- filter(k, rank == "S")
        add_raw("kraken2", mutate(sp, sample = s) |> select(sample, rank, taxid, name, percent, clade, taxon))
        top <- if (nrow(sp)) slice_max(sp, taxon, n = 1, with_ties = FALSE) else NULL
        tibble(sample = s,
               kraken2_top_species = if (is.null(top)) NA_character_ else top$name,
               kraken2_top_percent = if (is.null(top)) NA_real_ else top$percent,
               kraken2_unclassified_percent = {u <- filter(k, rank == "U"); if (nrow(u)) u$percent[1] else NA_real_})
    })
    add_wide(tab); report("kraken2", length(files), if (is.null(tab)) 0 else n_distinct(tab$sample),
                          "species rows kept", d)
}
if (want("bracken")) {
    d <- src_dir("bracken"); files <- if (is.na(d)) character() else Sys.glob(file.path(d, "*.tsv"))
    tab <- map_dfr(files, function(f) {
        b <- rd(f)
        if (!nrow(b) || !"fraction_total_reads" %in% names(b)) return(NULL)
        s <- sub("\\.tsv$", "", basename(f))
        add_raw("bracken", mutate(b, sample = s))
        top <- slice_max(b, fraction_total_reads, n = 1, with_ties = FALSE)
        tibble(sample = s, bracken_top_species = top$name, bracken_top_fraction = top$fraction_total_reads)
    })
    add_wide(tab); report("bracken", length(files), if (is.null(tab)) 0 else n_distinct(tab$sample), where = d)
}
if (want("quast")) {
    d <- src_dir("quast"); q <- NULL
    if (!is.na(opt$results) && file.exists(file.path(opt$results, "02_assembly", "quast_all.tsv"))) {
        q <- rd(file.path(opt$results, "02_assembly", "quast_all.tsv"))
    } else if (!is.na(d)) {                       # per-isolate report.tsv: metric in col 1, value in col 2
        files <- Sys.glob(file.path(d, "*.tsv"))
        q <- map_dfr(files, function(f) {
            m <- rd(f, col_names = FALSE)
            if (ncol(m) < 2) return(NULL)
            v <- setNames(as.list(m[[ncol(m)]]), m[[1]])
            as_tibble(v) |> mutate(sample = sub("\\.tsv$", "", basename(f)))
        })
    }
    q <- need(q, c("# contigs", "Total length", "N50"), "quast")
    if (!is.null(q)) {
        add_raw("quast", q)
        add_wide(tibble(sample = clean_id(q$sample),
                        contigs = suppressWarnings(as.numeric(q[["# contigs"]])),
                        assembly_length = suppressWarnings(as.numeric(q[["Total length"]])),
                        n50 = suppressWarnings(as.numeric(q[["N50"]]))))
    }
    report("quast", if (is.null(q)) 0 else nrow(q), if (is.null(q)) 0 else n_distinct(q$sample),
           "all metrics in tool table", d)
}
if (want("checkm2")) {
    d <- src_dir("checkm2")
    c2 <- NULL
    if (!is.na(d)) {
        f <- Sys.glob(file.path(d, "*quality_report.tsv")); if (!length(f)) f <- Sys.glob(file.path(d, "*.tsv"))
        if (length(f)) c2 <- rd(f[1])
    }
    c2 <- need(c2, c("Completeness", "Contamination"), "checkm2")
    if (!is.null(c2)) {
        id <- if ("sample" %in% names(c2)) c2$sample else c2$Name
        add_raw("checkm2", c2)
        add_wide(tibble(sample = clean_id(id),
                        completeness = suppressWarnings(as.numeric(c2$Completeness)),
                        contamination = suppressWarnings(as.numeric(c2$Contamination))))
    }
    report("checkm2", if (is.null(c2)) 0 else 1, if (is.null(c2)) 0 else nrow(c2), where = d)
}
if (want("kleborate")) {
    d <- src_dir("kleborate"); kb <- NULL
    if (!is.na(opt$results) && file.exists(file.path(opt$results, "03_typing", "kleborate_all.tsv"))) {
        kb <- rd(file.path(opt$results, "03_typing", "kleborate_all.tsv"))
    } else {
        kb <- read_each(d, "*.kleborate.txt", "\\.kleborate\\.txt")
        if (is.null(kb)) kb <- read_each(d, "*.txt", "\\.txt")
    }
    if (!is.null(kb) && nrow(kb)) {
        if (!"sample" %in% names(kb)) kb$sample <- kb[[1]]
        add_raw("kleborate", kb)
        w <- tibble(sample = clean_id(kb$sample),
                    kleborate_species = col_or_na(kb, "__species$"), ST = col_or_na(kb, "mlst__ST$"),
                    kleborate_K_locus = col_or_na(kb, "K_locus$"), kleborate_O_type = col_or_na(kb, "O_type$"),
                    resistance_score = col_or_na(kb, "resistance_score$"),
                    virulence_score = col_or_na(kb, "virulence_score$"))
        for (p in multi$kleborate_extra) {
            hit <- grep(p, names(kb), value = TRUE)
            if (length(hit)) w[[hit[1]]] <- as.character(kb[[hit[1]]])
        }
        add_wide(w)
    }
    report("kleborate", if (is.null(kb)) 0 else nrow(kb), if (is.null(kb)) 0 else n_distinct(clean_id(kb$sample)), where = d)
}
if (want("kaptive")) {
    d <- src_dir("kaptive"); seen <- 0; nf <- 0
    for (L in c("K", "O")) {
        kp <- NULL
        all_tab <- if (!is.na(opt$results)) file.path(opt$results, "03_typing", paste0("kaptive_", L, "_all.tsv")) else NA
        if (!is.na(all_tab) && file.exists(all_tab)) kp <- rd(all_tab)
        else kp <- read_each(d, paste0("*.kaptive_", L, ".tsv"), paste0("\\.kaptive_", L, "\\.tsv"))
        kp <- need(kp, c("Best match locus", "Best match type", "Match confidence"), paste0("kaptive ", L))
        if (is.null(kp)) next
        if (!"sample" %in% names(kp)) kp$sample <- kp[[1]]
        add_raw(paste0("kaptive_", L), kp)
        w <- tibble(sample = clean_id(kp$sample), a = kp[["Best match locus"]],
                    b = kp[["Best match type"]], c = kp[["Match confidence"]])
        names(w)[-1] <- paste0("kaptive_", L, "_", c("locus", "type", "confidence"))
        add_wide(w); seen <- max(seen, n_distinct(w$sample)); nf <- nf + nrow(kp)
    }
    report("kaptive", nf, seen, "K and O read separately", d)
}
amr_long <- NULL
if (want("amrfinderplus")) {
    d <- src_dir("amrfinderplus"); a <- NULL
    if (!is.na(opt$results) && file.exists(file.path(opt$results, "04_amr", "amrfinderplus_all.tsv"))) {
        a <- rd(file.path(opt$results, "04_amr", "amrfinderplus_all.tsv"))
    } else a <- read_each(d, "*.tsv", "\\.tsv")
    if (!is.null(a) && nrow(a)) {
        if (!"sample" %in% names(a) && "Name" %in% names(a)) a <- rename(a, sample = "Name")
        gcol <- intersect(c("Element symbol", "Gene symbol"), names(a))
        tcol <- intersect(c("Element type", "Type"), names(a))
        if (!length(gcol) || !"sample" %in% names(a)) {
            note("amrfinderplus: no gene or sample column - skipped")
        } else {
            add_raw("amrfinderplus", a)
            amr_long <- tibble(sample = clean_id(a$sample), gene = as.character(a[[gcol[1]]]),
                               element_type = if (length(tcol)) as.character(a[[tcol[1]]]) else NA_character_,
                               class = if ("Class" %in% names(a)) as.character(a$Class) else NA_character_,
                               subclass = if ("Subclass" %in% names(a)) as.character(a$Subclass) else NA_character_,
                               contig = if ("Contig id" %in% names(a)) as.character(a[["Contig id"]]) else NA_character_)
            w <- amr_long |>
                group_by(sample) |>
                summarise(amr_n_genes = sum(element_type == "AMR", na.rm = TRUE),
                          amr_genes   = paste(sort(unique(gene[element_type == "AMR"])), collapse = ";"),
                          amr_classes = paste(sort(unique(class[element_type == "AMR" & !is.na(class)])), collapse = ";"),
                          haystack    = paste(paste(class, subclass), collapse = " | "), .groups = "drop")
            for (nm in names(amr_spec)) {
                w[[nm]] <- as.integer(str_detect(w$haystack, regex(amr_spec[[nm]], ignore_case = TRUE)))
            }
            w$haystack <- NULL
            add_wide(w)
        }
    }
    report("amrfinderplus", if (is.null(a)) 0 else nrow(a),
           if (is.null(amr_long)) 0 else n_distinct(amr_long$sample),
           paste0(if (is.null(amr_long)) 0 else nrow(amr_long), " gene hits"), d)
}
plasmid_raw <- NULL
if (want("mobsuite")) {
    d <- src_dir("mobsuite")
    if (!is.na(opt$results) && file.exists(file.path(opt$results, "05_plasmids", "mobtyper_all.tsv"))) {
        plasmid_raw <- rd(file.path(opt$results, "05_plasmids", "mobtyper_all.tsv"))
    } else if (!is.na(d)) {
        files <- Sys.glob(file.path(d, "*", "mobtyper_results.txt"))
        plasmid_raw <- map_dfr(files, function(f) { x <- rd(f); if (nrow(x)) x$sample <- basename(dirname(f)); x })
    }
    if (!is.null(plasmid_raw) && nrow(plasmid_raw) && "sample" %in% names(plasmid_raw)) {
        plasmid_raw <- mutate(plasmid_raw, sample = clean_id(sample))
        add_raw("mobsuite", plasmid_raw)
        reps  <- if ("rep_type(s)" %in% names(plasmid_raw)) plasmid_raw[["rep_type(s)"]] else NA_character_
        sizes <- if ("size" %in% names(plasmid_raw)) suppressWarnings(as.numeric(plasmid_raw$size)) else NA_real_
        mob   <- if ("predicted_mobility" %in% names(plasmid_raw)) str_to_lower(plasmid_raw$predicted_mobility) else NA_character_
        clus  <- if ("primary_cluster_id" %in% names(plasmid_raw)) plasmid_raw$primary_cluster_id else NA_character_
        add_wide(tibble(sample = plasmid_raw$sample, rep = reps, size = sizes, mob = mob, cl = clus) |>
            group_by(sample) |>
            summarise(n_plasmids = n(), plasmid_total_bp = sum(size, na.rm = TRUE),
                      plasmid_rep_types = paste(sort(unique(str_trim(unlist(str_split(rep, ","))))), collapse = ";"),
                      plasmid_clusters  = paste(sort(unique(cl[!is.na(cl)])), collapse = ";"),
                      n_plasmids_conjugative     = sum(mob == "conjugative", na.rm = TRUE),
                      n_plasmids_mobilizable     = sum(mob == "mobilizable", na.rm = TRUE),
                      n_plasmids_non_mobilizable = sum(mob == "non-mobilizable", na.rm = TRUE),
                      .groups = "drop"))
    }
    report("mobsuite", if (is.null(plasmid_raw)) 0 else nrow(plasmid_raw),
           if (is.null(plasmid_raw)) 0 else n_distinct(plasmid_raw$sample),
           paste0(if (is.null(plasmid_raw)) 0 else nrow(plasmid_raw), " plasmids"), d)
}
for (tool in c("bakta", "prokka")) {
    if (!want(tool)) next
    d <- src_dir(tool); files <- if (is.na(d)) character() else Sys.glob(file.path(d, "*", "*.txt"))
    if (!length(files) && !is.na(d)) files <- Sys.glob(file.path(d, "*.txt"))
    tab <- map_dfr(files, function(f) {
        ln <- readLines(f, warn = FALSE); ln <- ln[str_detect(ln, ":")]
        if (!length(ln)) return(NULL)
        key <- str_to_lower(str_replace_all(str_trim(str_split_fixed(ln, ":", 2)[, 1]), " ", "_"))
        val <- suppressWarnings(as.numeric(str_trim(str_split_fixed(ln, ":", 2)[, 2])))
        s <- if (basename(dirname(f)) != basename(d)) basename(dirname(f)) else sub("\\.txt$", "", basename(f))
        keep <- key %in% c("cds", "cdss", "trna", "trnas", "rrna", "rrnas", "contigs", "count")
        if (!any(keep)) return(NULL)
        as_tibble(setNames(as.list(val[keep]), paste0(tool, "_", key[keep]))) |> mutate(sample = s)
    })
    if (!is.null(tab) && nrow(tab)) { add_wide(tab); add_raw(tool, tab) }
    report(tool, length(files), if (is.null(tab)) 0 else n_distinct(tab$sample), "feature counts only", d)
}
snp_long <- NULL
if (want("phylogeny")) {
    d <- src_dir("phylogeny")
    mats <- if (is.na(d)) character() else c(Sys.glob(file.path(d, "*", "*.snp_distances.tsv")),
                                             Sys.glob(file.path(d, "*.snp_distances.tsv")))
    near <- map_dfr(mats, function(f) {
        m <- rd(f); if (!nrow(m)) return(NULL)
        ids <- as.character(m[[1]]); v <- as.matrix(m[, -1, drop = FALSE])
        storage.mode(v) <- "numeric"; rownames(v) <- ids
        keep <- intersect(setdiff(ids, "Reference"), colnames(v))
        if (length(keep) < 2) return(NULL)
        v <- v[keep, keep, drop = FALSE]; diag(v) <- NA
        g <- sub("\\.snp_distances\\.tsv$", "", basename(f))
        pairs <- which(upper.tri(v), arr.ind = TRUE)
        snp_long <<- bind_rows(snp_long, tibble(group = g, sample1 = keep[pairs[, 1]],
                                                sample2 = keep[pairs[, 2]],
                                                snps = v[upper.tri(v)], method = "snp-dists"))
        tibble(sample = keep, phylo_group = g, phylo_group_n = length(keep),
               nearest_sample = colnames(v)[apply(v, 1, which.min)],
               nearest_snps = apply(v, 1, min, na.rm = TRUE),
               n_within = apply(v, 1, function(r) sum(r <= snp_threshold, na.rm = TRUE)))
    })
    if (!is.null(near) && nrow(near)) {
        names(near)[names(near) == "n_within"] <- paste0("n_within_", snp_threshold, "_snps")
        add_wide(near)
    }
    for (f in if (is.na(d)) character() else Sys.glob(file.path(d, "*", "*.ska_distances.tsv"))) {
        s <- rd(f)
        if (nrow(s) && all(c("Sample1", "Sample2") %in% names(s))) {
            snp_long <- bind_rows(snp_long, tibble(group = basename(dirname(f)), sample1 = s$Sample1,
                                                   sample2 = s$Sample2,
                                                   snps = if ("Mismatches" %in% names(s)) s$Mismatches else s$Distance,
                                                   method = "ska distance"))
        }
    }
    report("phylogeny", length(mats), if (is.null(near)) 0 else nrow(near),
           paste0(if (is.null(snp_long)) 0 else nrow(snp_long), " pairs"), d)
}

# ------------------------------------------------------------------------------------------------ join
expected <- character()
if (!is.na(opt$samplesheet)) {
    sep <- if (grepl("\\.csv$", opt$samplesheet)) "," else "\t"
    ss  <- suppressWarnings(read_delim(opt$samplesheet, delim = sep, show_col_types = FALSE, progress = FALSE))
    idc <- intersect(c("sample", "id", "Sample"), names(ss))
    idc <- if (length(idc)) idc[1] else names(ss)[1]
    expected <- unique(clean_id(ss[[idc]]))
}
ids <- unique(unlist(lapply(wide, function(d) d$sample)))
ids <- sort(ids[!is.na(ids)])
missing <- setdiff(expected, ids)
extra   <- if (length(expected)) setdiff(ids, expected) else character()
if (length(missing)) note(length(missing), " samplesheet isolate(s) have no results at all: ",
                          paste(head(missing, 10), collapse = ", "))
if (length(extra))   note(length(extra), " isolate(s) in the results are not in the samplesheet: ",
                          paste(head(extra, 10), collapse = ", "))
universe <- if (length(expected)) expected else ids
if (!length(universe)) note("no isolates found - check --results / --dir (see --list-sources)")

master <- tibble(sample = universe)
for (d in wide) master <- left_join(master, d, by = "sample")
if (!is.na(opt$metadata)) {
    sep <- if (grepl("\\.csv$", opt$metadata)) "," else "\t"
    md  <- suppressWarnings(read_delim(opt$metadata, delim = sep, show_col_types = FALSE, progress = FALSE))
    if (!opt$id_col %in% names(md)) stop("--id-col '", opt$id_col, "' is not a column of --metadata")
    md <- md |> rename(sample = all_of(opt$id_col)) |> mutate(sample = clean_id(sample)) |>
        distinct(sample, .keep_all = TRUE)
    master <- left_join(master, md, by = "sample")
    report("metadata", 1, sum(master$sample %in% md$sample), paste0("ID column '", opt$id_col, "'"))
}
if ("n_plasmids" %in% names(master)) master <- mutate(master, n_plasmids = replace_na(n_plasmids, 0))
ordered <- c(intersect(LEADING, names(master)), setdiff(names(master), LEADING))
master  <- master[, ordered]

written <- character()
put <- function(d, path, label) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    write_tsv(d, path, na = opt$na)
    written <<- c(written, path)
    say("  ", sprintf("%-30s %6d rows  ->  %s", label, nrow(d), path))
}
say("Isolates: ", nrow(master), "   Columns: ", ncol(master))
put(master, file.path(opt$outdir, paste0(opt$prefix, ".tsv")), "master_table (isolates)")

if (flags$tool_tables) {
    for (key in names(TOOL_TABLES)) {
        d <- raw[[key]]
        if (is.null(d) || !nrow(d)) next
        d <- relocate(d, sample)
        name <- TOOL_TABLES[[key]][1]; srcnm <- TOOL_TABLES[[key]][2]
        target <- file.path(opt$outdir, name)
        if (flags$beside_source) {
            sd <- src_dir(srcnm)
            if (!is.na(sd)) target <- file.path(dirname(sd), name)
            else note("--beside-source: no directory for ", srcnm, ", writing ", name, " to --outdir")
        }
        put(d, target, paste0(sub("\\.tsv$", "", name), " (tool table)"))
    }
}
if (flags$long) {
    if (!is.null(amr_long))    put(amr_long, file.path(opt$outdir, "amr_genes_long.tsv"), "amr_genes_long (gene hits)")
    if (!is.null(plasmid_raw)) put(plasmid_raw, file.path(opt$outdir, "plasmids_long.tsv"), "plasmids_long (plasmids)")
    if (!is.null(snp_long))    put(snp_long, file.path(opt$outdir, "snp_distances_long.tsv"), "snp_distances_long (pairs)")
}

if (flags$summaries) {
    n_iso     <- max(nrow(master), 1)
    flag_cols <- intersect(names(amr_spec), names(master))
    flag_cols <- intersect(c("esbl_or_ampc", "carbapenemase", "gentamicin_r"), flag_cols)
    pct_of    <- function(d, col) {
        v <- d[[col]][!is.na(d[[col]])]
        if (!length(v)) NA_real_ else round(100 * sum(as.integer(v)) / length(v), 1)
    }

    ov <- tibble(isolates = nrow(master))
    if ("ST" %in% names(master)) ov$with_ST <- sum(!is.na(master$ST))
    for (f in flag_cols) ov[[paste0(f, "_pct")]] <- pct_of(master, f)
    if ("n_plasmids" %in% names(master)) ov$median_plasmids <- median(master$n_plasmids, na.rm = TRUE)
    put(ov, file.path(opt$outdir, "summary_overall.tsv"), "summary_overall")

    for (col in str_trim(str_split(opt$group_by, ",")[[1]])) {
        if (!nzchar(col)) next
        if (!col %in% names(master)) { note("--group-by '", col, "' is not a column; skipped"); next }
        g <- master |>
            mutate(.g = replace_na(as.character(.data[[col]]), opt$na)) |>
            group_by(.g) |>
            summarise(n = n(),
                      across(all_of(flag_cols), ~ round(100 * mean(as.integer(.x), na.rm = TRUE), 1),
                             .names = "{.col}_pct"),
                      mean_plasmids = if ("n_plasmids" %in% names(master))
                          round(mean(n_plasmids, na.rm = TRUE), 1) else NA_real_,
                      .groups = "drop") |>
            arrange(desc(n)) |>
            rename(!!col := ".g")
        put(g, file.path(opt$outdir, paste0("summary_by_", col, ".tsv")), paste0("summary_by_", col))
    }

    if (!is.null(amr_long)) {
        put(amr_long |> filter(element_type == "AMR") |> distinct(sample, gene, class) |>
                count(gene, class, sort = TRUE, name = "isolates") |>
                mutate(percent = round(100 * isolates / n_iso, 1)),
            file.path(opt$outdir, "summary_amr_genes.tsv"), "summary_amr_genes")
    }
    spec <- c(kleborate_K_locus = "summary_k_locus.tsv", kleborate_O_type = "summary_o_type.tsv")
    for (cn in names(spec)) {
        if (!cn %in% names(master)) next
        put(master |> count(.data[[cn]], sort = TRUE, name = "isolates") |>
                mutate(percent = round(100 * isolates / n_iso, 1)),
            file.path(opt$outdir, spec[[cn]]), sub("\\.tsv$", "", spec[[cn]]))
    }
    if (!is.null(plasmid_raw) && "primary_cluster_id" %in% names(plasmid_raw)) {
        put(count(plasmid_raw, primary_cluster_id, sort = TRUE, name = "plasmids"),
            file.path(opt$outdir, "summary_plasmid_clusters.tsv"), "summary_plasmid_clusters")
    }
}

# ------------------------------------------------------------------------------------------------ summary
sum_path <- file.path(opt$outdir, "aggregate_summary.txt")
con <- file(sum_path, "w")
writeLines(c("aggregate_results.R summary", strrep("=", 70), "",
             paste0("Isolates in ", opt$prefix, ".tsv: ", nrow(master)), ""), con)
if (length(rows)) {
    rep_tbl <- bind_rows(rows)
    writeLines(sprintf("%-16s %6s %8s  %-26s %s", "source", "files", "samples", "note", "directory"), con)
    for (i in seq_len(nrow(rep_tbl))) {
        writeLines(sprintf("%-16s %6d %8d  %-26s %s", rep_tbl$source[i], rep_tbl$files[i],
                           rep_tbl$samples[i], rep_tbl$note[i], rep_tbl$directory[i]), con)
    }
}
writeLines(c("", "Tables written:", paste0("  ", written)), con)
writeLines(c("", if (length(problems)) c(paste0("Problems (", length(problems), "):"), paste0("  - ", problems))
             else "No problems reported."), con)
close(con)
say("  ", sprintf("%-30s %6s        ->  %s", "summary", "", sum_path))

if (flags$strict && length(missing)) {
    message("ERROR: --strict and ", length(missing), " isolate(s) missing")
    quit(status = 1)
}
