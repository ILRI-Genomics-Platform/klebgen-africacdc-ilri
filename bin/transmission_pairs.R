#!/usr/bin/env Rscript
#
# transmission_pairs.R - every pair of isolates in a phylogeny group, with SNP distance, ward and year,
# and a transmission classification (surveillance question 5).
#
#   bin/run_R.sh transmission_pairs.R \
#       --dist     results/06_phylogeny               # a folder (all groups) or one *.snp_distances.tsv
#       --metadata training_data/metadata_compact.tsv
#
# Options (defaults in brackets):
#   --dist PATH        snp-dists matrix, or a folder searched for */*.snp_distances.tsv      [required]
#   --metadata FILE    TSV or CSV, one row per isolate                                        [required]
#   --out FILE         output TSV                                                             [transmission_pairs.tsv]
#   --id-col NAME      sample ID column                                                       [sample]
#   --ward-col NAME    ward / unit column                                                     [ward]
#   --date-col NAME    collection date or year column (2017 or 2017-03-14)                   [collection_year]
#   --close N          at most N SNPs: genomically close                                     [10]
#   --near N           at most N SNPs: genomically near                                      [21]
#   --max-years N      at most N years apart counts as "close in time"                       [1]
#
# Classification of a pair:
#   probable transmission   close (<= --close SNPs), same ward, close in time
#   possible transmission   near (<= --near SNPs), same ward, close in time
#   linked, route unclear   close or near, but different ward, unknown ward, or years apart
#   unlinked                more than --near SNPs
# The thresholds are guides to agree on before looking at the data, not rules.

# ------------------------------------------------------------------------------------------------ options
args <- commandArgs(trailingOnly = TRUE)
opt <- list(dist = NA, metadata = NA, out = "transmission_pairs.tsv", id_col = "sample", ward_col = "ward",
            date_col = "collection_year", close = 10, near = 21, max_years = 1)
if (any(args %in% c("-h", "--help"))) {
    lines <- readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
    cat(sub("^# ?", "", lines[2:(which(lines == "")[1] - 1)]), sep = "\n")
    quit(status = 0)
}
if (length(args) %% 2 != 0) stop("Options come in pairs: --name value (see --help)")
for (k in seq(1, max(1, length(args)), by = 2)[seq_len(length(args) %/% 2)]) {
    key <- gsub("-", "_", sub("^--", "", args[k]))
    if (!key %in% names(opt)) stop("Unknown option: ", args[k], " (see --help)")
    opt[[key]] <- args[k + 1]
}
if (is.na(opt$dist) || is.na(opt$metadata)) stop("--dist and --metadata are required (see --help)")
for (k in c("close", "near", "max_years")) opt[[k]] <- as.numeric(opt[[k]])

# ------------------------------------------------------------------------------------------------ metadata
sep <- if (grepl("\\.csv$", opt$metadata, ignore.case = TRUE)) "," else "\t"
meta <- read.table(opt$metadata, sep = sep, header = TRUE, quote = "\"", comment.char = "",
                   colClasses = "character", check.names = FALSE)
for (col in c(opt$id_col, opt$ward_col, opt$date_col)) {
    if (!col %in% names(meta)) stop("Column '", col, "' not in ", opt$metadata, " (columns: ", paste(names(meta), collapse = ", "), ")")
}
# Year as a number from '2017', '2017.5' or '2017-03-14'; NA if missing
meta$.year <- suppressWarnings(as.numeric(substr(meta[[opt$date_col]], 1, 4)))
meta$.ward <- trimws(meta[[opt$ward_col]])
meta$.ward[is.na(meta$.ward) | meta$.ward == "" | toupper(meta$.ward) %in% c("UNKNOWN", "NA", "-")] <- NA

# ------------------------------------------------------------------------------------------------ distance matrices
files <- if (dir.exists(opt$dist)) {
    list.files(opt$dist, pattern = "\\.snp_distances\\.tsv$", recursive = TRUE, full.names = TRUE)
} else {
    opt$dist
}
if (length(files) == 0) stop("No *.snp_distances.tsv found in ", opt$dist)

pairs_of <- function(file) {
    m <- as.matrix(read.delim(file, row.names = 1, check.names = FALSE))   # first header cell is 'snp-dists x.y'
    keep <- setdiff(rownames(m), "Reference")                              # the mapping reference is not an isolate
    m <- m[keep, keep, drop = FALSE]
    if (nrow(m) < 2) return(NULL)
    idx <- which(upper.tri(m), arr.ind = TRUE)                             # each pair once
    data.frame(group = sub("\\.snp_distances\\.tsv$", "", basename(file)),
               isolate_1 = rownames(m)[idx[, 1]], isolate_2 = colnames(m)[idx[, 2]],
               SNPs = as.numeric(m[idx]), stringsAsFactors = FALSE)
}
pairs <- do.call(rbind, lapply(files, pairs_of))
if (is.null(pairs) || nrow(pairs) == 0) stop("No isolate pairs found")

# ------------------------------------------------------------------------------------------------ join metadata and classify
r1 <- match(pairs$isolate_1, meta[[opt$id_col]])
r2 <- match(pairs$isolate_2, meta[[opt$id_col]])
no_meta <- unique(c(pairs$isolate_1[is.na(r1)], pairs$isolate_2[is.na(r2)]))
if (length(no_meta)) cat("WARNING: no metadata for", paste(no_meta, collapse = ", "), "\n")

pairs$year_1 <- meta$.year[r1]
pairs$year_2 <- meta$.year[r2]
pairs$ward_1 <- meta$.ward[r1]
pairs$ward_2 <- meta$.ward[r2]
pairs$years_apart <- abs(pairs$year_1 - pairs$year_2)
pairs$same_ward <- !is.na(pairs$ward_1) & !is.na(pairs$ward_2) & pairs$ward_1 == pairs$ward_2
pairs$close_in_time <- !is.na(pairs$years_apart) & pairs$years_apart <= opt$max_years

pairs$genomic <- ifelse(pairs$SNPs <= opt$close, "close", ifelse(pairs$SNPs <= opt$near, "near", "distant"))
epi <- pairs$same_ward & pairs$close_in_time
pairs$classification <- ifelse(pairs$genomic == "distant", "unlinked",
                        ifelse(epi & pairs$genomic == "close", "probable transmission",
                        ifelse(epi & pairs$genomic == "near", "possible transmission",
                               "linked, route unclear")))

pairs <- pairs[order(pairs$group, pairs$SNPs, pairs$isolate_1, pairs$isolate_2), ]
write.table(pairs, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ------------------------------------------------------------------------------------------------ summary
cat(sprintf("%d pairs in %d group(s) -> %s\n", nrow(pairs), length(unique(pairs$group)), opt$out))
print(table(group = pairs$group, classification = pairs$classification))
linked <- pairs[pairs$classification != "unlinked", c("group", "isolate_1", "isolate_2", "SNPs", "ward_1", "ward_2",
                                                      "year_1", "year_2", "classification")]
if (nrow(linked)) {
    cat("\nPairs within", opt$near, "SNPs:\n")
    print(linked, row.names = FALSE)
} else {
    cat("\nNo pairs within", opt$near, "SNPs.\n")
}
