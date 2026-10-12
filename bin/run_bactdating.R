#!/usr/bin/env Rscript
#
# run_bactdating.R - time-scaled tree for one phylogeny group, from the pipeline's Gubbins output.
#
# Uses the recombination-filtered Gubbins tree (branch lengths in substitutions, plus the share of each branch
# not affected by recombination), so recombination is accounted for as in Didelot et al. 2018.
#
#   bin/run_R.sh run_bactdating.R \
#       --gubbins  results/06_phylogeny/gubbins/Kpneumoniae_ST39 \
#       --metadata training_data/metadata_compact.tsv \
#       --mrca     CAA65Z,CAA99X           # optional: date the common ancestor of these isolates
#
# Options (defaults in brackets):
#   --gubbins PREFIX     Gubbins output prefix: PREFIX.final_tree.tre, PREFIX.node_labelled.final_tree.tre and
#                        PREFIX.per_branch_statistics.csv must exist                       [required]
#   --metadata FILE      TSV or CSV with one row per isolate                               [required]
#   --id-col NAME        column with sample IDs (as in the tree)                          [sample]
#   --date-col NAME      column with collection dates: 2017, 2017-03 or 2017-03-14;
#                        a year or month alone is treated as a range, not a point         [collection_year]
#   --outdir DIR         output folder                                                    [bactdating_<group>]
#   --models LIST        comma-separated BactDating models; with more than one, each is
#                        run and the lowest DIC is reported (e.g. arc,poisson,relaxedgamma) [arc]
#   --iterations N       MCMC iterations per model                                         [1e6]
#   --mrca IDS           comma-separated sample IDs whose common ancestor to date (repeatable
#                        as IDS;IDS for several clusters)                                  [none]
#   --drop IDS           tips to remove before dating (comma-separated)                   [Reference]
#   --require-signal     stop if root-to-tip regression shows no temporal signal (p >= 0.05)
#   --seed N             random seed                                                      [1]
#
# Needs R packages ape, coda and BactDating (GitHub); missing ones are installed into your own R library.
# Run it through bin/run_R.sh, which starts a clean R (see that script):
#   bin/run_R.sh run_bactdating.R --gubbins ... --metadata ...

# Library paths (course library, then ~/R/<platform>-library/<version>) and package installation: bin/r_setup.R
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "r_setup.R"))
use_packages(cran = c("ape", "coda"), github = c(BactDating = "xavierdidelot/BactDating"))
suppressPackageStartupMessages({
    library(ape)
    library(BactDating)
})

# ------------------------------------------------------------------------------------------------ options
args <- commandArgs(trailingOnly = TRUE)
opt <- list(
    gubbins = NA, metadata = NA, id_col = "sample", date_col = "collection_year", outdir = NA,
    models = "arc", iterations = 1e6, mrca = NA, drop = "Reference", require_signal = FALSE, seed = 1
)
i <- 1
while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
        lines <- readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
        cat(sub("^# ?", "", lines[2:(which(lines == "")[1] - 1)]), sep = "\n")
        quit(status = 0)
    }
    if (a == "--require-signal") {
        opt$require_signal <- TRUE
        i <- i + 1
        next
    }
    if (!startsWith(a, "--") || i == length(args)) stop("Unexpected argument: ", a, " (see --help)")
    key <- gsub("-", "_", substring(a, 3))
    if (!key %in% names(opt)) stop("Unknown option: ", a, " (see --help)")
    opt[[key]] <- args[i + 1]
    i <- i + 2
}
if (is.na(opt$gubbins) || is.na(opt$metadata)) stop("--gubbins and --metadata are required (see --help)")
opt$iterations <- as.numeric(opt$iterations)
set.seed(as.integer(opt$seed))

group <- basename(opt$gubbins)
if (is.na(opt$outdir)) opt$outdir <- paste0("bactdating_", group)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(suffix) file.path(opt$outdir, paste0(group, suffix))

for (sfx in c(".final_tree.tre", ".per_branch_statistics.csv")) {
    if (!file.exists(paste0(opt$gubbins, sfx))) stop("Missing Gubbins file: ", opt$gubbins, sfx)
}

# ------------------------------------------------------------------------------------------------ tree
tree <- loadGubbins(opt$gubbins)                      # unrooted; branch lengths in substitutions; $unrec
drop <- intersect(strsplit(opt$drop, ",")[[1]], tree$tip.label)
if (length(drop)) tree <- drop.tip.useRec(tree, drop)   # keeps the recombination information in step
cat(sprintf("Group %s: %d isolates (removed: %s)\n", group, Ntip(tree), ifelse(length(drop), paste(drop, collapse = ", "), "none")))
if (Ntip(tree) < 3) stop("Fewer than 3 isolates: nothing to date")

# ------------------------------------------------------------------------------------------------ dates
sep <- if (grepl("\\.csv$", opt$metadata, ignore.case = TRUE)) "," else "\t"
meta <- read.table(opt$metadata, sep = sep, header = TRUE, quote = "\"", comment.char = "",
                   colClasses = "character", check.names = FALSE)
for (col in c(opt$id_col, opt$date_col)) {
    if (!col %in% names(meta)) stop("Column '", col, "' not in ", opt$metadata, " (columns: ", paste(names(meta), collapse = ", "), ")")
}

# Decimal-year range [lower, upper] for 2017, 2017.5, 2017-03 or 2017-03-14; NA if unknown
date_range <- function(x) {
    x <- trimws(x)
    if (is.na(x) || x == "" || toupper(x) %in% c("NA", "UNKNOWN", "-")) return(c(NA, NA))
    dec <- function(d) {
        y <- as.numeric(format(d, "%Y"))
        y + (as.numeric(format(d, "%j")) - 1) / ifelse((y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0, 366, 365)
    }
    if (grepl("^[0-9]{4}$", x)) return(as.numeric(x) + c(0, 0.999))
    if (grepl("^[0-9]{4}\\.[0-9]+$", x)) return(rep(as.numeric(x), 2))
    if (grepl("^[0-9]{4}-[0-9]{2}$", x)) {
        first <- as.Date(paste0(x, "-01"))
        last <- seq(first, by = "month", length.out = 2)[2] - 1
        return(c(dec(first), dec(last)))
    }
    if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)) return(rep(dec(as.Date(x)), 2))
    stop("Cannot read date '", x, "'; use 2017, 2017-03 or 2017-03-14")
}

rows <- match(tree$tip.label, meta[[opt$id_col]])
if (any(is.na(rows))) {
    cat("WARNING: no metadata for", paste(tree$tip.label[is.na(rows)], collapse = ", "), "- their dates will be estimated\n")
}
ranges <- t(sapply(seq_along(tree$tip.label), function(k) {
    if (is.na(rows[k])) c(NA, NA) else date_range(meta[[opt$date_col]][rows[k]])
}))
n_known <- sum(!is.na(ranges[, 1]))
if (n_known < 3) stop("Fewer than 3 isolates with a date: nothing to date")
mid <- rowMeans(ranges)                                # in tree tip order (no names)
if (diff(range(mid, na.rm = TRUE)) < 1) {
    cat("WARNING: all dates fall within one year; expect very wide intervals or no temporal signal\n")
}

# ------------------------------------------------------------------------------------------------ 1. temporal signal
rooted <- initRoot(tree, mid, useRec = TRUE)
pdf(out(".roottotip.pdf"), width = 9, height = 5)
rtt <- roottotip(rooted, mid, permTest = 10000)
invisible(dev.off())
signal <- !is.na(rtt$rate) && rtt$rate > 0 && rtt$pvalue < 0.05
cat(sprintf("Root-to-tip: rate %.2f substitutions/genome/year, p = %.4f (date permutation test) -> %s\n",
            rtt$rate, rtt$pvalue, ifelse(signal, "temporal signal", "NO clear temporal signal")))
if (!signal) {
    msg <- "No clear temporal signal: a dated tree from these data should not be reported as a date estimate."
    if (opt$require_signal) stop(msg, " (stopped because of --require-signal)")
    cat("WARNING:", msg, "Continuing so that you can inspect the output.\n")
}

# ------------------------------------------------------------------------------------------------ 2. BactDating
models <- trimws(strsplit(opt$models, ",")[[1]])
fits <- list()
for (m in models) {
    cat(sprintf("Running BactDating, model %s, %g iterations...\n", m, opt$iterations))
    fits[[m]] <- bactdate(tree, ranges, nbIts = opt$iterations, model = m, useRec = TRUE, showProgress = FALSE)
}
dic <- sapply(fits, function(f) f$dic)
best <- names(which.min(dic))
res <- fits[[best]]
if (length(models) > 1) {
    cat("DIC by model (lower is better; a difference under 5 is not meaningful):\n")
    print(round(sort(dic), 2))
}

# Posterior summaries (second half of the chain, as in BactDating's own print method)
post <- function(col) {
    v <- res$record[, col]
    v <- sort(v[(1 + floor(length(v) / 2)):length(v)])
    c(mean = mean(v), lower = v[max(1, floor(length(v) * 0.025))], upper = v[max(1, floor(length(v) * 0.975))])
}
root <- post(Ntip(res$tree) + 1)
mu <- post("mu")

ess <- c(mu = NA, alpha = NA)
if (requireNamespace("coda", quietly = TRUE)) {
    ess <- round(coda::effectiveSize(as.mcmc.resBactDating(res))[c("mu", "alpha")])
    if (any(ess < 200, na.rm = TRUE)) {
        cat("WARNING: effective sample size below 200 (mu ", ess["mu"], ", alpha ", ess["alpha"],
            "): run longer, e.g. --iterations ", format(opt$iterations * 10, scientific = TRUE), "\n", sep = "")
    }
}

# ------------------------------------------------------------------------------------------------ 3. common ancestors of chosen isolates
node_date <- function(node) {
    v <- res$record[, node]
    v <- sort(v[(1 + floor(length(v) / 2)):length(v)])
    c(mean = mean(v), lower = v[max(1, floor(length(v) * 0.025))], upper = v[max(1, floor(length(v) * 0.975))])
}
mrca_rows <- NULL
if (!is.na(opt$mrca)) {
    for (set in strsplit(opt$mrca, ";")[[1]]) {
        ids <- trimws(strsplit(set, ",")[[1]])
        missing <- setdiff(ids, res$tree$tip.label)
        if (length(missing)) {
            cat("WARNING: --mrca: not in the tree:", paste(missing, collapse = ", "), "\n")
            next
        }
        node <- getMRCA(res$tree, ids)
        d <- node_date(node)
        mrca_rows <- rbind(mrca_rows, data.frame(isolates = paste(ids, collapse = ","), mrca_date = round(d["mean"], 2),
                                                 lower_95 = round(d["lower"], 2), upper_95 = round(d["upper"], 2)))
        cat(sprintf("Common ancestor of %s: %.2f [%.2f; %.2f]\n", paste(ids, collapse = ", "), d["mean"], d["lower"], d["upper"]))
    }
}

# ------------------------------------------------------------------------------------------------ outputs
write.tree(res$tree, out(".dated_tree.nwk"))
saveRDS(list(result = res, all_models = fits, roottotip = rtt, tip_date_ranges = ranges), out(".bactdating.rds"))

pdf(out(".dated_tree.pdf"), width = 9, height = max(5, 0.25 * Ntip(res$tree)))
plot(res, "treeCI", show.tip.label = TRUE, cex = 0.7)
title(sprintf("%s: time-scaled tree (model %s); bars = 95%% credible intervals", group, best))
invisible(dev.off())
pdf(out(".mcmc_trace.pdf"), width = 9, height = 7)
plot(res, "trace")
invisible(dev.off())

summary_tab <- data.frame(
    group = group, isolates = Ntip(res$tree), isolates_with_date = n_known, removed = paste(drop, collapse = ","),
    model = best, iterations = opt$iterations,
    roottotip_rate = round(rtt$rate, 3), roottotip_pvalue = signif(rtt$pvalue, 3), temporal_signal = signal,
    root_date = round(root["mean"], 2), root_date_lower_95 = round(root["lower"], 2), root_date_upper_95 = round(root["upper"], 2),
    rate_subs_per_genome_per_year = signif(mu["mean"], 3), rate_lower_95 = signif(mu["lower"], 3), rate_upper_95 = signif(mu["upper"], 3),
    ess_mu = ess["mu"], ess_alpha = ess["alpha"], dic = round(res$dic, 2), row.names = NULL
)
write.table(summary_tab, out(".bactdating_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(mrca_rows)) write.table(mrca_rows, out(".mrca_dates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Root date %.2f [%.2f; %.2f]; rate %.2g [%.2g; %.2g] substitutions/genome/year (model %s)\n",
            root["mean"], root["lower"], root["upper"], mu["mean"], mu["lower"], mu["upper"], best))
cat("Outputs in", opt$outdir, ":", paste(basename(list.files(opt$outdir)), collapse = ", "), "\n")
