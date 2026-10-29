#!/usr/bin/env Rscript
#
# make_surveillance_figures.R - two figures that answer two standing questions in AMR surveillance.
#
#   A  Lineage composition through time   -> is the dominant clone turning over?
#   C  Within- vs between-outbreak SNP distances, per outbreak
#                                         -> do SNP distances actually resolve outbreaks here?
#   D  Empiric regimen coverage           -> would first-line therapy still work here?
#   E  Resistance markers by year         -> is resistance rising?
#   F  Drug classes per isolate           -> how multidrug-resistant are these isolates?
#   G  AMR gene emergence timeline        -> which genes are new, and which are established?
#
# A and C keep the lettering of the layout they were modelled on (there is no panel B); D-G are the
# public-health panels added on top. Each is skipped, with a message, when its inputs are missing, so
# a run without AST or without a year column still produces the rest.
#
# These complement bin/make_figures.R rather than repeating it: that script draws resistance
# prevalence, the gene-to-species-to-ward sunburst, the plasmid Sankey, annotated trees and the
# per-isolate gene heatmap. Nothing here duplicates those.
#
# Panel A is the question behind Heinz et al. 2024's "repeated replacement of dominant clones": a
# surveillance programme that reports only the current dominant ST cannot see replacement coming.
# Panel C is the Day 7 question: a SNP threshold is only useful if within-outbreak distances sit
# clearly below between-outbreak ones, and this panel shows whether they do, outbreak by outbreak.
#
#   bin/run_R.sh make_surveillance_figures.R --results results \
#       --metadata training_data/metadata_compact.tsv
#
# Options (defaults in brackets):
#   --results DIR       pipeline output folder (the run's --outdir)                     [required]
#   --master FILE       master table from bin/aggregate_results           [<results>/09_aggregate/
#                                                                           master_table.tsv]
#   --metadata FILE     TSV/CSV merged in when --master lacks the time or outbreak columns  [none]
#   --id-col NAME       sample ID column in --metadata                                    [sample]
#   --lineage-col NAME  column holding the lineage for panel A (ST, kleborate_K_locus ...)     [ST]
#   --time-col NAME     column holding the year                                 [collection_year]
#   --outbreak-cols LIST  columns pasted together to define an outbreak    [ward,collection_year]
#   --periods LIST      comma-separated cut years for panel A, e.g. 2009,2015; by default each
#                       year is its own period when there are 8 or fewer, else tertiles     [auto]
#   --top-n N           lineages coloured individually in panel A; the rest become "Other"     [10]
#   --min-group N       smallest outbreak group kept in panel C                                [3]
#   --snp-threshold N   cluster threshold drawn as a reference line in panel C                [20]
#   --xmax N            x limit for panel C; by default the 99th percentile of between-outbreak
#                       distances, so one very distant pair cannot squash the panel          [auto]
#   --ast FILE          TSV/CSV of S/I/R calls, one column per antibiotic. When given, panel D uses
#                       MEASURED susceptibility; without it, panel D falls back to genotype       [none]
#   --amr-long FILE     gene-level table for panel G    [<results>/09_aggregate/amr_genes_long.tsv]
#   --regimen-by NAME   split panel D's bars by this column (kleborate_species, ward ...)    [none]
#   --date-col NAME     column with a full date (YYYY-MM-DD or YYYY-MM), used when panels E and G
#                       are drawn at month resolution        [collection_date, then date, if present]
#   --time-resolution R  auto | year | month, for panels E and G. 'auto' picks month when the
#                       isolates span one calendar year or less and a date column parses     [auto]
#   --min-year-n N      time bins with fewer isolates than this are dropped from panel E       [5]
#   --genes-n N         genes shown in panel G, commonest first                                 [25]
#   --font NAME         typeface for every panel; falls back to the device default when it is not
#                       installed (NEONATAL_FONT sets it too)                       [Arial Narrow]
#   --outdir DIR        where the figures and tables are written           [<results>/09_figures]
#
# Figures use the themes from Wilke's "Fundamentals of Data Visualization" (bin/dviz_theme.R):
# a horizontal grid when the value is on y, a vertical grid when it is on x.
#
# Note on panel C: SNP distances only exist between isolates in the SAME phylogeny group, because
# each group is aligned to its own reference. "Between outbreaks" therefore means a different
# ward-year WITHIN one phylogeny group, not across species or STs. The counts printed by the script
# state how many groups contributed.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "r_setup.R"))
use_packages(cran = c("ggplot2", "tidyverse"))
use_packages(cran = "cowplot", required = FALSE)      # the dviz themes; they fall back without it
suppressPackageStartupMessages(library(tidyverse))
source(file.path(script_dir, "dviz_theme.R"))         # theme_dviz_open / _hgrid / _vgrid

# ------------------------------------------------------------------------------------------------ options
args <- commandArgs(trailingOnly = TRUE)
opt  <- list(results = NA, master = NA, metadata = NA, id_col = "sample", lineage_col = "ST",
             time_col = "collection_year", outbreak_cols = "ward,collection_year", periods = NA,
             top_n = "10", min_group = "3", snp_threshold = "20", xmax = NA, outdir = NA,
             ast = NA, amr_long = NA, regimen_by = NA, min_year_n = "5", genes_n = "25",
             font = "Arial Narrow", date_col = NA, time_resolution = "auto")

if (any(args %in% c("-h", "--help"))) {
    lines <- readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
    cat(sub("^# ?", "", lines[2:(which(lines == "")[1] - 1)]), sep = "\n")
    quit(status = 0)
}
if (length(args) %% 2 != 0) stop("Options come in pairs: --name value (see --help)")
for (k in seq_len(length(args) %/% 2)) {
    key <- gsub("-", "_", sub("^--", "", args[2 * k - 1]))
    if (!key %in% names(opt)) stop("Unknown option: ", args[2 * k - 1], " (see --help)")
    opt[[key]] <- args[2 * k]
}
if (is.na(opt$results)) stop("--results is required (see --help)")
if (!dir.exists(opt$results)) stop("--results is not a directory: ", opt$results)
if (is.na(opt$master)) opt$master <- file.path(opt$results, "09_aggregate", "master_table.tsv")
if (is.na(opt$outdir)) opt$outdir <- file.path(opt$results, "09_figures")
if (is.na(opt$amr_long)) opt$amr_long <- file.path(opt$results, "09_aggregate", "amr_genes_long.tsv")
dviz_use_font(opt$font)
min_year_n <- as.integer(opt$min_year_n)
genes_n    <- as.integer(opt$genes_n)
top_n     <- as.integer(opt$top_n)
min_group <- as.integer(opt$min_group)
snp_thr   <- as.numeric(opt$snp_threshold)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)
say <- function(...) cat(..., "\n", sep = "")
# grDevices::pdf() writes a single-byte encoding and cannot carry anything outside it: a ward or
# species name with an accent, or an en dash in the data, warns once per text element
# ("mbcsToSbcs: - substituted for ...") and is silently replaced. cairo_pdf handles UTF-8.
# Panel labels in this script are ASCII on purpose, so this only matters for text coming from data.
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
save_fig <- function(plot, name, width, height) {
    ggsave(file.path(opt$outdir, paste0(name, ".pdf")), plot, width = width, height = height,
           device = pdf_device, bg = "white", limitsize = FALSE)
    ggsave(file.path(opt$outdir, paste0(name, ".png")), plot, width = width, height = height,
           dpi = 300, bg = "white", limitsize = FALSE)
    say("  wrote ", file.path(opt$outdir, paste0(name, ".{pdf,png}")))
}
put <- function(d, name) {
    write_tsv(d, file.path(opt$outdir, name), na = "NA")
    say("  wrote ", file.path(opt$outdir, name), "  (", nrow(d), " rows)")
}

# ------------------------------------------------------------------------------------------------ look
# A qualitative palette, hand-coded: paletteer is not among the course's R packages.
PALETTE <- c("#5050FF", "#CE3D32", "#749B58", "#E4AF69", "#466983", "#BA6338", "#5DB1DD",
             "#802268", "#6BD76B", "#D595A7", "#924822", "#837B8D")
OTHER   <- "#C9C9C9"
BLUE <- "#0072B2"; VERM <- "#D55E00"
INK  <- "#1A1A1A"; INK2 <- "#4D4D4D"; INK3 <- "#808080"

# A stacked bar's in-bar label must read against whatever colour it sits on, so pick its ink from
# the fill's relative luminance rather than always using white.
ink_on <- function(hex) {
    m   <- grDevices::col2rgb(hex) / 255
    lum <- 0.2126 * m[1, ] + 0.7152 * m[2, ] + 0.0722 * m[3, ]
    unname(ifelse(lum > 0.55, "#111111", "white"))
}
# The dviz themes (bin/dviz_theme.R) carry the look; this adds only what is common to every panel.
# 'grid' picks the variant: a grid line is only useful along the axis the VALUE is on, so a plot with
# values on y gets "hgrid" and a horizontal plot with values on x gets "vgrid".
theme_fig <- function(grid = c("hgrid", "vgrid", "open"), base_size = 10) {
    grid <- match.arg(grid)
    base <- switch(grid,
                   hgrid = theme_dviz_hgrid(font_size = base_size),
                   vgrid = theme_dviz_vgrid(font_size = base_size),
                   open  = theme_dviz_open(font_size = base_size))
    base +
        theme(text                = element_text(family = dviz_font_family),
              plot.title          = element_text(family = dviz_font_family_bold, face = "bold",
                                                 size = base_size + 1, hjust = 0),
              plot.title.position = "plot",
              axis.title   = element_text(size = base_size - 1, colour = INK),
              axis.text    = element_text(size = base_size - 2, colour = INK2),
              legend.position = "bottom", legend.title = element_blank(),
              legend.text = element_text(size = base_size - 2, colour = INK2),
              legend.key.size = unit(9, "pt"), legend.margin = margin(t = -2),
              plot.margin = margin(6, 10, 6, 6))
}

# ------------------------------------------------------------------------------------------------ data
if (!file.exists(opt$master)) {
    stop("no master table at ", opt$master, "\n  run: bin/aggregate_results.py --results ", opt$results,
         "  (or pass --master)")
}
say("Reading ", opt$master)
master <- suppressWarnings(read_tsv(opt$master, show_col_types = FALSE, progress = FALSE))
if (!"sample" %in% names(master)) stop("the master table has no 'sample' column")

if (!is.na(opt$metadata)) {
    sep <- if (grepl("\\.csv$", opt$metadata)) "," else "\t"
    md  <- suppressWarnings(read_delim(opt$metadata, delim = sep, show_col_types = FALSE, progress = FALSE))
    if (!opt$id_col %in% names(md)) stop("--id-col '", opt$id_col, "' is not a column of --metadata")
    md <- md |> rename(sample = all_of(opt$id_col)) |> distinct(sample, .keep_all = TRUE)
    md <- md[, c("sample", setdiff(names(md), names(master)))]      # do not duplicate existing columns
    master <- left_join(master, md, by = "sample")
    say("  merged ", ncol(md) - 1, " column(s) from --metadata")
}

ob_cols <- str_trim(str_split(opt$outbreak_cols, ",")[[1]])
need_cols <- unique(c(opt$lineage_col, opt$time_col, ob_cols))
missing   <- setdiff(need_cols, names(master))
if (length(missing)) {
    stop("the master table has no column(s): ", paste(missing, collapse = ", "),
         "\n  available: ", paste(names(master), collapse = ", "))
}
master <- master |>
    mutate(.lin  = as.character(.data[[opt$lineage_col]]),
           .time = suppressWarnings(as.numeric(.data[[opt$time_col]])),
           .outbreak = apply(across(all_of(ob_cols)), 1, function(r)
               if (any(is.na(r))) NA_character_ else paste(r, collapse = " ")))
say("  ", nrow(master), " isolates; ", sum(!is.na(master$.lin)), " with ", opt$lineage_col,
    "; ", sum(!is.na(master$.time)), " with ", opt$time_col)

# ---- time resolution for panels E and G -------------------------------------------------------
# A collection whose isolates all come from one year gives a single point per gene in panel G and a
# single point per marker in panel E, which says nothing. When a real date is available the two
# panels drop to month resolution instead. Year-only values ("2017") cannot support that and are
# left as years. Both panels bin into a Date (.bin) so they share one x scale.
parse_partial_date <- function(x) {
    x   <- gsub("/", "-", trimws(as.character(x)))
    out <- as.Date(rep(NA_character_, length(x)))
    full <- grepl("^\\d{4}-\\d{1,2}-\\d{1,2}", x)
    out[full] <- as.Date(substr(x[full], 1, 10), format = "%Y-%m-%d")
    mon <- !full & grepl("^\\d{4}-\\d{1,2}$", x)
    out[mon] <- as.Date(paste0(x[mon], "-01"), format = "%Y-%m-%d")
    out
}
date_col <- if (!is.na(opt$date_col)) {
    if (!opt$date_col %in% names(master)) stop("--date-col '", opt$date_col, "' is not a column")
    opt$date_col
} else {
    cand <- intersect(c("collection_date", "date"), names(master))
    if (length(cand)) cand[1] else NA_character_
}
master$.date <- if (is.na(date_col)) {
    as.Date(rep(NA_character_, nrow(master)))
} else {
    parse_partial_date(master[[date_col]])
}
n_dated <- sum(!is.na(master$.date))
n_years <- n_distinct(master$.time[!is.na(master$.time)])

time_res <- opt$time_resolution
if (!time_res %in% c("auto", "year", "month")) stop("--time-resolution must be auto, year or month")
if (time_res == "auto") {
    span_m   <- if (n_dated >= 2) as.numeric(diff(range(master$.date, na.rm = TRUE))) / 30.44 else NA_real_
    time_res <- if (n_dated >= 2 && (n_years <= 1 || (!is.na(span_m) && span_m <= 24))) "month" else "year"
}
if (time_res == "month" && n_dated < 2) {
    warning("--time-resolution month needs a parseable date column (",
            if (is.na(date_col)) "none found" else date_col, "); using year", call. = FALSE)
    time_res <- "year"
}
master$.bin <- if (time_res == "month") {
    as.Date(format(master$.date, "%Y-%m-01"))
} else {
    as.Date(ifelse(is.na(master$.time), NA, paste0(master$.time, "-01-01")))
}
bin_span <- if (sum(!is.na(master$.bin)) >= 2)
    as.numeric(diff(range(master$.bin, na.rm = TRUE))) / 30.44 else 0
time_breaks <- if (time_res == "month") { if (bin_span <= 14) "1 month" else "3 months" } else "1 year"
time_labels <- if (time_res == "month") "%b %Y" else "%Y"
scale_time  <- function() scale_x_date(date_labels = time_labels, date_breaks = time_breaks)
say("  time resolution for panels E and G: ", time_res,
    if (time_res == "month") paste0(" (", n_dated, " isolates with a parseable date from '", date_col, "')")
    else if (!is.na(date_col)) paste0(" ('", date_col, "' present but ", n_years, " year(s) spanned)") else "")

# ================================================================================ Panel A: composition
dA <- master |> filter(!is.na(.lin), !is.na(.time))
pA <- NULL
if (nrow(dA) < 2) {
    warning("panel A: too few isolates with both ", opt$lineage_col, " and ", opt$time_col, call. = FALSE)
} else {
    years <- sort(unique(dA$.time))
    if (!is.na(opt$periods)) {
        cuts <- sort(as.numeric(str_trim(str_split(opt$periods, ",")[[1]])))
        labs <- c(paste0("<=", cuts[1]),
                  if (length(cuts) > 1) paste0(cuts[-length(cuts)] + 1, "-", cuts[-1]),
                  paste0(">", cuts[length(cuts)]))
        dA$period <- cut(dA$.time, breaks = c(-Inf, cuts, Inf), labels = labs)
    } else if (length(years) <= 8) {
        dA$period <- factor(dA$.time, levels = years)
        say("  panel A: one period per year (", length(years), " years)")
    } else {
        cuts <- unique(round(as.numeric(quantile(dA$.time, c(1/3, 2/3), na.rm = TRUE))))
        labs <- c(paste0("<=", cuts[1]),
                  if (length(cuts) > 1) paste0(cuts[-length(cuts)] + 1, "-", cuts[-1]),
                  paste0(">", cuts[length(cuts)]))
        dA$period <- cut(dA$.time, breaks = c(-Inf, cuts, Inf), labels = labs)
        say("  panel A: tertile periods at ", paste(cuts, collapse = ", "),
            " (", length(years), " distinct years)")
    }

    top <- dA |> count(.lin, sort = TRUE) |> slice_head(n = top_n) |> pull(.lin)
    dA  <- dA |> mutate(lineage = if_else(.lin %in% top, .lin, "Other"))
    cols <- setNames(PALETTE[seq_along(top)], top)
    if (any(dA$lineage == "Other")) cols <- c(cols, Other = OTHER)

    compA <- dA |>
        count(period, lineage, name = "n") |>
        complete(period, lineage = names(cols), fill = list(n = 0)) |>
        group_by(period) |>
        mutate(total = sum(n), prop = if_else(total > 0, n / total, 0)) |>
        ungroup()
    put(compA, "panelA_composition.tsv")

    # position_stack() puts the FIRST level on top, so order levels by ascending abundance and use
    # reverse = TRUE to stack them bottom-up; "Other" is forced to the very bottom.
    ab <- compA |> group_by(lineage) |> summarise(n = sum(n), .groups = "drop") |> arrange(n)
    lv <- c("Other", setdiff(ab$lineage, "Other"))
    lv <- lv[lv %in% names(cols)]
    compA <- compA |>
        mutate(lineage = factor(lineage, levels = lv),
               lab_ink = ink_on(cols[as.character(lineage)]))
    labsA <- compA |> distinct(period, total) |>
        mutate(lab = sprintf("%s\n(n=%d)", period, total)) |> arrange(period)

    pA <- ggplot(compA, aes(period, prop, fill = lineage)) +
        geom_col(width = 0.74, colour = "white", linewidth = 0.4,
                 position = position_stack(reverse = TRUE)) +
        geom_text(aes(label = if_else(prop >= 0.10, as.character(lineage), ""), colour = lab_ink),
                  position = position_stack(vjust = 0.5, reverse = TRUE),
                  fontface = "bold", size = 2.6, show.legend = FALSE) +
        scale_colour_identity() +
        scale_fill_manual(values = cols, breaks = c(setdiff(lv, "Other"), "Other"), drop = FALSE) +
        scale_x_discrete(labels = setNames(labsA$lab, labsA$period)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.02))) +
        guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
        labs(title = paste0("A   Temporal shifts in ", opt$lineage_col, " composition"),
             x = NULL, y = "Proportion of isolates") +
        theme_fig("hgrid")
    save_fig(pA, "panelA_composition", 5.6, 4.2)

    # the question the panel exists to answer, as a number
    shifts <- compA |>
        group_by(lineage) |>
        summarise(first = prop[which.min(as.integer(period))],
                  last  = prop[which.max(as.integer(period))], .groups = "drop") |>
        mutate(change_pct_points = round(100 * (last - first), 1)) |>
        arrange(desc(abs(change_pct_points)))
    put(shifts, "panelA_change.tsv")
    say("  largest shift: ", shifts$lineage[1], " ", sprintf("%+.1f", shifts$change_pct_points[1]),
        " percentage points between the first and last period")
}

# ================================================================== Panel C: within vs between outbreaks
mats <- c(Sys.glob(file.path(opt$results, "06_phylogeny", "*", "*.snp_distances.tsv")),
          Sys.glob(file.path(opt$results, "06_phylogeny", "*.snp_distances.tsv")))
pC <- NULL
if (!length(mats)) {
    warning("panel C: no *.snp_distances.tsv under ", file.path(opt$results, "06_phylogeny"), call. = FALSE)
} else {
    say("Reading ", length(mats), " SNP matri(ces)")
    pairs <- map_dfr(mats, function(f) {
        m <- suppressWarnings(read_tsv(f, show_col_types = FALSE, progress = FALSE))
        if (!nrow(m)) return(NULL)
        ids <- as.character(m[[1]])
        v   <- as.matrix(m[, -1, drop = FALSE])
        storage.mode(v) <- "numeric"
        rownames(v) <- ids
        keep <- intersect(setdiff(ids, "Reference"), colnames(v))
        if (length(keep) < 2) return(NULL)
        v   <- v[keep, keep, drop = FALSE]
        idx <- which(upper.tri(v), arr.ind = TRUE)
        tibble(phylo_group = sub("\\.snp_distances\\.tsv$", "", basename(f)),
               t1 = keep[idx[, 1]], t2 = keep[idx[, 2]], snps = v[upper.tri(v)])
    })
    grp <- setNames(master$.outbreak, master$sample)
    lin <- setNames(master$.lin, master$sample)
    sizes <- table(grp[!is.na(grp)])
    elig  <- names(sizes)[sizes >= min_group]
    say("  outbreak groups (", opt$outbreak_cols, ", n>=", min_group, "): ", length(elig),
        " covering ", sum(sizes[elig]), " isolates")

    pairs <- pairs |>
        filter(!is.na(snps)) |>
        mutate(g1 = grp[t1], g2 = grp[t2],
               comparison = case_when(!is.na(g1) & !is.na(g2) & g1 == g2 & g1 %in% elig ~ "Within outbreak",
                                      !is.na(g1) & !is.na(g2) & g1 != g2                ~ "Between outbreaks",
                                      TRUE ~ NA_character_)) |>
        filter(!is.na(comparison)) |>
        mutate(comparison = factor(comparison, levels = c("Within outbreak", "Between outbreaks")),
               group = if_else(comparison == "Within outbreak", g1, NA_character_))
    put(pairs |> select(phylo_group, comparison, group, snps), "panelC_distances.tsv")

    statsC <- pairs |> group_by(comparison) |>
        summarise(n_pairs = n(), median = median(snps), p25 = quantile(snps, .25),
                  p75 = quantile(snps, .75), .groups = "drop")
    put(statsC, "panelC_stats.tsv")
    print(statsC)

    if (!any(pairs$comparison == "Within outbreak")) {
        warning("panel C: no within-outbreak pairs; check --outbreak-cols and --min-group", call. = FALSE)
    } else {
        mb <- statsC$median[statsC$comparison == "Between outbreaks"]
        b25 <- statsC$p25[statsC$comparison == "Between outbreaks"]
        b75 <- statsC$p75[statsC$comparison == "Between outbreaks"]

        # A ward-year holding more than one lineage is not one outbreak; flagging them separates
        # genuine single-introduction groups from co-circulation.
        lin_by_group <- tibble(sample = names(grp), group = unname(grp)) |>
            filter(!is.na(group)) |>
            mutate(lineage = lin[sample]) |>
            group_by(group) |>
            summarise(n_isolates = n(), n_lineages = n_distinct(lineage[!is.na(lineage)]),
                      lineages = paste(sort(unique(lineage[!is.na(lineage)])), collapse = ";"),
                      .groups = "drop")

        obC <- pairs |> filter(comparison == "Within outbreak") |>
            group_by(group) |>
            summarise(n_pairs = n(), median_snps = median(snps), p25 = quantile(snps, .25),
                      p75 = quantile(snps, .75), min_snps = min(snps), max_snps = max(snps),
                      .groups = "drop") |>
            left_join(lin_by_group, by = "group") |>
            mutate(multi = n_lineages > 1) |>
            arrange(median_snps) |>
            mutate(label = factor(group, levels = group))
        put(obC, "panelC_outbreaks.tsv")

        med_single <- if (any(!obC$multi)) median(obC$median_snps[!obC$multi]) else NA_real_
        med_multi  <- if (any(obC$multi))  median(obC$median_snps[obC$multi])  else NA_real_
        xmax <- if (!is.na(opt$xmax)) as.numeric(opt$xmax)
                else max(c(obC$max_snps, as.numeric(quantile(pairs$snps[pairs$comparison == "Between outbreaks"],
                                                             0.99, na.rm = TRUE))), na.rm = TRUE) * 1.05

        pC <- ggplot(obC, aes(y = label, colour = multi)) +
            annotate("rect", xmin = b25, xmax = b75, ymin = -Inf, ymax = Inf, fill = VERM, alpha = 0.11) +
            geom_vline(xintercept = mb, colour = VERM, linewidth = 0.35, linetype = "42") +
            geom_vline(xintercept = snp_thr, colour = INK3, linewidth = 0.35, linetype = "22") +
            geom_linerange(aes(xmin = min_snps, xmax = max_snps), linewidth = 0.35, alpha = 0.45) +
            geom_linerange(aes(xmin = p25, xmax = p75), linewidth = 1.1) +
            geom_point(aes(x = median_snps), size = 1.8, shape = 21, fill = "white", stroke = 0.7) +
            annotate("text", x = mb, y = 0.6, hjust = -0.05, vjust = 1, size = 2.5, colour = VERM,
                     lineheight = 1.2,
                     label = sprintf("between-outbreak median %.0f (IQR shaded)", mb)) +
            annotate("text", x = snp_thr, y = length(levels(obC$label)) + 0.4, hjust = -0.05,
                     vjust = 0, size = 2.5, colour = INK3,
                     label = sprintf("threshold %g", snp_thr)) +
            scale_colour_manual(values = c(`FALSE` = BLUE, `TRUE` = VERM),
                                labels = c(`FALSE` = sprintf("One %s (n=%d), median %.0f SNPs",
                                                             opt$lineage_col, sum(!obC$multi), med_single),
                                           `TRUE` = sprintf("Two or more (n=%d), median %.0f SNPs",
                                                            sum(obC$multi), med_multi))) +
            scale_x_continuous(limits = c(-xmax * 0.02, xmax), expand = expansion(mult = c(0, 0))) +
            scale_y_discrete(expand = expansion(add = c(0.9, 1.8))) +
            guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1.3, shape = NA))) +
            labs(title = "C   Resolution within individual outbreaks",
                 x = "Pairwise core SNP differences", y = NULL) +
            theme_fig("vgrid") +
            theme(axis.text.y = element_text(size = 7.5, colour = INK), legend.margin = margin(t = 2))
        save_fig(pC, "panelC_outbreak_resolution", 6.2, max(3.2, 0.26 * nrow(obC) + 2.2))

        over <- sum(obC$median_snps > snp_thr)
        say("  ", over, " of ", nrow(obC), " outbreak group(s) have a median above the ",
            snp_thr, "-SNP threshold")
    }
}

# ======================================================================= Panel D: regimen coverage
# Which antibiotic columns stand for each regimen when --ast is given. Matching is case-insensitive
# and a regimen with no matching column is skipped.
REGIMEN_AST <- list("Ampicillin + gentamicin"        = c("gentamicin"),
                    "Third-generation cephalosporin" = c("ceftriaxone", "cefotaxime", "ceftazidime"),
                    "Carbapenem"                     = c("meropenem", "imipenem", "ertapenem"),
                    "Amikacin"                       = c("amikacin"))
# Without AST, the master table's genotype flags stand in. An isolate counts as covered when the
# flag for that regimen is NOT set. Ampicillin is intrinsically ineffective against the K. pneumoniae
# species complex, so the WHO first-line regimen rests on gentamicin alone and is scored on that flag.
REGIMEN_GENO <- list("Ampicillin + gentamicin"        = "gentamicin_r",
                     "Third-generation cephalosporin" = "esbl_or_ampc",
                     "Carbapenem"                     = "carbapenemase")

wilson <- function(k, n, z = 1.96) {
    if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
    p <- k / n; d <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / d
    half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    c(100 * (centre - half), 100 * (centre + half))
}

pD <- NULL
ast <- NULL
if (!is.na(opt$ast)) {
    sep <- if (grepl("\\.csv$", opt$ast)) "," else "\t"
    ast <- suppressWarnings(read_delim(opt$ast, delim = sep, show_col_types = FALSE, progress = FALSE))
    idc <- intersect(c(opt$id_col, "sample", "id"), names(ast))
    if (!length(idc)) { warning("--ast has no sample column; panel D falls back to genotype", call. = FALSE); ast <- NULL }
    else ast <- rename(ast, sample = all_of(idc[1]))
}
split_col <- if (!is.na(opt$regimen_by) && opt$regimen_by %in% names(master)) opt$regimen_by else NA
if (!is.na(opt$regimen_by) && is.na(split_col))
    warning("--regimen-by '", opt$regimen_by, "' is not a column; panel D is not split", call. = FALSE)

covD <- NULL
if (!is.null(ast)) {
    lc <- setNames(names(ast), tolower(names(ast)))
    d  <- left_join(master |> select(sample, any_of(na.omit(split_col))), ast, by = "sample")
    for (rn in names(REGIMEN_AST)) {
        cols <- unname(lc[intersect(names(lc), REGIMEN_AST[[rn]])])
        if (!length(cols)) next
        covered <- apply(d[, cols, drop = FALSE], 1,
                         function(r) any(toupper(trimws(r)) == "S", na.rm = TRUE))
        tested  <- apply(d[, cols, drop = FALSE], 1,
                         function(r) any(toupper(trimws(r)) %in% c("S", "I", "R"), na.rm = TRUE))
        covD <- bind_rows(covD, tibble(regimen = rn, sample = d$sample, covered = covered, tested = tested,
                                       split = if (is.na(split_col)) "All isolates" else as.character(d[[split_col]])))
    }
    src_lab <- "measured susceptibility (--ast)"
} else {
    for (rn in names(REGIMEN_GENO)) {
        fl <- REGIMEN_GENO[[rn]]
        if (!fl %in% names(master)) next
        covD <- bind_rows(covD, tibble(regimen = rn, sample = master$sample,
                                       covered = as.integer(master[[fl]]) == 0,
                                       tested  = !is.na(master[[fl]]),
                                       split = if (is.na(split_col)) "All isolates" else as.character(master[[split_col]])))
    }
    src_lab <- "predicted from resistance genes, NOT measured"
}
if (is.null(covD) || !nrow(covD)) {
    warning("panel D: no regimen could be scored (need --ast, or the genotype flag columns)", call. = FALSE)
} else {
    sumD <- covD |> filter(tested) |>
        group_by(regimen, split) |>
        summarise(n = n(), covered = sum(covered, na.rm = TRUE), .groups = "drop") |>
        mutate(percent = round(100 * covered / n, 1),
               lo = mapply(function(k, m) wilson(k, m)[1], covered, n),
               hi = mapply(function(k, m) wilson(k, m)[2], covered, n),
               regimen = factor(regimen, levels = rev(names(if (is.null(ast)) REGIMEN_GENO else REGIMEN_AST))))
    put(sumD, "panelD_regimen_coverage.tsv")

    pD <- ggplot(sumD, aes(percent, regimen, fill = split)) +
        geom_col(width = 0.6, position = position_dodge(width = 0.65)) +
        geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.4, colour = INK2,
                       position = position_dodge(width = 0.65)) +
        geom_text(aes(label = sprintf("%.0f%% (%d/%d)", percent, covered, n)),
                  position = position_dodge(width = 0.65), hjust = -0.08, size = 2.6, colour = INK2) +
        scale_fill_manual(values = setNames(PALETTE[seq_len(n_distinct(sumD$split))],
                                            sort(unique(sumD$split))), guide =
                              if (is.na(split_col)) "none" else guide_legend(nrow = 2)) +
        scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                           expand = expansion(mult = c(0, 0))) +
        labs(title = "D   Empiric regimen coverage",
             subtitle = paste0("Isolates each regimen would cover - ", src_lab),
             x = "Isolates covered (%), with 95% CI", y = NULL) +
        theme_fig("vgrid") +
        theme(plot.subtitle = element_text(size = 7.5, colour = INK2, margin = margin(b = 4)))
    save_fig(pD, "panelD_regimen_coverage", 6.4, 3.4)
    say("  panel D scored from ", src_lab)
}

# ========================================================================== Panel E: resistance by year
pE <- NULL
flagE <- intersect(c("esbl_or_ampc", "carbapenemase", "gentamicin_r"), names(master))
if (!length(flagE) || all(is.na(master$.time))) {
    warning("panel E: needs ", opt$time_col, " and the genotype flag columns", call. = FALSE)
} else {
    trend <- master |> filter(!is.na(.bin)) |>
        select(sample, .bin, all_of(flagE)) |>
        pivot_longer(all_of(flagE), names_to = "marker", values_to = "v") |>
        filter(!is.na(v)) |>
        group_by(period = .bin, marker) |>
        summarise(n = n(), k = sum(as.integer(v)), .groups = "drop") |>
        filter(n >= min_year_n) |>
        mutate(percent = 100 * k / n,
               lo = mapply(function(a, b) wilson(a, b)[1], k, n),
               hi = mapply(function(a, b) wilson(a, b)[2], k, n),
               marker = recode(marker, esbl_or_ampc = "ESBL / AmpC",
                               carbapenemase = "Carbapenemase", gentamicin_r = "Gentamicin"))
    if (!nrow(trend)) {
        warning("panel E: no ", time_res, " bin has at least --min-year-n (", min_year_n,
                ") isolates", call. = FALSE)
    } else {
        put(trend, "panelE_resistance_trend.tsv")
        pE <- ggplot(trend, aes(period, percent, colour = marker, group = marker)) +
            geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.35, alpha = 0.7) +
            geom_line(linewidth = 0.6) +
            geom_point(size = 1.9) +
            scale_colour_manual(values = c("ESBL / AmpC" = PALETTE[1], "Carbapenemase" = VERM,
                                           "Gentamicin" = PALETTE[3])) +
            scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
            scale_time() +
            labs(title = "E   Resistance markers through time",
                 subtitle = sprintf("by %s, bins with at least %d isolates; 95%% CI",
                                    time_res, min_year_n),
                 x = NULL, y = "Isolates carrying the marker (%)") +
            theme_fig("hgrid") +
            theme(plot.subtitle = element_text(size = 7.5, colour = INK2, margin = margin(b = 4)),
                  legend.margin = margin(t = 2))
        save_fig(pE, "panelE_resistance_trend", 6.4, 3.8)
    }
}

# ===================================================================== Panel F: drug classes per isolate
pF <- NULL
if (!"amr_classes" %in% names(master)) {
    warning("panel F: the master table has no amr_classes column", call. = FALSE)
} else {
    sp_col <- if ("kleborate_species" %in% names(master)) "kleborate_species" else NA
    burden <- master |>
        mutate(n_classes = if_else(is.na(amr_classes) | amr_classes == "", 0L,
                                   lengths(str_split(amr_classes, ";"))),
               species = if (is.na(sp_col)) "All isolates" else as.character(.data[[sp_col]])) |>
        mutate(species = replace_na(species, "unknown"))
    put(burden |> count(species, n_classes, name = "isolates"), "panelF_class_burden.tsv")
    med <- burden |> group_by(species) |> summarise(median_classes = median(n_classes), n = n(), .groups = "drop")
    pF <- ggplot(burden, aes(factor(n_classes), fill = species)) +
        geom_bar(width = 0.78, colour = "white", linewidth = 0.3) +
        scale_fill_manual(values = setNames(PALETTE[seq_len(n_distinct(burden$species))],
                                            sort(unique(burden$species)))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(title = "F   Acquired resistance classes per isolate",
             subtitle = paste(sprintf("%s: median %g (n=%d)", med$species, med$median_classes, med$n),
                              collapse = "   |   "),
             x = "Number of AMRFinderPlus drug classes", y = "Isolates") +
        theme_fig("hgrid") +
        theme(plot.subtitle = element_text(size = 7.5, colour = INK2, margin = margin(b = 4)))
    save_fig(pF, "panelF_class_burden", 6.0, 3.6)
}

# ====================================================================== Panel G: gene emergence timeline
pG <- NULL
if (!file.exists(opt$amr_long)) {
    say("  panel G: no ", opt$amr_long, " (run aggregate_results with --long); skipped")
} else if (all(is.na(master$.time))) {
    warning("panel G: needs ", opt$time_col, call. = FALSE)
} else {
    al <- suppressWarnings(read_tsv(opt$amr_long, show_col_types = FALSE, progress = FALSE))
    if (!all(c("sample", "gene") %in% names(al))) {
        warning("panel G: ", opt$amr_long, " has no sample/gene columns", call. = FALSE)
    } else {
        if ("element_type" %in% names(al)) al <- filter(al, element_type == "AMR")
        emg <- al |>
            left_join(master |> select(sample, .bin), by = "sample") |>
            filter(!is.na(.bin)) |>
            group_by(gene) |>
            summarise(first_seen = min(.bin), last_seen = max(.bin),
                      isolates = n_distinct(sample), .groups = "drop") |>
            arrange(desc(isolates)) |>
            slice_head(n = genes_n) |>
            arrange(first_seen, desc(isolates)) |>
            mutate(gene = factor(gene, levels = gene),
                   span_days = as.numeric(last_seen - first_seen))
        if (!nrow(emg)) {
            warning("panel G: no gene hits with a usable date", call. = FALSE)
        } else {
            put(emg, "panelG_gene_emergence.tsv")
            if (all(emg$span_days == 0)) {
                warning("panel G: every gene's first and last detection fall in the same ", time_res,
                        " bin, so the segments have no length. Pass --date-col, or ",
                        "--time-resolution month if a date column is present.", call. = FALSE)
            }
            pG <- ggplot(emg, aes(y = gene)) +
                geom_segment(aes(x = first_seen, xend = last_seen, yend = gene),
                             linewidth = 0.5, colour = INK3) +
                geom_point(aes(x = first_seen, size = isolates), colour = VERM) +
                geom_point(aes(x = last_seen), size = 1.1, colour = INK2) +
                scale_size_continuous(range = c(1.4, 4.4), name = "Isolates") +
                scale_time() +
                labs(title = "G   When each resistance gene was seen",
                     subtitle = sprintf("%d commonest genes, by %s; red = first detection, grey = last",
                                        nrow(emg), time_res),
                     x = NULL, y = NULL) +
                theme_fig("vgrid") +
                theme(plot.subtitle = element_text(size = 7.5, colour = INK2, margin = margin(b = 4)),
                      axis.text.y = element_text(size = 7, colour = INK),
                      legend.position = "right", legend.title = element_text(size = 7.5, colour = INK2),
                      legend.margin = margin(t = 2))
            save_fig(pG, "panelG_gene_emergence", 5.8, max(3.0, 0.22 * nrow(emg) + 1.8))
        }
    }
}

# ------------------------------------------------------------------------------------------------ compose
if (!is.null(pD) && !is.null(pE) && has("patchwork")) {
    save_fig(patchwork::wrap_plots(pD, pE, ncol = 2, widths = c(1, 1.05)), "Figure_amr_burden", 12, 4)
}
if (!is.null(pA) && !is.null(pC)) {
    if (has("patchwork")) {
        fig <- patchwork::wrap_plots(pA, pC, ncol = 2, widths = c(1, 1.1))
        save_fig(fig, "Figure_surveillance", 12, 5.6)
    } else {
        say("  patchwork is not installed: the panels are saved separately, not composed")
    }
}
say("\nDone: ", opt$outdir)
