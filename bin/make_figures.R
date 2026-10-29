#!/usr/bin/env Rscript
#
# make_figures.R - publication-style figures from the pipeline's results, modelled on Sands et al.,
# Nat Microbiol 2021 (BARNARDS; doi:10.1038/s41564-021-00870-7):
#
#   1_resistance_bars    % resistant per antibiotic, coloured by species (from --ast; otherwise % of isolates
#                        carrying genes of each AMRFinderPlus drug class)                     (paper Fig. 3a)
#   2_betalactamase_sunburst  carbapenemase and ESBL genes (inner ring) -> species -> ward  (paper Fig. 3b)
#   3_plasmid_sankey     plasmid replicon carrying each carbapenemase/ESBL gene -> species -> gene
#                        (MOB-suite contig report joined to AMRFinderPlus contig IDs)        (paper Fig. 4)
#   4_tree_<group>       circular tree per phylogeny group with rings for ST, ward, year, K locus and genes with
#                        resistance-conferring point mutations (ompK36, gyrA ...), tip labels coloured by ward and
#                        symbols for carbapenemase genes                                    (paper Figs 5-6)
#   5_amr_heatmap        presence/absence of every acquired gene and point mutation per isolate, with species,
#                        ST and ward annotations (ComplexHeatmap, as in ILRI's visualizeAMR.R)
#
#   bin/run_R.sh make_figures.R --results results --metadata training_data/metadata_compact.tsv \
#       --ast training_data/ast_compact.tsv
#
# Options (defaults in brackets):
#   --results DIR      pipeline output folder (--outdir of the run)                            [required]
#   --metadata FILE    TSV/CSV, one row per isolate                                            [required]
#   --ast FILE         TSV/CSV of S/I/R calls, one column per antibiotic (optional)            [none]
#   --tree FILES       comma-separated tree files; default: every 06_phylogeny/*/*.treefile    [auto]
#   --outgroup ID      root the trees on this isolate; default: midpoint-like root of ape      [none]
#   --outdir DIR       where to write PDF and PNG files                                         [figures]
#   --id-col, --ward-col, --date-col   metadata columns                       [sample, ward, collection_year]
#
# R packages: ggplot2, ape (always); ggtree (Bioconductor) and ggnewscale for trees; ggalluvial for the
# Sankey diagram. Missing packages are installed into your own R library (bin/r_setup.R); a figure whose
# packages are still missing is skipped with a message, and the others are still drawn. Run it through
# bin/run_R.sh, which starts a clean R.

# Library paths (course library, then ~/R/<platform>-library/<version>) and package installation: bin/r_setup.R
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "r_setup.R"))
use_packages(cran = "ggplot2")
# For the trees and the Sankey diagram; a figure whose packages cannot be installed is skipped
use_packages(cran = c("ape", "ggnewscale", "ggalluvial", "cowplot"), bioc = c("ggtree", "ComplexHeatmap"), required = FALSE)
suppressPackageStartupMessages(library(ggplot2))
# Wilke's dviz themes; NEONATAL_FONT overrides the typeface, which falls back when not installed
source(file.path(script_dir, "dviz_theme.R"))

# ------------------------------------------------------------------------------------------------ options
args <- commandArgs(trailingOnly = TRUE)
opt <- list(results = NA, metadata = NA, ast = NA, tree = NA, outgroup = NA, outdir = "figures",
            id_col = "sample", ward_col = "ward", date_col = "collection_year")
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
if (is.na(opt$results) || is.na(opt$metadata)) stop("--results and --metadata are required (see --help)")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)
# grDevices::pdf() cannot encode characters outside its single-byte encoding (an accented ward or
# species name warns "mbcsToSbcs" and is replaced); cairo_pdf handles UTF-8.
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
save_fig <- function(plot, name, width, height) {
    ggsave(file.path(opt$outdir, paste0(name, ".pdf")), plot, width = width, height = height,
           device = pdf_device, bg = "white", limitsize = FALSE)
    ggsave(file.path(opt$outdir, paste0(name, ".png")), plot, width = width, height = height,
           dpi = 300, bg = "white", limitsize = FALSE)
    cat("  wrote", file.path(opt$outdir, paste0(name, ".{pdf,png}")), "\n")
}
read_table <- function(f) {
    read.table(f, sep = if (grepl("\\.csv$", f, ignore.case = TRUE)) "," else "\t", header = TRUE,
               quote = "\"", comment.char = "", colClasses = "character", check.names = FALSE, fill = TRUE)
}
first_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit)) df[[hit[1]]] else rep(NA_character_, nrow(df))
}
short_species <- function(s) {
    w <- strsplit(ifelse(is.na(s), "unknown", s), " ")
    vapply(w, function(x) if (length(x) >= 2) paste0(substr(x[1], 1, 1), ". ", x[2]) else x[1], "")
}

# ------------------------------------------------------------------------------------------------ inputs
meta <- read_table(opt$metadata)
for (col in c(opt$id_col, opt$ward_col, opt$date_col)) {
    if (!col %in% names(meta)) stop("Column '", col, "' not in ", opt$metadata)
}
meta <- data.frame(sample = meta[[opt$id_col]], ward = meta[[opt$ward_col]],
                   year = substr(meta[[opt$date_col]], 1, 4), stringsAsFactors = FALSE)
meta$ward[is.na(meta$ward) | meta$ward == "" | toupper(meta$ward) %in% c("NA", "-")] <- "unknown"

# Kleborate: species, ST and K locus per sample
kfiles <- list.files(file.path(opt$results, "03_typing", "kleborate"), "\\.kleborate\\.txt$", full.names = TRUE)
if (!length(kfiles)) stop("No Kleborate results in ", file.path(opt$results, "03_typing", "kleborate"))
typing <- do.call(rbind, lapply(kfiles, function(f) {
    k <- read_table(f)
    if (!nrow(k)) return(NULL)
    data.frame(sample = sub("\\.kleborate\\.txt$", "", basename(f)),
               species = first_col(k, c("enterobacterales__species__species", "species"))[1],
               ST = first_col(k, c("klebsiella_pneumo_complex__mlst__ST", "escherichia__mlst_achtman__ST",
                                   "klebsiella_oxytoca_complex__mlst__ST", "ST"))[1],
               K_locus = first_col(k, c("klebsiella_pneumo_complex__kaptive__K_locus", "K_locus"))[1],
               stringsAsFactors = FALSE)
}))
typing$species_short <- short_species(typing$species)
info <- merge(typing, meta, by = "sample", all.x = TRUE)
info$ward[is.na(info$ward)] <- "unknown"
cat(sprintf("%d isolates with Kleborate results; %d with metadata\n", nrow(info), sum(!is.na(info$year))))

# AMRFinderPlus: one row per gene hit
afiles <- list.files(file.path(opt$results, "04_amr", "amrfinderplus"), "\\.tsv$", full.names = TRUE)
afiles <- afiles[!grepl("-mutations\\.tsv$", afiles)]
amr <- do.call(rbind, lapply(afiles, function(f) {
    a <- read_table(f)
    if (!nrow(a)) return(NULL)
    data.frame(sample = sub("\\.tsv$", "", basename(f)),
               gene = first_col(a, c("Element symbol", "Gene symbol")),
               name = first_col(a, c("Element name", "Sequence name")),
               contig = first_col(a, c("Contig id")),
               type = first_col(a, c("Type", "Element type")),
               subtype = first_col(a, c("Subtype", "Element subtype")),
               class = first_col(a, c("Class")), subclass = first_col(a, c("Subclass")),
               stringsAsFactors = FALSE)
}))
if (is.null(amr)) amr <- data.frame(sample = character(), gene = character(), name = character(), contig = character(),
                                     type = character(), subtype = character(), class = character(), subclass = character())
amr <- amr[is.na(amr$type) | amr$type == "AMR", ]
# Acquired genes vs point mutations (--plus reports both). Porin mutations (ompK36 ...) are labelled CARBAPENEM
# too, so carbapenemase / ESBL calls are made on acquired genes only.
amr$point <- amr$subtype %in% "POINT"
amr$category <- ifelse(amr$point, NA,
                ifelse(grepl("CARBAPENEM", amr$subclass), "Carbapenemase",
                ifelse(grepl("extended-spectrum", amr$name, ignore.case = TRUE) | grepl("^blaCTX-M", amr$gene),
                       "ESBL", NA)))
cat(sprintf("%d AMRFinderPlus gene hits (%d carbapenemase, %d ESBL)\n", nrow(amr),
            sum(amr$category %in% "Carbapenemase"), sum(amr$category %in% "ESBL")))

# ------------------------------------------------------------------------------------------------ 1. resistance bars
fig_resistance <- function() {
    if (!is.na(opt$ast)) {
        ast <- read_table(opt$ast)
        names(ast)[names(ast) == opt$id_col] <- "sample"
        drugs <- setdiff(names(ast), "sample")
        long <- data.frame(sample = rep(ast$sample, length(drugs)), drug = rep(drugs, each = nrow(ast)),
                           call = toupper(trimws(unlist(ast[drugs], use.names = FALSE))), stringsAsFactors = FALSE)
        long <- long[long$call %in% c("S", "I", "R"), ]
        long$resistant <- long$call == "R"
        ylab <- "Resistant (% of isolates tested)"
        sub <- "Phenotypic AST (I counted as not resistant)"
    } else {
        cls <- sort(unique(amr$class[!is.na(amr$class) & amr$class != ""]))
        long <- expand.grid(sample = info$sample, drug = cls, stringsAsFactors = FALSE)
        hits <- unique(paste(amr$sample, amr$class))
        long$resistant <- paste(long$sample, long$drug) %in% hits
        ylab <- "Isolates carrying genes (%)"
        sub <- "Genotypic: AMRFinderPlus drug class (no --ast given)"
    }
    long <- merge(long, info[c("sample", "species_short")], by = "sample")
    agg <- aggregate(resistant ~ drug + species_short, data = long, FUN = function(x) c(pct = 100 * mean(x), n = length(x)))
    agg <- data.frame(drug = agg$drug, species = agg$species_short, pct = agg$resistant[, "pct"], n = agg$resistant[, "n"])
    p <- ggplot(agg, aes(x = drug, y = pct, fill = species)) +
        geom_col(position = position_dodge(preserve = "single"), width = 0.8) +
        scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
        labs(x = NULL, y = ylab, fill = "Species", subtitle = sub) +
        theme_dviz_hgrid(font_size = 11) +
        theme(text = element_text(family = dviz_font_family),
              axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(face = "italic"))
    save_fig(p, "1_resistance_bars", width = max(7, 0.45 * length(unique(agg$drug)) + 3), height = 5)
}

# ------------------------------------------------------------------------------------------------ 2. sunburst
fig_sunburst <- function() {
    bl <- unique(amr[!is.na(amr$category), c("sample", "gene", "category")])
    if (!nrow(bl)) return(cat("  skipped: no carbapenemase or ESBL genes found\n"))
    bl <- merge(bl, info[c("sample", "species_short", "ward")], by = "sample")
    bl$gene <- paste0(bl$gene, ifelse(bl$category == "Carbapenemase", " (CP)", ""))
    keys <- c("gene", "species_short", "ward")
    rings <- do.call(rbind, lapply(seq_along(keys), function(level) {
        by <- bl[keys[seq_len(level)]]
        agg <- aggregate(list(n = rep(1, nrow(bl))), by = by, FUN = sum)
        agg <- agg[do.call(order, agg[keys[seq_len(level)]]), , drop = FALSE]
        data.frame(gene = agg$gene, label = agg[[keys[level]]], level = level,
                   xmax = cumsum(agg$n), xmin = cumsum(agg$n) - agg$n, n = agg$n, stringsAsFactors = FALSE)
    }))
    total <- sum(rings$n[rings$level == 1])
    rings$show <- rings$n / total >= 0.04
    p <- ggplot(rings) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = level, ymax = level + 1, fill = gene,
                      alpha = factor(level)), colour = "white", linewidth = 0.3) +
        geom_text(data = rings[rings$show, ], aes(x = (xmin + xmax) / 2, y = level + 0.5, label = label),
                  size = 2.6) +
        scale_alpha_manual(values = c(1, 0.7, 0.45), guide = "none") +
        coord_polar() + ylim(0, 4) +
        labs(fill = "Gene (CP = carbapenemase)",
             subtitle = "Inner ring: gene; middle: species; outer: ward (isolates with each gene)") +
        theme_void(base_size = 11, base_family = dviz_font_family) +
        theme(text = element_text(family = dviz_font_family))
    save_fig(p, "2_betalactamase_sunburst", width = 8, height = 7)
}

# ------------------------------------------------------------------------------------------------ 3. Sankey
fig_sankey <- function() {
    if (!has("ggalluvial")) return(cat("  skipped: install the ggalluvial package\n"))
    suppressPackageStartupMessages(library(ggalluvial))
    cfiles <- list.files(file.path(opt$results, "05_plasmids", "mobsuite"), "^contig_report\\.txt$",
                         recursive = TRUE, full.names = TRUE)
    if (!length(cfiles)) return(cat("  skipped: no MOB-suite contig_report.txt files\n"))
    contigs <- do.call(rbind, lapply(cfiles, function(f) {
        cr <- read_table(f)
        rel <- sub(".*/mobsuite/", "", f)
        data.frame(sample = strsplit(rel, "/")[[1]][1], contig = first_col(cr, c("contig_id")),
                   molecule = first_col(cr, c("molecule_type")), rep = first_col(cr, c("rep_type(s)", "rep_type")),
                   stringsAsFactors = FALSE)
    }))
    contigs$contig <- sub(" .*", "", contigs$contig)
    bl <- unique(amr[!is.na(amr$category), c("sample", "gene", "contig", "category")])
    if (!nrow(bl)) return(cat("  skipped: no carbapenemase or ESBL genes found\n"))
    bl$contig <- sub(" .*", "", bl$contig)
    bl <- merge(bl, contigs, by = c("sample", "contig"), all.x = TRUE)
    bl$carrier <- ifelse(is.na(bl$molecule), "not resolved",
                  ifelse(bl$molecule == "chromosome", "chromosome",
                  ifelse(is.na(bl$rep) | bl$rep %in% c("", "-"), "plasmid (no replicon)", bl$rep)))
    bl <- merge(bl, info[c("sample", "species_short")], by = "sample")
    if (all(bl$carrier == "not resolved")) {
        cat("  note: no gene contig matched a MOB-suite contig ID; check that both use the same assembly headers\n")
    }
    sk <- aggregate(list(n = rep(1, nrow(bl))), by = bl[c("carrier", "species_short", "gene")], FUN = sum)
    p <- ggplot(sk, aes(axis1 = carrier, axis2 = species_short, axis3 = gene, y = n)) +
        geom_alluvium(aes(fill = carrier), width = 1 / 5, alpha = 0.8) +
        geom_stratum(width = 1 / 5, fill = "grey95", colour = "grey40") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.8) +
        scale_x_discrete(limits = c("Replicon / location", "Species", "Gene"), expand = c(0.15, 0.05)) +
        labs(y = "Gene copies (isolate-gene pairs)", fill = "Replicon / location",
             subtitle = "Carbapenemase and ESBL genes: where they sit, and in which species") +
        theme_dviz_open(font_size = 11) +
        theme(text = element_text(family = dviz_font_family), panel.grid = element_blank())
    save_fig(p, "3_plasmid_sankey", width = 10, height = max(5, 0.25 * nrow(sk) + 3))
}

# ------------------------------------------------------------------------------------------------ 4. trees
fig_trees <- function() {
    if (!has("ggtree") || !has("ggnewscale") || !has("ape")) {
        return(cat("  skipped: install ape, ggtree (Bioconductor) and ggnewscale\n"))
    }
    suppressPackageStartupMessages(library(ggtree))
    trees <- if (!is.na(opt$tree)) strsplit(opt$tree, ",")[[1]] else
        list.files(file.path(opt$results, "06_phylogeny"), "\\.treefile$", recursive = TRUE, full.names = TRUE)
    if (!length(trees)) return(cat("  skipped: no .treefile in 06_phylogeny\n"))
    cpg <- amr[amr$category %in% "Carbapenemase", ]
    cp <- if (nrow(cpg)) {
        aggregate(list(carbapenemase = cpg$gene), by = list(sample = cpg$sample),
                  FUN = function(g) paste(sort(unique(g)), collapse = "+"))
    } else {
        data.frame(sample = character(), carbapenemase = character())
    }
    for (tf in trees) {
        tr <- ape::read.tree(tf)
        tr <- ape::drop.tip(tr, intersect("Reference", tr$tip.label))
        if (ape::Ntip(tr) < 3) { cat("  skipped", tf, ": fewer than 3 isolates\n"); next }
        if (!is.na(opt$outgroup) && opt$outgroup %in% tr$tip.label) tr <- ape::root(tr, opt$outgroup, resolve.root = TRUE)
        group <- sub("\\.treefile$", "", basename(tf))
        tips <- data.frame(label = tr$tip.label, stringsAsFactors = FALSE)
        tips <- merge(tips, info[c("sample", "ward", "year", "ST", "K_locus")], by.x = "label", by.y = "sample", all.x = TRUE)
        tips <- merge(tips, cp, by.x = "label", by.y = "sample", all.x = TRUE)
        depth <- max(ape::node.depth.edgelength(tr))
        p <- ggtree(tr, layout = "circular") %<+% tips +
            geom_tiplab(aes(colour = ward), size = 2.4, offset = 0.02 * depth) +
            geom_tippoint(aes(shape = carbapenemase), size = 2.5, na.rm = TRUE) +
            labs(colour = "Ward (tip label)", shape = "Carbapenemase")
        rings <- data.frame(ST = tips$ST, Ward = tips$ward, Year = tips$year, K_locus = tips$K_locus,
                            row.names = tips$label, stringsAsFactors = FALSE)
        if (length(unique(na.omit(rings$ST))) < 2) rings$ST <- NULL       # one ST per group: no need for a ring
        # One ring per gene with resistance-conferring point mutations in these isolates (as in visualizeAMR.R)
        pm <- amr[amr$point & amr$sample %in% tips$label, ]
        if (nrow(pm)) {
            pm$locus <- sub("_.*", "", pm$gene)
            for (locus in sort(unique(pm$locus))) {
                hits <- pm[pm$locus == locus, ]
                muts <- tapply(sub(paste0("^", locus, "_"), "", hits$gene), hits$sample,
                               function(m) paste(sort(unique(m)), collapse = ", "))
                rings[[locus]] <- unname(muts[rownames(rings)])
            }
        }
        start <- 0.35 * depth    # leave room for the tip labels
        for (i in seq_along(rings)) {
            if (i > 1) p <- p + ggnewscale::new_scale_fill()
            p <- gheatmap(p, rings[i], offset = start + (i - 1) * 0.12 * depth, width = 0.1,
                          colnames = FALSE, color = "white") +
                scale_fill_discrete(name = names(rings)[i], na.value = "grey90")
        }
        p <- p + labs(title = group, subtitle = paste("Rings, inside to outside:", paste(names(rings), collapse = ", "),
                                                    "(gene names: resistance-conferring point mutations)")) +
            theme(legend.position = "right")
        save_fig(p, paste0("4_tree_", group), width = 11, height = 9)
    }
}

# ------------------------------------------------------------------------------------------------ 5. gene heatmap
fig_heatmap <- function() {
    if (!has("ComplexHeatmap")) return(cat("  skipped: install ComplexHeatmap (Bioconductor)\n"))
    if (!nrow(amr)) return(cat("  skipped: no AMRFinderPlus hits\n"))
    feat <- unique(amr[c("sample", "gene", "point")])
    feats <- unique(feat[c("gene", "point")])
    feats <- feats[order(feats$point, feats$gene), ]
    samples <- info$sample[order(info$species_short, suppressWarnings(as.numeric(sub("\\D.*", "", sub("^ST", "", info$ST)))), info$sample)]
    mat <- matrix("absent", nrow(feats), length(samples), dimnames = list(feats$gene, samples))
    mat[cbind(match(feat$gene, rownames(mat)), match(feat$sample, colnames(mat)))] <- "present"
    anno <- data.frame(Species = info$species_short, ST = info$ST, Ward = info$ward)[match(samples, info$sample), ]
    anno[is.na(anno)] <- "unknown"
    top <- ComplexHeatmap::HeatmapAnnotation(df = anno, which = "column", gap = grid::unit(1, "mm"),
                                             simple_anno_size = grid::unit(0.5, "cm"))
    ht <- ComplexHeatmap::Heatmap(mat, name = "AMR determinant", col = c(absent = "grey92", present = "#E41A1C"),
                                  cluster_rows = FALSE, cluster_columns = FALSE,
                                  row_split = factor(ifelse(feats$point, "Point mutations", "Acquired genes"),
                                                     levels = c("Acquired genes", "Point mutations")),
                                  rect_gp = grid::gpar(col = "white", lwd = 1), top_annotation = top,
                                  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9),
                                  column_title = "AMRFinderPlus: acquired genes and point mutations per isolate")
    width <- 0.22 * ncol(mat) + 5
    height <- 0.2 * nrow(mat) + 4
    pdf(file.path(opt$outdir, "5_amr_heatmap.pdf"), width = width, height = height)
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    invisible(dev.off())
    png(file.path(opt$outdir, "5_amr_heatmap.png"), width = width, height = height, units = "in", res = 300)
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    invisible(dev.off())
    cat("  wrote", file.path(opt$outdir, "5_amr_heatmap.{pdf,png}"), "\n")
}

# ------------------------------------------------------------------------------------------------ run each figure
for (f in c("fig_resistance", "fig_sunburst", "fig_sankey", "fig_trees", "fig_heatmap")) {
    cat(f, "\n")
    tryCatch(get(f)(), error = function(e) cat("  FAILED:", conditionMessage(e), "\n"))
}
cat("Figures in", opt$outdir, "\n")
