#!/usr/bin/env Rscript
#
# install_r_packages.R - install every R package used by the course's R scripts, once.
#
#   Every user, once, into their own library (~/R/<platform>-library/<version>):
#       bin/run_R.sh install_r_packages.R
#   Another folder:  bin/run_R.sh install_r_packages.R --lib /path/to/folder   (then set NEONATAL_R_LIB to it)
#
# Needs internet access (CRAN, Bioconductor, GitHub): run it on a login node, not a compute node.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "r_setup.R"))

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) == 2 && args[1] == "--lib") args[2] else user_lib
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(lib)            # only this library and R's built-in packages

packages <- list(
    cran   = c("ggplot2", "ape", "coda", "ggnewscale", "ggalluvial"),
    bioc   = c("ggtree", "ComplexHeatmap"),
    github = c(BactDating = "xavierdidelot/BactDating")
)
cat("Installing into", lib, "\n")
installed_version <- function(p) {
    v <- tryCatch(as.character(packageVersion(p, lib.loc = c(lib, .Library))), error = function(e) NA)
    if (is.na(v)) "MISSING" else v
}

# 1. CRAN packages
install_missing(cran = packages$cran, lib = lib)

# 2. On R 4.4, Bioconductor 3.20's ggtree was released before ggplot2 4.0 and does not load with it
#    ("object 'check_linewidth' not found"): keep ggplot2 at 3.5.2, the last 3.x release. Newer ggtree (R 4.5) is fine.
if (getRversion() < "4.5.0") {
    options(neonatal.pinned = "ggplot2")
    if (installed_version("ggplot2") != "MISSING" && package_version(installed_version("ggplot2")) >= "4.0.0") {
        cat("R", as.character(getRversion()), ": installing ggplot2 3.5.2 (ggtree for R 4.4 does not work with ggplot2 4.x)\n")
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib, repos = cran_repo)
        remotes::install_version("ggplot2", version = "3.5.2", lib = lib, repos = cran_repo, upgrade = "never")
    }
}

# 3. Bioconductor and GitHub packages
install_missing(bioc = packages$bioc, github = packages$github, lib = lib)

all_pkgs <- c(packages$cran, packages$bioc, names(packages$github))
versions <- vapply(all_pkgs, installed_version, character(1))
for (p in all_pkgs) cat(sprintf("  %-15s %s\n", p, versions[p]))
ok <- versions != "MISSING"
if (!all(ok)) quit(status = 1)
cat("All packages installed.\n")
