# r_setup.R - library paths and package loading shared by the course's R scripts (sourced, not run).
#
# Every user installs the packages into their own library and R uses only that library plus R's built-in
# packages, so packages from other libraries (e.g. an older shared one) cannot get mixed in:
#     ~/R/<platform>-library/<major.minor>, e.g. ~/R/x86_64-pc-linux-gnu-library/4.4
#     (NEONATAL_R_LIB=/other/folder to use another folder)
# A missing package is installed there, unless NEONATAL_R_NO_INSTALL=1 is set (compute nodes without internet).

r_version <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
user_lib  <- Sys.getenv("NEONATAL_R_LIB",
                        file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), r_version))

# R's own 'methods' package must load, or every other package fails and looks "missing"
if (!suppressWarnings(requireNamespace("methods", quietly = TRUE))) {
    stop("R cannot load its built-in 'methods' package on ", Sys.info()[["nodename"]], " (R ", R.home(), ").\n",
         "  1. A stray 'methods' folder in your library? ls -ld ", file.path(user_lib, "methods"),
         "  -> move it away\n",
         "  2. Otherwise R/", r_version, " is broken on this node: run on another node, and tell the HPC administrator.",
         call. = FALSE)
}

if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(user_lib)       # = your library, then R's built-in library (.Library); R_LIBS and others are dropped

cran_repo <- "https://cloud.r-project.org"

# Dependencies (Depends/Imports/LinkingTo, recursively) of `pkgs` that are missing, or older than in the
# repository, in the copy R would load (the first one on .libPaths()). install.packages() leaves an installed
# dependency alone even when it is too old, so a new package installed into your library can still load an old
# dependency from the shared course library and fail with "lazy loading failed"; these are installed too.
outdated_deps <- function(pkgs, db) {
    deps <- unique(unlist(tools::package_dependencies(pkgs, db = db, which = c("Depends", "Imports", "LinkingTo"),
                                                      recursive = TRUE)))
    # pinned packages (e.g. ggplot2 3.5.2 on R 4.4, see install_r_packages.R) are never "upgraded" here
    deps <- setdiff(deps, c(rownames(installed.packages(priority = "base")), pkgs, getOption("neonatal.pinned")))
    deps <- deps[deps %in% rownames(db)]
    if (!length(deps)) return(character())
    inst <- installed.packages()                          # rows follow .libPaths(): first row = copy R loads
    inst <- inst[!duplicated(inst[, "Package"]), , drop = FALSE]
    have <- inst[match(deps, inst[, "Package"]), "Version"]
    # compareVersion() < 0: the copy R would load is older than the repository's
    old <- is.na(have) | mapply(function(h, v) !is.na(h) && utils::compareVersion(h, v) < 0, have, db[deps, "Version"])
    deps[old]
}

# Install missing packages into `lib`. cran / bioc: package names; github: c(Package = "owner/repo")
install_missing <- function(cran = character(), bioc = character(), github = character(), lib = user_lib) {
    missing <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
    if (file.access(lib, 2) != 0) stop("Cannot write to the R library ", lib)
    # Leftover lock folders from interrupted installs block new ones
    locks <- Sys.glob(file.path(lib, "00LOCK*"))
    if (length(locks)) { cat("Removing stale lock folder(s):", basename(locks), "\n"); unlink(locks, recursive = TRUE) }
    # Install into `lib` and load from it first, so freshly installed dependencies are the ones used
    old_paths <- .libPaths()
    .libPaths(unique(c(lib, old_paths)))
    on.exit(.libPaths(old_paths), add = TRUE)

    todo <- missing(cran)
    if (length(todo)) {
        db <- available.packages(repos = cran_repo)
        deps <- outdated_deps(todo, db)
        if (length(deps)) cat("Also installing newer dependencies into", lib, ":", paste(deps, collapse = ", "), "\n")
        install.packages(c(deps, todo), lib = lib, repos = cran_repo)
    }
    todo <- missing(bioc)
    if (length(todo)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib, repos = cran_repo)
        db <- available.packages(repos = BiocManager::repositories())
        deps <- outdated_deps(todo, db)
        if (length(deps)) cat("Also installing newer dependencies into", lib, ":", paste(deps, collapse = ", "), "\n")
        BiocManager::install(c(deps, todo), lib = lib, update = FALSE, ask = FALSE)
    }
    todo <- github[missing(names(github))]
    if (length(todo)) {
        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib, repos = cran_repo)
        for (repo in todo) remotes::install_github(repo, lib = lib, upgrade = "never")
    }
    invisible(TRUE)
}

# Make packages available: install what is missing (unless NEONATAL_R_NO_INSTALL is set), then report.
# required = TRUE stops with instructions if a package is still missing; otherwise returns TRUE/FALSE per package.
use_packages <- function(cran = character(), bioc = character(), github = character(), required = TRUE) {
    all_pkgs <- c(cran, bioc, names(github))
    have <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
    if (!all(have) && Sys.getenv("NEONATAL_R_NO_INSTALL") == "") {
        cat("Installing missing R package(s) into", user_lib, ":", paste(all_pkgs[!have], collapse = ", "), "\n")
        tryCatch(install_missing(cran[cran %in% all_pkgs[!have]], bioc[bioc %in% all_pkgs[!have]],
                                 github[names(github) %in% all_pkgs[!have]]),
                 error = function(e) cat("  installation failed:", conditionMessage(e), "\n"))
        have <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
    }
    if (required && !all(have)) {
        stop("Missing R package(s): ", paste(all_pkgs[!have], collapse = ", "), "\n",
             "  Install them once on a node with internet access:  bin/run_R.sh install_r_packages.R\n",
             "  Libraries searched: ", paste(.libPaths(), collapse = ", "), call. = FALSE)
    }
    invisible(have)
}
