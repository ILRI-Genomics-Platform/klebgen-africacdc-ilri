# dviz_theme.R - the themes from Wilke, "Fundamentals of Data Visualization" (the dataviz package),
# shared by the course's figure scripts (sourced, not run).
#
#   theme_dviz_open()    axis lines, no grid          - scatter, categorical comparisons
#   theme_dviz_hgrid()   horizontal grid only         - plots whose VALUE is on y (bars, lines)
#   theme_dviz_vgrid()   vertical grid only           - plots whose VALUE is on x (horizontal bars)
#
# Two things this file handles so the figure scripts do not have to:
#
# 1. cowplot. The three themes are thin wrappers around cowplot's. cowplot is installed by
#    bin/install_r_packages.R, but when it is missing these fall back to close equivalents built on
#    ggplot2's own themes, so a figure script still runs (with a message) instead of failing.
#
# 2. The font. "Arial Narrow" is a Windows/macOS font and is usually ABSENT on a Linux cluster; an
#    unavailable family makes ggplot2 warn once per text element and then substitute anyway. The
#    family is therefore checked once here and replaced with the device default when it is not
#    installed. Override with the NEONATAL_FONT environment variable or dviz_use_font("Family").
#    Good Linux substitutes, if you want one: "Nimbus Sans Narrow", "DejaVu Sans Condensed",
#    "Liberation Sans Narrow".

dviz_font_family                <- "Arial Narrow"
dviz_font_family_bold           <- "Arial Narrow"
dviz_font_family_condensed      <- "Arial Narrow"
dviz_font_family_bold_condensed <- "Arial Narrow"

# TRUE when the family is installed, or when we have no way to check (then ggplot2's own warning
# is the signal). "" is ggplot2's device default and always counts as available.
dviz_font_available <- function(family) {
    if (!nzchar(family)) return(TRUE)
    if (requireNamespace("systemfonts", quietly = TRUE)) {
        f <- try(systemfonts::system_fonts(), silent = TRUE)
        if (!inherits(f, "try-error") && "family" %in% names(f)) {
            return(any(tolower(f$family) == tolower(family)))
        }
    }
    TRUE
}

# Set the family used by every theme below. Returns the family actually in force.
dviz_use_font <- function(family = "Arial Narrow", quiet = FALSE) {
    if (nzchar(family) && !dviz_font_available(family)) {
        if (!quiet) {
            message("font '", family, "' is not installed; using the device default instead ",
                    "(set NEONATAL_FONT, or pass --font, to choose another)")
        }
        family <- ""
    }
    dviz_font_family                <<- family
    dviz_font_family_bold           <<- family
    dviz_font_family_condensed      <<- family
    dviz_font_family_bold_condensed <<- family
    invisible(family)
}
dviz_use_font(Sys.getenv("NEONATAL_FONT", "Arial Narrow"), quiet = TRUE)

.dviz_have_cowplot <- requireNamespace("cowplot", quietly = TRUE)
if (!.dviz_have_cowplot) {
    message("cowplot is not installed: theme_dviz_* fall back to ggplot2 equivalents ",
            "(bin/run_R.sh install_r_packages.R installs it)")
}

theme_dviz_open <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                            rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14) {
    half_line <- font_size / 2

    base <- if (.dviz_have_cowplot) {
        cowplot::theme_half_open(font_size = font_size, font_family = font_family,
                                 line_size = line_size, rel_small = rel_small,
                                 rel_tiny = rel_tiny, rel_large = rel_large)
    } else {
        ggplot2::theme_classic(base_size = font_size, base_family = font_family,
                               base_line_size = line_size) %+replace%
            ggplot2::theme(
                axis.text  = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                legend.text = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                plot.title = ggplot2::element_text(size = ggplot2::rel(rel_large), hjust = 0))
    }
    base %+replace%
        ggplot2::theme(
            text = ggplot2::element_text(family = dviz_font_family),
            plot.margin = ggplot2::margin(half_line / 2, 1.5, half_line / 2, 1.5),
            complete = TRUE
        )
}

theme_dviz_hgrid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
                             colour = "grey90") {
    half_line <- font_size / 2

    base <- if (.dviz_have_cowplot) {
        cowplot::theme_minimal_hgrid(font_size = font_size, font_family = font_family,
                                     line_size = line_size, rel_small = rel_small,
                                     rel_tiny = rel_tiny, rel_large = rel_large, colour = colour)
    } else {
        ggplot2::theme_minimal(base_size = font_size, base_family = font_family) %+replace%
            ggplot2::theme(
                panel.grid         = ggplot2::element_blank(),
                panel.grid.major.y = ggplot2::element_line(colour = colour, linewidth = line_size),
                axis.line.x        = ggplot2::element_line(colour = "black", linewidth = line_size),
                axis.ticks.x       = ggplot2::element_line(colour = "black", linewidth = line_size),
                axis.text          = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                legend.text        = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                plot.title         = ggplot2::element_text(size = ggplot2::rel(rel_large), hjust = 0))
    }
    base %+replace%
        ggplot2::theme(
            text = ggplot2::element_text(family = dviz_font_family),
            plot.margin = ggplot2::margin(half_line / 2, 1.5, half_line / 2, 1.5),
            complete = TRUE
        )
}

theme_dviz_vgrid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
                             colour = "grey90") {
    half_line <- font_size / 2

    base <- if (.dviz_have_cowplot) {
        cowplot::theme_minimal_vgrid(font_size = font_size, font_family = font_family,
                                     line_size = line_size, rel_small = rel_small,
                                     rel_tiny = rel_tiny, rel_large = rel_large, colour = colour)
    } else {
        ggplot2::theme_minimal(base_size = font_size, base_family = font_family) %+replace%
            ggplot2::theme(
                panel.grid         = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_line(colour = colour, linewidth = line_size),
                axis.line.y        = ggplot2::element_line(colour = "black", linewidth = line_size),
                axis.ticks.y       = ggplot2::element_line(colour = "black", linewidth = line_size),
                axis.text          = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                legend.text        = ggplot2::element_text(size = ggplot2::rel(rel_small)),
                plot.title         = ggplot2::element_text(size = ggplot2::rel(rel_large), hjust = 0))
    }
    base %+replace%
        ggplot2::theme(
            text = ggplot2::element_text(family = dviz_font_family),
            plot.margin = ggplot2::margin(half_line / 2, 1.5, half_line / 2, 1.5),
            complete = TRUE
        )
}
