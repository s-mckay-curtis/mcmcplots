.img.to.base64 <- function(png_path) {
    raw_bytes <- readBin(png_path, "raw", file.info(png_path)$size)
    b64_chars <- c(LETTERS, letters, 0:9, "+", "/")
    n <- length(raw_bytes)
    pad <- (3 - (n %% 3)) %% 3
    if (pad > 0) raw_bytes <- c(raw_bytes, as.raw(rep(0, pad)))
    m <- matrix(as.integer(raw_bytes), nrow = 3)
    idx1 <- bitwShiftR(m[1, ], 2) + 1
    idx2 <- bitwOr(bitwShiftL(bitwAnd(m[1, ], 3), 4), bitwShiftR(m[2, ], 4)) + 1
    idx3 <- bitwOr(bitwShiftL(bitwAnd(m[2, ], 15), 2), bitwShiftR(m[3, ], 6)) + 1
    idx4 <- bitwAnd(m[3, ], 63) + 1
    res <- c(rbind(b64_chars[idx1], b64_chars[idx2], b64_chars[idx3], b64_chars[idx4]))
    if (pad >= 1) res[length(res)] <- "="
    if (pad == 2) res[length(res) - 1] <- "="
    paste0("data:image/png;base64,", paste0(res, collapse = ""))
}

mcmcplot <- function(mcmcout, parms = NULL, regex = NULL, random = NULL,
                     leaf.marker = "[\\[_]", dir = tempdir(), filename = "MCMCoutput",
                     extension = "html", title = NULL, heading = title,
                     col = NULL, lty = 1, xlim = NULL, ylim = NULL,
                     style = c("clean", "plain", "gray"), greek = FALSE,
                     browse = TRUE, retina = TRUE, res = NULL,
                     embed.img = FALSE) {
    ## This must come before mcmcout is evaluated in any other expression
    if (is.null(title))
        title <- paste("MCMC Plots: ", deparse(substitute(mcmcout)), sep = "")
    if (is.null(heading))
        heading <- title

    style <- match.arg(style)

    ## Scale and resolution for retina/high-DPI
    scale <- if (isTRUE(retina)) 2 else 1
    if (is.null(res)) {
        res.plot <- if (isTRUE(retina)) 150 else 96
    } else {
        res.plot <- res
    }

    ## Turn off graphics device if interrupted in the middle of plotting
    current.devices <- dev.list()
    on.exit(sapply(dev.list(), function(dev) if (!(dev %in% current.devices)) dev.off(dev)))

    ## Convert input mcmcout to mcmc.list object
    mcmcout <- convert.mcmc.list(mcmcout)
    nchains <- length(mcmcout)
    if (is.null(col)) {
        col <- mcmcplotsPalette(nchains)
    }
    css.file <- system.file("MCMCoutput.css", package = "mcmcplots")
    if (!dir.exists(dir)) {
        stop("Directory in argument 'dir' must exist.")
    }
    htmlfile <- .html.begin(dir, filename, extension, title = title, cssfile = css.file, embed.css = TRUE)

    ## Select parameters for plotting
    if (is.null(varnames(mcmcout))) {
        warning("Argument 'mcmcout' did not have valid variable names, so names have been created for you.")
        varnames(mcmcout) <- varnames(mcmcout, allow.null = FALSE)
    }
    parnames <- parms2plot(varnames(mcmcout), parms, regex, random, leaf.marker, do.unlist = FALSE)
    if (length(parnames) == 0)
        stop("No parameters matched arguments 'parms' or 'regex'.")
    np <- length(unlist(parnames))

    cat('\n<header class="mcmc-header"><h1>', heading, '</h1></header>\n', sep = "", file = htmlfile, append = TRUE)
    cat('<div id="outer">\n', file = htmlfile, append = TRUE)
    cat('<div id="toc">\n', file = htmlfile, append = TRUE)
    cat('<h2>Parameters</h2>\n', file = htmlfile, append = TRUE)
    cat('<input type="text" id="param_search" placeholder="Filter parameters..." onkeyup="filterPlots()" onsearch="filterPlots()">\n', file = htmlfile, append = TRUE)
    cat('<ul id="toc_items">\n', file = htmlfile, append = TRUE)
    for (group.name in names(parnames)) {
        cat(sprintf('<li class="toc_item" data-group="%s"><a href="#group-%s">%s</a></li>\n', group.name, group.name, group.name), file = htmlfile, append = TRUE)
    }
    cat('</ul>\n</div>\n', file = htmlfile, append = TRUE)

    cat('<div class="main">\n', file = htmlfile, append = TRUE)
    htmlwidth <- 800
    htmlheight <- 600
    for (group.name in names(parnames)) {
        cat(sprintf('<section class="group-section" id="group-%s">\n', group.name), file = htmlfile, append = TRUE)
        cat(sprintf('<h2>Plots for %s</h2>\n', group.name), file = htmlfile, append = TRUE)
        for (p in parnames[[group.name]]) {
            pctdone <- round(100 * match(p, unlist(parnames)) / np)
            cat("\r", rep(" ", getOption("width")), sep = "")
            cat("\rPreparing plots for ", group.name, ".  ", pctdone, "% complete.", sep = "")
            gname <- paste(p, ".png", sep = "")
            gpath <- file.path(dir, gname)
            if (capabilities("cairo")) {
                png(gpath, width = htmlwidth * scale, height = htmlheight * scale,
                    res = res.plot, type = "cairo", antialias = "subpixel")
            } else {
                png(gpath, width = htmlwidth * scale, height = htmlheight * scale,
                    res = res.plot)
            }
            plot_err <- tryCatch({
                mcmcplot1(mcmcout[, p, drop = FALSE], col = col, lty = lty, xlim = xlim, ylim = ylim, style = style, greek = greek)
            }, error = function(e) {e})
            dev.off()
            if (inherits(plot_err, "error")) {
                cat(sprintf('<div class="plot-card" data-param="%s"><h3>%s</h3><p class="plot_err">%s. %s</p></div>\n', p, p, p, plot_err),
                    file = htmlfile, append = TRUE)
            } else {
                img_src <- if (isTRUE(embed.img)) .img.to.base64(gpath) else gname
                cat(sprintf('<div class="plot-card" data-param="%s"><h3>%s</h3>\n', p, p), file = htmlfile, append = TRUE)
                .html.img(file = htmlfile, class = "mcmcplot", src = img_src,
                          width = htmlwidth, height = htmlheight, alt = p)
                cat('</div>\n', file = htmlfile, append = TRUE)
            }
        }
        cat('</section>\n', file = htmlfile, append = TRUE)
    }
    cat("\r", rep(" ", getOption("width")), "\r", sep = "")
    cat('\n</div>\n</div>\n', file = htmlfile, append = TRUE)
    .html.end(htmlfile)
    full.name.path <- paste("file://", normalizePath(htmlfile), sep = "")
    if (browse) browseURL(full.name.path)
    invisible(full.name.path)
}
