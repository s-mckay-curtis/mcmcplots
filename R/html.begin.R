.html.begin <- function(outdir = tempdir(), filename = "index", extension = "html", title, cssfile = NULL, embed.css = TRUE) {
    doctype <- '<!DOCTYPE html>\n<html lang="en">\n'
    css_content <- NULL
    csslink <- NULL

    if (!is.null(cssfile) && length(cssfile) > 0 && cssfile != "") {
        # Check if local file path
        raw_css_path <- gsub("^file://(/+)?", "", cssfile)
        if (embed.css && file.exists(raw_css_path)) {
            css_text <- paste(readLines(raw_css_path, warn = FALSE), collapse = "\n")
            css_content <- paste0("<style>\n", css_text, "\n</style>")
        } else {
            csslink <- paste0('<link rel="stylesheet" type="text/css" href="', cssfile, '">')
        }
    }

    file <- file.path(outdir, paste(filename, extension, sep = "."))
    out <- c(
        doctype,
        '<head>',
        '  <meta charset="utf-8">',
        '  <meta name="viewport" content="width=device-width, initial-scale=1.0">',
        paste0('  <title>', title, '</title>'),
        if (!is.null(csslink)) paste0('  ', csslink) else NULL,
        if (!is.null(css_content)) css_content else NULL,
        '</head>\n<body>'
    )
    out <- paste(out[!sapply(out, is.null)], collapse = "\n")
    cat(out, "\n", file = file, append = FALSE)
    invisible(file)
}
