mcmcplotsPalette <- function(n, type = c("colorblind", "viridis", "rainbow", "sequential", "grayscale"), seq = NULL) {
    type <- match.arg(type)
    if (type == "colorblind") {
        # Chromatic Okabe-Ito sequence (accessible and high-contrast)
        pal <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#E6AB02", "#334155")
        if (n == 1) return(pal[1])
        return(rep_len(pal, n))
    }
    if (type == "viridis") {
        return(hcl.colors(n, palette = "Viridis"))
    }
    if (type == "rainbow") {
        if (n == 1)
            return(rainbow_hcl(1, start = 240, l = 50, c = 100))
        return(rainbow_hcl(n, start = 0, end = 240, c = 100))
    }
    if (type == "sequential") {
        return(sequential_hcl(n))
    }
    if (type == "grayscale") {
        return(gray((1:n / (n + 1))))
    }
}
