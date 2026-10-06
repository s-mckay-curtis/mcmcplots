mcmcplot1 <- function(x, col=mcmcplotsPalette(n), lty=1, xlim=NULL, ylim=NULL, style=c("clean", "plain", "gray"), greek = FALSE){
    x <- convert.mcmc.list(x)
    style <- match.arg(style)
    n <- length(x)
    parname <- varnames(x)
    label <- parname
    if (greek) {
      label <- .to.greek(label)
    }
    if (style == "clean") {
        opar <- par(mar = c(3.8, 4.4, 1.8, 1.2), oma = c(0, 0, 2.5, 0), mgp = c(2.6, 0.6, 0),
                    col.lab = "#334155", font.lab = 1, cex.lab = 0.95)
    } else {
        opar <- par(mar=c(5, 4, 2, 1) + 0.2, oma=c(0, 0, 2, 0) + 0.1)
    }
    on.exit(par(opar))
    layout(matrix(c(1, 2, 1, 3, 4, 4), 3, 2, byrow=TRUE))
    denoverplot1(x, col=col, lty=lty, xlim=xlim, ylim=ylim, style=style, xlab=label, ylab="Density")
    autplot1(x, style=style, col=col[1])
    rmeanplot1(x, col=col, lty=lty, style=style)
    traplot1(x, col=col, lty=lty, style=style, ylab=label, xlab="Iteration")
    title_col <- if (style == "clean") "#0F172A" else "black"
    if (greek) {
        title(parse(text=paste("paste('Diagnostics for ', ", label, ")")), cex.main=1.4, outer=TRUE, col.main=title_col)
    } else {
        title(paste("Diagnostics for ", parname, sep=""), cex.main=1.4, outer=TRUE, col.main=title_col)
    }
}
