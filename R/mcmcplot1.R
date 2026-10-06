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
        opar <- par(mar = c(4.2, 5.4, 2.0, 1.5), oma = c(0, 0, 3.2, 0), mgp = c(2.8, 0.7, 0),
                    col.lab = "#1E293B", font.lab = 2, cex.lab = 1.15)
    } else {
        opar <- par(mar = c(4.8, 5.4, 2.0, 1.2), oma = c(0, 0, 3.0, 0), mgp = c(2.8, 0.7, 0))
    }
    on.exit(par(opar))
    layout(matrix(c(1, 2, 1, 3, 4, 4), 3, 2, byrow=TRUE))
    denoverplot1(x, col=col, lty=lty, xlim=xlim, ylim=ylim, style=style, xlab=label, ylab="Density")
    autplot1(x, style=style, col=col[1])
    rmeanplot1(x, col=col, lty=lty, style=style)
    traplot1(x, col=col, lty=lty, style=style, ylab=label, xlab="Iteration")
    title_col <- if (style == "clean") "#0F172A" else "black"
    if (greek) {
        title(parse(text=paste("paste('Diagnostics for ', ", label, ")")), cex.main=1.6, outer=TRUE, col.main=title_col, font.main=2)
    } else {
        title(paste("Diagnostics for ", parname, sep=""), cex.main=1.6, outer=TRUE, col.main=title_col, font.main=2)
    }
}
