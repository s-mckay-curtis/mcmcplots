denoverplot1 <- function(..., ci = NULL, col=NULL, lty=1, xlim=NULL, ylim=NULL, xlab = "", ylab = "Density", main = NULL, style=c("clean", "plain", "gray"), gpar=NULL){
    dat <- list(...)
    style <- match.arg(style)
    if (length(dat)==1 && is.list(dat[[1]])) dat <- dat[[1]]
    n <- length(dat)
    if (n==1 && is.list(dat[[1]])) {
        dat <- dat[[1]]
        n <- length(dat)
    }
    if (is.null(col)){
        col <- mcmcplotsPalette(n)
    }
    denout <- lapply(dat, density, bw="SJ")
    xx <- sapply(denout, function(den) den$x)
    yy <- sapply(denout, function(den) den$y)
    if (style=="plain")
        do.call("matplot", c(list(x=xx, y=yy, col=col, lty=lty, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, type="l"), gpar))
    if (style=="clean"){
        do.call("matplot", c(list(x=xx, y=yy, type="n", bty="n", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)))
        .cleanpr()
        clean_gpar <- gpar
        if (is.null(clean_gpar$lwd)) clean_gpar$lwd <- 1.8
        do.call("matlines", c(list(x=xx, y=yy, col=col, lty=lty), clean_gpar))
    }
    if (style=="gray"){
        do.call("matplot", c(list(x=xx, y=yy, type="n", bty="n", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)))
        .graypr()
        do.call("matlines", c(list(x=xx, y=yy, col=col, lty=lty), gpar))
    }
    if (!is.null(ci)){
        lb <- sapply(dat, quantile, (1 - ci)/2)
        ub <- sapply(dat, quantile, ci + (1 - ci)/2)
        ci_lty <- if (style == "clean") 2 else lty
        ci_lwd <- if (style == "clean") 1.2 else 1
        do.call("abline", c(list(v=lb, col=col, lty=ci_lty, lwd=ci_lwd), gpar))
        do.call("abline", c(list(v=ub, col=col, lty=ci_lty, lwd=ci_lwd), gpar))
    }
}
