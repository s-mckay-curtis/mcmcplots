traplot1 <- function(x, col=NULL, lty=1, style=c("clean", "plain", "gray"), ...){
    style <- match.arg(style)
    nchains <- nchain(x)
    if (is.null(col)){
        col <- mcmcplotsPalette(nchains)
    }
    xx <- as.vector(time(x))
    yy <- do.call("cbind", as.list(x))
    if (style=="plain"){
        matplot(xx, yy, type="l", col=col, lty=lty, ...)
    }
    if (style=="clean"){
        matplot(xx, yy, type="n", xaxt="n", yaxt="n", bty="n", ...)
        .cleanpr()
        matlines(xx, yy, col=col, lty=lty, lwd=1.2)
    }
    if (style=="gray"){
        matplot(xx, yy, type="n",  xaxt="n", yaxt="n", bty="n", ...)
        .graypr()
        matlines(xx, yy, col=col, lty=lty)
    }
}
