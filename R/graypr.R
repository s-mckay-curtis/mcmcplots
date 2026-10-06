.cleanpr <- function(x.axis = TRUE, y.axis = TRUE,
                     x.major = TRUE, y.major = TRUE,
                     x.minor = FALSE, y.minor = FALSE,
                     bg = "#FFFFFF",
                     grid.col = "#F1F5F9",
                     axis.col = "#CBD5E1",
                     text.col = "#475569",
                     lwd.grid = 1,
                     lwd.axis = 1,
                     tcl = -0.35,
                     cex.axis = 0.9) {
    u <- par("usr")
    rect(u[1], u[3], u[2], u[4], border = NA, col = bg)
    x.ticks <- axTicks(1)
    y.ticks <- axTicks(2)
    if (x.major) abline(v = x.ticks, col = grid.col, lwd = lwd.grid)
    if (y.major) abline(h = y.ticks, col = grid.col, lwd = lwd.grid)
    if (x.minor) {
        x.sep <- diff(x.ticks)[1] / 2
        abline(v = c(min(x.ticks) - x.sep, x.ticks + x.sep), col = grid.col, lwd = lwd.grid * 0.7, lty = 3)
    }
    if (y.minor) {
        y.sep <- diff(y.ticks)[1] / 2
        abline(h = c(min(y.ticks) - y.sep, y.ticks + y.sep), col = grid.col, lwd = lwd.grid * 0.7, lty = 3)
    }
    if (x.axis) axis(1, at = x.ticks, col = axis.col, col.ticks = axis.col,
                     col.axis = text.col, tcl = tcl, cex.axis = cex.axis,
                     mgp = c(2.2, 0.45, 0))
    if (y.axis) axis(2, at = y.ticks, col = axis.col, col.ticks = axis.col,
                     col.axis = text.col, tcl = tcl, cex.axis = cex.axis,
                     mgp = c(2.2, 0.45, 0))
    box(col = axis.col, lwd = lwd.axis)
}

.graypr <- function(x.axis=TRUE, y.axis=TRUE, x.major=TRUE, y.major=TRUE, x.minor=TRUE, y.minor=TRUE, x.malty=1, y.malty=1, x.milty=1, y.milty=1){
    if (x.axis)
        axis(1, lwd=0, lwd.ticks=1)
    if (y.axis)
        axis(2, lwd=0, lwd.ticks=1)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border=NA, col=gray(0.85))
    x.ticks <- axTicks(1)
    y.ticks <- axTicks(2)
    if (x.major){
        abline(v=x.ticks, col=gray(0.90), lty=x.malty, lwd=2)
    }
    if (y.major){
        abline(h=y.ticks, col=gray(0.90), lty=y.malty, lwd=2)
    }
    if (x.minor){
        x.sep <- diff(x.ticks)[1]/2
        x.minorgrid <- c(min(x.ticks)-x.sep, x.ticks+x.sep)
        abline(v=x.minorgrid, col=gray(0.90), lty=x.milty)
    }
    if (y.minor){
        y.sep <- diff(y.ticks)[1]/2
        y.minorgrid <- c(min(y.ticks)-y.sep, y.ticks+y.sep)
        abline(h=y.minorgrid, col=gray(0.90), lty=y.milty)
    }
}
