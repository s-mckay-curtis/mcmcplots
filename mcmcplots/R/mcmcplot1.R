mcmcplot1 <- function(x, col=mcmcplotsPalette(n), lty=1, xlim=NULL, ylim=NULL, style=c("gray", "plain"), greek = FALSE){
    x <- convert.mcmc.list(x)
    style <- match.arg(style)
    n <- length(x)
    parname <- varnames(x)
    label <- parname
    if (greek==TRUE){
        greek.alphabet <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega")
        has.greek <- any(sapply(greek.alphabet, function(x) length(grep(x, label))>0))
        if (has.greek){
            label <- parse(text=gsub(",", "*','*", label))
        }
    }
    ## layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE))
    opar <- par(mar=c(5, 4, 2, 1) + 0.2, oma=c(0, 0, 2, 0) + 0.1)
    on.exit(par(opar))
    layout(matrix(c(1, 2, 1, 3, 4, 4), 3, 2, byrow=TRUE))
    denoverplot1(x, col=col, lty=lty, xlim=xlim, ylim=ylim, style=style, xlab=label, ylab="Density")
    autplot1(x, style=style)
    autplot1(x, style=style, partial=TRUE)
    traplot1(x, col=col, lty=lty, style=style, ylab=label, xlab="Iteration")
    if (greek) {
        title(parse(text=paste("paste('Diagnostics for ', ", parname, ")")), cex.main=1.5, outer=TRUE)
    } else {
        title(paste("Diagnostic Plots for ", parname, sep=""), cex.main=1.5, outer=TRUE)
    }
}
