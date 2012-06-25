denplot <- function(mcmcout, parms=NULL, regex=NULL, random=NULL, ci = NULL, xlim=NULL, ylim=NULL, auto.layout=TRUE, mar=c(2.0, 2.0, 1.5, 0.25)+0.1, col=NULL, lty=1, xlab="", ylab="", plot.title=NULL, main=NULL, greek = FALSE, collapse=FALSE, style=c("gray", "plain"), ...){
    gpar.args <- list(...)
    style <- match.arg(style)
    mcmcout <- convert.mcmc.list(mcmcout)
    if (collapse) mcmcout <- as.mcmc.list(as.mcmc(as.matrix(mcmcout)))
    nchains <- length(mcmcout)
    if (is.null(col)){
        col <- mcmcplotsPalette(nchains)
    }
    col <- rep(col, length=nchains)
    lty <- rep(lty, length=nchains)

    if (is.null(varnames(mcmcout))){
        warning("Argument 'mcmcout' did not have valid variable names, so names have been created for you.")
        varnames(mcmcout) <- varnames(mcmcout, allow.null=FALSE)
    }
    parnames <- parms2plot(varnames(mcmcout), parms, regex, random)
    if (length(parnames)==0)
        stop("No parameters matched arguments 'parms' or 'regex'.")
    p <- length(parnames)
    ## parnames <- varnames(mcmcout)
    ## if (is.null(parnames))
    ##     stop("Argument 'mcmcout' must have valid variable names, chump!")
    ## parnames <- parms2plot(parnames, parms, regex, random)
    ## p <- length(parnames)

    if (auto.layout){
        op <- mult.fig(p, main=plot.title, mar=mar)$old.par
        on.exit(par(op))
    }
    if (is.null(main)){
        main <- parnames
        names(parnames) <- main
    }
    if (greek==TRUE){
        greek.alphabet <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega")
        has.greek <- sapply(main, function(l) any(sapply(greek.alphabet, function(x) length(grep(x, l))>0)))
        main[has.greek] <- parse(text=gsub(",", "*','*", main[has.greek]))
    }
    main <- rep(main, length.out=length(parnames))
    names(main) <- parnames
    for (pmtr in parnames){
        denoverplot1(mcmcout[, pmtr, drop=FALSE], ci=ci, col=col, lty=lty, xlim=xlim, ylim=ylim, style=style, xlab=xlab, ylab=ylab, main=main[pmtr], gpar=gpar.args)
    }
    return(invisible(parnames))
}
