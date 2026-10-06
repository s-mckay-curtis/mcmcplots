.array_to_mcmc_list <- function(arr) {
    d <- dim(arr)
    dn <- dimnames(arr)
    n_iter <- d[1]
    n_chains <- d[2]
    n_vars <- d[3]
    var_names <- if (length(dn) >= 3 && !is.null(dn[[3]])) dn[[3]] else paste0("var", seq_len(n_vars))

    chain_list <- lapply(seq_len(n_chains), function(i) {
        mat <- matrix(arr[, i, ], nrow = n_iter, ncol = n_vars)
        colnames(mat) <- var_names
        coda::mcmc(mat)
    })
    coda::as.mcmc.list(chain_list)
}

convert.mcmc.list <- function(x) {
    if (is.mcmc.list(x)) {
        return(x)
    }

    ## CmdStanMCMC objects from cmdstanr
    if (inherits(x, "CmdStanMCMC")) {
        x <- x$draws()
    }

    ## posterior::draws objects (draws_array, draws_df, draws_matrix, draws_list)
    if (inherits(x, "draws")) {
        if (requireNamespace("posterior", quietly = TRUE)) {
            x <- posterior::as_draws_array(x)
        }
    }

    ## rstan::stanfit objects
    if (inherits(x, "stanfit")) {
        x <- as.array(x)
    }

    ## brms::brmsfit objects
    if (inherits(x, "brmsfit")) {
        if (requireNamespace("posterior", quietly = TRUE)) {
            x <- posterior::as_draws_array(x)
        } else if (!is.null(x$fit)) {
            return(convert.mcmc.list(x$fit))
        }
    }

    ## 3D arrays: iterations x chains x variables
    if (is.array(x) && length(dim(x)) == 3) {
        return(.array_to_mcmc_list(x))
    }

    ## Standard coda / list / bugs / rjags conversions
    if (!is.mcmc(x)) {
        if ("list" %in% class(x)) {
            x <- lapply(x, as.mcmc)
        } else {
            x <- as.mcmc(x)
        }
    }
    if (!is.mcmc.list(x)) {
        x <- coda::mcmc.list(x)
    }
    return(x)
}
