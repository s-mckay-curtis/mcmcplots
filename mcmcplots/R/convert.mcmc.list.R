convert.mcmc.list <- function(x){
    if (!is.mcmc.list(x)){
        if (!is.mcmc(x)){
            x <- lapply(x, as.mcmc)
        }
        x <- mcmc.list(x)
    }
    return(x)
}
