pkgname <- "mcmcplots"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('mcmcplots')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("as.mcmc.bugs")
### * as.mcmc.bugs

flush(stderr()); flush(stdout())

### Name: as.mcmc.bugs
### Title: Convert a bugs Object to an mcmc or mcmc.list Object
### Aliases: as.mcmc.bugs
### Keywords: manip

### ** Examples

## Not run: 
##D ## Data object "schools.sim" generated from the examples
##D ## in the bugs function of the R2WinBUGS package.
##D outmcmc <- as.mcmc(schools.sim)
##D 
##D ## Gelman Rubin diagnostics
##D coda:::gelman.diag(outmcmc)
##D coda:::mcmc.plot(outmcmc)
## End(Not run)



cleanEx()
nameEx("as.mcmc.rjags")
### * as.mcmc.rjags

flush(stderr()); flush(stdout())

### Name: as.mcmc.rjags
### Title: Convert an rjags Object to an mcmc or mcmc.list Object.
### Aliases: as.mcmc.rjags
### Keywords: manip

### ** Examples

## None ##



cleanEx()
nameEx("autplot1")
### * autplot1

flush(stderr()); flush(stdout())

### Name: autplot1
### Title: Autocorrelation Plot of MCMC Output
### Aliases: autplot1
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

autplot1(fakemcmc[, "alpha[1]", drop=FALSE])
autplot1(fakemcmc[, "alpha[1]", drop=FALSE], chain=2, style="plain")
autplot1(fakemcmc[, "alpha[1]", drop=FALSE], partial=TRUE)




cleanEx()
nameEx("caterplot")
### * caterplot

flush(stderr()); flush(stdout())

### Name: caterplot
### Title: Caterpillar Plots of MCMC Output
### Aliases: caterplot
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[1,", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## caterplot plots of the fake MCMC output
par(mfrow=c(2,2))
caterplot(fakemcmc, "alpha", collapse=FALSE)
caterplot(fakemcmc, "gamma", collapse=FALSE)
caterplot(fakemcmc, "alpha", labels.loc="axis", greek=TRUE, col="blue")
caterplot(fakemcmc, "gamma", labels.loc="above", greek=TRUE, col="red")

caterplot(fakemcmc, "alpha", collapse=FALSE, denstrip=TRUE)
caterplot(fakemcmc, "gamma", collapse=FALSE, denstrip=TRUE)
caterplot(fakemcmc, "alpha", labels.loc="axis", col="blue", denstrip=TRUE)
caterplot(fakemcmc, "gamma", labels.loc="above", col="red", denstrip=TRUE)

caterplot(fakemcmc, "alpha", collapse=FALSE, style="plain")
caterplot(fakemcmc, "gamma", collapse=FALSE, style="plain")
caterplot(fakemcmc, "alpha", labels.loc="axis")
caterplot(fakemcmc, "gamma", labels.loc="above")

caterplot(fakemcmc, "alpha", horizontal=FALSE)
caterplot(fakemcmc, horizontal=FALSE)
caterpoints(rnorm(10, 21, 2), horizontal=FALSE, pch="x", col="red")
caterplot(fakemcmc, horizontal=FALSE, denstrip=TRUE, col="blue", pch=NA)
caterplot(fakemcmc, horizontal=FALSE, col="red", pch=19, add=TRUE)
caterplot(fakemcmc, denstrip=TRUE, col="blue", pch=NA)
caterplot(fakemcmc, col="purple", pch=19, add=TRUE)

## What happens with NULL varnames?
varnames(fakemcmc) <- NULL
caterplot(fakemcmc)
caterplot(fakemcmc, collapse=FALSE)

## Not run: 
##D ## caterplot works on bugs objects too:
##D library(R2WinBUGS)
##D example("openbugs", "R2WinBUGS")
##D ## from the help file for openbugs:
##D schools.sim <- bugs(data, inits, parameters, model.file,
##D                     n.chains = 3, n.iter = 5000,
##D                     program = "openbugs", working.directory = NULL)
##D caterplot(schools.sim, "theta")
## End(Not run)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("caterpoints")
### * caterpoints

flush(stderr()); flush(stdout())

### Name: caterpoints
### Title: Points on a "caterplot"
### Aliases: caterpoints
### Keywords: aplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## caterplot plots of the fake MCMC output
par(mfrow=c(2,2))
caterplot(fakemcmc, "alpha", collapse=FALSE)
caterpoints(runif(5, -0.1, 0.1), pch="x", col="red")
caterplot(fakemcmc, "alpha", horizontal=FALSE)
caterpoints(runif(5, -0.1, 0.1), horizontal=FALSE, pch="x", col="red")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("corplot")
### * corplot

flush(stderr()); flush(stdout())

### Name: corplot
### Title: Plot a Correlation Matrix
### Aliases: corplot
### Keywords: hplot

### ** Examples

Rho <- matrix(c(
 1.00,  0.35, -0.65, -0.66,  0.46,  0.42,
 0.35,  1.00, -0.69, -0.64,  0.40, -0.06,
-0.65, -0.69,  1.00,  0.70, -0.57, -0.11,
-0.66, -0.64,  0.70,  1.00, -0.15, -0.10,
 0.46,  0.40, -0.57, -0.15,  1.00,  0.18,
 0.42, -0.06, -0.11, -0.10,  0.18,  1.00), 6, 6)
dimnames(Rho) <- list(paste("rho[", 1:6, "]", sep=""), paste("rho[", 1:6, "]", sep=""))
corplot(Rho)
corplot(Rho, greek=TRUE)



cleanEx()
nameEx("denoverplot")
### * denoverplot

flush(stderr()); flush(stdout())

### Name: denoverplot
### Title: Overlaying Densities for Parameters from two MCMC Simulations.
### Aliases: denoverplot
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))
fakemcmc2 <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## Plot the fake MCMC output
denoverplot(fakemcmc, fakemcmc2)
denoverplot(fakemcmc, fakemcmc2, style="plain")
denoverplot(fakemcmc, fakemcmc2, plot.title="Comparison of densities of fake data")
denoverplot(fakemcmc, fakemcmc2, plot.title="Comparison of densities of fake data", greek=TRUE)




cleanEx()
nameEx("denoverplot1")
### * denoverplot1

flush(stderr()); flush(stdout())

### Name: denoverplot1
### Title: Plot Overlaying Densities
### Aliases: denoverplot1
### Keywords: hplot

### ** Examples

denoverplot1(rnorm(1000), rnorm(1000))
denoverplot1(rnorm(1000, 0.0, 1.0), rnorm(1000, 0.1, 1.0), style="plain")
denoverplot1(list(rgamma(1000, 1, 1), rgamma(1000, 1, 1)))



cleanEx()
nameEx("denplot")
### * denplot

flush(stderr()); flush(stdout())

### Name: denplot
### Title: Density Plots for MCMC Parameters on a Single Plot
### Aliases: denplot
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## Plot densities of the fake MCMC output
denplot(fakemcmc)
denplot(fakemcmc, style="plain")
denplot(fakemcmc, collapse=TRUE, greek=TRUE)
denplot(fakemcmc, xlim=range(unlist(fakemcmc)), plot.title="Density plots of fake data.  Yawn.")
denplot(fakemcmc, "gamma")
denplot(fakemcmc, "gamma", "alpha\\[[12]]$") # all gamma and alpha[1] and alpha[2]

## What happens with NULL varnames?
varnames(fakemcmc) <- NULL
denplot(fakemcmc)

## Not run: 
##D ## denplot works on bugs objects too
##D library(R2WinBUGS)
##D example("openbugs", "R2WinBUGS")
##D ## from the help file for openbugs:
##D schools.sim <- bugs(data, inits, parameters, model.file,
##D                     n.chains = 3, n.iter = 5000,
##D                     program = "openbugs", working.directory = NULL)
##D denplot(schools.sim, "theta")
## End(Not run)



cleanEx()
nameEx("graypr")
### * graypr

flush(stderr()); flush(stdout())

### Name: .graypr
### Title: Create a Gray Plotting Region
### Aliases: .graypr
### Keywords: aplot

### ** Examples

## None.



cleanEx()
nameEx("mcmcplot")
### * mcmcplot

flush(stderr()); flush(stdout())

### Name: mcmcplot
### Title: Diagnostics Plots for MCMC in HTML format
### Aliases: mcmcplot
### Keywords: hplot

### ** Examples

## Not run: 
##D ## Create fake MCMC output
##D nc <- 10; nr <- 1000
##D pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
##D means <- rpois(10, 20)
##D fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))
##D 
##D ## Use mcmcplot to plot
##D ## the fake MCMC output
##D mcmcplot(fakemcmc)
##D mcmcplot(fakemcmc, greek=TRUE)
##D mcmcplot(fakemcmc, xlim=range(fakemcmc)) # put the densities on the same scale
##D mcmcplot(fakemcmc, "gamma")
##D mcmcplot(fakemcmc, regex="alpha\\[[12]", style="plain")
##D mcmcplot(fakemcmc, "gamma", regex="alpha\\[[12]")
##D mcmcplot(fakemcmc, random=2)
##D mcmcplot(fakemcmc, random=c(2, 3))
##D 
##D ## What happens with NULL varnames?
##D varnames(fakemcmc) <- NULL
##D mcmcplot(fakemcmc)
##D 
##D ## mcmcplot works on bugs objects too
##D library(R2WinBUGS)
##D example("openbugs", "R2WinBUGS")
##D ## from the help file for openbugs:
##D schools.sim <- bugs(data, inits, parameters, model.file,
##D                     n.chains = 3, n.iter = 5000,
##D                     program = "openbugs", working.directory = NULL)
##D mcmcplot(schools.sim)
## End(Not run)



cleanEx()
nameEx("mcmcplot1")
### * mcmcplot1

flush(stderr()); flush(stdout())

### Name: mcmcplot1
### Title: MCMC Diagnostics Plots for one Model Parameter
### Aliases: mcmcplot1
### Keywords: hplot

### ** Examples

## Create fake MCMC output
fakemcmc <- as.mcmc.list(mcmc(sapply(1:5, function(dum) rnorm(1000))))
varnames(fakemcmc) <- c("gamma[1,1]", "gamma[1,2]", "gamma[1,3]", "sigma[1]", "sigma[2]")

mcmcplot1(fakemcmc[, "sigma[1]", drop=FALSE])
mcmcplot1(fakemcmc[, "gamma[1,3]", drop=FALSE], style="plain")



cleanEx()
nameEx("mcmcplots-package")
### * mcmcplots-package

flush(stderr()); flush(stdout())

### Name: mcmcplots-package
### Title: Plots for MCMC Output
### Aliases: mcmcplots-package mcmcplots
### Keywords: package

### ** Examples

## Not run: 
##D ## mcmcplots functions work on bugs objects too
##D library(R2WinBUGS)
##D example("openbugs", "R2WinBUGS")
##D ## from the help file for openbugs:
##D schools.sim <- bugs(data, inits, parameters, model.file,
##D                     n.chains = 3, n.iter = 5000,
##D                     program = "openbugs", working.directory = NULL)
##D caterplot(schools.sim, "theta")
##D traplot(schools.sim, "theta")
##D denplot(schools.sim, "theta")
##D mcmcplot(schools.sim)
## End(Not run)

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## Use mcmcplot to plot
## the fake MCMC output
## Not run: 
##D mcmcplot(fakemcmc)
##D mcmcplot(fakemcmc, "gamma")
##D mcmcplot(fakemcmc, regex="alpha\\[[12]")
##D mcmcplot(fakemcmc, "gamma", "alpha\\[[12]")
##D mcmcplot(fakemcmc, random=2)
##D mcmcplot(fakemcmc, random=c(2, 3))
## End(Not run)

## Use traplot to create
## trace plots of fake MCMC data
traplot(fakemcmc)
traplot(fakemcmc, "gamma")
traplot(fakemcmc, "gamma", "alpha\\[[12]]$") # all gamma and alpha[1] and alpha[2]

## Use denplot to create
## density plots of fake MCMC data
denplot(fakemcmc)
denplot(fakemcmc, "gamma")
denplot(fakemcmc, "gamma", "alpha\\[[12]]$") # all gamma and alpha[1] and alpha[2]

## Use caterplot to create
## caterpillar plots of fake MCMC data
## caterplot plots of the fake MCMC output
par(mfrow=c(2,2))
caterplot(fakemcmc, "alpha", collapse=FALSE)
caterplot(fakemcmc, "gamma", collapse=FALSE)
caterplot(fakemcmc, "alpha", labels.loc="axis", col="blue")
caterplot(fakemcmc, "gamma", labels.loc="above", col="red")

## Use denoverplot to create overlaying density plots
## of all parameters in fake MCMC data
fakemcmc2 <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))
denoverplot(fakemcmc, fakemcmc2)

## Use corplot to create a "heat plot" of a
## correlation matrix of the fake MCMC draws
corplot(cor(as.matrix(fakemcmc)), cex.axis=0.75)  ## not exciting
Rho1 <- outer(1:10, 1:10, function(i, j) 0.5^(abs(i-j)))
Rho2 <- outer(1:5, 1:5, function(i, j) 0.25^(i!=j))
dat1 <- t(apply(matrix(rnorm(10*1000), 1000, 10), 1, function(z, Rho1) crossprod(Rho1, z), Rho1))
dat2 <- t(apply(matrix(rnorm(5*1000), 1000, 5), 1, function(z, Rho2) crossprod(Rho2,z), Rho2))
colnames(dat1) <- paste("theta[", 1:10, "]", sep="")
colnames(dat2) <- paste("alpha[", 1:5, "]", sep="")
dat <- cbind(dat1, dat2)
parcorplot(dat, "theta", col=gray(31:0/31), cex.axis=0.75)  ## just theta parameters
parcorplot(dat, col=heat.colors(31), cex.axis=0.75)
parcorplot(dat, col=topo.colors(31), cex.axis=0.75)
parcorplot(dat, col=terrain.colors(31), cex.axis=0.75)
parcorplot(dat, col=cm.colors(31), cex.axis=0.75)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mcmcplotsPalette")
### * mcmcplotsPalette

flush(stderr()); flush(stdout())

### Name: mcmcplotsPalette
### Title: Color Palette for the mcmcplots Package
### Aliases: mcmcplotsPalette
### Keywords: color

### ** Examples

colorpie <- function(n, type="rainbow") pie(rep(1, n), col=mcmcplotsPalette(n, type=type))
colorpie(1)
colorpie(8)
colorpie(4, type="sequential")
colorpie(4, type="grayscale")

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

denplot(fakemcmc)
denplot(fakemcmc, style="plain", col=mcmcplotsPalette(3, type="sequential"))
denplot(fakemcmc, style="plain", col=mcmcplotsPalette(3, type="grayscale"))




cleanEx()
nameEx("parcorplot")
### * parcorplot

flush(stderr()); flush(stdout())

### Name: parcorplot
### Title: Correlation Plot for Posterior Draws of Model Parameters
### Aliases: parcorplot
### Keywords: hplot

### ** Examples

Rho1 <- outer(1:10, 1:10, function(i, j) 0.5^(abs(i-j)))
Rho2 <- outer(1:5, 1:5, function(i, j) 0.25^(i!=j))
dat1 <- t(apply(matrix(rnorm(10*1000), 1000, 10), 1, function(z, Rho1) t(Rho1)%*%z, Rho1))
dat2 <- t(apply(matrix(rnorm(5*1000), 1000, 5), 1, function(z, Rho2) t(Rho2)%*%z, Rho2))
colnames(dat1) <- paste("theta[", 1:10, "]", sep="")
colnames(dat2) <- paste("alpha[", 1:5, "]", sep="")
dat <- cbind(dat1, dat2)
parcorplot(dat, "theta", col=gray(31:0/31), cex.axis=0.75)
parcorplot(dat, col=heat.colors(31), cex.axis=0.75)
parcorplot(dat, col=topo.colors(31), cex.axis=0.75)
parcorplot(dat, col=terrain.colors(31), cex.axis=0.75)
parcorplot(dat, col=cm.colors(31), cex.axis=0.75)



cleanEx()
nameEx("parms2plot")
### * parms2plot

flush(stderr()); flush(stdout())

### Name: parms2plot
### Title: Matches groups of parameters to plot in MCMC diagnostics plots.
### Aliases: parms2plot
### Keywords: utilities

### ** Examples

prm <- c(paste("gamma[", 1:30, "]", sep=""),paste("alpha[", 1:20, "]", sep=""))

parms2plot(prm, NULL, NULL, NULL)      # returns all
parms2plot(prm, NULL, NULL, 5)         # returns 5 randomly from each group
parms2plot(prm, NULL, NULL, c(5, 10))  # 5 from gamma, 10 from alpha
parms2plot(prm, NULL, NULL, c(10, NA)) # 10 from gamma, all from alpha
parms2plot(prm, "alpha", NULL, NULL)   # all alphas
parms2plot(prm, "gamma", NULL, NULL)   # all gamma
parms2plot(prm, NULL, "alpha\\[1[[:digit:]]\\]$", NULL)   # alpha[10]-alpha[19]
parms2plot(prm, "gamma", "alpha\\[1[[:digit:]]\\]$", NULL)  # all gamma and alpha[10]-alpha[19]



cleanEx()
nameEx("traplot")
### * traplot

flush(stderr()); flush(stdout())

### Name: traplot
### Title: Traceplots of Multiple Parameters.
### Aliases: traplot
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## Plot traces of the fake MCMC output
traplot(fakemcmc)
traplot(fakemcmc, style="plain")
traplot(fakemcmc, "gamma", greek=TRUE)

## What happens with NULL varnames?
varnames(fakemcmc) <- NULL
traplot(fakemcmc)

## Not run: 
##D ## traplot works on bugs objects too
##D library(R2WinBUGS)
##D example("openbugs", "R2WinBUGS")
##D ## from the help file for openbugs:
##D schools.sim <- bugs(data, inits, parameters, model.file,
##D                     n.chains = 3, n.iter = 5000,
##D                     program = "openbugs", working.directory = NULL)
##D traplot(schools.sim, "theta")
## End(Not run)




cleanEx()
nameEx("traplot1")
### * traplot1

flush(stderr()); flush(stdout())

### Name: traplot1
### Title: Trace Plot for a Single Parameter in MCMC Output
### Aliases: traplot1
### Keywords: hplot

### ** Examples

## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means,each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

traplot(fakemcmc[, "alpha[5]", drop=FALSE])
traplot(fakemcmc[, "alpha[5]", drop=FALSE], style="plain")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
