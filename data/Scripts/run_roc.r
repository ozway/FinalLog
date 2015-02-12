#Rprof(filename = "BeforeChange.rprof")

## check if all packages are installed to run this script
check.pack <- c( "cubfits" %in% rownames(installed.packages()), "psych" %in% rownames(installed.packages()), 
                "VGAM" %in% rownames(installed.packages()), "coda" %in% rownames(installed.packages()),
                "getopt" %in% rownames(installed.packages()) )#, "Rmpi" %in% rownames(installed.packages()))
requ.pack <- c("cubfits", "psych", "VGAM", "coda", "getopt") #, "Rmpi")
if(sum(check.pack) != length(check.pack))
{
  cat("Missing package(s): ");cat(requ.pack[!check.pack]); cat("\n")
  stop("Install missing package(s)")
}


#suppressMessages(library(cubfits, quietly = TRUE))
suppressMessages(library(cubfits, lib.loc="~/cubfitsBuild/", quietly = TRUE))
suppressMessages(library(psych, quietly = TRUE))
##suppressMessages(library(Rmisc, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))

##load helper functions used to process current project.
## not general enough to include in package.
source("run_utility.r")
source("config.r")


##Process input from command line
# using getopt library

##define matrix for options and flags that are passed to "Rscript run_roc.r ..." command 
spec = matrix(c(
  'sdlog' , 's', 2, "character", # one or more sd logs
  'cubmethod', 'c', 2, "character", # cub method to use
  'fnin' , 'f', 1, "character", # genome (fasta)
  'fnpin' , 'p', 2, "character", # x obs
  'fnout' , 'o', 1, "character", # output folder
  'fname' , 'n', 1, "character", # filename
  'pinit' , 'i', 2, "character" # p matrix initial values
), byrow=TRUE, ncol=4);

cat("===================================================================================\n")
cat("================================ START HEADER =====================================\n")
cat("===================================================================================\n")
args <- commandArgs(trailingOnly = TRUE)
cat("Function call:\nRscript run_roc.r ");cat(args);cat("\n\n")

opt <- getopt(spec, args)
## print config file to log
print.config(config)
cat("===================================================================================\n")
cat("================================= END HEADER ======================================\n")
cat("===================================================================================\n\n\n")

##Assigning options to variables
##Question: why do they have different names?
debug.code <- FALSE
if(debug.code){
  ## hard coded variables to be used when debugging.
  ## to be removed later?
  cat("=================== WARNING ===================\n")
  cat("======== Script running in debug mode! ========\n")
  cat("===============================================\n")
  ## I/O variables
  
  sdlog.phi.init <- c(0.5,1,2,4) # has to be non 0 for cubappr
  fn.in <- "../ecoli/data/ecoli_K12_MG1655_genome_filtered.fasta"
  fn.phi.in <- "../ecoli/data/ecoli_X_obs.csv"
  fname <- "test"
  out.folder <- "../results/test/"
  fn.phi.out <- paste(out.folder, fname, ".phi", sep="")
  fn.out <- paste(out.folder, fname, ".dat", sep="")
  cubmethods <- "cubappr"
  
  ##Hard coding of initial values for pMat of each chain
  ## c(sd_epsilon, meanlog(phi), sdlog(phi), K_x scaling constant)
  p.init <- list(NULL)
  length(p.init) <- config$n.chains
  if(cubmethods == "cubfits"){
    p.init[[1]] <- c(0.1, -0.125, 0.5, 1.0)
    p.init[[2]] <- c(0.5, -0.5, 1.0, 7.0)
    p.init[[3]] <- c(1.5, -2.0, 2.0, 5.0)
    p.init[[4]] <- c(1.0, -8.0, 4.0, 10.0)
  } else if(cubmethods == "cubappr") {
    p.init[[1]] <- c(-0.125, 0.5)
    p.init[[2]] <- c(-0.5, 1.0)
    p.init[[3]] <- c(-2.0, 2.0)
    p.init[[4]] <- c(-8.0, 4.0)
  }  
}else{    
  sdlog.phi.init <- opt$sdlog
  sdlog.phi.init <- as.double(unlist(strsplit(sdlog.phi.init, " ")))
  ## set method for either with Xobs, cubfits, or without Xobs, cubappr,
  ## NOTE: cross validation approach, cubpred, is not supported in this script
  cubmethods <- opt$cubmethod
  fn.in <- opt$fnin
  fn.phi.in <- opt$fnpin
  fname <- opt$fname
  out.folder <- opt$fnout
  fn.phi.out <- paste(out.folder, fname, ".phi", sep="")
  fn.out <- paste(out.folder, fname, ".dat", sep="")
  
  if(!is.null(opt$pinit))
  {
    p.init <- as.list(read.table(opt$pinit, header=F, sep=","))
  }else{ # start from default, no initial values set
    p.init <- list(NULL)
    length(p.init) <- config$n.chains
  }
}

cat(paste("started at:", Sys.time(), "\n"))
cat(paste("using", cubmethods, "\n"))

cat(paste("reading sequences from file", fn.in, "\n"))
## readGenome function in cubfits package
seq.string <- readGenome(fn.in, config$rm.short, config$rm.first.aa)


## if using Xobs data 
if(cubmethods == "cubfits")
{
  cat(paste("reading gene expression measurements (Xobs) from file\n", fn.phi.in, "\nand compare to ORF list from FASTA file\n", fn.in, "\n"))
  ret.list <- read.empirical.data(fn.phi.in, seq.string, config$selected.env, 0)
  phi.obs <- ret.list$empirical
  phi.obs <- phi.obs[order(names(phi.obs))]
  seq.string <- ret.list$genome
  rm("ret.list")
  
  ## setting the sdlog.phi.init to sd(log(Xobs)) if sdlog.phi.init == 0. 
  for(i in 1:config$n.chains)
  {  
    sdlog.phi.init[i] <- ifelse(sdlog.phi.init[i] == 0, sdlog.phi.init[i] <- sd(log(phi.obs)), sdlog.phi.init[i]) 
  }
}

cat("generating list of codon position in ORFs for each AA...\ngenerating list of number of AA occurences per ORF...\ngenerating list of codon counts per ORF...\n")
data <- generate.data(seq.string, config$aa)


## 
cat("generate initial phi using ")
init.phi <- list()
length(init.phi) <- config$n.chains
if(config$use.scuo)
{
  cat("SCUO with sd(ln(phi)) values\n")
  scuo <- gen.scuo(seq.string, config$aa)
  scuo <- calc_scuo_values(scuo)$SCUO
  for(i in 1:config$n.chains)
  {
    cat( paste("\t-", sdlog.phi.init[i],"\n") )
    randscuo <- scuo.random(scuo, meanlog = -sdlog.phi.init[i]^2.0 / 2.0, sdlog = sdlog.phi.init[i])
    randscuo <- randscuo / mean(randscuo)
    names(randscuo) <- names(seq.string)
    init.phi[[i]] <- randscuo
  }
}else{
  cat("Xobs data. Note: All chains starting with same initial phi values\n")
  if(cubmethods == "cubfits")
  {
    for(i in 1:config$n.chains){ init.phi[[i]] <- phi.obs }
  }else{
    ret.list <- read.empirical.data(fn.phi.in, seq.string, config$selected.env)
    phi.obs <- ret.list$empirical
    phi.obs <- phi.obs[order(names(phi.obs))]
    for(i in 1:config$n.chains){ init.phi[[i]] <- phi.obs }
  }
}

##define objects to hold initial values and output from MCMC chains
results <- list()
length(results) <- config$n.chains

funct <- ifelse(config$n.chains > 1, "cubmultichain", "cubsinglechain") 
cat(paste("running", cubmethods, "using", funct, "\n"))
seeds <- round(runif(config$n.chains, 1, 100000))
cat("\t with seeds: ");cat(seeds);cat("\n")
runtime.info <- system.time(
  {
    .CF.CT$parallel <- config$parallel
    .CF.CONF$compute.logL <- T
    .CF.CT$prior.dist <- "normal"
    iterations <- config$n.samples*config$chain.thin # between convergence checks
    if(cubmethods == "cubfits"){
      .CF.CT$type.p <- "lognormal_bias"
      .CF.CONF$scale.phi.Obs <- F
      .CF.CONF$estimate.bias.Phi <- T
      if(config$n.chains < 2)
      {
        results <- cubsinglechain(cubmethods, frac1=config$frac1, frac2=config$frac2, 
                               reset.qr=config$reset.qr, seed=seeds[1], teston="sphi",
                               min=config$min.samples, max=config$max.samples, conv.thin=config$conv.thin, eps=config$eps, 
                               reu13.df.obs=data$reu13.df, phi.Obs=phi.obs, y=data$y, n=data$n, phi.Init=init.phi[[1]],
                               nIter=iterations, p.Init=p.init[[1]], iterThin=config$chain.thin,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }else{
        results <- cubmultichain(cubmethods, reset.qr=config$reset.qr, seeds=seeds, teston="sphi", 
                               swap=config$swap, swapAt=config$swapAt, min=config$min.samples, max=config$max.samples, 
                               nchains=config$n.chains, conv.thin=config$conv.thin, eps=config$eps,
                               ncores=config$n.cores, reu13.df.obs=data$reu13.df, phi.Obs=phi.obs, y=data$y, n=data$n, phi.Init=init.phi,
                               nIter=iterations, p.Init=p.init, iterThin=config$chain.thin,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }
    }else if(cubmethods == "cubappr") {
      if(config$n.chains < 2)
      {
        results <- cubsinglechain(cubmethods, frac1=config$frac1, frac2=config$frac2, 
                               reset.qr=config$reset.qr, seed=seeds[1], teston="sphi",
                               min=config$min.samples, max=config$max.samples, conv.thin=config$conv.thin, eps=config$eps, 
                               reu13.df.obs=data$reu13.df, y=data$y, n=data$n, phi.pred.Init=init.phi[[1]],
                               nIter=iterations, p.Init=p.init[[1]], iterThin=config$chain.thin,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }else{
        results <- cubmultichain(cubmethods, reset.qr=config$reset.qr, seeds=seeds, teston="sphi", 
                               swap=config$swap, swapAt=config$swapAt, min=config$min.samples, max=config$max.samples, 
                               nchains=config$n.chains, conv.thin=config$conv.thin, eps=config$eps,
                               #monitor=function(x, i){cat(paste("I am monitor for chain", i, "\n"))},
                               ncores=config$n.cores, reu13.df.obs=data$reu13.df, y=data$y, n=data$n, phi.pred.Init=init.phi,
                               nIter=iterations, p.Init=p.init, iterThin=config$chain.thin,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }
    }
  }
) 

cat(paste("Elapsed time for", config$n.chains, "chains doing", config$n.iter, "iterations on", config$n.cores, "cores was", round(runtime.info["elapsed"]/60, digits=2), "min\n" ))
seq.string.names <- names(seq.string)
rm("seq.string")



## process results
cat("process results...\n")
phi.pred <- list()
length(phi.pred) <- config$n.chains
b.mat <- list()
length(b.mat) <- config$n.chains
if(config$n.chains > 1)
{
  for(i in 1:config$n.chains)
  {
    if(cubmethods == "cubfits"){
      phi.pred[[i]] <- do.call("cbind", results$chains[[i]]$phi.Mat)### matrix n.G by nIter+1
    }
    if(cubmethods == "cubappr") {
      phi.pred[[i]] <- do.call("cbind", results$chains[[i]]$phi.pred.Mat) ### matrix n.G by nIter+1
    }
    b.mat[[i]] <- do.call("cbind", results$chains[[i]]$b.Mat) ### matrix n.G by nIter+1
  }
}else{
  if(cubmethods == "cubfits"){
    phi.pred[[1]] <- do.call("cbind", results$chains$phi.Mat)### matrix n.G by nIter+1
  }
  if(cubmethods == "cubappr") {
    phi.pred[[1]] <- do.call("cbind", results$chains$phi.pred.Mat) ### matrix n.G by nIter+1
  }
  b.mat[[1]] <- do.call("cbind", results$chains$b.Mat) ### matrix n.G by nIter+1  
}

mean.phis <- list()
median.phis <- list()
geo.mean.phis <- list()
harm.mean.phis <- list()
sd.phi <- list()
mean.b.mat <- list()
sd.b.mat <- list()
interval <- (dim(phi.pred[[i]])[2]-config$use.n.samples):dim(phi.pred[[i]])[2]
for(i in 1:config$n.chains)
{
  mean.phis[[i]] <- rowMeans(phi.pred[[i]][, interval]) ### mean of the last X iterations
  median.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, median)
  geo.mean.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, geometric.mean)
  harm.mean.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, harmonic.mean)
  sd.phi[[i]] <- apply(phi.pred[[i]][, interval], 1, sd)
  mean.b.mat[[i]] <- rowMeans(b.mat[[i]][, interval])
  sd.b.mat[[i]] <- apply(b.mat[[i]][, interval], 1, sd)
}
## save results and input
cat("saving results...\n")
mean.phis <- do.call("cbind", mean.phis)
median.phis <- do.call("cbind", median.phis)
geo.mean.phis <- do.call("cbind", geo.mean.phis)
harm.mean.phis <- do.call("cbind", harm.mean.phis)
sd.phi <- do.call("cbind", sd.phi)
mean.b.mat <- do.call("cbind", mean.b.mat)
sd.b.mat <- do.call("cbind", sd.b.mat)

bmat.names <- mapBMatNames(rownames(mean.b.mat), config$aa)
if(config$n.chains > 1) # only need that if there is more than one run
{
  for(i in 1:length(results$chains))
  {
    dir.create(file.path(paste(out.folder, "chain_", i, sep="")), showWarnings = FALSE)
    means <- cbind(seq.string.names, mean.phis[,i], median.phis[,i], geo.mean.phis[,i], harm.mean.phis[,i], sd.phi[,i])
    colnames(means) <- c("ORF_Info", "Phi_Post_Arith_Mean", "Phi_Post_Median", "Phi_Post_Geom_Mean", "Phi_Post_Harm_Mean", "Phi_Post_SD") 
    write.csv(means, file = paste(out.folder, "chain_", i, "/", fname, ".phi", sep=""), row.names=F, quote=F)
    
    mean.b <- cbind(bmat.names, mean.b.mat[,i], sd.b.mat[,1])
    colnames(mean.b) <- c("Parameter", "Value", "SD")
    write.csv(mean.b, file = paste(out.folder, "chain_", i, "/", fname, ".bmat", sep=""), row.names=F, quote=F)
  } 
  mean.phis <- rowMeans(mean.phis)
  median.phis <- apply(median.phis, 1, median)
  geo.mean.phis <- apply(geo.mean.phis, 1, geometric.mean)
  harm.mean.phis <- apply(harm.mean.phis, 1, harmonic.mean)
  mean.b.mat <- rowMeans(mean.b.mat)
  sd.b.mat <- rowMeans(sd.b.mat)
  sd.phi <- rowMeans(sd.phi)
}
means <- cbind(seq.string.names, mean.phis, median.phis, geo.mean.phis, harm.mean.phis, sd.phi)
colnames(means) <- c("ORF_Info", "Phi_Post_Arith_Mean", "Phi_Post_Median", "Phi_Post_Geom_Mean", "Phi_Post_Harm_Mean", "Phi_Post_SD") 
write.csv(means, file = fn.phi.out, row.names=F, quote=F)


if(config$n.chains > 1){
#mean.b.mat <- cbind(bmat.names, mean.b.mat, sd.b.mat)
mean.b.mat <- cbind(names(results$chains[[1]]$b.Mat[[1]]), mean.b.mat, sd.b.mat)
}
else{
mean.b.mat <- cbind(names(results$chains$b.Mat[[1]]), mean.b.mat, sd.b.mat)
}


colnames(mean.b.mat) <- c("Parameter", "Value", "SD")
write.csv(mean.b.mat, file = paste(out.folder, fname, ".bmat", sep=""), row.names=F, quote=F)

if(config$n.chains > 1)
{
  for(i in 1:config$n.chains)
  {
    chain <- results$chains[[i]]
    convergence <- results$convergence
    list.save <- c("data", "mean.phis", "seeds", "chain", "seq.string.names", "config", "means", "sdlog.phi.init", "convergence")    
    save(list = list.save, file = paste(out.folder, "chain_", i, "/", fname, ".dat", sep=""))
  }
}else{
  chain <- results$chains
  convergence <- results$convergence
  list.save <- c("data", "mean.phis", "seeds", "chain", "seq.string.names", "config", "means", "sdlog.phi.init", "convergence")    
  save(list = list.save, file = paste(out.folder, fname, ".dat", sep=""))
}
cat(paste("finished at:", Sys.time(), "\n"))
rm("phi.pred")


#Rprof(NULL)
