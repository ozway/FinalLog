print.config <- function(config)
{
  cat("MCMC parameters:\n")
  cat(paste("Number of samples between checks:", config$n.samples, "\n"))
  cat(paste("Min samples:", config$min.samples, "\n"))
  cat(paste("Max samples:", config$max.samples, "\n"))
  cat(paste("Reset QR until:", config$reset.qr, "samples\n"))
  cat(paste("Thining: store every", config$chain.thin, "iteration as sample\n"))
  cat(paste("Swap", config$swap * 100, "% of b matrix\n"))
  cat(paste("Swap b matrix if L1-norm is <", config$swapAt, "\n"))
  
  cat("\n")
  
  cat("Simulation parameters\n")
  cat(paste("Number of Cores", config$n.cores, "\n"))
  cat(paste("Number of Chains:", config$n.chains, "\n"))
  cat(paste("Parallel mode within chain:", config$parallel, "\n"))
  cat(paste("Samples used:", config$use.n.samples, "\n"))
  cat(paste("First", config$rm.first.aa, " AAs removed from sequence\n"))
  cat(paste("Sequences with less than", config$rm.short, "AAs ignored\n"))
  cat("List of AAs taken into account:\n\t");cat(config$aa);cat("\n")
  
  cat("\n") 
  
  cat("Convergence criteria\n")
  if(config$n.chains < 2)
  {
    cat("Convergence test: Geweke\n")
    cat(paste("Convergence criterium: Z Score <", config$eps, "\n"))
    cat(paste("% chain used at the begining:", config$frac1, "\n"))
    cat(paste("% chain used at the end:", config$frac2, "\n"))
  }else{
    cat("Convergence test: Gelman & Rubin\n")
    cat(paste("Convergence criterium: Gelman Score <", config$eps, "\n"))
  }
  cat(paste("Use every", config$conv.thin, "sample for convergence test\n"))
  
  cat("\n")
}

NUMBER <- 6000

config <- list(
  delta.a_12 = 0,
  a_2 = 1,

  n.samples = NUMBER,  # number of samples between convergence checks.
    #Set equal to max/min samples to avoid resetting the scale
  use.n.samples = NUMBER/5, #sample size for testing for convergence. 
#  n.samples = 1000,  # number of samples between convergence checks.
#  use.n.samples = 1000, #sample size for testing for convergence. 
                     #If convergence threshold, set in eps, is reached the sample size of our posteriors equals this value. 
		                 #IMPORTANT: Make sure this one is smaller than the thined chain, otherwise saving will crash!
  n.chains = 1, # num chains
  n.cores = 1, # total num of cpu cores (should be about 5*n.chains when using parallel method other then "lapply")
  selected.env = 1, # deprecated (us if more than one dataset is stored in csv, e.g. different conditions)

  min.samples=NUMBER, # minimum samples each chain has to do, convergence criterium is ignored until min.iter is reached 
  max.samples= NUMBER, # maximum samples for each chain. MCMC will be stoped if the chain reached max.iter iterations (convergence criterium is ignored)

  reset.qr=0, # stop resetting qr matrix when checking for convergence after this many samples (after thining)
  conv.thin=1,  # thining for convergence test (recommend 1 if chain.thin != 1, otherwise double thining)
  chain.thin=10, # thining of the chain during runtime. This is done before gathering convergence test sample. See note for conv.thin.
  rm.first.aa=0, # remove first rm.first.aa AAs (after the first codon which is expected to be the start codon)
  rm.short=0, # # ignore sequences with length < rm.short AAs after the first rm.first.aa AAs are removed
  parallel="lapply", # parallel method within chain 
                     # lapply = no parallelization within chain)
                     # mclapply = parallelization within chain.
		                 # Other options are also possible.
  eps=0.15,  # Convergence threshold.
            # Multichain MCMC uses with Gelman test where threshold is |Gelman Score - 1| < eps 
            # Single chain MCMC uses Geweke test where Geweke score < eps.
  gf=0.4, ## growthfactor for convergence test window (n.samples + gf*currSamples)
  frac1=0.1, # Used with single chain MCMC. Part of Geweke test. Value is proportion of chain at the beginning of the chain (use.n.iter window) used to calculate mean for convergence test
  frac2=0.5, # Used with single chain MCMC. Part of Geweke test. Value is proportion of chain at the end of the chain (use.n.iter window) used to calculate mean for convergence test
  swap=0.0, # proportion of b matrix (deltat, logmu) swapped. Which parameters are swapped is random
  swapAt=0.0, # manhattan distance between two consequtive convergence test below which b matrix is swapped
  
  #aa = c("D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "Z"),
  #aa = c("A", "V"), 
  aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z"), # AAs to take into account
  
  use.scuo = T # false means empirical data is used as initial conditions
)

rm(NUMBER)
