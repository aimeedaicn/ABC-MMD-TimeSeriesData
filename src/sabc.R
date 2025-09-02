library(winference)

sabc <- function(args, ...)
#'@description  Implement the ABC-SMC with a given data discrepancy metric.
{
  param_algo <- list(nthetas = args$nthetas, nmoves = 1, proposal = randomwalk_proposal(), 
                   minimum_diversity = 0.3, R = 2, maxtrials = 1e5)
  smcresults <- wsmc(args$discrepancy, args, param_algo,
       parallel = FALSE, ...)
  return(smcresults)
}


sabc_get_last_samples <- function(sabc.results)
#'@description  Extract the samples from the last SMC step.
{
  samples.df <- wsmc_to_dataframe(sabc.results)
  samples.df <- dplyr::filter(samples.df, step == max(step))
  return(samples.df) 
}
