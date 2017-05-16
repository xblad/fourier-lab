rm(list=ls())
source('Spread-option-pricing-functions.R')

underlyings = list()
# underlyings templates (for SV and GBM)
underlyings[[1]] <- create.underlyingAsset(S_0 = 100, sigma = 1.0, delta = 0.05)
underlyings[[2]] <- create.underlyingAsset(S_0 = 96, sigma = 0.5, delta = 0.05)
underlyings[[3]] <- create.underlyingAsset();
underlyings[[3]]$copyFromUnderlyingAsset(underlyings[[1]]);
underlyings[[3]]$setSigma(0.2);
underlyings[[4]] <- create.underlyingAsset();
underlyings[[4]]$copyFromUnderlyingAsset(underlyings[[2]]);
underlyings[[4]]$setSigma(0.1);


volatility <- create.stochasticVolatility(nu_0 = 0.04, sigma = 0.05,
                                          kappa = 1.0, mu = 0.04)

corrs <- create.correlations()
corrs$setCorrelation(underlyings[[1]], underlyings[[2]], 0.5)
corrs$setCorrelation(underlyings[[1]], volatility, -0.5)
corrs$setCorrelation(underlyings[[2]], volatility, 0.25)


## PARAMS SETTINGS ##
# time to expiry
T_t = 1.0
# interest rate and dividend rates
r = 0.1
# the strikes range (from, to, step)
K = seq(-40,40,1)
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$SV
# Monte carlo params
n_sim = 10^4
sim_timesteps = 1 # irrelevant in case of MonteCarloGBM3 or MonteCarloGBM4
## END ~~ PARAMS SETTINGS ~~ END ##

if (modelType == modelNames$SV) {
  underlying1 = underlyings[[1]]
  underlying2 = underlyings[[2]]
} else if (modelType == modelNames$GBM) {
  underlying1 = underlyings[[3]]
  underlying2 = underlyings[[4]]
}

#===============================================================#
#               MONTE CARLO SIMULATION                          #
#===============================================================#
t1 = proc.time()

if (modelType == modelNames$GBM) {
  mc_sprds <- monteCarloGBM4(underlying1 = underlying1, underlying2 = underlying2,
    corrs = corrs, r = r, T_t = T_t)#,sim_timesteps = sim_timesteps)
} else if (modelType == modelNames$SV) {
  mc_sprds = NULL
  for (ii in 1:100) {
  mc_sprds_temp <- monteCarloSV4(underlying1 = underlying1, underlying2 = underlying2,
    corrs = corrs, volatility = volatility, r = r, T_t = T_t)
  mc_sprds = cbind(mc_sprds, mc_sprds_temp)
  }
} else {
  stop("Unknown modelType!")
}

cat("\nCalls:\n")
mc_calls <- optionPrices(S_T=mc_sprds,K=K,r=r,T_t=T_t,call=T)
cat(sprintf("K = %.1f: %.6f (%.3f)\n",K,mc_calls$prices,mc_calls$se     ))
cat("\nPuts:\n")
mc_puts <- optionPrices(S_T=mc_sprds,K=K,r=r,T_t=T_t,call=F)
cat(sprintf("K = %.1f: %.6f (%.3f)\n",K,mc_puts$prices,mc_puts$se))

t2 = proc.time()
cat("\nTime elapsed:\n")
cat(unname(t2[3] - t1[3]))
