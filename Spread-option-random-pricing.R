rm(list=ls())
source('Spread-option-pricing-functions.R')
## FIXED PARAMS ##
# time to expiry
T_t = 1.0
# number of discretization, steps etc.
N = 256
# trancation error param
u_min = 40
# decay terms for square-intergability
alpha_1 = -3
alpha_2 = 1
# interest rate
r = 0.1
# Monte carlo params
n_sim = 10^6
sim_timesteps = 2000 # irrelevant in case of MonteCarloGBM3 or MonteCarloGBM4
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$GBM
## END ~~ FIXED PARAMS ~~ END ##
set.seed(42);

rOptNum = 1
SpreadOptionPrices = data.frame(FFT = seq_len(rOptNum)*NA,
  MC = seq_len(rOptNum)*NA)
timeConsum = SpreadOptionPrices
for (j in 1:rOptNum) {
  underlyings = list()
  # underlyings templates (for SV and GBM)
  underlyings[[1]] <- create.underlyingAsset(S_0 = 100,
    sigma = runif(1,0.01,0.9), delta = 0.05)
  underlyings[[2]] <- create.underlyingAsset(S_0 = runif(1,60,140),
    sigma = runif(1,0.01,0.9), delta = 0.05)

  # Filler condition sigma
  volatility <- create.stochasticVolatility(nu_0 = 0.04,
    sigma = runif(1,0.01,floor(sqrt(2*1.0*0.04)*100)/100), kappa = 1.0, mu = 0.04)

  corrs <- create.correlations()
  corrs$setCorrelation(underlyings[[1]], underlyings[[2]], runif(1,-0.9,0.9))
  corrs$setCorrelation(underlyings[[1]], volatility, runif(1,-0.9,0.9))
  corrs$setCorrelation(underlyings[[2]], volatility, runif(1,-0.9,0.9))

  K = runif(1,-2,2) * underlyings[[2]]$getInitPrice()

  optType = sample(c("call","put"),1)

  print(underlyings[[1]]$getParams())
  print(underlyings[[2]]$getParams())
  print(volatility$getParams())
  print(corrs$getCorrelationMatrix())
  cat("K =",K)
  cat("\nType =",optType)

  charFunction = getCharFunction(modelType, initConds = FALSE)
  underlying1 = underlyings[[1]]
  underlying2 = underlyings[[2]]
  underlying1.negStrike = underlyings[[2]]
  underlying2.negStrike = underlyings[[1]]

  #===============================================================#
  #               HURD ZHOU FOURIER APPROACH                      #
  #===============================================================#

  t1 = proc.time()[3]

  SpreadOptionPrices$FFT[j] = getSpreadOptionPricesHurdZhou(K = K, r = r, T_t = T_t)

  t2 = proc.time()[3]

  timeConsum$FFT[j] = t2-t1

  #===============================================================#
  #               MONTE CARLO SIMULATION                          #
  #===============================================================#
  t1 = proc.time()[3]

  if (modelType == modelNames$GBM) {
    mc_sprds <- monteCarloGBM4(underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, r = r, T_t = T_t)
  } else if (modelType == modelNames$SV) {
    mc_sprds <- monteCarloSV(underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, volatility = volatility, r = r, T_t = T_t)
  } else {
    stop("Unknown modelType!")
  }

  SpreadOptionPrices$MC[j] = optionPrices(S_T=mc_sprds,K=K,r=r,T_t=T_t,
    call=(optType == "call"))

  t2 = proc.time()[3]
  timeConsum$MC[j] = t2-t1
}
