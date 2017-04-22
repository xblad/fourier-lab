rm(list=ls()[!ls() %in% "times"])
source('Spread-option-pricing-functions.R')
## test params
# NOTE: params are flipped in comparision with DEMPSTER HONG
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
# number of discretization, steps etc.
N = 256
# trancation error param
u_min = 20
# decay terms for square-intergability
alpha_1 = -3
alpha_2 = 1
# interest rate
r = 0.1
# the strikes range (from, to, step)
K = seq(-40,40,1)
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$SV
## END ~~ PARAMS SETTINGS ~~ END ##

#===============================================================#
#               HURD ZHOU FOURIER APPROACH                      #
#===============================================================#
charFunction = getCharFunction(modelType, initConds = FALSE)
if (modelType == modelNames$SV) {
  underlying1 = underlyings[[1]]
  underlying2 = underlyings[[2]]
  underlying1.negStrike = underlyings[[2]]
  underlying2.negStrike = underlyings[[1]]
} else if (modelType == modelNames$GBM) {
  underlying1 = underlyings[[3]]
  underlying2 = underlyings[[4]]
  underlying1.negStrike = underlyings[[4]]
  underlying2.negStrike = underlyings[[3]]
}

t1 = proc.time()

SpreadOptionPrices = getSpreadOptionPricesHurdZhou(K = K, r = r, T_t = T_t)

t2 = proc.time()

print(SpreadOptionPrices)
#cat(sprintf("\n%.2f: %.4f",K, SpreadOptionPrices))
cat("\n")

matplot(SpreadOptionPrices,type=c("l"),pch=1,x = K,ylab = "Spread option value")
legend("top",legend = c("Call","Put"), col = 1:2, pch=1)

# results (time consumtion)
times = if ("times" %in% ls()) {
  times
} else {
  data.frame(N = 2^(3:12), time = 0)
}
times[times$N == N,"time"] = unname(t2[3] - t1[3])
cat("\nTime elapsed:\n")
print(times)
