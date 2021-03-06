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
# number of discretization, steps etc.
N = 512
# trancation error param
u_min = 158
# decay terms for square-intergability
alpha_1 = -3
alpha_2 = 1
# interest rate and dividend rates
r = 0.1
# upper bound params epsilon
epsilon = 1e-1
# the strikes range (from, to, step)
K = seq(1,4,0.1)
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$SV
## END ~~ PARAMS SETTINGS ~~ END ##

# NOTE: same templates, but reversed order
if (modelType == modelNames$SV) {
  underlying1 = underlyings[[2]]
  underlying2 = underlyings[[1]]
} else if (modelType == modelNames$GBM) {
  underlying1 = underlyings[[4]]
  underlying2 = underlyings[[3]]
}

# last strike
K_last = underlying2$getInitPrice() - underlying1$getInitPrice()

setLambdasAndDeltas(K = 1, u_minimum = u_min,
  underlying1 = underlying1, underlying2 = underlying2)

#===============================================================#
#               BASIC FOURIER APPROACH                          #
#===============================================================#
t1 = proc.time()
inverse = F
# 1. preparation of input
# 2. transformation of input
# 3. calculate Pi components
if (modelType == modelNames$SV) {
  cat("Preparation of input for Fourier transform ...\n")
  fourierInputXChi1SV = fourierInputX(charFuncChi1, charFuncSV)
  fourierInputXChi2SV = fourierInputX(charFuncChi2, charFuncSV)

  cat("Calculation of Fourier transform output ...\n")
  fourierOutputDempsterHongChi1SV = fftw_c2c_2d(fourierInputXChi1SV,inverse = inverse)
  fourierOutputDempsterHongChi2SV = fftw_c2c_2d(fourierInputXChi2SV,inverse = inverse)

  cat("Calculation of Pi components ...\n")
  componentPi1.Under.SV = componentPi1.Under(modelType = modelType)
  componentPi2.Under.SV = componentPi2.Under(modelType = modelType)
  componentPi1.Over.SV = componentPi1.Over(modelType = modelType)
  componentPi2.Over.SV = componentPi2.Over(modelType = modelType)

} else if (modelType == modelNames$GBM) {
  cat("Preparation of input for Fourier transform ...\n")
  fourierInputXChi1GBM = fourierInputX(charFuncChi1, charFuncGBM)
  fourierInputXChi2GBM = fourierInputX(charFuncChi2, charFuncGBM)

  cat("Calculation of Fourier transform output ...\n")
  fourierOutputDempsterHongChi1GBM = fftw_c2c_2d(fourierInputXChi1GBM,inverse = inverse)
  fourierOutputDempsterHongChi2GBM = fftw_c2c_2d(fourierInputXChi2GBM,inverse = inverse)

  cat("Calculation of Pi components ...\n")
  componentPi1.Under.GBM = componentPi1.Under(modelType = modelType)
  componentPi2.Under.GBM = componentPi2.Under(modelType = modelType)
  componentPi1.Over.GBM = componentPi1.Over(modelType = modelType)
  componentPi2.Over.GBM = componentPi2.Over(modelType = modelType)

} else {
  stop("Undefined modelType")
}

spreadsFFT <- getSpreadOptionPriceRange(K = K, modelType = modelType, r = r, T_t = T_t)
print(spreadsFFT)

t2 = proc.time()

# results (time consumtion)
times = if ("times" %in% ls()) {
  times
} else {
  data.frame(N = 2^(9:12), time = 0)
}
times[times$N == N,"time"] = unname(t2[3] - t1[3])
cat("\nTime elapsed:\n")
print(times)
cat("\n\n")
