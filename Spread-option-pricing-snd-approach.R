rm(list=ls())
source('Spread-option-pricing-functions.R')
## test params
# NOTE: params are flipped in comparision with DEMPSTER HONG
underlyings = list()
# underlyings templates (for SV and GBM)
underlyings[[1]] <- create.underlyingAsset(S_0 = 100, sigma = 1.0, delta = 0.05)
underlyings[[2]] <- create.underlyingAsset(S_0 = 96, sigma = 0.5, delta = 0.05)
underlyings[[3]] = underlyings[[1]]; underlyings[[3]]$setSigma(0.2)
underlyings[[4]] = underlyings[[2]]; underlyings[[4]]$setSigma(0.1)


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
u_min = 40
# decay terms for square-intergability
alpha_1 = -3
alpha_2 = 1
# interest rate and dividend rates
r = 0.1
# the strikes range (from, to, step)
K = seq(-4,4,0.5)
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$SV
## END ~~ PARAMS SETTINGS ~~ END ##

charFunction = getCharFunction(modelType, initConds = FALSE)
if (modelType == modelNames$SV) {
  underlying1 = underlyings[[1]]
  underlying2 = underlyings[[2]]
} else if (modelType == modelNames$GBM) {
  underlying1 = underlyings[[3]]
  underlying2 = underlyings[[4]]
}

## DEDICATED METHOD OF PARAMS SETUP ##
# actual value of underlyings and SV
S_1_0 = 100
S_2_0 = 96
nu_0 = 0.04
delta_1 = 0.05
delta_2 = 0.05
# volatilities of underlyings and SV
sigma_1 = 1.0
sigma_2 = 0.5
sigma_nu = 0.05
# correlations
ro = 0.5
ro_1 = -0.5
ro_2 = 0.25
# SV params
kappa = 1.0
mu = 0.04

# logarithmisation
s_1_0 = log(S_1_0)
s_2_0 = log(S_2_0)

# helper variables
sigma_12 = sigma_21 = sigma_1*sigma_2*ro
sigma_1nu = sigma_nu1 = sigma_1*sigma_nu*ro_1
sigma_2nu = sigma_nu2 = sigma_2*sigma_nu*ro_2
S_cov = matrix(
  c(sigma_1^2, sigma_12,
    sigma_21, sigma_2^2), nrow = 2)

sigma_nu_covs = t(c(sigma_1nu, sigma_2nu)) # row vector


#===============================================================#
#               HURD ZHOU FOURIER APPROACH                      #
#===============================================================#
t1 = proc.time()

# 1. calculation of input
# 2. transformation of input
# 3. calculate Spread option value
if (modelType == modelNames$SV) {
  sigma_1 = 1.0
  sigma_2 = 0.5

  SpreadOptionPrice = K * NA
  for (kk in seq_along(K)) {
    if (K[kk] == 0) K[kk] = 1e-5 # small K instead zero

    setLambdasAndDeltas(K = abs(K[kk]), u_min)
    if (K[kk] > 0) {
      fourierInputHurdZhouSV = fourierInputHurdZhou(charFunc_T = charFuncSV_T,
                                                    negativeStrike = FALSE)

      fourierOutputHurdZhouSV = fftw_c2c_2d(fourierInputHurdZhouSV,inverse = T)
    } else if (K[kk] < 0) {
      fourierInputHurdZhouSV.negStrike = fourierInputHurdZhou(charFunc_T = charFuncSV_T,
                                                              negativeStrike = TRUE)
      fourierOutputHurdZhouSV.negStrike = fftw_c2c_2d(fourierInputHurdZhouSV.negStrike,inverse = T)
    }

    SpreadOptionPrice[kk] = abs(K[kk]) *
      abs(SprHurdZhou(K = K[kk], modelType = modelType))
  }
    cat(sprintf("\n%.2f: %.4f",K, SpreadOptionPrice))
    cat("\n")

} else if (modelType == modelNames$GBM) {
  sigma_1 = 0.2
  sigma_2 = 0.1

  K = seq(0.5,4,0.5)
  SpreadOptionPrice = K * NA
  for (kk in seq_along(K)) {
    setLambdasAndDeltas(K = K[kk], u_min)

    fourierInputHurdZhouGBM = fourierInputHurdZhou(charFunc_T = charFuncGBM_T)

    fourierOutputHurdZhouGBM = fftw_c2c_2d(fourierInputHurdZhouGBM,inverse = T)

    SpreadOptionPrice[kk] = K[kk] *
      abs(SprHurdZhou(K = K[kk], modelType = modelType))
  }
  cat(sprintf("\n%.2f: %.4f",K, SpreadOptionPrice))
  cat("\n")

} else {
  stop("Undefined modelType")
}

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

# m = n = seq_len(N)-1
# combs = expand.grid(n = n, m = m) # all possible combinations of n and m
# m = combs$m; n = combs$n
# u_1 = v_1(m)+i*alpha_1
# u_2 = v_2(n)+i*alpha_2
#
# t11 = proc.time()
# charFuncResult = charFunc_T(u_1, u_2)
# t12 = proc.time()
#
# print(unname(t12[3] - t11[3]))
#
# t11 = proc.time()
# complexGammaPayoffCore = gammaz(i*(u_1+u_2)-1)
# t12 = proc.time()
#
# print(unname(t12[3] - t11[3]))
