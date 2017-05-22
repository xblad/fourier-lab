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
# model for joint characteristic function (GBM, SV etc.)
modelType = modelNames$SV
# default strike
K = 5
# random variables params
optType = c('Call','Put')[1] # Call or put
# underlyings, volatility and correlations
underlyings = list()
# underlyings templates (for SV and GBM)
underlyings[[1]] <- create.underlyingAsset(S_0 = 100, sigma = 1.0, delta = 0.05)
underlyings[[2]] <- create.underlyingAsset(S_0 = 96, sigma = 0.5, delta = 0.05)


charFunction = getCharFunction(modelType, initConds = FALSE)
if (modelType == modelNames$GBM) {
        underlyings[[1]]$setSigma(0.2);
        underlyings[[2]]$setSigma(0.1);
}

volatility <- create.stochasticVolatility(nu_0 = 0.04,
                                          sigma = 0.05, kappa = 1.0, mu = 0.04)

corrs <- create.correlations()
corrs$setCorrelation(underlyings[[1]], underlyings[[2]], 0.5)
corrs$setCorrelation(underlyings[[1]], volatility, -0.5)
corrs$setCorrelation(underlyings[[2]], volatility, 0.25)

## END ~~ FIXED PARAMS ~~ END ##

# fine of interval
num = 21

# variables
S_2_0s = seq(60,140,length.out = num);
Ks = seq(-40,40,length.out = num)
Ts = seq(0,3,length.out = num)
rhos = seq(-1,1,length.out = num)
mus = seq(0,1,length.out = num)
sigma_2s = seq(0,1,length.out = num)
params = expand.grid(sigma_2s = sigma_2s, mus = mus)#Ks = Ks, rhos = rhos)#Ts = Ts)#S_2_0s = S_2_0s)#
SpreadOptionPrices = data.frame(sigma_2s = params$sigma_2s,#rhos = params$rhos,#Ts = params$Ts,#S_2_0s = params$S_2_0s,#
                                mus = params$mus,#Ks = params$Ks,
                                FFT_c = NA,
                                FFT_p = NA)
timeConsum = SpreadOptionPrices
for (j in 1:nrow(params)) {
        progress(j,nrow(params))
        #underlyings[[2]]$setInitPrice(params$S_2_0s[j]);
        #K = params$Ks[j]
        #T_t = params$Ts[j]
        #corrs$setCorrelation(underlyings[[1]], underlyings[[2]], params$rhos[j])
        volatility$setMu(params$mus[j])
        underlyings[[2]]$setSigma(params$sigma_2s[j]);
        
        if (F) { # print or not
                print(underlyings[[1]]$getParams())
                print(underlyings[[2]]$getParams())
                print(volatility$getParams())
                print(corrs$getCorrelationMatrix())
                cat("K =",K)
                cat("\n")
        }
        
        underlying1 = underlyings[[1]]
        underlying2 = underlyings[[2]]
        underlying1.negStrike = underlyings[[2]]
        underlying2.negStrike = underlyings[[1]]
        
        #===============================================================#
        #               HURD ZHOU FOURIER APPROACH                      #
        #===============================================================#
        
        t1 = proc.time()[3]
        
        sop = getSpreadOptionPricesHurdZhou(K = K, r = r, T_t = T_t)
        SpreadOptionPrices$FFT_c[j] = sop[,"Call"]
        SpreadOptionPrices$FFT_p[j] = sop[,"Put"]
        
        t2 = proc.time()[3]
        
        timeConsum$FFT_c[j] = t2-t1
        
}
