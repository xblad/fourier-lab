library(fftwtools)
library(magrittr)
library(pracma)
library(digest)
library(dplyr)
options(stringsAsFactors = FALSE)
#===============================================================#
#                    HELPER FUNCTIONS                           #
#===============================================================#
# progress function from klmr/rcane on github.com
progress <- function (x, max = 100, stepTime=0,totalTime=0) {
    percent <- x / max * 100
    cat(sprintf('\r[%-50s] %d%% (%.3fs,%.2fs)',
                paste(rep('=', percent / 2), collapse = ''),
                floor(percent),stepTime,totalTime))
    if (x == max)
        cat('\n')
}

copyToClipboard <- function(dataFrame, col.names = T) {
  write.table(dataFrame, file = paste("clipboard",2^15,sep="-"), sep = "\t", row.names = F, col.names = col.names)
}

copyMatrixToClipboard <- function(mat) {
	write.table(mat, file = paste("clipboard",2^15,sep="-"), sep = "\t")
}

#===============================================================#
#               PREPARATION FUNCTIONS                           #
#===============================================================#
# imaginary unit
i = complex(imaginary = 1)

modelNames = data.frame(
    SV = "SV",
    GBM = "GBM"
  )

create.underlyingAsset <- function(S_0 = NA, sigma = NA, delta = NA) {
  hash <- substring(digest(rnorm(10)), 0, 6)
  list(
    getHash = function() return(hash),

    setInitPrice = function(S_0) S_0 <<- S_0,
    setSigma = function(sigma) sigma <<- sigma,
    setDividendRate = function(delta) delta <<- delta,

    getInitPrice = function() return(S_0),
    getSigma = function() return(sigma),
    getDividendRate = function() return(delta),

    setParams = function(S_0, sigma, delta) {
      S_0 <<- S_0; sigma <<- sigma; delta <<- delta;
    },
    copyFromUnderlyingAsset = function(underlyingAsset, withHash = TRUE) {
      S_0 <<- underlyingAsset$getInitPrice()
      sigma <<- underlyingAsset$getSigma()
      delta <<- underlyingAsset$getDividendRate()
      if (withHash) hash <<- underlyingAsset$getHash()
    },
    getParams = function()
      return(data.frame(
        S_0 = S_0,
        sigma = sigma,
        delta = delta
      )),
    showParams = function() {
	     cat(sprintf("Underlying has following parameters:
       Initial price:\t%.2f;
       Sigma:\t\t%.2f;
       Dividend rate:\t%.2f.", S_0, sigma, delta), "\n\n")
    }
  )
}

create.stochasticVolatility <- function(nu_0, sigma, kappa, mu) {
  hash <- substring(digest(rnorm(10)), 0, 6)
  list(
    getHash = function() return(hash),
    setInitVolatility = function(nu_0) nu_0 <<- nu_0,
    setSigma = function(sigma) sigma <<- sigma,
    setKappa = function(kappa) kappa <<- kappa,
    setMu = function(mu) mu <<- mu,

    getInitVolatility = function() return(nu_0),
    getSigma = function() return(sigma),
    getKappa = function() return(kappa),
    getMu = function() return(mu),

    setParams = function(nu_0, sigma, kappa, mu) {
      nu_0 <<- nu_0; sigma <<- sigma; kappa <<- kappa; mu <<- mu
    },
    getParams = function()
      return(data.frame(
        nu_0 = nu_0,
        sigma = sigma,
        kappa = kappa,
        mu = mu
      )),
    showParams = function() {
	     cat(sprintf("stochastic volatility has following parameters:
       Initial value:\t%.2f;
       Sigma:\t\t%.2f;
       Kappa:\t\t%.2f;
       Mu:\t\t%.2f.", nu_0, sigma, kappa, mu), "\n\n")
    }
  )
}

create.correlations <- function() {
  correlations <- list()
  procNames <- data.frame(hash = character(), name = character())

  list(
    setCorrelation = function(process_1, process_2, corr) {
      if (!process_1$getHash() %in% procNames$hash) {
        procNames <<- rbind(procNames,
          data.frame(hash = process_1$getHash(),
                     name = deparse(substitute(process_1))))
      }
      if (!process_2$getHash() %in% procNames$hash) {
        procNames <<- rbind(procNames,
          data.frame(hash = process_2$getHash(),
                     name = deparse(substitute(process_2))))
      }
      correlations[[process_1$getHash()]][[process_2$getHash()]]$corr <<- corr
      correlations[[process_2$getHash()]][[process_1$getHash()]]$corr <<- corr
    },
    getCorrelation = function(process_1, process_2) {
      if (process_1$getHash() == process_2$getHash()) return(1)
      return(correlations[[process_1$getHash()]][[process_2$getHash()]]$corr)
    },
    getCorrelationMatrix = function() {
      uniqueNames = sort(unique(procNames$name))
      if (length(uniqueNames) == 0) return(NULL)
      grid = expand.grid(row = uniqueNames, col = uniqueNames, corr = NA)
      for (ii in 1:nrow(grid)) {
        if (grid$row[ii] == grid$col[ii]) {
         grid$corr[ii] = 1
        } else {
          corr = correlations[[
          procNames$hash[procNames$name == grid$row[ii]]
          ]][[
          procNames$hash[procNames$name ==
          grid$col[ii]]
          ]]
          if (!is.null(corr)) {
            grid$corr[ii] = corr
          }
        }
      }
      return(matrix(grid$corr, nrow = length(uniqueNames),
                    dimnames = list(uniqueNames,uniqueNames)))
    }
  )
}

getLogPricesIncrements <- function(K, u_minimum = u_min,
  underlying1, underlying2) {
  if (K < 0) stop("Function is undefined for non-positive strikes!")
  S_0 = c(underlying1$getInitPrice(), underlying2$getInitPrice())
  lambdas = NA * S_0
  for (m in seq_along(S_0)) {
    # only positive j make sense, so that we start by N/2,
    # than first u_test will be zero
    for (j in (N/2):(N-1)) {
      u_test = pi*(j-N/2)/log(S_0[m]/K)
      if (u_test > u_minimum) {
        lambdas[m] = pi/u_test
        break
      } else if (j == (N-1)) {
        lambdas[m] = pi/u_test
        warning(
          sprintf("Minimum bound is too high for K=%s,
          last value have been used!",K)
        )
      }
    }
  }

  return(lambdas)
}

setLambdasAndDeltas <- function(K, u_minimum = u_min,
  identicalIncrements = FALSE, underlying1, underlying2) {

  lambdas = getLogPricesIncrements(K = K, u_minimum = u_min,
    underlying1 = underlying1, underlying2 = underlying2)
  lambda_1 <<- lambdas[1]
  lambda_2 <<- lambdas[2]
  if (identicalIncrements) lambda_2 <<- lambda_1
  Delta_1 <<- 2*pi/N/lambda_1
  Delta_2 <<- 2*pi/N/lambda_2
}

#===============================================================#
#               DISCRETIZATION FUNCTIONS                        #
#===============================================================#

# v_1,m, v_2,n
v_1 <- function(m) (m-N/2)*Delta_1
v_2 <- function(n) (n-N/2)*Delta_2

# k_1,p, k_2,q, p, q
k_1 <- function(p) (p-N/2)*lambda_1
k_2 <- function(...,type=c("normal","under","over")) {
  # based on type, we expect neither q or p
  type = head(type,1)
  dots = list(...)
  if (type == "normal") {
    q = dots$q
    return((q-N/2)*lambda_2)
  } else if (type == "under") {
    q = seq_len(N)-1
    k_2_q = k_2(q = q)
    p = dots$p
    x = p * NaN
    for (pp in seq_along(p)) {
      x[pp] = min(subset(k_2_q, exp(k_2_q)-exp(k_1(p[pp]+1)) >= K_last)
                ,last(k_2_q))
    }
    return(x)
  } else if (type == "over") {
    q = seq_len(N)-1
    k_2_q = k_2(q = q)
    p = dots$p
    x = p * NaN
    for (pp in seq_along(p)) {

      x[pp] = max(subset(k_2_q, exp(k_2_q)-exp(k_1(p[pp])) < K_last)
                ,first(k_2_q))
    }
    return(x)
  }
}
pp <- function(k_1_p) {
  as.numeric(round(k_1_p/lambda_1 + N/2,0))
}
qq <- function(k_2_q) {
  as.numeric(round(k_2_q/lambda_2 + N/2,0))
}

#===============================================================#
#               CHARACTERISTIC FUNCTIONS                        #
#===============================================================#

# phi_T(u_1, u_2)
charFunc_T <- function(u_1, u_2, type = modelNames$SV,
  underlying1, underlying2, corrs, volatility = NA, r, T_t) {
  # general params and GBM
  und1 = underlying1$getParams()
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  rho = corrs$getCorrelation(underlying1, underlying2)

  # zeta definition
  zeta = -1/2 * (
    (sigma_1^2*u_1^2 + sigma_2^2*u_2^2 + 2*rho*sigma_1*sigma_2*u_1*u_2) +
    i*(sigma_1^2*u_1 + sigma_2^2*u_2)
  );
  if (type == modelNames$GBM) {
    return(
      exp(
        zeta*T_t +
        i*u_1*(r - delta_1)*T_t + i*u_2*(r - delta_2)*T_t
      )
    )
  } else if (type == modelNames$SV) {
    # params SV
    vol = volatility$getParams()
    nu_0 = vol$nu_0
    sigma_nu = vol$sigma
    kappa = vol$kappa
    mu = vol$mu

    rho_1 = corrs$getCorrelation(underlying1, volatility)
    rho_2 = corrs$getCorrelation(underlying2, volatility)

    # gamma definition
    gamma = kappa - i*(rho_1*sigma_1*u_1 + rho_2*sigma_2*u_2)*sigma_nu;
    # theta definition
    theta = sqrt(gamma^2 - 2*sigma_nu^2*zeta);

    return(
      exp(
        ( 2*zeta*(1-exp(-theta*T_t)) /
          (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) ) * nu_0 +
        i*u_1*(r - delta_1)*T_t + i*u_2*(r - delta_2)*T_t -
        kappa*mu/sigma_nu^2 *
        (2 * log( (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) /(2*theta) ) +
        (theta - gamma)*T_t)
      )
    )
  } else {
    stop("Unknown type of characteristic function!")
  }
}

# phi(u_1, u_2)
charFunc <- function(u_1, u_2, type = modelNames$SV,
  underlying1, underlying2, corrs, volatility = NA, r, T_t) {
  s_1_0 = log(underlying1$getInitPrice())
  s_2_0 = log(underlying2$getInitPrice())

  return(
    exp(i*u_1*s_1_0 + i*u_2*s_2_0) *
    charFunc_T(u_1 = u_1, u_2 = u_2, type = type,
      underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, volatility = volatility, r = r, T_t = T_t)
  )
}

# phi_sv(u_1, u_2)
charFuncSV <- function(u_1, u_2,
  underlying1, underlying2, corrs, volatility, r, T_t) {
    charFunc(u_1 = u_1, u_2 = u_2, type = modelNames$SV,
      underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, volatility = volatility, r = r, T_t = T_t)
}

# phi_gbm(u_1, u_2)
charFuncGBM <- function(u_1, u_2,
  underlying1, underlying2, corrs, volatility = NA, r, T_t) {
  charFunc(u_1 = u_1, u_2 = u_2, type = modelNames$GBM,
    underlying1 = underlying1, underlying2 = underlying2,
    volatility = volatility,
    corrs = corrs, r = r, T_t = T_t)
}

# phi_T,sv(u_1, u_2)
charFuncSV_T <- function(u_1, u_2,
  underlying1, underlying2, corrs, volatility, r, T_t) {
    charFunc_T(u_1 = u_1, u_2 = u_2, type = modelNames$SV,
      underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, volatility = volatility, r = r, T_t = T_t)
}

# phi_T,gbm(u_1, u_2)
charFuncGBM_T <- function(u_1, u_2,
  underlying1, underlying2, corrs, volatility = NA, r, T_t) {
  charFunc_T(u_1 = u_1, u_2 = u_2, type = modelNames$GBM,
    underlying1 = underlying1, underlying2 = underlying2,
    volatility = volatility,
    corrs = corrs, r = r, T_t = T_t)
}

getCharFunction <- function(modelType, initConds = FALSE) {
  # function returns characteristic function depend on selected model
  # and whather with initial conditions (phi) or without (phi_T)
    if (modelType == modelNames$SV) {
        if (initConds) {charFuncSV}
        else {charFuncSV_T}
    } else if (modelType == modelNames$GBM) {
        if (initConds) {charFuncGBM}
        else {charFuncGBM_T}
    }
}

#===============================================================#
#               DEMPSTER-HONG APPROACH                          #
#===============================================================#

# chi_1(v_1, v_2)
charFuncChi1 <- function(v_1, v_2, charFunc = charFuncSV) {
  return(
    (
      charFunc(v_1-alpha_1*i, v_2-(alpha_2+1)*i,
        underlying1 = underlying1, underlying2 = underlying2,
        corrs = corrs, volatility = volatility, r = r, T_t = T_t) -
      charFunc(v_1-(alpha_1+1)*i, v_2-alpha_2*i,
        underlying1 = underlying1, underlying2 = underlying2,
        corrs = corrs, volatility = volatility, r = r, T_t = T_t)
    ) /
    ((alpha_1+i*v_1)*(alpha_2+i*v_2))
  )
}

# chi_2(v_1, v_2)
charFuncChi2 <- function(v_1, v_2, charFunc = charFuncSV) {
  return(
    (charFunc(v_1-alpha_1*i, v_2-alpha_2*i,
      underlying1 = underlying1, underlying2 = underlying2,
      corrs = corrs, volatility = volatility, r = r, T_t = T_t)) /
    ((alpha_1+i*v_1)*(alpha_2+i*v_2))
  )
}

## Fourier transformation
# chi(m, n)
fourierInputX <- function(charFuncChi, charFunc = charFuncSV) {
  m = n = seq_len(N)-1
  combs = expand.grid(n = n, m = m) # all possible combinations of n and m
  m = combs$m; n = combs$n
  X = (-1)^(m+n) * charFuncChi(v_1(m), v_2(n), charFunc)
  X = matrix(X, nrow = N, byrow = TRUE)

  return(X)
}

# Gamma(p, q) in DEMPSTER HONG
fourierOutputDempsterHongSelect <- function(p, q, Chi = 1,
  modelType = modelNames$SV) {
  if (Chi == 1) {
    if (modelType == modelNames$SV) {
      fourierOutputDempsterHong = fourierOutputDempsterHongChi1SV
    } else if (modelType == modelNames$GBM) {
      fourierOutputDempsterHong = fourierOutputDempsterHongChi1GBM
    }
  }
  else if (Chi == 2) {
    if (modelType == modelNames$SV) {
      fourierOutputDempsterHong = fourierOutputDempsterHongChi2SV
    } else if (modelType == modelNames$GBM) {
      fourierOutputDempsterHong = fourierOutputDempsterHongChi2GBM
    }
  }
  # q can be +-Inf, we need to replace it by max/min q
  q = pmax(pmin(q, N-1),0)
  p = pmax(pmin(p, N-1),0) # p also can't be out of bounds
  combs = expand.grid(q = q, p = p) # all possible combinations of p and q
  signCombs = (-1)^(combs$p + combs$q)
  signMatrix = matrix(signCombs, nrow = length(p), byrow = TRUE)

  return(signMatrix * fourierOutputDempsterHong[p+1,q+1])
}

# Pi(k_1,p, k_2,q)
componentPi <- function(k_1, k_2, Chi, modelType = modelNames$SV) {
  p = pp(k_1); q = qq(k_2)
  combs = expand.grid(k_2 = k_2, k_1 = k_1)
  expCombs = exp(-alpha_1*combs$k_1 -alpha_2*combs$k_2)
  expMatrix = matrix(expCombs, nrow = length(p), byrow = TRUE)
  return(
    expMatrix / (2*pi)^2 *
    Delta_2 * Delta_1 *
    fourierOutputDempsterHongSelect(p, q, Chi, modelType)
  )
}

# sum(Pi(k_1,p, k_2(p)) - Pi(k_1,p, k_2(p))), sum of differences
componentPi.SumOfDiffs <- function(Chi, type, modelType = modelNames$SV) {
  p = seq_len(N)-1
  Pi = sum(
    diag(componentPi(k_1 = k_1(p = p),
                     k_2 = k_2(p = p, type = type),
                     Chi = Chi, modelType = modelType)
    ) -
    diag(componentPi(k_1 = k_1(p = p + 1),
                     k_2 = k_2(p = p, type = type),
                     Chi = Chi, modelType = modelType)
    )
  )
  return(Pi)
}

# wrapper for Pi_1(k_1,p, k_2,q), underlined
componentPi1.Under <- function(modelType = modelNames$SV) {
  return(
    componentPi.SumOfDiffs(
      Chi = 1, type = "under",
      modelType = modelType)
  )
}

# wrapper for Pi_2(k_1,p, k_2,q), underlined
componentPi2.Under <- function(modelType = modelNames$SV) {
  return(
    componentPi.SumOfDiffs(
      Chi = 2, type = "under",
      modelType = modelType)
  )
}

# wrapper for Pi_1(k_1,p, k_2,q), overlined
componentPi1.Over <- function(modelType = modelNames$SV) {
  return(
    componentPi.SumOfDiffs(
      Chi = 1, type = "over",
      modelType = modelType)
  )
}

# wrapper for Pi_2(k_1,p, k_2,q), overlined
componentPi2.Over <- function(modelType = modelNames$SV) {
  return(
    componentPi.SumOfDiffs(
      Chi = 2, type = "over",
      modelType = modelType)
  )
}

# spread option price range for V(K)
getSpreadOptionPriceRange <- function(K, modelType = modelNames$SV, r, T_t) {
  priceRange = matrix(NA, nrow = length(K), ncol = 2,
    dimnames = list(K, c("lower","upper")))
  if (modelType == modelNames$SV) {
    priceRange[,1] = componentPi1.Under.SV -
                    K * componentPi2.Under.SV
    priceRange[,2] = componentPi1.Over.SV -
                    (K - epsilon) * componentPi2.Over.SV -
                    epsilon * componentPi2.Under.SV
  } else if (modelType == modelNames$GBM) {
    priceRange[,1] = componentPi1.Under.GBM -
                    K * componentPi2.Under.GBM
    priceRange[,2] = componentPi1.Over.GBM -
                    (K - epsilon) * componentPi2.Over.GBM -
                    epsilon * componentPi2.Under.GBM
  } else {
    stop("Undefined model type!")
  }

  return(Re(exp(-r*T_t)*priceRange))
}

#===============================================================#
#                    HURD-ZHOU APPROACH                         #
#===============================================================#

# P^hat(u)
payoffCoreFuncComplexGamma <- function(u_1, u_2) {
  return( gammaz(i*(u_1+u_2)-1) * gammaz(-i*u_2) / gammaz(i*u_1+1) )
}

# HurdZhou inverse fourier transform input
fourierInputHurdZhou <- function(charFunc_T = charFuncSV_T,
  underlying1, underlying2, corrs, volatility = NA, r, T_t) {
  if (alpha_2 <= 0) {
    stop("Second decay term (alpha_2) should be higher than zero!")
  }
  if (alpha_1+alpha_2 >= -1) {
    stop("Sum of decay terms should be lower than minus one!")
  }
  m = n = seq_len(N)-1
  combs = expand.grid(n = n, m = m) # all possible combinations of n and m
  m = combs$m; n = combs$n

  u_1 = v_1(m)+i*alpha_1
  u_2 = v_2(n)+i*alpha_2

  charFuncResult = charFunc_T(u_1 = u_1, u_2 = u_2,
    underlying1 = underlying1, underlying2 = underlying2, corrs = corrs,
    volatility = volatility, r = r, T_t = T_t)
  complexGammaPayoffCore = payoffCoreFuncComplexGamma(u_1, u_2)
  H = (-1)^(m+n) *
      charFuncResult *
      complexGammaPayoffCore
  H = matrix(H, nrow = N, byrow = TRUE)

  return(H)
}

# Spr (spread value) function from HurdZhou paper
SprHurdZhou <- function(K, fourierOutputHurdZhou, modelType = modelNames$SV,
  underlying1 = underlying1, underlying2 = underlying2) {

  if (K == 0) K = 1e-5 # small K instead zero

  k_1 = log(underlying1$getInitPrice()) - log(abs(K))
  k_2 = log(underlying2$getInitPrice()) - log(abs(K))
  p = pp(k_1); q = qq(k_2)
  if (any(!c(p,q) %in% (seq_len(N)-1))) {
    stop("p or q in Spr function is out of bound!")
  }

  return(K * Re(
    (-1)^(p + q) * exp(-r*T_t) *
    (Delta_1 * Delta_2)/(2*pi)^2 * # Delta_x*N/(2*pi) = 1/lambda_x
    exp(-alpha_1*k_1 -alpha_2*k_2) *
    fourierOutputHurdZhou[p+1,q+1]
  ))
}

getSpreadOptionPricesHurdZhou <- function(K, r, T_t, smallStrike = 1e-3) {
  # We are going to calculate value for K = 0 as average of value for
  # extremely small negative strike and extremely small positive strike
  # i.e. (V[0-] + V[0+])/ 2, where V[K] stands for value of option with given
  # strike K
  zeroCloseK = smallStrike # extremely small K instead of zero
  if (0 %in% K) K = c(K,c(-1,1)*zeroCloseK)

  SpreadOptionPrices = matrix(rep(K,2) * NA,ncol = 2,
    dimnames=list(K,c("Call","Put")))
  # Calculation of Calls (including negative strikes)
  t_begin = proc.time()[3]
  for (kk in seq_along(K)) {
    t_start = proc.time()[3]
    if (K[kk] == 0) next

    if (K[kk] > 0) {
      Strike = K[kk]
      underlyingX = underlying1
      underlyingY = underlying2
    } else if (K[kk] < 0) {
      # switch underlying assets
      Strike = -1*K[kk]
      underlyingX = underlying1.negStrike
      underlyingY = underlying2.negStrike
    }

    undX = underlyingX$getParams()
    undY = underlyingY$getParams()
    S_1_0 = undX$S_0; delta_1 = undX$delta
    S_2_0 = undY$S_0; delta_2 = undY$delta

    CallPutParityDiff = Strike*exp(-r*T_t) -
                        ( S_1_0*exp(-delta_1*T_t) - S_2_0*exp(-delta_2*T_t) )

    setLambdasAndDeltas(K = Strike, u_minimum = u_min,
      underlying1 = underlyingX, underlying2 = underlyingY)
    fourierInput = fourierInputHurdZhou(charFunc_T = charFunction,
      underlying1 = underlyingX, underlying2 = underlyingY, corrs = corrs,
      volatility = volatility, r = r, T_t = T_t)

    fourierOutput = fftw_c2c_2d(fourierInput,inverse = T)

    SpreadOptionPrices[kk,1] = SprHurdZhou(K = Strike,
      fourierOutputHurdZhou = fourierOutput, modelType = modelType,
      underlying1 = underlyingX, underlying2 = underlyingY)

    if (K[kk] < 0)
      SpreadOptionPrices[kk,1] = SpreadOptionPrices[kk,1] + CallPutParityDiff

    t_finish = proc.time()[3]
    progress(kk,length(K),t_finish-t_start,t_finish-t_begin)
  }
  # zero strike handling
  if (all((c(-1,1,0)*zeroCloseK) %in% K)) {
    SpreadOptionPrices["0",1] = mean(SpreadOptionPrices[paste(-zeroCloseK),1],
                                 SpreadOptionPrices[paste(zeroCloseK),1])
    K = K[!K %in% (c(-1,1)*zeroCloseK)]
    SpreadOptionPrices = SpreadOptionPrices[paste(K),,drop=F]
  }
  # Puts calculated using Call-put parity
  und1 = underlying1$getParams()
  und2 = underlying2$getParams()
  S_1_0 = und1$S_0; delta_1 = und1$delta
  S_2_0 = und2$S_0; delta_2 = und2$delta

  CallPutParityDiff = K*exp(-r*T_t) -
                      ( S_1_0*exp(-delta_1*T_t) - S_2_0*exp(-delta_2*T_t) )

  SpreadOptionPrices[,2] = SpreadOptionPrices[,1] + CallPutParityDiff

  return(SpreadOptionPrices)
}

#===============================================================#
#               MONTE CARLO SIMULATION                          #
#===============================================================#

optionPayoff <- function(S_T, K, call=TRUE) {
  return(
    pmax(ifelse(call,1,-1) * (S_T-K),0)
    )
}

optionPrices <- function(S_T, K, r, T_t, call=TRUE, se=TRUE) {
  opt = list()
  opt$prices = sapply(K,function(K)
    mean(exp(-r*T_t)*optionPayoff(S_T,K,call)))
  names(opt$prices) <- K
  if (se) {
    opt$se = sd(optionPayoff(S_T,1,call))/sqrt(length(S_T))
    return(opt)
  }
  return(opt$prices)
}

monteCarloSV <- function(underlying1, underlying2, volatility, corrs, r, T_t) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  vol = volatility$getParams()
  nu_0 = vol$nu_0
  sigma_nu = vol$sigma
  kappa = vol$kappa
  mu = vol$mu

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()), nrow = 3)

  sprds = lsims*NA
  t_begin = proc.time()[3]
  lapply((2*lsims-1), function(ss){
    t_start = proc.time()[3]
    V = matrix(rep(tsteps, 6)*NA, nrow = 6,
    dimnames = list(c("s_1", "s_2", "nu", "s_1_alt", "s_2_alt", "nu_alt")))
    dW = matrix(rnorm(sim_timesteps*3, sd = sqrt(dt)), nrow = 3)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2","dW_nu"))

    V[,1] = c(s_1_0,s_2_0,nu_0)
    for (ts in tsteps[-1]) {
      V[c("s_1","s_1_alt"), ts] = V[c("s_1","s_1_alt"), ts-1] +
        (r-delta_1-1/2*sigma_1^2*pmax(V[c("nu","nu_alt"), ts-1],0))*dt +
        c(1,-1)*sigma_1*sqrt(pmax(V[c("nu","nu_alt"), ts-1],0))*dW["dW_1",ts]
      V[c("s_2","s_2_alt"), ts] = V[c("s_2","s_2_alt"), ts-1] +
        (r-delta_2-1/2*sigma_2^2*pmax(V[c("nu","nu_alt"), ts-1],0))*dt +
        c(1,-1)*sigma_2*sqrt(pmax(V[c("nu","nu_alt"), ts-1],0))*dW["dW_2",ts]
      V[c("nu","nu_alt"), ts] = V[c("nu","nu_alt"), ts-1] +
        kappa*(mu-pmax(V[c("nu","nu_alt"), ts-1],0))*dt +
        c(1,-1)*sigma_nu*sqrt(pmax(V[c("nu","nu_alt"),ts-1],0))*dW["dW_nu",ts]
    }
    sprds[ss] <<- exp(V["s_1", sim_timesteps]) - exp(V["s_2", sim_timesteps])
    sprds[ss+1] <<- exp(V["s_1_alt", sim_timesteps]) -
                    exp(V["s_2_alt", sim_timesteps])

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  })
  return(sprds)
}

monteCarloSV2 <- function(underlying1, underlying2, volatility, corrs, r, T_t) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  vol = volatility$getParams()
  nu_0 = vol$nu_0
  sigma_nu = vol$sigma
  kappa = vol$kappa
  mu = vol$mu

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()), nrow = 3)

  sprds = lsims*NA
  t_begin = proc.time()[3]
  lapply((2*lsims-1), function(ss){
    t_start = proc.time()[3]
    V = matrix(rep(tsteps, 6)*NA, nrow = 6,
    dimnames = list(c("s_1", "s_2", "nu", "s_1_alt", "s_2_alt", "nu_alt")))
    dW = matrix(rnorm(sim_timesteps*3, sd = sqrt(dt)), nrow = 3)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2","dW_nu"))

    V[,1] = c(s_1_0,s_2_0,nu_0)
    for (ts in tsteps[-1]) {
      V["s_1", ts] = V["s_1", ts-1] +
        (r-delta_1-1/2*sigma_1^2*max(V["nu", ts-1],0))*dt +
        sigma_1*sqrt(max(V["nu", ts-1],0))*dW["dW_1",ts]
      V["s_2", ts] = V["s_2", ts-1] +
        (r-delta_2-1/2*sigma_2^2*max(V["nu", ts-1],0))*dt +
        sigma_2*sqrt(max(V["nu", ts-1],0))*dW["dW_2",ts]
      V["nu", ts] = V["nu", ts-1] +
        kappa*(mu-max(V["nu", ts-1],0))*dt +
        sigma_nu*sqrt(max(V["nu",ts-1],0))*dW["dW_nu",ts]
      V["s_1_alt", ts] = V["s_1_alt", ts-1] +
        (r-delta_1-1/2*sigma_1^2*max(V["nu_alt", ts-1],0))*dt +
        sigma_1*sqrt(max(V["nu_alt", ts-1],0))*dW["dW_1",ts]
      V["s_2_alt", ts] = V["s_2_alt", ts-1] +
        (r-delta_2-1/2*sigma_2^2*max(V["nu_alt", ts-1],0))*dt +
        sigma_2*sqrt(max(V["nu_alt", ts-1],0))*dW["dW_2",ts]
      V["nu_alt", ts] = V["nu_alt", ts-1] +
        kappa*(mu-max(V["nu_alt", ts-1],0))*dt +
        sigma_nu*sqrt(max(V["nu_alt",ts-1],0))*dW["dW_nu",ts]
    }
    sprds[ss] <<- exp(V["s_1", sim_timesteps]) - exp(V["s_2", sim_timesteps])
    sprds[ss+1] <<- exp(V["s_1_alt", sim_timesteps]) -
                    exp(V["s_2_alt", sim_timesteps])

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  })
  return(sprds)
}

monteCarloSV3 <- function(underlying1, underlying2, volatility, corrs, r, T_t) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  rmdelta_1 = r - und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  rmdelta_2 = r - und2$delta

  vol = volatility$getParams()
  nu_0 = vol$nu_0
  sigma_nu = vol$sigma
  kappa = vol$kappa
  mu = vol$mu

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()), nrow = 3)
  tcholW_cor = t(chol(W_cor))
  V_empty = matrix(rep(tsteps, 6)*NA, nrow = 6,
    dimnames = list(c("s_1", "s_2", "nu", "s_1_alt", "s_2_alt", "nu_alt")))
  V_empty[,1] = c(s_1_0,s_2_0,nu_0)

  sprds = lsims*NA
  t_begin = proc.time()[3]

  lapply((2*lsims-1), function(ss){
    t_start = proc.time()[3]
    V = V_empty

    dW = matrix(rnorm(sim_timesteps*3, sd = sqrt(dt)), nrow = 3)
    dW = tcholW_cor %*% dW
    dimnames(dW) = list(c("dW_1","dW_2","dW_nu"))

    for (ts in tsteps[-1]) {
      V_nu_p = max(V["nu", ts-1],0)
      V_nu_alt_p = max(V["nu_alt", ts-1],0)
      sV_nu_p = sqrt(V_nu_p)
      sV_nu_alt_p = sqrt(V_nu_alt_p)

      V["s_1", ts] = V["s_1", ts-1] + (rmdelta_1-1/2*sigma_1^2*V_nu_p)*dt +
        sigma_1*sV_nu_p*dW["dW_1",ts]
      V["s_2", ts] = V["s_2", ts-1] + (rmdelta_2-1/2*sigma_2^2*V_nu_p)*dt +
        sigma_2*sV_nu_p*dW["dW_2",ts]
      V["nu", ts] = V["nu", ts-1] + kappa*(mu-V_nu_p)*dt +
        sigma_nu*sV_nu_p*dW["dW_nu",ts]
      V["s_1_alt", ts] = V["s_1_alt", ts-1] +
        (rmdelta_1-1/2*sigma_1^2*V_nu_alt_p)*dt +
        sigma_1*sV_nu_alt_p*dW["dW_1",ts]
      V["s_2_alt", ts] = V["s_2_alt", ts-1] +
        (rmdelta_2-1/2*sigma_2^2*V_nu_alt_p)*dt +
        sigma_2*sV_nu_alt_p*dW["dW_2",ts]
      V["nu_alt", ts] = V["nu_alt", ts-1] +
        kappa*(mu-V_nu_alt_p)*dt +
        sigma_nu*sV_nu_alt_p*dW["dW_nu",ts]
    }
    sprds[ss] <<- exp(V["s_1", sim_timesteps]) - exp(V["s_2", sim_timesteps])
    sprds[ss+1] <<- exp(V["s_1_alt", sim_timesteps]) -
                    exp(V["s_2_alt", sim_timesteps])

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  })
  return(sprds)
}

monteCarloSV4 <- function(underlying1, underlying2, volatility, corrs, r, T_t) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  rmdelta_1 = r - und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  rmdelta_2 = r - und2$delta

  vol = volatility$getParams()
  nu_0 = vol$nu_0
  sigma_nu = vol$sigma
  kappa = vol$kappa
  mu = vol$mu

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()), nrow = 3)
  tcholW_cor = t(chol(W_cor))
  V_empty = matrix(rep(tsteps, 6)*NA, nrow = 6,
    dimnames = list(c("s_1", "s_2", "nu", "s_1_alt", "s_2_alt", "nu_alt")))
  V_empty[,1] = c(s_1_0,s_2_0,nu_0)

  s_1_empty = V_empty["s_1",]
  s_2_empty = V_empty["s_2",]
  nu_empty = V_empty["nu",]
  s_1_alt_empty = V_empty["s_1_alt",]
  s_2_alt_empty = V_empty["s_2_alt",]
  nu_alt_empty = V_empty["nu_alt",]

  sprds = lsims*NA
  t_begin = proc.time()[3]

  for (ss in (2*lsims-1)) {
    t_start = proc.time()[3]
    s_1 = s_1_empty; s_2 = s_2_empty
    nu = nu_empty; nu_alt = nu_alt_empty
    s_1_alt = s_1_alt_empty; s_2_alt = s_2_alt_empty

    dW = matrix(rnorm(sim_timesteps*3, sd = sqrt(dt)), nrow = 3)
    dW = tcholW_cor %*% dW
    dimnames(dW) = list(c("dW_1","dW_2","dW_nu"))

    for (ts in tsteps[-1]) {
      nu_p = max(nu[ts-1],0)
      nu_alt_p = max(nu_alt[ts-1],0)
      sqrt_nu_p = sqrt(nu_p)
      sqrt_nu_alt_p = sqrt(nu_alt_p)

      s_1[ts] = s_1[ts-1] + (rmdelta_1-1/2*sigma_1^2*nu_p)*dt +
        sigma_1*sqrt_nu_p*dW["dW_1",ts]
      s_2[ts] = s_2[ts-1] + (rmdelta_2-1/2*sigma_2^2*nu_p)*dt +
        sigma_2*sqrt_nu_p*dW["dW_2",ts]
      nu[ts] = nu[ts-1] + kappa*(mu-nu_p)*dt +
        sigma_nu*sqrt_nu_p*dW["dW_nu",ts]
      s_1_alt[ts] = s_1_alt[ts-1] +
        (rmdelta_1-1/2*sigma_1^2*nu_alt_p)*dt +
        sigma_1*sqrt_nu_alt_p*dW["dW_1",ts]
      s_2_alt[ts] = s_2_alt[ts-1] +
        (rmdelta_2-1/2*sigma_2^2*nu_alt_p)*dt +
        sigma_2*sqrt_nu_alt_p*dW["dW_2",ts]
      nu_alt[ts] = nu_alt[ts-1] +
        kappa*(mu-nu_alt_p)*dt +
        sigma_nu*sqrt_nu_alt_p*dW["dW_nu",ts]
    }
    sprds[ss] <- exp(s_1[sim_timesteps]) - exp(s_2[sim_timesteps])
    sprds[ss+1] <- exp(s_1_alt[sim_timesteps]) -
                    exp(s_2_alt[sim_timesteps])

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  }
  return(sprds)
}

monteCarloGBM <- function(underlying1, underlying2, corrs, r, T_t,
  sim_timesteps) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()[1:2,1:2]), nrow = 2)

  sprds = lsims*NA
  t_begin = proc.time()[3]
  for (ss in 2*lsims-1) {
    t_start = proc.time()[3]
    V = matrix(rep(tsteps, 4)*NA, nrow = 4,
    dimnames = list(c("s_1", "s_2", "s_1_alt", "s_2_alt")))
    dW = matrix(rnorm(sim_timesteps*2, sd = sqrt(dt)), nrow = 2)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2"))

    V[,1] = c(s_1_0,s_2_0)
    for (ts in tsteps[-1]) {
      V[c("s_1","s_1_alt"), ts] = V[c("s_1","s_1_alt"), ts-1] +
        (r-delta_1-1/2*sigma_1^2)*dt + c(1,-1)*sigma_1*dW["dW_1",ts]
      V[c("s_2","s_2_alt"), ts] = V[c("s_2","s_2_alt"), ts-1] +
        (r-delta_2-1/2*sigma_2^2)*dt + c(1,-1)*sigma_2*dW["dW_2",ts]
    }
    sprds[ss] = exp(V["s_1", sim_timesteps]) - exp(V["s_2", sim_timesteps])
    sprds[ss+1] = exp(V["s_1_alt", sim_timesteps]) -
                  exp(V["s_2_alt", sim_timesteps])

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  }
  return(sprds)
}

monteCarloGBM2 <- function(underlying1, underlying2, corrs, r, T_t,
  sim_timesteps) {
  dt = T_t/sim_timesteps
  sim_timesteps = sim_timesteps + 1
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  S_1_0 = und1$S_0
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  S_2_0 = und2$S_0
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()[1:2,1:2]), nrow = 2)

  sprds = lsims*NA
  t_begin = proc.time()[3]
  for (ss in 2*lsims-1) {
    t_start = proc.time()[3]
    V = matrix(rep(tsteps, 4)*NA, nrow = 4,
               dimnames = list(c("S_1", "S_2", "S_1_alt", "S_2_alt")))
    dW = matrix(rnorm(sim_timesteps*2, sd = sqrt(dt)), nrow = 2)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2"))

    V[,1] = c(S_1_0,S_2_0)
    for (ts in tsteps[-1]) {
      V[c("S_1","S_1_alt"), ts] = V[c("S_1","S_1_alt"), ts-1]*(1 +
        (r-delta_1)*dt + c(1,-1)*sigma_1*dW["dW_1",ts])
      V[c("S_2","S_2_alt"), ts] = V[c("S_2","S_2_alt"), ts-1]*(1 +
        (r-delta_2)*dt + c(1,-1)*sigma_2*dW["dW_2",ts])
    }
    sprds[ss] = V["S_1", sim_timesteps] - V["S_2", sim_timesteps]
    sprds[ss+1] = V["S_1_alt", sim_timesteps] - V["S_2_alt", sim_timesteps]

    t_finish = proc.time()[3]
    progress((ss+1)/2, n_sim, t_finish-t_start, t_finish-t_begin)
  }
  return(sprds)
}

monteCarloGBM3 <- function(underlying1, underlying2, corrs, r, T_t) {
  sim_timesteps = 1
  return(monteCarloGBM(underlying1, underlying2, corrs, r, T_t, sim_timesteps))
}

#USE THIS
monteCarloGBM4 <- function(underlying1, underlying2, corrs, r, T_t) {
  # optimized function, based on analytical formula
  # (terminal S_T without steps)
  lsims = seq_len(n_sim)

  und1 = underlying1$getParams()
  s_1_0 = log(und1$S_0)
  sigma_1 = und1$sigma
  delta_1 = und1$delta

  und2 = underlying2$getParams()
  s_2_0 = log(und2$S_0)
  sigma_2 = und2$sigma
  delta_2 = und2$delta

  W_cor = matrix(unlist(corrs$getCorrelationMatrix()[1:2,1:2]), nrow = 2)

  sprds = seq_len(n_sim*2)*NA
  V = matrix(rep(lsims, 4)*NA, nrow = 4,
            dimnames = list(c("s_1", "s_2", "s_1_alt", "s_2_alt")))
  dW = matrix(rnorm(n_sim*2, sd = sqrt(T_t)), nrow = 2)
  dW = t(chol(W_cor)) %*% dW
  dimnames(dW) = list(c("dW_1","dW_2"))
  V["s_1",] = s_1_0 + (r-delta_1-1/2*sigma_1^2)*T_t + sigma_1*dW["dW_1",]
  V["s_2",] = s_2_0 + (r-delta_2-1/2*sigma_2^2)*T_t + sigma_2*dW["dW_2",]
  V["s_1_alt",] = s_1_0 + (r-delta_1-1/2*sigma_1^2)*T_t - sigma_1*dW["dW_1",]
  V["s_2_alt",] = s_2_0 + (r-delta_2-1/2*sigma_2^2)*T_t - sigma_2*dW["dW_2",]

  sprds = c(
    exp(V["s_1",]) - exp(V["s_2",]),
    exp(V["s_1_alt",]) - exp(V["s_2_alt",])
  )

  return(sprds)
}

#===============================================================#
#                         TESTING                               #
#===============================================================#
K_LimKoef = 0.5
sigma_LimFrom = 0.01
sigma_LimTo = 0.9
S_0_LimFrom = 60
S_0_LimTo = 140
getPosStrikes = function(n, basis, LimKoef = K_LimKoef)
  return(runif(n, 1e-3, LimKoef) * basis)

getSigmas = function(n, LimFrom = sigma_LimFrom, LimTo = sigma_LimTo)
  return(runif(n,LimFrom,LimTo))

getInitPrices = function(n, LimFrom = S_0_LimFrom, LimTo = S_0_LimTo)
  return(runif(n,LimFrom,LimTo))
