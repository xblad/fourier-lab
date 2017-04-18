library(fftwtools)
library(pracma)
library(digest)
library(dplyr)
options(stringsAsFactors = FALSE)

#===============================================================#
#               PREPARATION FUNCTIONS                           #
#===============================================================#
# imaginary unit
i = complex(imaginary = 1)

modelNames = data.frame(
    SV = "SV",
    GBM = "GBM"
  )

create.underlyingAsset <- function(S_0, sigma, delta) {
  hash <- substring(digest(rnorm(10)), 0, 6)
  list(
    getHash = function() return(hash),

    setUnderlyingPrice = function(S_0) S_0 <<- S_0,
    setSigma = function(sigma) sigma <<- sigma,
    setDividendRate = function(delta) delta <<- delta,

    getUnderlyingPrice = function() return(S_0),
    getSigma = function() return(sigma),
    getDividendRate = function() return(delta),

    setParams = function(S_0, sigma, delta) {
      S_0 <<- S_0; sigma <<- sigma; delta <<- delta;
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

getLogPricesIncrements <- function(K, u_minimum = u_min) {
  if (K < 0) stop("Function is undefined for non-positive strikes!")
  S_0 = c(S_1_0, S_2_0)
  lambdas = NA * S_0
  for (m in seq_along(S_0)) {
    # only positive j make sense, so that we start by N/2,
    # than first u_test will be zero
    for (j in (N/2):(N-1)) {
      u_test = pi*(j-N/2)/log(S_0[m]/K)
      if (u_test > u_minimum) {
        lambdas[m] = pi/u_test
        break
      }
    }
  }

  if (any(is.na(lambdas)))
    stop("Impossible to find value for defined minimum bound!")

  return(lambdas)
}

setLambdasAndDeltas <- function(K, u_minimum = u_min, identicalIncrements = FALSE) {
  lambdas = getLogPricesIncrements(K = K, u_min)
  lambda_1 <<- lambdas[1]
  lambda_2 <<- lambdas[2]
  if (identicalIncrements) lambda_2 <<- lambda_1
  Delta_1 <<- 2*pi/N/lambda_1
  Delta_2 <<- 2*pi/N/lambda_2
}

#===============================================================#
#               FIRST FOURIER APPROACH                          #
#===============================================================#


## Characteristic functions
# phi_T(u_1, u_2)
charFunc_T <- function(u_1, u_2, type = modelNames$SV) {
  type = head(type,1)

  # zeta definition
  zeta = -1/2 * (
    (sigma_1^2*u_1^2 + sigma_2^2*u_2^2 + 2*ro*sigma_1*sigma_2*u_1*u_2) +
    i*(sigma_1^2*u_1 + sigma_2^2*u_2)
  );
  if (type == modelNames$GBM) {
    return(
      exp(
        zeta*T_t +
        # NOTE: strange... somewhere with i, somewhere without i...
        i*u_1*(r - delta_1)*T_t + i*u_2*(r - delta_2)*T_t
      )
    )
  } else if (type == modelNames$SV) {
    # gamma definition
    gamma = kappa - i*(ro_1*sigma_1*u_1 + ro_2*sigma_2*u_2)*sigma_nu;
    # theta definition
    theta = sqrt(gamma^2 - 2*sigma_nu^2*zeta);

    return(
      exp(
        ( 2*zeta*(1-exp(-theta*T_t)) / (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) ) * nu_0 +
        # NOTE: strange... somewhere with i, somewhere without i...
        i*u_1*(r - delta_1)*T_t + i*u_2*(r - delta_2)*T_t -
        kappa*mu/sigma_nu^2 * (2 * log( (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) / (2*theta) ) + (theta - gamma)*T_t)
      )
    )
  } else {
    stop("Unknown type of characteristic function!")
  }
}

# phi(u_1, u_2)
charFunc <- function(u_1, u_2, type = modelNames$SV) {
  return(
    exp(i*u_1*s_1_0 + i*u_2*s_2_0) * charFunc_T(u_1, u_2, type)
  )
}

# charFunc_V2 <- function(u_1, u_2, type = c("SV","GBM")) {
#   type = head(type,1)
#
#   s_0 = t(c(s_1_0, s_2_0)) # row vector
#   delta = t(c(delta_1, delta_2)) # row vector
#   u = matrix(c(u_1,u_2),ncol=2) # row vector
#
#   # initial condition
#   # NOTE: probably zero condition is irrelevant for phi_T (only for phi)???
#   init_cond = i * u %*% t(s_0)
#
#   # zeta definition
#   zeta = -1/2 * ( diag(u %*% S_cov %*% t(u)) + i*(u %*% diag(S_cov)) );
#   if (type == "GBM") {
#     return(
#       exp(
#         # zero condition is irrelevant for phi_T (only for phi)
#         #init_cond +
#         zeta*T_t +
#         u_1*(r - delta_1)*T_t + u_2*(r - delta_2)*T_t
#       )
#     )
#   } else if (type == "SV") {
#     # gamma definition
#     gamma = kappa - i * t(sigma_nu_covs %*% t(u));
#     # theta definition
#     theta = sqrt(gamma^2 - 2*sigma_nu^2*zeta);
#
#     return(
#       exp(
#         # NOTE: probably zero condition is irrelevant for phi_T (only for phi)???
#         init_cond +
#         ( 2*zeta*(1-exp(-theta*T_t)) / (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) ) * nu_0 +
#         # NOTE: strange... somewhere with i, somewhere without i...
#         i * u %*% t(r - delta) * T_t -
#         kappa*mu/sigma_nu^2 * (2 * log( (2*theta - (theta-gamma)*(1-exp(-theta*T_t))) / (2*theta) ) + (theta - gamma)*T_t)
#       )
#     )
#   } else {
#     stop("Unknown type of characteristic function!")
#   }
# }

# phi_sv(u_1, u_2)
charFuncSV <- function(u_1, u_2) charFunc(u_1, u_2, type = modelNames$SV)

# phi_gbm(u_1, u_2)
charFuncGBM <- function(u_1, u_2) charFunc(u_1, u_2, type = modelNames$GBM)

# phi_T,sv(u_1, u_2)
charFuncSV_T <- function(u_1, u_2) charFunc_T(u_1, u_2, type = modelNames$SV)

# phi_T,gbm(u_1, u_2)
charFuncGBM_T <- function(u_1, u_2) charFunc_T(u_1, u_2, type = modelNames$GBM)

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

# chi_1(v_1, v_2)
charFuncChi1 <- function(v_1, v_2, charFunc = charFuncSV) {
  return(
    (
      charFunc(v_1-alpha_1*i, v_2-(alpha_2+1)*i) -
      charFunc(v_1-(alpha_1+1)*i, v_2-alpha_2*i)
    ) /
    ((alpha_1+i*v_1)*(alpha_2+i*v_2))
  )
}

# chi_2(v_1, v_2)
charFuncChi2 <- function(v_1, v_2, charFunc = charFuncSV) {
  return(
    (charFunc(v_1-alpha_1*i, v_2-alpha_2*i)) /
    ((alpha_1+i*v_1)*(alpha_2+i*v_2))
  )
}

## Fourier transformation
# v_1,m, v_2,n
v_1 <- function(m) (m-N/2)*Delta_1
v_2 <- function(n) (n-N/2)*Delta_2

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
fourierOutputDempsterHongSelect <- function(p, q, Chi = 1, modelType = modelNames$SV) {
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
getSpreadOptionPriceRange <- function(K, modelType = modelNames$SV) {
  priceRange = matrix(NA, nrow = length(K), ncol = 2, dimnames = list(K, c("lower","upper")))
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

  return(abs(exp(-r*T_t)*priceRange))
}

#===============================================================#
#               SECOND FOURIER APPROACH                         #
#===============================================================#

# P^hat(u)
payoffCoreFuncComplexGamma <- function(u_1, u_2, negativeStrike = FALSE) {
  if (!negativeStrike) {
    return( gammaz(i*(u_1+u_2)-1) * gammaz(-i*u_2) / gammaz(i*u_1+1) )
  } else {
   return( -1 * payoffCoreFuncComplexGamma(u_2, u_1) )
  }
  #return(1/((1-i*u_1)*u_2*(u_1+u_2+i)))
}

# HurdZhou inverse fourier transform input
fourierInputHurdZhou <- function(charFunc_T = charFuncSV_T, negativeStrike = FALSE) {
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
  if (negativeStrike) {
    u_1 = v_1(m)+i*alpha_2; u_2 = v_2(n)+i*alpha_1
  }
  charFuncResult = charFunc_T(u_1, u_2)
  complexGammaPayoffCore = payoffCoreFuncComplexGamma(u_1, u_2, negativeStrike)
  # NOTE: with or without zero condition??? Now without
  H = (-1)^(m+n) *
      # NOTE: in Hurd Zhou exp() function in characteristic function
      # for SV is missed!
      charFuncResult *
      complexGammaPayoffCore
  H = matrix(H, nrow = N, byrow = TRUE)

  return(H)
}

# HurdZhou inverse fourier transform output (must be calculated)
fourierOutputHurdZhouSelect <- function(modelType = modelNames$SV, K = 1) {
  if (K > 0) {
    if (modelType == modelNames$SV) {
      fourierOutputHurdZhou = fourierOutputHurdZhouSV
    } else if (modelType == modelNames$GBM) {
      fourierOutputHurdZhou = fourierOutputHurdZhouGBM
    }
  } else if (K < 0) {
    if (modelType == modelNames$SV) {
      fourierOutputHurdZhou = fourierOutputHurdZhouSV.negStrike
    } else if (modelType == modelNames$GBM) {
      fourierOutputHurdZhou = fourierOutputHurdZhouGBM.negStrike
    }
  }

  return(fourierOutputHurdZhou)
}

# Spr (spread value) function from HurdZhou paper
SprHurdZhou <- function(K, modelType = modelNames$SV) {

  if (K == 0) K = 1e-5 # small K instead zero

  k_1 = s_1_0 - log(abs(K))
  k_2 = s_2_0 - log(abs(K))
  p = pp(k_1); q = qq(k_2)
  if (any(!c(p,q) %in% (seq_len(N)-1))) {
    stop("p or q in Spr function is out of bound!")
  }

  return(
    (-1)^(p + q) * exp(-r*T_t) *
    (Delta_1 * Delta_2)/(2*pi)^2 * # Delta_x*N/(2*pi) = 1/lambda_x
    exp(-alpha_1*k_1 -alpha_2*k_2) *
    fourierOutputHurdZhouSelect(modelType, K)[p+1,q+1]
  )
}

# Spr2 (spread value) function from HurdZhou paper (test for interpolation)
SprHurdZhou2 <- function(p, q, modelType = modelNames$SV) {

  if (any(!c(p,q) %in% (seq_len(N)-1))) {
    stop("p or q in Spr function is out of bound!")
  }
  k_1 = k_1(p=p)
  k_2 = k_2(q=q)

  return(
    (-1)^(p + q) * exp(-r*T_t) *
    (Delta_1 * Delta_2)/(2*pi)^2 * # Delta_x*N/(2*pi) = 1/lambda_x
    exp(-alpha_1*k_1 -alpha_2*k_2) *
    fourierOutputHurdZhouSelect(modelType)[p+1,q+1]
  )
}

#===============================================================#
#               MONTE CARLO SIMULATION                          #
#===============================================================#

monteCarloSV <- function() {
  dt = T_t/sim_timesteps
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  W_cor = matrix(
    c(1,    ro,   ro_1,
      ro,   1,    ro_2,
      ro_1, ro_2,   1), nrow = 3)

  sims = list()
  sprds = lsims*NA
  for (ss in 2*lsims-1) {
    V = matrix(rep(tsteps, 6)*NA, nrow = 6,
    dimnames = list(c("s_1", "s_2", "nu", "s_1_alt", "s_2_alt", "nu_alt")))
    dW = matrix(rnorm(sim_timesteps*3, sd = sqrt(dt)), nrow = 3)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2","dW_nu"))

    V[,1] = c(s_1_0,s_2_0,nu_0)
    for (ts in tsteps[-1]) {
      V[c("s_1","s_1_alt"), ts] = V[c("s_1","s_1_alt"), ts-1] + (r-delta_1-1/2*sigma_1^2*V[c("nu","nu_alt"),ts-1])*dt + c(1,-1)*sigma_1*sqrt(V[c("nu","nu_alt"),ts-1])*dW["dW_1",ts]
      V[c("s_2","s_2_alt"), ts] = V[c("s_2","s_2_alt"), ts-1] + (r-delta_2-1/2*sigma_2^2*V[c("nu","nu_alt"),ts-1])*dt + c(1,-1)*sigma_2*sqrt(V[c("nu","nu_alt"),ts-1])*dW["dW_2",ts]
      V[c("nu","nu_alt"), ts] = V[c("nu","nu_alt"), ts-1] + kappa*(mu-pmax(V[c("nu","nu_alt"), ts-1],0))*dt + c(1,-1)*sigma_nu*sqrt(pmax(V[c("nu","nu_alt"),ts-1],0))*dW["dW_nu",ts]
      V[c("nu","nu_alt"), ts] = pmax(V[c("nu","nu_alt"), ts],0)
    }
    sims[[(ss+1)/2]] = V
    sprds[ss] = exp(V["s_2", sim_timesteps]) - exp(V["s_1", sim_timesteps])
    sprds[ss+1] = exp(V["s_2_alt", sim_timesteps]) - exp(V["s_1_alt", sim_timesteps])
  }
  return(sprds)
}

monteCarloGBM <- function() {
  dt = T_t/sim_timesteps
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  W_cor = matrix(
    c(1,   ro,
      ro,   1), nrow = 2)

  sims = list()
  sprds = lsims*NA
  for (ss in 2*lsims-1) {
    V = matrix(rep(tsteps, 4)*NA, nrow = 4,
    dimnames = list(c("s_1", "s_2", "s_1_alt", "s_2_alt")))
    dW = matrix(rnorm(sim_timesteps*2, sd = sqrt(dt)), nrow = 2)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2"))

    V[,1] = c(s_1_0,s_2_0)
    for (ts in tsteps[-1]) {
      V[c("s_1","s_1_alt"), ts] = V[c("s_1","s_1_alt"), ts-1] + (r-delta_1-1/2*sigma_1^2)*dt + c(1,-1)*sigma_1*dW["dW_1",ts]
      V[c("s_2","s_2_alt"), ts] = V[c("s_2","s_2_alt"), ts-1] + (r-delta_2-1/2*sigma_2^2)*dt + c(1,-1)*sigma_2*dW["dW_2",ts]
    }
    sims[[(ss+1)/2]] = V
    sprds[ss] = exp(V["s_2", sim_timesteps]) - exp(V["s_1", sim_timesteps])
    sprds[ss+1] = exp(V["s_2_alt", sim_timesteps]) - exp(V["s_1_alt", sim_timesteps])
  }
  return(sprds)
}

monteCarloGBM2 <- function() {
  dt = T_t/sim_timesteps
  tsteps = seq_len(sim_timesteps)
  lsims = seq_len(n_sim)

  W_cor = matrix(
    c(1,   ro,
      ro,   1), nrow = 2)

  sims = list()
  sprds = lsims*NA
  for (ss in 2*lsims-1) {
    V = matrix(rep(tsteps, 4)*NA, nrow = 4,
               dimnames = list(c("S_1", "S_2", "S_1_alt", "S_2_alt")))
    dW = matrix(rnorm(sim_timesteps*2, sd = sqrt(dt)), nrow = 2)
    dW = t(chol(W_cor)) %*% dW
    dimnames(dW) = list(c("dW_1","dW_2"))

    V[,1] = c(S_1_0,S_2_0)
    for (ts in tsteps[-1]) {
      V[c("S_1","S_1_alt"), ts] = V[c("S_1","S_1_alt"), ts-1]*(1 + (r-delta_1)*dt + c(1,-1)*sigma_1*dW["dW_1",ts])
      V[c("S_2","S_2_alt"), ts] = V[c("S_2","S_2_alt"), ts-1]*(1 + (r-delta_2)*dt + c(1,-1)*sigma_2*dW["dW_2",ts])
    }
    sims[[(ss+1)/2]] = V
    sprds[ss] = V["S_2", sim_timesteps] - V["S_1", sim_timesteps]
    sprds[ss+1] = V["S_2_alt", sim_timesteps] - V["S_1_alt", sim_timesteps]
  }
  return(sprds)
}
