library(Matrix)
library(glmnet)
library(foreach)
library(doParallel)

estimate3 <- function(Data, I, statedim, delta, powers) {
  
  # Define functions under the names given in the paper
  Lambda <- function(u) 1 / (1 + exp(-u))
  dLambda <- function(u) exp(-u) * (Lambda(u)^2)
  lambda <- function(u) dLambda(u) / (Lambda(u) * (1 - Lambda(u)))
  H <- function(p) -log(1 - p)
  H_p <- function(p) 1 / (1 - p)
  
  rule_of_thumb <- 0 # Set this to 1 to use rule-of-thumb penalty parameters
  # instead of cross-validation
  
  # Rename columns in Data to match up with the paper
  y <- 1 - Data[, 1]        # Y_2 in the document (y(i) = 1 if continue with engine and y(i) = 0 if replace)
  x <- Data[, 2:(1 + statedim)]
  y_n <- 1 - Data[, (2 + statedim)]
  x_n <- Data[, (3 + statedim):ncol(Data)]
  N <- nrow(Data)
  
  # Generate basis functions corresponding to the choice given in the input 'powers'
  if (powers == 1) {
    X <- x
    X_n <- x_n
  }
  if (powers == 1.5) {
    X <- cbind(x, x^2)
    X_n <- cbind(x_n, x_n^2)
  }
  if (powers >= 2) {
    res <- fullinteracts(x, powers)
    powers <- res$inters
    X <- x2fx(x, powers)
    X_n <- x2fx(x_n, powers)
  }
  
  rule_of_thumb_penalty <- 2 * 0.5 * qnorm(1 - 0.1 / (2 * ncol(X)))
  
  # We will save nonparametric regression functions here
  g_1_hat <- array(NA, dim = c(N, ncol(I), ncol(I)))
  Hg_1_hat <- matrix(NA, nrow = N, ncol = ncol(I))
  g_1_tilde <- matrix(NA, nrow = N, ncol = ncol(I))
  g_1_sub <- matrix(NA, nrow = N, ncol = 1)
  g_3_hat <- matrix(NA, nrow = N, ncol = 1)
  
  for (i in 1:ncol(I)) {
    # Estimate gamma_1 and H() of gamma_1 with two-way sample splitting
    for (j in 1:ncol(I)) {
      comp <- as.logical(1 - I[, i] - I[, j])
      if (rule_of_thumb) {
        lambda <- rule_of_thumb_penalty / sum(comp, 1)
        cvfit <- cv.glmnet(as.matrix(X[comp, ]), y[comp], nfolds = 10)
      } else {
        # Perform lasso logistic regression with cross-validation
        cvfit <- cv.glmnet(as.matrix(X[comp, ]), y[comp], family = "binomial", nfolds = 3)
      }
      
      lambdas <- cvfit$lambda
      # # Access the deviance values for each lambda
      # deviance <- cvfit$cvm
      # # Access the index of the lambda with minimum deviance
      # idx <- which.min(deviance)      # # Access the lambda values
      # lambdas <- cvfit$lambda
      idx <- cvfit$index[1]
      # Access the coefficients for the lambda with minimum deviance
      params <- coef(cvfit, s = lambdas[idx])
      coef <- c(params[1], params[-1])
      # Compute predicted log-odds
      log_odds <- coef[1] + X_n %*% coef[-1]
      # Convert log-odds to probabilities using the logistic function
      g_1_hat[,i,j] <- 1 / (1 + exp(-log_odds))
      # Apply lower and upper bounds to predicted probabilities
      g_1_hat[,i,j] <- pmax(pmin(g_1_hat[,i,j], 0.9999), 0.0001)
      
      Hg_1_hat[I[, j], i] <- H(g_1_hat[I[, j], i, j])
    }
    
    comp <- as.logical(rep(1, N) - I[, i])
    
    # Estimate gamma_1 and H() of gamma_1 with one-way sample splitting
    if (rule_of_thumb) {
      cvfit <- glmnet::glmnet(X[comp, ], y[comp], family = "binomial", alpha = 1, lambda = rule_of_thumb_penalty / sum(comp))
    } else {
      cvfit <- cv.glmnet(as.matrix(X[comp, ]), y[comp], family = "binomial", nfolds = 3)
    }
    # # Access the lambda values
    lambdas <- cvfit$lambda
    # Access the index of the lambda with minimum deviance
    # idx <- which.min(deviance)
    idx <- cvfit$index[1]
    # Access the coefficients for the lambda with minimum deviance
    params <- coef(cvfit, s = lambdas[idx])
    coef <- c(params[1], params[-1])
    # Compute predicted log-odds
    log_odds <- coef[1] + X_n %*% coef[-1]
    # Convert log-odds to probabilities using the logistic function
    g_1_tilde[,i] <- 1 / (1 + exp(-log_odds))
    # Apply lower and upper bounds to predicted probabilities
    g_1_tilde[,i] <- pmax(pmin(g_1_tilde[,i], 0.9999), 0.0001)
    g_1_sub[I[, i],] <- g_1_tilde[I[, i], i]
    # Estimate gamma_3 with one-way sample splitting
    g_3_hat[I[, i], ] <- mean(H(g_1_tilde[which(comp), i]) * (1 - y[comp])) / mean(1 - y[comp])
  }
  
  g_2_hat <- rep(NA, N)
  for (i in 1:ncol(I)) {
    comp <- as.logical(rep(1, N) - I[, i])
    if (rule_of_thumb) {
      lambda <- rule_of_thumb_penalty / sum(comp) 
      cvfit <- cv.glmnet(X_sub, y_sub, alpha = 1, lambda = lambda, nfolds = 10)
    } else {
      cvfit <- cv.glmnet(x = y[comp] * X[comp, ], y = y[comp] * Hg_1_hat[comp, i], alpha = 1, family = "gaussian", type.measure = "mse", nfolds = 3)
    }
    lambdas <- cvfit$lambda
    # Access the index of the lambda with minimum deviance
    # idx <- which.min(deviance)
    idx <- cvfit$index[1]
    # Access the coefficients for the lambda with minimum deviance
    params <- coef(cvfit, s = lambdas[idx])
    coef <- c(params[1], params[-1])
    # coef <- coef(fit, s = fit$lambcvfit$indexda[idx])
    g_2_hat[I[, i]] <- cbind(1, X[I[, i], ]) %*% coef
  }
  
  Vdiff <- delta * (g_2_hat - g_3_hat)
  d <- cbind(sqrt(x[, 1]), -1)
  
  estimate <- matrix(NA, nrow = 2, ncol = ncol(I))
  v_hat <- rep(NA, N)
  
  # Loop through each column of I
  for (i in 1:ncol(I)) {
    # Obtain the logical vector comp
    comp <- as.logical(1 - I[, i])
    # Create a function 'mom' using the formula provided
    mom <- function(theta) {
      d[comp, ] * t(matrix(rep(lambda((d[comp, ]%*%theta) + Vdiff[comp, ]) * (y[comp] - Lambda((d[comp, ]%*%theta) + Vdiff[comp, ])), each = 2), nrow = 2))
    }
    
    # Create a function 'G' using the formula provided
    G <- function(theta) {
      # mean_mom <- colMeans(mom(theta))
      # sum(mean_mom^2)
      sum(colMeans(mom(theta))^2)
    }
    # Estimate the parameters using 'optim'
    # estimate[, i] <- optim(c(0, 0), G)$par
    estimate[, i] <- optim(par = c(0, 0), fn = G, method = "BFGS")$par
    # Update v_hat using the estimated parameters
    v_hat[I[, i]] <- (d[I[, i], ] %*% estimate[, i]) + Vdiff[I[, i], ]
  }

  Lambda_hat <- Lambda(v_hat)
  dLambda_hat <- dLambda(v_hat)
  lambda_hat <- lambda(v_hat)
  lamb_part <- (lambda_hat * dLambda_hat * y) / Lambda_hat
  
  alpha2 <- -delta * d * lamb_part
  phi2 <- alpha2 * as.vector(H(g_1_sub) - g_2_hat)
    # t(matrix(rep((H(g_1_sub) - g_2_hat), each = 2), nrow = 2))
  
  A <- matrix(NA, nrow = N, ncol = 2)
  P1 <- rep(NA, N)
  
  # phi3 <- matrix(NA, nrow = N, ncol = 1)
  b_rho3 <- matrix(NA, nrow = N, ncol = 1)
  
  for (i in 1:ncol(I)) {
    comp <- as.logical(1 - I[, i])
    if (rule_of_thumb) {
      lambda <- rule_of_thumb_penalty / sum(comp) 
      cvfit <- cv.glmnet(X_sub, y_sub, alpha = 1, lambda = lambda, nfolds = 10)
    } else {
      cvfit <- cv.glmnet(x = X_n[comp, ], y = 1 - y[comp], alpha = 1, family = "gaussian", type.measure = "mse", nfolds = 3)
    }
    # idx <- cvfit$lambda.min
    lambdas <- cvfit$lambda
    # Access the index of the lambda with minimum deviance
    # idx <- which.min(deviance)
    idx <- cvfit$index[1]
    # Access the coefficients for the lambda with minimum deviance
    params <- coef(cvfit, s = lambdas[idx])
    coef <- c(params[1], params[-1])
    # coef <- coef(fit, s = fit$lambda[idx])
    b_rho3[I[, i]] <- cbind(1, X_n[I[, i], ]) %*% coef
    A[I[, i], ] <- t(replicate(sum(I[, i]), colMeans(d[comp, ] * lambda_hat[comp] * dLambda_hat[comp])))
    P1[I[, i]] <- mean(1 - y[comp])
  }
  
  alpha3 <- A * as.vector(b_rho3) * as.vector(H_p(g_1_sub)) / P1
  phi3 <- alpha3 * as.vector(y_n - g_1_sub) + (A * (1 - y) / P1) * as.vector((H(g_1_sub) - g_3_hat))
  
  zeta <- matrix(NA, nrow = N, ncol = 2)
  
  # phi1 <- matrix(NA, nrow = N, ncol = 1)
  for (j in 1:2) {
    for (i in 1:ncol(I)) {
      comp <- as.logical(1 - I[, i])
      if (rule_of_thumb) {
        lambda <- rule_of_thumb_penalty / sum(comp) 
        cvfit <- cv.glmnet(X_sub, y_sub, alpha = 1, lambda = lambda, nfolds = 10)
      } else {
        cvfit <- cv.glmnet(x = X_n[comp, ], y = -delta * d[comp, j] * lamb_part[comp], alpha = 1, family = "gaussian", type.measure = "mse", nfolds = 3)
      }
      lambdas <- cvfit$lambda
      # Access the index of the lambda with minimum deviance
      # idx <- which.min(deviance)
      idx <- cvfit$index[1]
      # Access the coefficients for the lambda with minimum deviance
      params <- coef(cvfit, s = lambdas[idx])
      coef <- c(params[1], params[-1])
      # coef <- coef(fit, s = fit$lambda[idx])
      zeta[I[, i], j] <- as.vector(cbind(1, X_n[I[, i], ]) %*% coef)
    }
  }
  
  alpha1 <- zeta * as.vector(H_p(g_1_sub))
  phi1 <- alpha1 * as.vector(y_n - g_1_sub)
  
  # Uncorrected estimates
  # mom <- function(theta) d * matrix(rep(lambda((d %*% theta) + Vdiff) * (y - Lambda((d %*% theta) + Vdiff)), each = 2), ncol = 2)
  # G <- function(theta) sum(colMeans(mom(theta))^2)
  # est1 <- optim(c(0, 0), G)$par
  # Define the functions mom and G
  mom <- function(theta) {
    # lambda <- lambda((d %*% theta) + Vdiff)
    # y_minus_lambda <- y - Lambda((d %*% theta) + Vdiff)
    # d * matrix(rep(lambda * y_minus_lambda, each = 2), ncol = 2)
    d * t(matrix(rep(lambda((d%*%theta) + Vdiff) * (y - Lambda((d%*%theta) + Vdiff)), each = 2), nrow = 2))
  }
  G <- function(theta) {
    sum(colMeans(mom(theta))^2)
  }
  # Perform optimization using optim
  est1 <- optim(par = c(0, 0), fn = G, method = "BFGS")$par
  # est1
  
  # Estimates with correction term phi_2
  # mom_2 <- function(theta) mom(theta) + phi2
  # G <- function(theta) sum(colMeans(mom_2(theta))^2)
  mom_2 <- function(theta) {
    mom(theta)+phi2
  }
  G <- function(theta) {
    sum(colMeans(mom_2(theta))^2)
  }
  est2 <- optim(par = c(0, 0), fn = G, method = "BFGS")$par
  
  # Estimates with correction terms phi_2 and phi_3
  # mom_23 <- function(theta) mom(theta) + phi2 + phi3
  # G <- function(theta) sum(colMeans(mom_23(theta))^2)
  mom_23 <- function(theta) {
    mom(theta) + phi2 + phi3
  }
  G <- function(theta) {
    sum(colMeans(mom_23(theta))^2)
  }
  est3 <- optim(par = c(0, 0), fn = G, method = "BFGS")$par
  
  # Fully corrected estimates with phi_1, phi_2, and phi_3
  # mom_123 <- function(theta) mom(theta) + phi1 + phi2 + phi3
  # G <- function(theta) sum(colMeans(mom_123(theta))^2)
  mom_123 <- function(theta) {
    mom(theta) + phi1 + phi2 + phi3
  }
  G <- function(theta) {
    sum(colMeans(mom_123(theta))^2)
  }
  est4 <- optim(par = c(0, 0), fn = G, method = "BFGS")$par
  
  # Calculate standard errors using the standard formula with fully adjusted moments
  diff <- 0.00001  # First difference for numerical derivative calculation
  
  dG <- matrix(NA, nrow = 2, ncol = 2)
  mM <- function(theta) colMeans(mom(theta))
  dG[1, ] <- (mM(est1 + c(diff, 0)) - mM(est1)) / diff
  dG[2, ] <- (mM(est1 + c(0, diff)) - mM(est1)) / diff
  Om <- cov(mom(est1))
  covar <- solve(dG) %*% Om %*% solve(t(dG)) / N
  sd1 <- sqrt(diag(covar))
  
  dG <- matrix(NA, nrow = 2, ncol = 2)
  mM <- function(theta) colMeans(mom_2(theta))
  dG[1, ] <- (mM(est2 + c(diff, 0)) - mM(est2)) / diff
  dG[2, ] <- (mM(est2 + c(0, diff)) - mM(est2)) / diff
  Om <- cov(mom_123(est2))
  covar <- solve(dG) %*% Om %*% solve(t(dG)) / N
  sd2 <- sqrt(diag(covar))
  
  dG <- matrix(NA, nrow = 2, ncol = 2)
  mM <- function(theta) colMeans(mom_23(theta))
  dG[1, ] <- (mM(est3 + c(diff, 0)) - mM(est3)) / diff
  dG[2, ] <- (mM(est3 + c(0, diff)) - mM(est3)) / diff
  Om <- cov(mom_123(est3))
  covar <- solve(dG) %*% Om %*% solve(t(dG)) / N
  sd3 <- sqrt(diag(covar))
  
  dG <- matrix(NA, nrow = 2, ncol = 2)
  mM <- function(theta) colMeans(mom_123(theta))
  dG[1, ] <- (mM(est4 + c(diff, 0)) - mM(est4)) / diff
  dG[2, ] <- (mM(est4 + c(0, diff)) - mM(est4)) / diff
  Om <- cov(mom_123(est4))
  covar <- solve(dG) %*% Om %*% solve(t(dG)) / N
  sd4 <- sqrt(diag(covar))
  
  # Return the results
  estimates <- cbind(est1, est2, est3, est4)
  standard_errors <- cbind(sd1, sd2, sd3, sd4)
  
  return(list(estimates = estimates, standard_errors = standard_errors))
  # return(list(estimates = est1, est2, est3, est4, standard_errors = sd1, sd2, sd3, sd4))
  
}

fullinteracts <- function(x, maxpower) {
  k <- ncol(x)
  Ik <- rbind(matrix(0, 1, k), diag(k))
  tmp <- Ik
  
  for (i in 1:(maxpower-1)) {
    tmp <- do.call(rbind, replicate(nrow(tmp), Ik, simplify = FALSE)) + tmp[rep(1:nrow(tmp), each = k + 1), ]
  }
  
  inters <- unique(tmp)
  # Sort the unique rows
  inters <- inters[do.call(order, as.data.frame(inters)), ]
  
  # Initialize an empty matrix D
  allinters <- matrix(NA, nrow = nrow(x), ncol = nrow(inters))
  
  # Loop over each row of x
  for (i in 1:nrow(x)) {
    # Loop over each row of inters
    for (j in 1:nrow(inters)) {
      # Compute the element-wise power
      allinters[i, j] <- prod(x[i, ] ^ inters[j, ])
    }
  }
  
  varx <- apply(allinters, 2, var)
  interacts <- allinters[, varx > 0]
  inters <- inters[varx > 0, ]
  
  # Find the unique rows
  transposed_interacts <- t(interacts)
  # Find the unique rows
  unique_rows <- unique(t(interacts))
  # Sort the unique rows
  a <- unique_rows[do.call(order, as.data.frame(unique_rows)), ]
  
  b <- match(apply(sorted_unique_rows, 1, paste, collapse = ","), apply(transposed_interacts, 1, paste, collapse = ","))
  
  interacts <- t(a)
  
  inters <- inters[b, ]
  
  return(list(interacts = interacts, inters = inters))
}

periods <- 10
sample_size <- 100

powers <- 1 # order of polynomial terms for the basis functions
# 1 - only linear terms
# 1.5 - linear and squared terms
# 2 - linear, squared terms, and interactions
# >=2 - all powers and interactions up to that order

if (powers == 1) {
  suffix <- ''
} else if (powers == 1.5) {
  suffix <- '_squares'
} else if (powers == 2) {
  suffix <- '_quad'
} else if (powers > 2) {
  suffix <- paste0("_", powers)
}

setwd("C:\\Users\\mengs\\Documents\\GitHub\\value_function_HD\\Rfiles")

data_path <- paste0("MCdata_t=", periods, "_n=", sample_size, suffix, ".RData")
save_path_est <- paste0("MCestimates_t=", periods, "_n=", sample_size, suffix, ".RData")
save_path_sd <- paste0("MCstandard_errors_t=", periods, "_n=", sample_size, suffix, ".RData")

load(data_path) # load simulated data from just_data.R

data <- mat_data$data

draws <- dim(data)[4] # find number of Monte Carlo draws
periods <- dim(data)[2] # number of individual engines in each Monte Carlo sample
sample_size <- dim(data)[1] # number of individual engines in each Monte Carlo sample

number_multi_index_states <- ((dim(data)[3] - 2) / 2) # number of multi-index (iid) states
delta <- 0.9 # discount factor

splitnum <- 5 # number of subsamples to use for sample splitting

splitsize <- floor(sample_size / splitnum)

if ((sample_size / splitnum) - splitsize > 0) {
  data <- data[1:(splitsize * splitnum), , , ]
}

# flatten the data
Data <- drop(aperm(array(aperm(data, c(2, 1, 3, 4)), c(dim(data)[1] * dim(data)[2], 1, dim(data)[3], dim(data)[4])), c(1, 3, 4, 2)))

# figure out all the split indices
I <- as.logical(kronecker(diag(splitnum), matrix(1, nrow = periods * splitsize, ncol = 1)))
I <- matrix(I, nrow = periods * splitsize * splitnum, ncol = splitnum)

# estimate parameters for each Monte Carlo draw
est <- array(NA, dim = c(2, 4, draws)) # we will save the coefficient estimates here
sd <- array(NA, dim = c(2, 4, draws)) # we will save the estimated standard errors here

# note: estimation code returns two matrices. The first 'estimates' or
# 'esttmp' below, is a 2-by-4 matrix whose columns correspond to estimates
# of the two parameters. The first column has no bias correction, the final
# column is fully bias corrected. The 2nd and 3rd columns are partially
# adjusted using only \phi_2 and using \phi_2 and \phi_3 respectively.
# 'standard_errors' or 'sdtmp' below is a 2-by-4 matrix whose entries are
# estimated standard errors for the corresponding elements of 'estimates'.

# draws = 50

for (i in 1:draws) {

  result <- estimate3(Data[, , i], I, number_multi_index_states, delta, powers)

  est[, , i] <- result$estimates
  sd[, , i] <- result$standard_errors
}

# # Set the number of cores for parallel processing
# num_cores <- detectCores()
# registerDoParallel(num_cores)
# 
# # Define a function to estimate for each iteration
# estimate_function <- function(i) {
#   result <- estimate3(Data[, , i], I, number_multi_index_states, delta, powers)
#   list(estimates = result$estimates, standard_errors = result$standard_errors)
# }
# 
# # Apply the function in parallel using foreach
# results <- foreach(i = 1:draws, .combine = rbind) %dopar% {
#   estimate_function(i)
# }
# 
# # Unpack the results
# est <- array(unlist(lapply(results, function(x) x$estimates)), dim = dim(est))
# sd <- array(unlist(lapply(results, function(x) x$standard_errors)), dim = dim(sd))
# 
# # Stop parallel processing
# stopImplicitCluster()

# save all the estimates from the various bootstrap draws. We save the
# estimated coefficients and the standard errors in separate files
save(est, file = save_path_est)
save(sd, file = save_path_sd)
