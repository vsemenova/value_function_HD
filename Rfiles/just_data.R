DPsolve_cont3 <- function(params, s_grid, z_grid, tn) {
  
  z_dimension <- length(params$multi_index)  # number of iid state variables that affect transitions
  
  size_s <- length(s_grid)
  size_z <- length(z_grid)
  
  # Create a mesh grid using expand.grid()
  grid <- expand.grid(x = s_grid, y = z_grid)
  
  # Extract the grid points as separate vectors
  S_grid <- array(grid$x, c(size_s, size_z))
  Z_grid <- array(grid$y, c(size_s, size_z))
  
  EV1 <- matrix(0, nrow = size_z, ncol = size_s)  # initial guess of the expected value function
  step_change <- 1
  
  # Mileage transitions (renewal)
  replacement_transitions <- 1 + rnorm(prod(c(size_z, size_s, tn)), mean = params$state_transition[1], sd = params$state_transition[2])^2
  replacement_transitions <- array(replacement_transitions, dim = c(size_z, size_s, tn))
  
  # Mileage transitions (non-renewal)
  z_grid <- matrix(z_grid)  # Assuming z_grid is a vector, convert it to a matrix if necessary
  maintenance_transitions <- rep(s_grid, each = size_z) + rnorm(prod(c(size_z, size_s, tn)), mean = params$state_transition[1] + t(z_grid), sd = params$state_transition[2])^2
  maintenance_transitions <- array(maintenance_transitions, dim = c(size_z, size_s, tn))
  
  # Zchi
  Zchi <- rchisq(n = tn * ceiling(z_dimension/2), df = 1)
  Zchi <- matrix(Zchi, nrow = ceiling(z_dimension/2), ncol = tn)
  Zchi <- array(Zchi, dim = c(1, ceiling(z_dimension/2), tn))
  
  # Zdisc
  Zdisc <- rbinom(n = floor(z_dimension/2) * tn, size = 1, prob = 0.5)
  Zdisc <- matrix(Zdisc, nrow = floor(z_dimension/2), ncol = tn)
  Zdisc <- array(Zdisc, dim = c(1, floor(z_dimension/2), tn))
  
  # Ztemp
  # Calculate the number of slices to include from Zchi based on the z_dimension
  slices <- 1:(nrow(Zchi) - (z_dimension %% 2))
  # Extract the slices from Zchi and concatenate with Zdisc along the second dimension
  Ztemp <- abind(Zchi[, slices, , drop = FALSE], Zdisc, along = 2)
  Ztemp <- array(Ztemp, dim = c(2*floor(z_dimension/2), tn))
  
  # z_draws_temp
  slices <- (dim(Zchi)[1] - (z_dimension %% 2) + 1):dim(Zchi)[1]
  Zchi_sub <- array(Zchi[slices, 1, ], dim = c(1, tn))
  z_draws_temp <- abind(Ztemp, Zchi_sub, along = 1)
  
  # # Concatenate Ztemp and Zchi_sub along the first dimension
  # z_draws_temp <- abind(Ztemp, Zchi[, slices, , drop = FALSE], along = 2)
  result <- params$multi_index %*% z_draws_temp
  permuted_array <- array(result, dim = c(1, dim(result))) 
  z_draws <- array(permuted_array, dim = c(size_z, size_s, tn))
  
  EV <- matrix(0, nrow = size_z, ncol = size_s)  # placeholder for expected value
  
  # Reshape the data into long format
  S_vec <- as.vector(S_grid)
  Z_vec <- as.vector(Z_grid)
  EV1_vec <- as.vector(EV1)
  
  replacement_transitions_vec <- as.vector(replacement_transitions)
  z_draws_vec <- as.vector(z_draws)
  
  replacement_EV <- array(NA, dim = c(size_z, size_s, tn))
  maintenance_EV <- array(NA, dim = c(size_z, size_s, tn))
  
  while (step_change > 0.0000001) {
    
    for (i in 1:tn) {
      replacement_result <- interp2(s_grid, z_grid, EV1, as.vector(replacement_transitions[,,i]), as.vector(z_draws[,,i]), method = "linear")
      replacement_EV[ , ,i] <- array(replacement_result, c(size_z, size_s))
      
      maintenance_result <- interp2(s_grid, z_grid, EV1, as.vector(maintenance_transitions[,,i]), as.vector(z_draws[,,i]), method = "linear")
      maintenance_EV[ , ,i] <- array(maintenance_result, c(size_z, size_s))
    }
    
    EV <- log(exp(params.maintenance_factor * sqrt(S_grid) + params.beta * mean(maintenance_EV, 3)) +
                exp(params.replacement_cost + params.beta * mean(replacement_EV, 3)))
    
    step_change <- abs(EV - EV1)
    step_change <- max(step_change)
    
    EV1 <- EV
      
  }
  
  EV1_vec <- as.vector(EV1)
  
  for (i in 1:tn) {
    replacement_result <- interp(S_grid, Z_grid, EV1, as.vector(replacement_transitions[,,i]), as.vector(z_draws[,,i]), linear = FALSE, extrap = TRUE)
    replacement_EV[ , ,i] <- array(diag(replacement_result), c(size_z, size_s))
    
    maintenance_result <- interp(S_grid, S_grid, EV1, as.vector(maintenance_transitions[,,i]), as.vector(z_draws[,,i]), linear = FALSE, extrap = TRUE)
    maintenance_EV[ , ,i] <- array(diag(maintenance_result), c(size_z, size_s))
  }
  
  grid <- meshgrid(s_grid, z_grid)
  S_grid <- grid$S_grid
  Z_grid <- grid$Z_grid
  
  # Evaluate choice probabilities on the grid of states
  CCP_grid <- 1 / (1 + exp(params$maintenance_factor * sqrt(S_grid) + 
                             params$beta * mean(maintenance_EV, dims = 3) - 
                             (params$replacement_cost + params$beta * mean(replacement_EV, dims = 3))))
  
  # Define CCP function using interpolation
  CCP_function <- function(s, z) {
    for (i in 1:tn) {
      interp_result <- interp2(s_grid, z_grid, , z = CCP_grid, xi = s, yi = z, method = "linear")
      replacement_EV[ , ,i] <- array(replacement_result, c(size_z, size_s))
      return(array(interp_result, c(size_z, size_s)))
    }
    
    interp2(x = S_grid, y = Z_grid, z = CCP_grid, xi = s, yi = z, method = "linear", extrap = TRUE)
  }
  
  return(CCP_function)
}

MCdata3 <- function(CCP_function, params, sample_size, periods) {
  # The arguments:
  # CCP function - function that returns probability of replacement at each state.
  # params - list containing parameter values.
  # sample_size - number of 'buses' for which to draw data.
  # periods - number of periods for which each bus is observed.
  
  burn_in <- 10  # burn-in period so that the mileage for each bus does not start at 1
  z_dimension <- length(params$multi_index)  # number of iid state variables that affect transitions
  data <- array(0, dim = c(sample_size, periods, 4 + 2 * z_dimension))  # where we will store the data
  
  # Generate transition state variables (iid chi-squared or Bernoulli)
  Zchi <- rchisq(n = sample_size*(periods+burn_in+1) * ceiling(z_dimension/2), df = 1)
  Zchi <- matrix(Zchi, nrow = ceiling(z_dimension/2), ncol = sample_size*(periods+burn_in+1))
  Zchi <- array(Zchi, dim = c(1, ceiling(z_dimension/2), sample_size*(periods+burn_in+1)))
  
  # Zdisc
  Zdisc <- rbinom(n = floor(z_dimension/2) * sample_size*(periods+burn_in+1), size = 1, prob = 0.5)
  Zdisc <- matrix(Zdisc, nrow = floor(z_dimension/2), ncol = sample_size*(periods+burn_in+1))
  Zdisc <- array(Zdisc, dim = c(1, floor(z_dimension/2), sample_size*(periods+burn_in+1)))
  
  # Ztemp
  # Calculate the number of slices to include from Zchi based on the z_dimension
  slices <- 1:(nrow(Zchi) - (z_dimension %% 2))
  # Extract the slices from Zchi and concatenate with Zdisc along the second dimension
  Ztemp <- abind(Zchi[, slices, , drop = FALSE], Zdisc, along = 2)
  Ztemp <- array(Ztemp, dim = c(2*floor(z_dimension/2), sample_size*(periods+burn_in+1)))

  # z_draws_temp
  slices <- (dim(Zchi)[1] - (z_dimension %% 2) + 1):dim(Zchi)[1]
  Zchi_sub <- array(Zchi[slices, 1, ], dim = c(1, sample_size*(periods+burn_in+1)))
  Z_draws <- abind(Ztemp, Zchi_sub, along = 1)
  
  # Compute z_draws by multiplying Z_draws with params.multi_index
  z_draws <- t(Z_draws) %*% params$multi_index
  
  # Define the indices to save each variable
  z_block <- 3:(3 + z_dimension - 1)
  choicen_block <- 3 + z_dimension
  sn_block <- 3 + z_dimension + 1
  zn_block <- 3 + z_dimension + 2:(3 + 2 * z_dimension + 1)
  
  durations <- vector()
  
  # Generate the data
  for (i in 1:sample_size) {
    s <- 1
    last_replacement <- 1
    for (t in 1:(periods + burn_in)) {
      index <- (i - 1) * (periods + burn_in + 1) + t
      
      # Choices drawn using precalculated CCP function
      choice <- rbinom(n = 1, size = 1, prob = CCP_function(s, z_draws[index]))
      
      if (t > burn_in) {
        # Save this period's data
        data[i, t - burn_in, 1] <- choice
        data[i, t - burn_in, 2] <- s
        data[i, t - burn_in, z_block] <- Z_draws[index, ]
        
        data[i, t - burn_in, zn_block] <- Z_draws[index + 1, ]
        if (t > burn_in + 1) {
          data[i, t - burn_in - 1, choicen_block] <- choice
        }
        if (choice == 1) {
          durations <- c(durations, t - max(last_replacement, burn_in))
          last_replacement <- t
        }
      }
      
      if (choice == 1) {
        s <- 1
      }
      
      # Simulate the next period mileage
      s <- s + rnorm(n = 1, mean = params$state_transition[1] + z_draws[index], sd = params$state_transition[2])^2
      if (t > burn_in) {
        data[i, t - burn_in, sn_block] <- s
      }
    }
    
    choice <- rbinom(n = 1, size = 1, prob = CCP_function(s, z_draws[index]))
    data[i, t - burn_in, choicen_block] <- choice
  }
  
  return(list(data = data, durations = durations))
}

# library(matrixStats)
library(MASS)
library(abind)
library(akima)
library(stats)

# Set the current working directory

set.seed(2)  # set the seed

sample_size <- 100  # number of individual engines in each Monte Carlo sample
constant_heterogeneity <- 2  # is the heterogeneity in the utility function constant or iid?
number_multi_index_states <- 5  # number of states that affect mileage transitions through a multi-index
params <- list(
  maintenance_factor = -0.5,  # per-period maintenance cost coefficient
  replacement_cost = -1,  # replacement cost
  state_transition = c(0.2, 1),  # the constants in the mileage transition density
  beta = 0.9,  # discount factor
  multi_index = 0.1 * (1:number_multi_index_states)^(-2)  # multi-index coefficients for the states that affect mileage transitions
)

draws <- 500  # number of Monte Carlo draws for generating data
periods <- 10  # number of individual engines in each Monte Carlo sample

save_path <- paste0('MCdata_t=', periods, '_n=', sample_size, '.RData')

s_grid <- seq(0, 100, length.out = 100)  # mileage grid points for iterating the Bellman equation
z_grid <- seq(0, 100, length.out = 100)  # multi-index grid points for iterating the Bellman equation

tn <- 10000  # number of draws used to form Monte Carlo expectation in DP solution and data generation

cat('Solving the dynamic programming problem...\n')
CCP_function <- DPsolve_cont3(params, s_grid, z_grid, tn)
cat('Solved!\n')

data <- array(NA, dim = c(sample_size, periods, 4 + 2 * number_multi_index_states, draws))

for (b in 1:draws) {
  datatmp <- MCdata3(CCP_function, params, sample_size, periods)
  data[, , , b] <- datatmp
}

save(data, file = save_path)
