# # generate table with results with quadratic and interactions
# rm(list = ls())
# 
# meanbias <- matrix(NA, nrow = 8, ncol = 4)
# meansd <- matrix(NA, nrow = 8, ncol = 4)
# truesd <- matrix(NA, nrow = 8, ncol = 4)
# medianbias <- matrix(NA, nrow = 8, ncol = 4)
# mediansd <- matrix(NA, nrow = 8, ncol = 4)
# coverageprob <- matrix(NA, nrow = 8, ncol = 4)
# 
# load('MCestimates_t=10_n=100_quad.Rdata')
# load('MCstandard_errors_t=10_n=100_quad.Rdata')
# 
# meanbias[1:2, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[1:2, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[1:2, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[1:2, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[1:2, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[1:2, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 
# 
# load('MCestimates_t=10_n=300_quad.Rdata')
# load('MCstandard_errors_t=10_n=300_quad.Rdata')
# 
# meanbias[3:4, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[3:4, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[3:4, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[3:4, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[3:4, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[3:4, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 
# 
# load('MCestimates_t=10_n=1000_quad.Rdata')
# load('MCstandard_errors_t=10_n=1000_quad.Rdata')
# 
# meanbias[5:6, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[5:6, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[5:6, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[5:6, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[5:6, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[5:6, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 
# 
# load('MCestimates_t=10_n=10000_quad.Rdata')
# load('MCstandard_errors_t=10_n=10000_quad.Rdata')
# 
# meanbias[7:8, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[7:8, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[7:8, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[7:8, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[7:8, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[7:8, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 
# 
# rownames <- c('continuation (t=10, n=100)', 'replacement (t=10, n=100)', 'continuation (t=10, n=300)', 'replacement (t=10, n=300)', 'continuation (t=10, n=1000)', 'replacement (t=10, n=1000)', 'continuation (t=10, n=10000)', 'replacement (t=10, n=10000)')
# columnnames <- c('absolute mean bias (not corrected)', 'absolute mean bias (corrected)', 'median standard error', 'standard deviation (not corrected)', 'standard deviation (corrected)', '95% confidence interval coverage (not corrected)', '95% confidence interval coverage (corrected)')
# 
# T <- data.frame(meanbias[, 1], meanbias[, 4], mediansd[, 4], truesd[, 1], truesd[, 4], coverageprob[, 1], coverageprob[, 4])
# colnames(T) <- columnnames
# rownames(T) <- rownames
# 
# write.csv(T, 'some_results_quad.csv', row.names = TRUE)

# generate table with results with linear terms
rm(list = ls())

meanbias <- matrix(NA, nrow = 8, ncol = 4)
meansd <- matrix(NA, nrow = 8, ncol = 4)
truesd <- matrix(NA, nrow = 8, ncol = 4)
medianbias <- matrix(NA, nrow = 8, ncol = 4)
mediansd <- matrix(NA, nrow = 8, ncol = 4)
coverageprob <- matrix(NA, nrow = 8, ncol = 4)

load('MCestimates_t=10_n=100.Rdata')
load('MCstandard_errors_t=10_n=100.Rdata')

meanbias[1:2, ] <- abs(apply(est[,,1:50], c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
meansd[1:2, ] <- apply(sd[,,1:50], c(1, 2), mean, na.rm = TRUE)
truesd[1:2, ] <- sqrt(apply(est[,,1:50], c(1, 2), function(x) var(x, na.rm = TRUE)))
medianbias[1:2, ] <- abs(apply(est[,,1:50], c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
mediansd[1:2, ] <- apply(sd[,,1:50], c(1, 2), median, na.rm = TRUE)
coverageprob[1:2, ] <- apply(abs(est[,,1:50] - c(-0.5, -1)) / sd[,,1:50] <= 1.96, c(1, 2), mean, na.rm = TRUE) 

load('MCestimates_t=10_n=300.Rdata')
load('MCstandard_errors_t=10_n=300.Rdata')

meanbias[3:4, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
meansd[3:4, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
truesd[3:4, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
medianbias[3:4, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
mediansd[3:4, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
coverageprob[3:4, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 

load('MCestimates_t=10_n=1000.Rdata')
load('MCstandard_errors_t=10_n=1000.Rdata')

meanbias[5:6, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
meansd[5:6, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
truesd[5:6, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
medianbias[5:6, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
mediansd[5:6, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
coverageprob[5:6, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE)

load('MCestimates_t=10_n=10000.Rdata')
load('MCstandard_errors_t=10_n=10000.Rdata')

meanbias[7:8, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
meansd[7:8, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
truesd[7:8, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
medianbias[7:8, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
mediansd[7:8, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
coverageprob[7:8, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 

rownames <- c('continuation (t=10, n=100)', 'replacement (t=10, n=100)', 'continuation (t=10, n=300)', 'replacement (t=10, n=300)', 'continuation (t=10, n=1000)', 'replacement (t=10, n=1000)', 'continuation (t=10, n=10000)', 'replacement (t=10, n=10000)')
columnnames <- c('absolute mean bias (not corrected)', 'absolute mean bias (corrected)', 'median standard error', 'standard deviation (not corrected)', 'standard deviation (corrected)', '95% confidence interval coverage (not corrected)', '95% confidence interval coverage (corrected)')

T <- data.frame(meanbias[, 1], meanbias[, 4], mediansd[, 4], truesd[, 1], truesd[, 4], coverageprob[, 1], coverageprob[, 4])
colnames(T) <- columnnames
rownames(T) <- rownames

write.csv(T, 'some_results.csv', row.names = TRUE)

# # generate table with results with linear and squared terms
# rm(list = ls())
# 
# meanbias <- matrix(NA, nrow = 8, ncol = 4)
# meansd <- matrix(NA, nrow = 8, ncol = 4)
# truesd <- matrix(NA, nrow = 8, ncol = 4)
# medianbias <- matrix(NA, nrow = 8, ncol = 4)
# mediansd <- matrix(NA, nrow = 8, ncol = 4)
# coverageprob <- matrix(NA, nrow = 8, ncol = 4)
# 
# load('MCestimates_t=10_n=100_squares.Rdata')
# load('MCstandard_errors_t=10_n=100_squares.Rdata')
# 
# meanbias[1:2, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[1:2, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[1:2, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[1:2, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[1:2, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[1:2, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE)
# 
# load('MCestimates_t=10_n=300_squares.Rdata')
# load('MCstandard_errors_t=10_n=300_squares.Rdata')
# 
# meanbias[3:4, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[3:4, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[3:4, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[3:4, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[3:4, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[3:4, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE) 
# 
# load('MCestimates_t=10_n=1000_squares.Rdata')
# load('MCstandard_errors_t=10_n=1000_squares.Rdata')
# 
# meanbias[5:6, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[5:6, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[5:6, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[5:6, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[5:6, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[5:6, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE)
# 
# load('MCestimates_t=10_n=10000_squares.Rdata')
# load('MCstandard_errors_t=10_n=10000_squares.Rdata')
# 
# meanbias[7:8, ] <- abs(apply(est, c(1, 2), mean, na.rm = TRUE) - c(-0.5, -1))
# meansd[7:8, ] <- apply(sd, c(1, 2), mean, na.rm = TRUE)
# truesd[7:8, ] <- sqrt(apply(est, c(1, 2), function(x) var(x, na.rm = TRUE)))
# medianbias[7:8, ] <- abs(apply(est, c(1, 2), median, na.rm = TRUE) - c(-0.5, -1))
# mediansd[7:8, ] <- apply(sd, c(1, 2), median, na.rm = TRUE)
# coverageprob[7:8, ] <- apply(abs(est - c(-0.5, -1)) / sd <= 1.96, c(1, 2), mean, na.rm = TRUE)
# 
# rownames <- c('continuation (t=10, n=100)', 'replacement (t=10, n=100)', 'continuation (t=10, n=300)', 'replacement (t=10, n=300)', 'continuation (t=10, n=1000)', 'replacement (t=10, n=1000)', 'continuation (t=10, n=10000)', 'replacement (t=10, n=10000)')
# columnnames <- c('absolute mean bias (not corrected)', 'absolute mean bias (corrected)', 'median standard error', 'standard deviation (not corrected)', 'standard deviation (corrected)', '95% confidence interval coverage (not corrected)', '95% confidence interval coverage (corrected)')
# 
# T <- data.frame(meanbias[, 1], meanbias[, 4], mediansd[, 4], truesd[, 1], truesd[, 4], coverageprob[, 1], coverageprob[, 4])
# colnames(T) <- columnnames
# rownames(T) <- rownames
# 
# write.csv(T, 'some_results_squares.csv', row.names = TRUE)
