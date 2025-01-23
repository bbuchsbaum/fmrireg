# library(devtools)
# load_all()
# 
# library(CVXR)
# ## k is the number of condition
# ## d is dimension of basis
# ## n is the number of rows
# ## k + d free parameters (1 for each basis, one amplitude h for each condition)
# 
# onsets <- c(5, 25, 32, 44, 46, 49, 62, 76, 82, 91)
# cond <- rep(letters[1:2], 5)
# des <- data.frame(onsets=onsets, cond=cond, block=rep(1,length(onsets)))
# sframe <- sampling_frame(120, TR=1)
# ev <- event_model(onsets ~ hrf(cond, basis="spmg3"), block = ~ block, data=des, sampling_frame=sframe)
# 
# reg <- regressor(onsets, hrf=HRF_SPMG1, duration=2, amplitude=runif(1:length(onsets))*2)
# ytruth <- evaluate(reg, 0:119)
# y <- ytruth + rnorm(1:length(ytruth), sd=.2)
# X <- as.matrix(design_matrix(ev))
# 
# 
# library(mgcv)
# ht <- hrf_time
# ev <- event_model(onsets ~ trialwise(basis=ht), block = ~ block, data=des, sampling_frame=sframe)
# X <- as.matrix(design_matrix(ev))
# 
# gam.1 <- bam(ytruth ~ s(X))
# 
# p <- ncol(X)
# n <- nrow(X)
# beta <- CVXR::Variable(p) 
# beta0 <- Variable(1)
# loss <- sum((y - beta0 + X %*% beta)^2) / (2 * n)
# obj <- loss
# prob <- Problem(Minimize(obj))
# solve(prob)
# 
# beta2 <- Variable(2)
# h <- Variable(3)
# loss2 <- sum((y - beta0 + (X %*% (h %*% t(beta2))))^2) / (2 * n)
# 
# obj <- loss2
# prob <- Problem(Minimize(obj))
# solve(prob, gp=TRUE)
# 
# 
# fn <- function(b) {
#   o <- sum((y - X %*% c(b[1], kronecker(b[2:3], b[4:6])))^2) / (2 * length(y))
#   print(o)
#   o
# }




beta_vals <- matrix(0, nrow = ncol(X), ncol = length(y))
lambda_vals <- 10^seq(-2, log10(50), length.out = TRIALS)
beta <- Variable(p)  
loss <- sum((y - X %*% beta)^2) / (2 * n)



