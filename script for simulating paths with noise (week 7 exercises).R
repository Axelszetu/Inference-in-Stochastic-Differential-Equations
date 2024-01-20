#Script for investigating stability of 1d dynamical systems with noise
dynamics_1 <- function(x){
  4*(x^2) - 16
}

noisy_path_maker <- function(dynamics, x0, tmax, deltat, sigma){
  t <- seq(from = 0, to = tmax, by = deltat)
  dW <- rnorm(n = length(t)-1, sd = sigma * sqrt(deltat))
  x <- numeric(length = length(t))
  x[1] <- x0
  for (i in (2:length(t))){
    x[i] <- x[i-1] + dynamics(x[i-1])*deltat + dW[i-1]
  }
  out <- data.frame('x' = x, 't' = t)
  plot(x = out$t, y = out$x)
}

noisy_path_maker(dynamics = dynamics_1, x0 = 0, tmax = 5, deltat = 0.01, sigma = 0.1)
noisy_path_maker(dynamics = dynamics_1, x0 = 2, tmax = 5, deltat = 0.01, sigma = 0.1)
  
