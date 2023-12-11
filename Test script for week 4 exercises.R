#Test script for functions in Week4 exercises
path <- make_path(tmax = 10, x1_0 = 600, x2_0 = 30, x3_0 = 10^5, rho = 0.108, delta = 0.5, eta = 9.5*10^(-6), lambda = 36, N1 = 1000, C = 3, sigma1 = 0.1, sigma2 = 0.1, sigma3 = 0.1, deltat = 0.001)
drift <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  y1 <- lambda - rho*x1 - eta*x1*x3
  y2 <- eta*x1*x3 - delta*x2
  y3 <- N1*delta*x2 - C*x3
  y <- c(y1, y2, y3)
  y
}
drift(c(1,2,3))
