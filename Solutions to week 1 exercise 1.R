#Solutions to week 1 exercise 1
library(ggplot2)

dx <- function(x){
  x*(1-x^2)
}

analytic_solver <- function(x0, t){
  sign(x0)*exp(t)*(x0^(-2) + exp(2*t) -1)^(-1/2)
}

solution_curve <- function(x0){
  force(x0)
  f <- function(t){
    sign(x0)*exp(t)*(x0^(-2) + exp(2*t) -1)^(-1/2)
  }
}

curve_2 <- solution_curve(2)

double_well_plot <- ggplot() +
                    xlim(0,5) +
                    geom_function(fun = solution_curve(2), color = "red") +
                    geom_function(fun = solution_curve(0.2), color = "blue") +
                    geom_function(fun = solution_curve(-2), color = "purple") +
                    geom_function(fun = solution_curve(-0.2), color = "green")
double_well_plot

double_well_plot_500 <- ggplot() +
  xlim(0,500) +
  geom_function(fun = solution_curve(2), color = "red") +
  geom_function(fun = solution_curve(0.2), color = "blue") +
  geom_function(fun = solution_curve(-2), color = "purple") +
  geom_function(fun = solution_curve(-0.2), color = "green")
double_well_plot_500

#It seems that when we reach t=350, the last term becomes so close to zero
#that it numerically kills the process. One solution could be rescaling of t.

double_well_potential <- ggplot() + xlim(-1.5,1.5) + geom_function(fun = function(x) (x^4)/4 - (x^2)/2)
double_well_potential

#The naming makes sense, since there are two fixed points that attract solutions.

#Implementing Euler scheme for variable parameters
Euler_solver <- function(x0, tmax, deltat){
  #browser()
  t <- seq(from = 0, to = tmax, by = deltat)
  length <- length(t)
  x <- numeric(length)
  x[1] <- x0
  for (i in (2:length)){
    x[i] <- x[i-1] + dx(x[i-1])*deltat
  }
  path <- data.frame("t" = t, "X" = x)
  path
}
Euler_path <- Euler_solver(2, 10, 0.05)

Euler_path_plot <- ggplot(data = Euler_path, mapping = aes(x = t, y = X)) + geom_line()
Euler_path_plot

#Implementing Heun scheme for variable parameters
Heun_solver <- function(x0, tmax, deltat){
  #browser()
  t <- seq(from = 0, to = tmax, by = deltat)
  length <- length(t)
  x <- numeric(length)
  x[1] <- x0
  for (i in (2:length)){
    x_tilde <- x[i-1] + dx(x[i-1])*deltat
    x[i] <- x[i-1] + (dx(x[i-1]) + dx(x_tilde))*deltat*0.5
  }
  path <- data.frame("t" = t, "X" = x)
  path
}
Heun_path <- Euler_solver(2, 10, 0.05)
Heun_path_plot <- ggplot(data = Heun_path, mapping = aes(x = t, y = X)) + geom_line()
Heun_path_plot

#Implementing RK4 scheme for variable parameters
RK4_solver <- function(x0, tmax, deltat){
  t <- seq(from = 0, to = tmax, by = deltat)
  length <- length(t)
  x <- numeric(length)
  x[1] <- x0
  for (i in 2:length){
    k1 <- dx(x[i-1])
    k2 <- dx(x[i-1] + deltat*k1*0.5)
    k3 <- dx(x[i-1] + deltat*k2*0.5)
    k4 <- dx(x[i-1] + deltat*k3)
    x[i] <- x[i-1] + deltat*(k1+2*k2+2*k3+k4)/6
  }
  path <- data.frame("t" = t, "X" = x)
  path
}
RK4_path <- RK4_solver(2, 10, 0.05)
RK4_path_plot <- ggplot(data = RK4_path, mapping = aes(x = t, y = X)) + geom_line()
RK4_path_plot

#Now, we have the slave work of doing making plots for different values of dt, for each method.
#If I am clever I might make a function computing the error and making the plot for a given method, and vector of dt values.

abs_error_computer <- function(x0, tmax, deltat, method){
  t <- seq(from = 0, to = tmax, deltat)
  path_sim <- method(x0, tmax, deltat)
  analytic_sol <- function(t){
    analytic_solver(x0 = x0, t = t)
  }
  path_an <- analytic_sol(t)
  abs_error <- abs(path_sim[,2] - path_an)
  plt_data <- data.frame('t' = t, 'error' = abs_error)
  ggplot(data = plt_data, mapping = aes(x = t, y = error)) + geom_line() + labs(title = paste("Method =", deparse(substitute(method)), "deltat =", deltat))
}
abs_error_computer(x0 = 2, tmax = 10, deltat = 0.05, method = Heun_solver)
#As suspected, the absolute error goes to 0 quite quickly, as we would expect since the system has an attracting stationary point.
#We need to shorten the time horizon to get something interesting: We need the error of the numerical method to propagate faster than the system converges.

abs_error_computer(x0 = 2, tmax = 1, deltat = 0.005, method = Euler_solver)

#Now, lets make a function that takes a vector of deltats and makes the plot for each of them.

timestep_comparison <- function(x0, tmax, deltats, method){
  individual_computer <- function(deltat){
    abs_error_computer(x0 = x0, tmax = tmax, deltat = deltat, method = method)
  }
  out <- lapply(X = deltats, FUN = individual_computer)
  out
}

deltats <- c(0.05, 0.005, 0.0005)
timestep_comparison(x0 = 2, tmax = 1, deltats = deltats, method = Euler_solver)
timestep_comparison(x0 = 2, tmax = 1, deltats = deltats, method = Heun_solver)
timestep_comparison(x0 = 2, tmax = 1, deltats = deltats, method = RK4_solver)

#From the plots we conclude that the Euler solver is worse than the Heun solver, which again is worse than the RK4 solver.

#Now, we implement a Strang splitting scheme
Strang_solver <- function(x0, tmax, deltat){
  t <- seq(from = 0, to = tmax, by = deltat)
  length <- length(t)
  x <- numeric(length)
  x[1] <- x0
  for (i in 2:length){
    x[i] <- sign(x[i-1])*((x[i-1]^(-2) + deltat) * exp(-2*deltat) + deltat)^(-(1/2))
  }
  path <- data.frame("t" = t, "X" = x)
  path
}

timestep_comparison(x0 = 2, tmax = 1, deltats = deltats, method = Strang_solver)
#It seems to work, and better than the Heun solver,yet worse than the RK$ solver in this case.
#Predrags solution outlines some of the peculiar properties of this model and how it interacts specifically with the Strang and LT algorithms.