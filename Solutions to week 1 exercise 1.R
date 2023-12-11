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
