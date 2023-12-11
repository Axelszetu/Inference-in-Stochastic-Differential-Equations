#Solution to week 1 exercise 3
#1
increments <- rbinom(n = 100, size = 1, prob = 0.5)*2 - 1
W <- cumsum(increments)
S0 <- 0
S <- S0 + W
k <- (1:100)

A <- numeric(100)
A[1] <- -W[1]*S[1]
for (i in (2:100)){
  A[i] <- A[i-1] + (W[i-1] - W[i])*S[i]
}

compensator <- W^2/2 - k/2
intW <- cumsum(W)

profit_path <- data.frame('A' = A, 'k' = k, 'W' = W, 'l' = compensator, 'int' = intW)
profit_plot <- ggplot(data = profit_path, aes(y = A, x = k)) +
  geom_line(mapping = aes(x = k, y = A), color = 'red') +
  geom_line(mapping = aes(x = k, y = W), color = 'blue') +
  geom_line(mapping = aes(x = k, y = l), color = 'purple') +
  geom_line(mapping = aes(x = k, y = W - l), color = 'green')
profit_plot
