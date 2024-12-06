#### sheet 7
# problem 9: two-sample testing
n <- 16
p <- 2
S <- matrix(c(63.7, -3, -3, 33.0), ncol=2)
S_inv <- solve(S)
d_mean <- c(0.13, 3.81)
T2 <- (n-2)*n/((n-1)*2)*t(d_mean) %*% S_inv %*% d_mean
T2 > qf(0.95, p, n-p)
# cannot reject H0

# problem 10: MANOVA
n <- 5
g <- 3
p <- 3
X1_mean <- c(0.28, -0.12, -0.08)
X2_mean <- c(0.08, -0.44, 0.12)
X3_mean <- c(-0.16, -0.8, -0.04)
x_mean <- c(0.067, -0.453, 0)

X1_outer <- X1_mean %*% t(X1_mean)
X2_outer <- X2_mean %*% t(X2_mean)
X3_outer <- X3_mean %*% t(X3_mean)
Xmean_outer <- x_mean %*% t(x_mean)
B <- n*(X1_outer + X2_outer + X3_outer-g*Xmean_outer)
Total_T <- matrix(c(4.853, 0.093, 0.44, 0.093, 2.437, 0.64, 0.44, 0.64, 5.28), ncol=3)
W <- Total_T-B
wilks_l <- det(W)/det(Total_T)
F_stat <- (1-sqrt(wilks_l))/sqrt(wilks_l)*(g*(n-1)-p+1)/p
qf(0.95, 2*p, 2*(g*(n-1)-p+1))
# cannot reject H0
