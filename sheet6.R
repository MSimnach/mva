### Sheet 6 -- Problem 8
# Setup
n <- 35
p <- 2
alpha <- 0.05
x_f <- c(22.860, 24.297)
mu_0 <- c(20, 25)
Sigma_f <- matrix(c(20, 0, 0, 25), ncol = 2)
S_f <- matrix(c(17.683, 20.290, 20.290, 24.407), ncol = 2)

# a)
T2 <- n*t(xf-mu_0)%*%solve(Sigma_f)%*%(xf-mu_0)
chisq_2 <- qchisq(0.95,2)

# b)
TS2 <- 33/68* n*t(xf-mu_0)%*%solve(S_f)%*%(xf-mu_0)
F_2_33 <- qf(0.95, 2, 33)

# c)
# Univariate confidence intervals (t-test)
ci_x1_uni <- c(x_f[1] - qt(1-alpha/2, n-1)*sqrt(S_f[1]/n),
               x_f[1] + qt(1-alpha/2, n-1)*sqrt(S_f[1]/n))
ci_x2_uni <- c(x_f[2] - qt(1-alpha/2, n-1)*sqrt(S_f[4]/n),
               x_f[2] + qt(1-alpha/2, n-1)*sqrt(S_f[4]/n))
ci_x1_uni
ci_x2_uni

# Univariate confidence intervals with Bonferroni adjustment
ci_x1_bon <- c(x_f[1] - qt(1-alpha/4, n-1)*sqrt(S_f[1]/n),
               x_f[1] + qt(1-alpha/4, n-1)*sqrt(S_f[1]/n))
ci_x2_bon <- c(x_f[2] - qt(1-alpha/4, n-1)*sqrt(S_f[4]/n),
               x_f[2] + qt(1-alpha/4, n-1)*sqrt(S_f[4]/n))
ci_x1_bon
ci_x2_bon
#Bonferroni-square
x1 <- c(ci_x1_bon[1], ci_x1_bon[1], ci_x1_bon[2], ci_x1_bon[2], ci_x1_bon[1])
x2 <- c(ci_x2_bon[1], ci_x2_bon[2], ci_x2_bon[2], ci_x2_bon[1], ci_x2_bon[1])

# Bivariate confidence region based on Mahalanobis distance
d2_fun <-function(mu1, mu2, S, x_bar){
  S_inv <- solve(S)
  res <- S_inv[1,1]*(x_bar[1] - mu1)^2 +
    (S_inv[1, 2]+ S_inv[2, 1])*(x_bar[2]-mu2)*(x_bar[1] - mu1) +
    S_inv[2, 2]*(x_bar[2] - mu2)^2
  return(res)
}

# Add a constant for plotting regions
c <- 2
# Create a grid for the confidence ellipsoid to evaluate it with d2_fun
x1_seq <- seq(ci_x1_bon[1] - c, ci_x1_bon[2] + c, length = 100)
x2_seq <- seq(ci_x2_bon[1] - c, ci_x2_bon[2] + c, length = 100)

# Use diagonal Sigma from a) and S from b)
z_Sigma <- outer(x1_seq, x2_seq, d2_fun, S = Sigma_f, x_bar = x_f)
z_S <- outer(x1_seq, x2_seq, d2_fun, S = S_f, x_bar = x_f)

# Evaluate corresponding quantile functions
ellipse_Sigma <- contourLines(x1_seq, x2_seq, z_Sigma,
                              levels = 1/n*qchisq(1-alpha, p))
ellipse_S <- contourLines(x1_seq, x2_seq, z_S,
                          levels = (n-1)*p/((n-p)*n)*qf(1-alpha, p, n-p))
plot(x1, x2, # Bonferroni rectangle
     type = "l", xlab = "x1", ylab = "x2",
     xlim = c(ci_x1_bon[1] - c, ci_x1_bon[2] + c),
     ylim = c(ci_x2_bon[1] - c, ci_x2_bon[2] + c))
points(x_f[1], x_f[2], pch = 18) # sample mean
lines(ellipse_Sigma[[1]]$x, ellipse_Sigma[[1]]$y, col = "blue") # a)
lines(ellipse_S[[1]]$x, ellipse_S[[1]]$y, col = "red") # b)
points(20, 25) # mu_0

# d)
Sf <- matrix(c(17.683, 20.290, 20.290, 24.407), ncol = 2)
Sm <- matrix(c(18.479, 19.095, 19.095, 19.273), ncol = 2)
Sinv <- solve((34*Sf +13*Sm)/47)
Sinv
xm <- c(21.821, 22.843)
xf <- c(22.86, -24.397)
x <- (xf-xm)
T2 <- 35*14/49*t(x) %*% Sinv %*% x
46/94*T2
qf(0.95, 2, 46)
