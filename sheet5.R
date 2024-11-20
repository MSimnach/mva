#### Sheet 5 -- Problem 6
load("income_share.Rdata")

## a) marginal and joint density plots
suppressWarnings(library(ggplot2))
library(gridExtra)
p1 <- ggplot(data = dat, aes(x = income, y = share)) +
  geom_point()
p2 <- ggplot(data = dat, aes(x = share)) +
  geom_density()
p3 <- ggplot(data = dat, aes(x = income)) +
  geom_density()
p4 <- ggplot(data = dat, aes(x = income, y = share)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                  colour = "white")
grid.arrange(p1, p2, p3, p4, ncol = 2)

# measure of dependency
energy::dcor(dat$income, dat$share)

## b) inadequacy of mvn
library(mvtnorm)
mvn_mean <- c(mean(dat$income), mean(dat$share))
mvn_cor <- cor(dat$income, dat$share)
mvn_cor
mvn_Sigma <- matrix(c(var(dat$income), mvn_cor*sd(dat$income)*sd(dat$share),
                      mvn_cor*sd(dat$income)*sd(dat$share), var(dat$share)),
                    ncol = 2)
sim_dat_mvn <- data.frame(rmvnorm(300, mean = mvn_mean, mvn_Sigma))
p1 + geom_point(data = sim_dat_mvn, mapping = aes(x = X1, y = X2),
                color = "red")
p5 <- ggplot(data = sim_dat_mvn, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                  colour = "white")
grid.arrange(p4,p5)

## c) method of moments to find appropriate marginal distributions
# Estimate Gamma distribution
income_mean <- mean(dat$income)
income_var <- var(dat$income)
income_rate <- income_mean / income_var
income_shape <- income_mean^2 /income_var
income_shape2 <- income_mean^2/income_var
income_scale <- income_var/income_mean
# Estimate Beta distribution
share_mean <- mean(dat$share)
share_var <- var(dat$share)
share_alpha <- share_mean*(share_mean*(1-share_mean)/share_var - 1)
share_beta <- (1-share_mean) * (share_mean*(1-share_mean)/share_var - 1)
p5 <- p2 + geom_function(fun = dbeta,
                         args = list(shape1 = share_alpha, shape2 = share_beta),
                         color = "red")
p6 <- p3 + geom_function(fun = dgamma,
                         args = list(shape = income_shape, rate = income_rate),
                         color = "red")
grid.arrange(p1, p5, p6, p4, ncol = 2)

# Copulas are a helpful approach to modeling the dependency structure for non-standard
# random vectors. Copulas are multivariate cumulative distribution functions defined on the
# unit cube for random vectors with uniform marginal distributions.

## d) uniform distributed pseudo-observations for income and share
# Pseudo observations
ps_obs <- data.frame(inc_ps = pgamma(dat[, 1], shape = income_shape,
                                     rate = income_rate),
                     sha_ps = pbeta(dat[, 2], shape1 = share_alpha,
                                    shape2 = share_beta))
p7 <- ggplot(ps_obs, aes(x = sha_ps))+
  geom_histogram(aes(y = ..density..), binwidth = 0.1, boundary = 0, color = "black") +
  stat_function(fun = function(x) 1, color = "red", linetype = "dashed") +
  labs(
    title = "Standardized Histogram with Uniform Density Reference",
    x = "Transformed Variable",
    y = "Density"
  ) +
  theme_minimal()
p8 <- ggplot(ps_obs, aes(x = inc_ps)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.1, boundary = 0, color = "black") +
  stat_function(fun = function(x) 1, color = "red", linetype = "dashed") +
  labs(
    title = "Standardized Histogram with Uniform Density Reference",
    x = "Transformed Variable",
    y = "Density"
  ) +
  theme_minimal()
p9 <- ggplot(ps_obs, aes(x = inc_ps, y = sha_ps)) +
  geom_point()
grid.arrange(p1, p7, p8, p9, ncol = 2)

### Alternative: nonparametric estimation of marginals
# Empirical CDF transformation
income_ecdf <- ecdf(dat$income)  # ECDF for income
share_ecdf <- ecdf(dat$share)   # ECDF for share

# Transform the data to uniform [0, 1] using the ECDF
income_uniform <- income_ecdf(dat$income)
share_uniform <- share_ecdf(dat$share)

# Combine into a data frame
transformed_data <- data.frame(
  income = dat$income,
  share = dat$share,
  income_uniform = income_uniform,
  share_uniform = share_uniform
)

# Histogram for transformed income
p10 <- ggplot(transformed_data, aes(x = income_uniform)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.1, boundary = 0, fill = "skyblue", color = "black") +
  labs(
    title = "Nonparametric Marginal Transformation: Income",
    x = "Transformed Income (Uniform [0, 1])",
    y = "Density"
  ) +
  theme_minimal()

# Histogram for transformed share
p11 <- ggplot(transformed_data, aes(x = share_uniform)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.1, boundary = 0,fill = "lightgreen", color = "black") +
  labs(
    title = "Nonparametric Marginal Transformation: Share",
    x = "Transformed Share (Uniform [0, 1])",
    y = "Density"
  ) +
  theme_minimal()
p12 <- ggplot(transformed_data, aes(x = income_uniform, y = share_uniform)) +
  geom_point()
grid.arrange(p1, p10, p11, p12, ncol = 2)

## e) plot independence copula and Gaussian copula for rho=0.3,0.9
par(mfrow = c(2, 3))
library(copula)
cop_ind <- indepCopula()
cop_gauss03 <- normalCopula(param = 0.3)
cop_gauss09 <- normalCopula(param = 0.9)
persp(cop_ind, pCopula)
persp(cop_gauss03, pCopula)
persp(cop_gauss09, pCopula)
contour(cop_ind, pCopula)
contour(cop_gauss03, pCopula)
contour(cop_gauss09, pCopula)

# Density of copula
persp(cop_ind, dCopula)
persp(cop_gauss03, dCopula)
persp(cop_gauss09, dCopula)
contour(cop_ind, dCopula)
contour(cop_gauss03, dCopula)
contour(cop_gauss09, dCopula)

## f) Use these copulas to construct joint distributions for income and share using the estimated
# marginal distributions. Use these joint distributions to generate new observations and
# compare them to the original data.
mv_dens_ind <- mvdc(copula = indepCopula(),
                    margins = c("gamma", "beta"),
                    paramMargins = list(list(shape = income_shape,
                                             rate = income_rate),
                                        list(shape1 = share_alpha,
                                             shape2 = share_beta)))
mv_dens_gauss03 <- mvdc(copula = normalCopula(param = 0.3),
                        margins = c("gamma", "beta"),
                        paramMargins = list(list(shape = income_shape,
                                                 rate = income_rate),
                                            list(shape1 = share_alpha,
                                                 shape2 = share_beta)))
mv_dens_gauss09 <- mvdc(copula = normalCopula(param = 0.9),
                        margins = c("gamma", "beta"),
                        paramMargins = list(list(shape = income_shape,
                                                 rate = income_rate),
                                            list(shape1 = share_alpha,
                                                 shape2 = share_beta)))
p13 <- ggplot() +
  geom_point(data = dat, mapping = aes(x = income, y = share)) +
  geom_point(data = data.frame(rMvdc(n = 300, mv_dens_ind)),
             mapping = aes(x = X1, y = X2), color = "red")
p14 <- ggplot() +
  geom_point(data = dat, mapping = aes(x = income, y = share)) +
  geom_point(data = data.frame(rMvdc(n = 300, mv_dens_gauss03)),
             mapping = aes(x = X1, y = X2), color = "red")
p15 <- ggplot() +
  geom_point(data = dat, mapping = aes(x = income, y = share)) +
  geom_point(data = data.frame(rMvdc(n = 300, mv_dens_gauss09)),
             mapping = aes(x = X1, y = X2), color = "red")
grid.arrange(p13, p14, p15, nrow = 1)

## g) Fit a Gaussian copula to the data using the function fitCopula() of the copula package.
# Fit copula
fit <- fitCopula(normalCopula(), as.matrix(ps_obs))
coef(fit)
## rho.1
## 0.8656992
# Joint density
mv_dens_est <- mvdc(copula = normalCopula(param = coef(fit)),
                    margins = c("gamma", "beta"),
                    paramMargins = list(list(shape = income_shape,
                                             rate = income_rate),
                                        list(shape1 = share_alpha,
                                             shape2 = share_beta)))
sim_dat <- data.frame(rMvdc(n = 300, mv_dens_est))
p13 <- p1 +geom_point(data = sim_dat, mapping = aes(x = X1, y = X2), color = "red")
# Joint density with all parameters estimated jointly
fit_mvdc <- fitMvdc(data = as.matrix(dat), mvdc = mv_dens_est,
                    start = c(unlist(mv_dens_est@paramMargins),
                              mv_dens_est@copula@parameters))
p14 <- p1 +
  geom_point(data = data.frame(rMvdc(n = 300, fit_mvdc@mvdc)),
             mapping = aes(x = X1, y = X2), color = "red")
grid.arrange(p13, p14, nrow = 1)
c(income_shape, income_rate, share_alpha, share_beta, coef(fit))
## rho.1
## 1.5857942 0.5050517 0.9047526 5.3570057 0.8656992
fit_mvdc@estimate
## shape rate shape1 shape2 rho.1
## 1.5838678 0.5107294 0.9251847 5.5585465 0.8634597

## h) Use the estimated joint distribution to find the expected share of income attributed to
# property assets for incomes between 4000 and 5000€.
sample_cond <- data.frame(rMvdc(n = 10000, mv_dens_est))
sample_cond$range <- sample_cond[, 1] > 4 & sample_cond[, 1] < 5
ggplot(sample_cond, aes(x = X1, y = X2, color = range)) +
  geom_point() +
  geom_vline(xintercept = c(4, 5), color = "blue") +
  scale_color_manual(values = c("black", "red"), guide = "none")
summary(sample_cond[sample_cond$range, 2])
