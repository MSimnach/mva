library(mvtnorm)
#### Problem 3
# Function for 3D plot of bivariate normal
binorm_plot <- function(sigma1, sigma2, rho, type = c("contour", "3d"),
                        phi = NULL, theta = NULL){
  if (is.null(phi) | is.null(theta)) type <- "contour"
  # Set mu and covariance matrix
  mu <- rep(0, 2)
  sigma <- matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2,
                    sigma2^2), ncol = 2, nrow = 2)
  # Set grid to compute density
  x1 <- x2 <- seq(-5, 5, 0.25)
  grid <- expand.grid(x1, x2)
  # Compute
  fvals <- matrix(dmvnorm(grid, mean = mu, sigma = sigma),
                  ncol = length(x1), nrow = length(x2))
  if (type == "contour") {
    # Contour plot
    contour(x = x1, y = x2, z = fvals, xlab = "x_1", ylab = "x_2")
    title(main = paste("sigma1=", as.character(sigma1), ", sigma2=",
                       as.character(sigma2), ", rho=", as.character(rho)))
  } else {
    # 3D plot of density
    persp(x = x1, y = x2, z = fvals, ticktype = "detailed", xlab = "x_1",
          ylab = "x_2", phi = phi, theta = theta)
    title(main = paste("sigma1=", as.character(sigma1), ", sigma2=",
                       as.character(sigma2), ", rho=", as.character(rho)))
  }
}

par(mfrow=c(1,2))
binorm_plot(sigma1 = 1, sigma2 = 1, rho = 0.9)
binorm_plot(sigma1 = 1, sigma2 = 1, rho = 0.9, type = "3d", phi = 35,
            theta = 120)

library(manipulate)
manipulate(binorm_plot(sigma1, sigma2, rho, type = "3d", phi, theta),
           sigma1 = slider(0.5, 4, step = 0.5),
           sigma2 = slider(0.5, 4, step = 0.5),
           rho = slider(-.9, .9, step = 0.1),
           phi = slider(0, 90, step = 10),
           theta = slider(0, 90, step = 10))
manipulate(binorm_plot(sigma1, sigma2, rho),
           sigma1 = slider(0.5, 4, step = 0.5),
           sigma2 = slider(0.5, 4, step = 0.5),
           rho = slider(-.9, .9, step = 0.1))


### Problem 4
suppressWarnings(library(ggplot2))
library(gridExtra)
# Function for 3D plot of bivariate normal (joint, marginal, conditional)
binorm_jmc_plot <- function(x, sigma_y, sigma_x, rho){
  # Set mu and covariance matrix
  mu <- rep(0, 2)
  Sigma <- matrix(c(sigma_y^2, rho*sigma_y*sigma_x, rho*sigma_y*sigma_x,
                    sigma_x^2), ncol = 2, nrow = 2)
  # Set grid to compute densities
  x1 <- x2 <- seq(-5, 5, 0.1)
  joint <- data.frame(expand.grid(x1, x2))
  joint$values <- dmvnorm(joint, mean = mu, sigma = Sigma)
  margin <- data.frame(x = x2, values = dnorm(x2, mean = 0, sd = sigma_x))
  mu_c <- rho*sigma_y/sigma_x*x
  var_c <- sigma_y^2*(1-rho^2)
  condi <- data.frame(x = x1, values_m = dnorm(x1, mean = 0, sd = sigma_y),
                      values_c = dnorm(x1, mean = mu_c, sd = sqrt(var_c)))
  p1 <- ggplot(joint, aes(x = Var2, y = Var1, z = values)) +
    geom_raster(aes(fill = values)) +
    geom_contour(colour = "white") +
    geom_vline(xintercept = x, colour = "red") +
    geom_abline(intercept = 0, slope = rho*sigma_y/sigma_x, colour = "orange") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "x", y = "y")
  p2 <- ggplot(margin, aes(x = x, y = values)) +
    geom_line() +
    geom_vline(xintercept = x, colour = "red") +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1)) +
    xlab("x") + ylab("f(x)")
  p3 <- ggplot(condi, aes(x = x)) +
    geom_line(aes(y = values_m)) +
    geom_line(aes(y = values_c), lty = 2) +
    geom_vline(xintercept = mu_c, colour = "orange") +
    geom_vline(xintercept = mu_c, colour = "red", lty = "ff") +
    theme_bw() +
    xlab("y") + ylab("f(y), f(y|x)") +
    scale_y_continuous(limits = c(0, 1)) +
    coord_flip()
  grid.arrange(p1, p2, p3, layout_matrix = matrix(c(1, 2, 3, NA), nrow = 2))
}
binorm_jmc_plot(x = -2, sigma_y = 2.5, sigma_x = 2.5, rho = 0.6)
binorm_jmc_plot(x = 3.1, sigma_y = 2.5, sigma_x = 2.5, rho = 0.6)
binorm_jmc_plot(x = 3.1, sigma_y = 2.5, sigma_x = 2.5, rho = 0)
binorm_jmc_plot(x = 0, sigma_y = 0.8, sigma_x = 1, rho = 0.8)
