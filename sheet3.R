library(mvtnorm)
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
