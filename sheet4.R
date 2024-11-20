n<-1000
set.seed(12)
x<-rnorm(n)
z<-x^2+rnorm(n)
plot(x,z)
cor(x,z)

# Create distance matrix
A <- matrix(NA, nrow = n, ncol = n)
for (j in seq_along(x)) { A[, j] <- abs(x- x[j]) }
A[1:6, 1:6]

H <- diag(n)- matrix(1, ncol = n, nrow = n)/n
A_tilde <- H %*% A %*% H
# A_tilde is double centered
head(rowSums(A_tilde))
head(colSums(A_tilde))
sum(A_tilde)

# Double centering and distance correlation implemented in package 'energy'
library(energy)
dimnames(A_tilde) <- list(seq_len(n), seq_len(n))
all.equal(A_tilde, Dcenter(x))
?Dcenter

# Distance correlation by multiplication of matrix elements
B_tilde <- Dcenter(z)
D_cov2_xz <- (1/n)^2 * sum(A_tilde * B_tilde)
D_var2_x <- (1/n)^2 * sum(A_tilde * A_tilde)
D_var2_z <- (1/n)^2 * sum(B_tilde * B_tilde)
sqrt(D_cov2_xz) / sqrt(sqrt(D_var2_x) * sqrt(D_var2_z))

dcor(x, z)
cor(x, z)
