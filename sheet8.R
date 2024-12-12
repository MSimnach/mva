### Sheet 8
# Problem 18 Correspondence Analysis
library(ade4)
data(housetasks)
# a) chisq test
ht <- as.matrix(housetasks)
xl <- rowSums(ht)
xj <- colSums(ht)
x <- sum(ht)
E <- xl%*%t(xj)/x
C <- (ht - E)/sqrt(E)
t <- sum(C^2)
df <- (nrow(ht) - 1)*(ncol(ht) - 1)
chisq <- qchisq(0.95, df)
t>chisq
chisq.test(ht)
# NOTE
lambda_sqrt <- svd(C)$d
lambda <- lambda_sqrt^2
sum(lambda)

# b) Correspondence Analysis
cumsum(lambda_sqrt^2)/sum(lambda_sqrt^2)
gamma2 <- svd(C)$u[, 1:2]
delta2 <- svd(C)$v[, 1:2]
A <- diag(c(xl^(-1/2)))
B <- diag(c(xj^(-1/2)))
r12 <- A%*%C%*%delta2
s12 <- B%*%t(C)%*%gamma2

library(FactoMineR)
house_ca <- CA(housetasks)
house_ca$col$coord
house_ca$row$coord
library(ggplot2)
dat <- data.frame(rbind(r12, s12),
                  Factor = c(rep("Row", nrow(r12)), rep("Column", nrow(s12))),
                  name = c(rownames(ht), colnames(ht)))
ggplot(data = dat, aes(x = X1, y = X2, color = Factor, label = name)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_text() +
  labs(x = paste0("Dim1 (",
                  round(lambda_sqrt[1]^2/sum(lambda_sqrt^2), 3)*100,
                  "%)"),
       y = paste0("Dim2 (",
                  round(lambda_sqrt[2]^2/sum(lambda_sqrt^2), 3)*100,
                  "%)")) +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  theme_bw()
