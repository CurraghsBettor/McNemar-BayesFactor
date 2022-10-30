## by normal law
rm(list = ls())

Ztest <- function(data, alternative = c("less", "alternative", "two.sided")) {
  if(is.matrix(data)) {
    if(nrow(data) < 2 || ncol(data) != nrow(data)) {
      stop(" 'data' has to be a 2x2 matrix, and the number of rows have to match the number of columns")
    } else if (nrow(data) == 3 && ncol(data) == nrow) {
      stop("'data' has to be a 2x2 matrix, in the case of 3x3 matrix use the McNemar-Bowker exact test")
    } else if (any(data < 0) || anyNA(data)) {
      stop("all entries of 'data' have to be nonnegative and finite")
    } else {
      b <- data[1,2]
      n <- sum(data[1,2], data[2,1])
      Z <- (b - (n/2))/(sqrt(n)/2)
      palpha <- 0.05
      Zalpha <- qnorm(1 - palpha, mean = 0, sd = 1)
      ## palpha <- 1 - pnorm(1.64, mean = 0, sd = 1)
      if (alternative == "greater") {
        pvaleur <- integrate(dnorm, lower = Z, upper = Inf)$value
      } else if (alternative == "less") {
        pvaleur <- integrate(dnorm, lower = -Inf, upper = Z)$value
      } else if (alternative == "two.sided") {
        pvaleur <- integrate(dnorm, lower = abs(Z), upper = Inf)$value
        pvaleur <- 2*pvaleur
      }
      ## pvaleur <- 1 - pnorm(Z, mean = 0, sd = 1)
      data <- deparse(substitute(data))
      method <- "McNemar's z-test"
      names(Z) <- "McNemar's z score"
      names(pvaleur) <- "p-value"
      structure(list(statistic = Z, p.value = pvaleur, method = method, alternative = alternative, data.name = data), class = "htest")
    }
  } else {
    stop("'data' has to be converted into matrix")
  }
}

## Example
dataset <- matrix(c(88, 37, 22, 53), 
                     nrow = 2, byrow = T); print(dataset)

## McNemar chi-square test
mcnemar.test(dataset, correct = F)

## McNemar exact test
library(exact2x2)
mcnemar.exact(dataset)

## exact binomial test = McNemar exact test
binom.test(dataset[1,2], sum(dataset[1,2], dataset[2,1]), p = 0.5)

## Ztest

Ztest(dataset, alternative = "two.sided")
## extract the statistical value
stats <- Ztest(dataset, alternative ="greater")$statistic; print(stats)
## extract the p-value 
p <- Ztest(dataset, alternative = "greater")$p.value
2*p
test <- Ztest(dataset, alternative = "greater"); print(test)

