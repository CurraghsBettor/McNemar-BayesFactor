rm(list=ls())

set.seed(12345)

BFMcNemar <- function(data, hypothesis = c("null", "alternative")) {
  if(is.matrix(data)) {
    if(nrow(data) < 2 || ncol(data) != nrow(data)) {
      stop(" 'data' has to be a 2x2 matrix, and the number of rows has to match the number of columns")
    } else if (nrow(data) == 3 && ncol(data) == nrow) {
      stop("'data' has to be a 2x2 matrix")
    } else if (any(data < 0) || anyNA(data)) {
      stop("all entries of 'data' have to be nonnegative and finite")
    } else {
      n <- sum(data[1,2], data[2,1])
      BF01num <- gamma(sum(n, 2))
      BF01denom <- 2^n*gamma(sum(data[1,2], 1))*gamma(sum(data[2,1], 1))
      if (hypothesis != "alternative") {
        BF01 <- BF01num/BF01denom
      } else {
        BF01 <- BF01num/BF01denom
        BF10 <- 1/BF01
      }
      data <- deparse(substitute(data))
      method <- "McNemar-Bayes Factor"
      if (hypothesis != "alternative") {
        names(BF01) <- "Bayes Factor Null: pi12/pi12+pi21 = 0.5 vs Alternative: pi12/pi12+pi21 != 0.5"
        Bayes <- list(data = data, method = method, BF01 =  BF01)
        return(Bayes)
      } else {
        names(BF10) <- "Bayes Factor Alternative: pi12/pi12+pi21 != 0.5 vs Null: pi12/pi12+pi21 = 0.5"
        Bayes <- list(data = data, method = method, BF10 = BF10)
        return(Bayes)
      }
    }
  } else {
    stop("'data' has to be converted into matrix")
  }
}


BFMNdirichlet <- function(data, a12, a21, hypothesis = c("null", "alternative")) {
  if(is.matrix(data)) {
    if(nrow(data) < 2 || ncol(data) != nrow(data)) {
      stop(" 'data' has to be a 2x2 matrix, and the number of rows has to match the number of columns")
    } else if (nrow(data) == 3 && ncol(data) == nrow) {
      stop("'data' has to be a 2x2 matrix")
    } else if (any(data < 0) || anyNA(data)) {
      stop("all entries of 'data' have to be nonnegative and finite")
    } else {
      n <- sum(data[1,2], data[2,1])
      BFnum <- gamma(a12)*gamma(a21)*gamma(sum(n, a12, a21))
      BFdenom <- 2^n*gamma(sum(a12, a21))*gamma(sum(data[1,2], a12))*gamma(sum(data[2,1], a21))
      if (hypothesis != "alternative") {
        BFHsHa <- BFnum/BFdenom
      } else {
        BFHsHa <- BFnum/BFdenom
        BFHaHs <- 1/BFHsHa
      }
      data <- deparse(substitute(data))
      method <- "McNemar-Bayes Factor"
      if (hypothesis != "alternative") {
        names(BFHsHa) <- "Bayes Factor Null - Hs: pi12 = pi21 vs Alternative - Ha: pi12 != pi21"
        Bayes <- list(data = data, method = method, BFHsHa =  BFHsHa)
        return(Bayes)
      } else {
        names(BFHaHs) <- "Bayes Factor Alternative - Ha: pi12 != pi21 vs Null - Hs: pi12 = pi21"
        Bayes <- list(data = data, method = method, BFHaHs = BFHaHs)
        return(Bayes)
      }
    }
  } else {
    stop("'data' has to be converted into matrix")
  }
}

## example
dataset <- matrix(c(88, 22, 37, 53),
                  nrow = 2, byrow = T); print(dataset)

a <- BFMcNemar(dataset, "alternative")$BF10; a ## att

## example from Kateri et al. (2001)
ex1 <- matrix(c(2, 3, 13, 2),
              nrow = 2, byrow = T); print(ex1)
ex2 <- matrix(c(42, 3, 13, 42),
              nrow = 2, byrow = T); print(ex2)

# Dirichlet prior parameters (alpha12, a21, a1) = (1,1,1)
BFMNdirichlet(ex1, 1, 1, hypothesis = "null")
BFMNdirichlet(ex2, 1, 1, hypothesis = "null")

# Dirichlet prior parameters (alpha12, a21, a1) = (0.5,0.5,0.5)
BFMNdirichlet(ex1, 0.5, 0.5, hypothesis = "null")
BFMNdirichlet(ex2, 0.5, 0.5, hypothesis = "null")

# Dirichlet prior parameters (alpha12, a21, a1) = (0.25,0.25,0.25)
BFMNdirichlet(ex1, 0.25, 0.25, hypothesis = "null")
BFMNdirichlet(ex2, 0.25, 0.25, hypothesis = "null")