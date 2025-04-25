# Runge Kutta based Search Mechanism
RungeKutta <- function(XB, XW, DelX)
{
  dim <- length(XB)
  C <- sample(1:2, 1) * (1 - runif(1))
  r1 <- runif(dim)
  r2 <- runif(dim)
  
  K1 <- 0.5 * (runif(dim) * XW - C * XB)
  K2 <- 0.5 * (runif(dim) * (XW + r2 * K1 * DelX / 2) - (C * XB + r1 * K1 * DelX / 2))
  K3 <- 0.5 * (runif(dim) * (XW + r2 * K2 * DelX / 2) - (C * XB + r1 * K2 * DelX / 2))
  K4 <- 0.5 * (runif(dim) * (XW + r2 * K3 * DelX) - (C * XB + r1 * K3 * DelX))
  
  XRK <- (K1 + 2 * K2 + 2 * K3 + K4)
  SM <- 1 / 6 * XRK
  return(SM)
}

# Random Initialization
initialization <- function(nP, dim, ub, lb)
{
  if (length(ub) == 1) {
    X <- matrix(runif(nP * dim), nrow = nP) * (ub - lb) + lb
  } else {
    X <- matrix(0, nP, dim)
    for (i in 1:dim) {
      X[, i] <- runif(nP) * (ub[i] - lb[i]) + lb[i]
    }
  }
  return(X)
}

# Uniform Random
Unifrnd <- function(a, b, c, dim)
{
  a2 <- a / 2
  b2 <- b / 2
  mu <- a2 + b2
  sig <- b2 - a2
  z <- mu + sig * (2 * matrix(runif(c * dim), c, dim) - 1)
  return(z)
}

# Random indices
RndX <- function(nP, i)
{
  Qi <- sample(setdiff(1:nP, i))
  return(Qi[1:3])
}

# Main RUN optimizer
RUN <- function(nP, MaxIt, lb, ub, fobj)
{
  Cost <- rep(0, nP)
  dim <- length(lb)
  X <- initialization(nP, dim, ub, lb)
  Convergence_curve <- numeric(MaxIt)
  
  for (i in 1:nP) Cost[i] <- fobj(X[i, ])
  Best_Cost <- min(Cost)
  Best_X <- X[which.min(Cost), ]
  Convergence_curve[1] <- Best_Cost
  
  for (it in 2:MaxIt)
  {
    f <- 20 * exp(-12 * (it / MaxIt))
    Xavg <- colMeans(X)
    SF <- 2 * (0.5 - runif(nP)) * f
    
    for (i in 1:nP)
    {
      ind_l <- which.min(Cost)
      lBest <- X[ind_l, ]
      idx <- RndX(nP, i)
      A <- idx[1]; B <- idx[2]; C <- idx[3]
      ind1 <- which.min(Cost[c(A, B, C)])
      ind1 <- c(A, B, C)[ind1]
      
      gama <- runif(dim) * (X[i, ] - runif(dim) * (ub - lb)) * exp(-4 * it / MaxIt)
      Stp <- runif(dim) * ((Best_X - runif(dim) * Xavg) + gama)
      DelX <- 2 * runif(dim) * abs(Stp)
      
      if (Cost[i] < Cost[ind1])
      {
        Xb <- X[i, ]; Xw <- X[ind1, ]
      } else {
        Xb <- X[ind1, ]; Xw <- X[i, ]
      }
      
      SM <- RungeKutta(Xb, Xw, DelX)
      L <- runif(dim) < 0.5
      Xc <- ifelse(L, X[i, ], X[A, ])
      Xm <- ifelse(L, Best_X, lBest)
      
      vec <- c(1, -1)
      r <- sample(vec, dim, replace = TRUE)
      g <- 2 * runif(1)
      mu <- 0.5 + 0.1 * rnorm(dim)
      
      if (runif(1) < 0.5)
      {
        Xnew <- (Xc + r * SF[i] * g * Xc) + SF[i] * SM + mu * (Xm - Xc)
      } else {
        Xnew <- (Xm + r * SF[i] * g * Xm) + SF[i] * SM + mu * (X[A, ] - X[B, ])
      }
      
      Xnew <- pmin(pmax(Xnew, lb), ub)
      CostNew <- fobj(Xnew)
      
      if (CostNew < Cost[i])
      {
        X[i, ] <- Xnew
        Cost[i] <- CostNew
      }
      
      # Enhanced solution quality (ESQ)
      if (runif(1) < 0.5)
      {
        EXP <- exp(-5 * runif(1) * it / MaxIt)
        r_val <- sample(-1:1, 1)
        u <- 2 * runif(dim)
        w <- Unifrnd(0, 2, 1, dim) * EXP
        
        idx <- RndX(nP, i)
        Xavg_esq <- colMeans(X[idx, ])
        beta <- runif(dim)
        Xnew1 <- beta * Best_X + (1 - beta) * Xavg_esq
        
        Xnew2 <- numeric(dim)
        for (j in 1:dim)
        {
          if (w[j] < 1)
          {
            Xnew2[j] <- Xnew1[j] + r_val * w[j] * abs((Xnew1[j] - Xavg_esq[j]) + rnorm(1))
          } else {
            Xnew2[j] <- (Xnew1[j] - Xavg_esq[j]) + r_val * w[j] * abs((u[j] * Xnew1[j] - Xavg_esq[j]) + rnorm(1))
          }
        }
        
        Xnew2 <- pmin(pmax(Xnew2, lb), ub)
        CostNew <- fobj(Xnew2)
        
        if (CostNew < Cost[i])
        {
          X[i, ] <- Xnew2
          Cost[i] <- CostNew
        } else if (runif(1) < w[sample(1:dim, 1)]) {
          SM <- RungeKutta(X[i, ], Xnew2, DelX)
          Xnew <- (Xnew2 - runif(dim) * Xnew2) + SF[i] * (SM + (2 * runif(dim) * Best_X - Xnew2))
          Xnew <- pmin(pmax(Xnew, lb), ub)
          CostNew <- fobj(Xnew)
          if (CostNew < Cost[i])
          {
            X[i, ] <- Xnew
            Cost[i] <- CostNew
          }
        }
      }
      
      if (Cost[i] < Best_Cost)
      {
        Best_X <- X[i, ]
        Best_Cost <- Cost[i]
      }
    }
    
    Convergence_curve[it] <- Best_Cost
    cat(sprintf("It %d | Best Cost: %.4e\n", it, Best_Cost))
  }
  return(list(Best_Cost = Best_Cost, Best_X = Best_X, Convergence_curve = Convergence_curve))
}
