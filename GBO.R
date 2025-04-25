GradientSearchRule<-function(ro1, Best_X, Worst_X, X, Xr1, DM, eps, Xm, Flag)
{
    nV <- length(X)
    Delta <- 2 * runif(1) * abs(Xm - X)
    Step <- ((Best_X - Xr1) + Delta) / 2
    DelX <- runif(nV) * abs(Step)
    GSR <- rnorm(1) * ro1 * (2 * DelX * X) / (Best_X - Worst_X + eps)
    if (Flag == 1) {
        Xs <- X - GSR + DM
    } else {
        Xs <- Best_X - GSR + DM
    }
    yp <- runif(1) * (0.5 * (Xs + X) + runif(nV) * DelX)
    yq <- runif(1) * (0.5 * (Xs + X) - runif(nV) * DelX)
    GSR <- rnorm(1) * ro1 * (2 * DelX * X) / (yp - yq + eps)
    return(GSR)
}

initialization<-function(nP, dim, ub, lb)
{
    X <- matrix(0, nrow = nP, ncol = dim)
    for (i in 1:dim) {
        X[, i] <- runif(nP, lb[i], ub[i])
    }
    return(X)
}

GBO<-function(nP, MaxIt, lb, ub, dim, fobj)
{
    # Initialization
    pr <- 0.5
    lb <- rep(lb, dim)
    ub <- rep(ub, dim)
    X <- initialization(nP, dim, ub, lb)
    Cost <- apply(X, 1, fobj)
    
    Best_Cost <- min(Cost)
    Best_X <- X[which.min(Cost), ]
    Worst_Cost <- max(Cost)
    Worst_X <- X[which.max(Cost), ]
    
    Convergence_curve <- numeric(MaxIt)
    Best_Solutions <- matrix(0, nrow = MaxIt, ncol = dim)
    
    for (it in 1:MaxIt) {
        beta <- 0.2 + (1.2 - 0.2) * (1 - (it / MaxIt)^3)^2
        alpha <- abs(beta * sin((3 * pi / 2 + sin(3 * pi / 2 * beta))))
        
        for (i in 1:nP) {
            A1 <- sample(1:nP, 4, replace = TRUE)
            r1 <- A1[1]; r2 <- A1[2]; r3 <- A1[3]; r4 <- A1[4]
            
            Xm <- (X[r1, ] + X[r2, ] + X[r3, ] + X[r4, ]) / 4
            ro <- alpha * (2 * runif(1) - 1)
            ro1 <- alpha * (2 * runif(1) - 1)
            eps <- 5e-3 * runif(1)
            
            DM <- runif(1) * ro * (Best_X - X[r1, ])
            GSR <- GradientSearchRule(ro1, Best_X, Worst_X, X[i, ], X[r1, ], DM, eps, Xm, 1)
            DM <- runif(1) * ro * (Best_X - X[r1, ])
            X1 <- X[i, ] - GSR + DM
            
            DM <- runif(1) * ro * (X[r1, ] - X[r2, ])
            GSR <- GradientSearchRule(ro1, Best_X, Worst_X, X[i, ], X[r1, ], DM, eps, Xm, 2)
            DM <- runif(1) * ro * (X[r1, ] - X[r2, ])
            X2 <- Best_X - GSR + DM
            
            Xnew <- numeric(dim)
            for (j in 1:dim) {
                ro <- alpha * (2 * runif(1) - 1)
                X3 <- X[i, j] - ro * (X2[j] - X1[j])
                ra <- runif(1); rb <- runif(1)
                Xnew[j] <- ra * (rb * X1[j] + (1 - rb) * X2[j]) + (1 - ra) * X3
            }
            
            # Local escaping operator
            if (runif(1) < pr) {
                k <- sample(1:nP, 1)
                f1 <- -1 + 2 * runif(1)
                f2 <- -1 + 2 * runif(1)
                ro <- alpha * (2 * runif(1) - 1)
                Xk <- runif(dim, min = lb, max = ub)
                
                L1 <- runif(1) < 0.5
                u1 <- if (L1) 2 * runif(1) else 1
                u2 <- if (L1) runif(1) else 1
                u3 <- if (L1) runif(1) else 1
                L2 <- runif(1) < 0.5
                
                Xp <- if (L2) Xk else X[k, ]
                if (u1 < 0.5) {
                    Xnew <- Xnew + f1 * (u1 * Best_X - u2 * Xp) + f2 * ro * (u3 * (X2 - X1) + u2 * (X[r1, ] - X[r2, ])) / 2
                } else {
                    Xnew <- Best_X + f1 * (u1 * Best_X - u2 * Xp) + f2 * ro * (u3 * (X2 - X1) + u2 * (X[r1, ] - X[r2, ])) / 2
                }
            }
            
            # Boundary control
            Xnew <- pmin(pmax(Xnew, lb), ub)
            Xnew_Cost <- fobj(Xnew)
            
            if (Xnew_Cost < Cost[i]) {
                X[i, ] <- Xnew
                Cost[i] <- Xnew_Cost
                if (Xnew_Cost < Best_Cost) {
                    Best_X <- Xnew
                    Best_Cost <- Xnew_Cost
                }
            }
            
            if (Cost[i] > Worst_Cost) {
                Worst_X <- X[i, ]
                Worst_Cost <- Cost[i]
            }
        }
        
        Convergence_curve[it] <- Best_Cost
        Best_Solutions[it, ] <- Best_X
        cat(sprintf("Iteration %d: Best Fitness = %e\n", it, Best_Cost))
    }
    
    return(list(
        Best_Cost = Best_Cost,
        Best_X = Best_X,
        Convergence_curve = Convergence_curve,
        Best_Solutions = Best_Solutions
    ))
}
