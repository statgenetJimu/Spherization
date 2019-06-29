library(mvtnorm)
library(TruncatedNormal)

# Calculates marginal counts of a contingency table "A".
calc.marg<-function (A){
    d <- dim(A)
    ret <- list()
    for (i in 1:length(d)) {
        ret[[i]] <- apply(A, i, sum)
    }
    ret
}

# Calculates expected table of a contingency table "A".
make.exp.table <- function(A){
  n <- sum(A)
  marg <- calc.marg(A)
  tmp <- marg[[1]]
  for(i in 2:length(marg)){
    tmp <- t(marg[[i]]/n) %x% tmp
  }
  tmp
}

# Calculates a k × k simplex-rotation matrix for an integer k.
# Example: For a 2×3 table, k = 2,3. For each k, this function returns two rotation matrices with dimensions, 2×2 and 3×3, respectively.
make.simplex <- function(k){
  ret <- matrix(0,k,k)
  for(i in 1:(k-1)){
    for(j in 1:k){
      if(j < i){
      }else if(j==i){
        ret[i,j] <-  sqrt((k-i)/(k-i+1))
      }else{
        ret[i,j] <- -sqrt(1/((k-i)*(k-i+1)))
      }
    }
  }
  ret[k,] <- sqrt(1/k)
  ret
}

# Calculates the Kronecker product of multiple simplex-rotation matrices.
# r: a vector of numbers of the levels of the variables of a contingency table.
# Example: For a 2×3 table, r = c(2,3), this function returns a 6×6 matrix.
make.simplex.multi <- function(r){
  X <- make.simplex(r[1])
  k <- length(r)
  if(k > 1){
    for(i in 2:k){
      X <- make.simplex(r[i]) %x% X
    }
  }
  X
}

# Calculates the index vector to define the positions of zero and non-zero elements produced from the rotation.
# Example: For a 2×3 table, d = c(2,3).
arrive.index <- function(d){
  n <- prod(d-1)
  x <- numeric(n)
  ind <- 1
    for(i in 1:(d[2]-1)){
      for(j in 1:(d[1]-1)){
      x[ind] <- (i-1)*d[1]+j
      ind <- ind+1
    }
  }
  x
}

# Spherization of table.
# Returns the matrices that will transfer test vectors to the df-dimensional vectors. 
calc.rotate <- function(table){
  dim <- dim(table)
  exp.table <- make.exp.table(table)
  e.vec <- c(exp.table)
  
  R <- make.simplex.multi(dim)
  Einv <- diag(1/e.vec)
  d.vec <- c(table)-e.vec
  
  Z <- diag(prod(dim))[arrive.index(dim),]
  XtX <- Z %*% R %*% Einv %*% t(R) %*% t(Z)
  
  eigenout <- eigen(XtX)
  X <- diag((eigenout[[1]])^0.5) %*% solve(eigenout[[2]])
  
  P <- X %*% Z %*% R

  Pinv <- t(R) %*% t(Z) %*% solve(X)
  return(list(P=P,Pinv=Pinv))
}

# Calculates the df-dimensional test vectors.
# tests: a list of tables representing test models.
# rotation: the spherization information from the “calc.rotate” function.
make.test.vecs <- function(rotation,tests){
  test.vecs <- tests %*% rotation
  L.test.vecs <- sqrt(apply(test.vecs^2,1,sum))
  test.vecs / L.test.vecs
}


pmway.table.null.botev <- function(stat,table,tests,lower.tail,one.side,n){
  d <- dim(table)
  df <- prod(d-1)
  rotate <- calc.rotate(table)$Pinv
  test.vecs <- make.test.vecs(rotate,tests)
  
  
  sig <- test.vecs %*% t(test.vecs)
  u <- rep(stat,nrow(sig))
  pr <- 1 - mvNcdf(-u,u,sig,n)$prob
  
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr)/2
    }else{
      ret <- pr
    }
  }else{
    if(one.side){
      ret <- (1 - pr)/2
    }else{
      ret <- 1 - pr
    }
  }
  ret
}

# Calculates the p-values for the Sph-Btv method.
# stat: the χ^2-value(s) of the table from the proportion trend test.
# table: the contingency table.
# lower.tail: logical; if TRUE, probabilities are P[X≤x], otherwise, P[X>x].
# one.side: logical; if TRUE, the test vectors indicate one-sided test, otherwise, two-sided.
# n: the number of Monte Carlo samples.
pmway.table.null.multi.botev <- function(stat,table,tests,lower.tail,one.side,n){
  d <- dim(table)
  df <- prod(d-1)
  rotate <- calc.rotate(table)$Pinv
  test.vecs <- make.test.vecs(rotate,tests)
  
  
  sig <- test.vecs %*% t(test.vecs)
  
  pr <- numeric(length(stat))
  for(i in 1:length(stat)){
    u <- rep(stat[i],nrow(sig))
    pr[i] <- 1 - mvNcdf(-u,u,sig,n)$prob
  }
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr)/2
    }else{
      ret <- pr
    }
  }else{
    if(one.side){
      ret <- (1 - pr)/2
    }else{
      ret <- 1 - pr
    }
  }
  ret
}


pmway.table.botev <- function(stat,table,tests,lower.tail,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  
  d <- dim(table)
  df <- prod(d-1)
  PPinv <- calc.rotate(table)
  rotate <- PPinv$Pinv
  test.vecs <- make.test.vecs(rotate,tests)
  
 alt.vec <- PPinv$P %*% c(table-make.exp.table(table))
  
  sig <- test.vecs %*% t(test.vecs)
  u <- rep(stat,nrow(sig))
  pr.null <- 1 - mvNcdf(-u,u,sig,n)$prob
  
  pr.alt <- 1 - mvNcdf(-u - alt.vec, u - alt.vec, sig, n)$prob
  
  pr <- c(pr.null,pr.alt)
  
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr)/2
    }else{
      ret <- pr
    }
  }else{
    if(one.side){
      ret <- (1 - pr)/2
    }else{
      ret <- 1 - pr
    }
  }
  ret
}


qmway.table.botev <- function(p,table,tests,lower.tail,range,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(range))range <- c(1e-4,100)
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  
  lx <- range[1]
  ux <- range[2]
  
  if(lower.tail){
    lower <- pmway.table.null.botev(lx,table,tests,lower.tail,one.side,n)
    while(lower < p){
      lx <- lx*0.5
      lower <- pmway.table.null.botev(lx,table,tests,lower.tail,one.side,n)
    }
    upper <- pmway.table.null.botev(ux,table,tests,lower.tail,one.side,n)
    while(upper > p){
      ux <- ux*2
      upper <- pmway.table.null.botev(ux,table,tests,lower.tail,one.side,n)
    }
    for(it in 1:10){
      hx <- (l+u)/2
      half <- pmway.table.null.botev(hx,table,tests,lower.tail,one.side,n)
      if(half < p){
        ux <- hx
      }
      else{
        lx <- hx
      }
    }
  }
  else{
    lower <- pmway.table.null.botev(lx,table,tests,lower.tail,one.side,n)
    while(lower > p){
      lx <- lx*0.5
      lower <- pmway.table.null.botev(lx,table,tests,lower.tail,one.side,n)
    }
    upper <- pmway.table.null.botev(ux,table,tests,lower.tail,one.side,n)
    while(upper < p){
      ux <- ux*2
      upper <- pmway.table.null.botev(ux,table,tests,lower.tail,one.side,n)
    }
    
    for(it in 1:10){
      hx <- (lx+ux)/2
      half <- pmway.table.null.botev(hx,table,tests,lower.tail,one.side,n)
      if(half < p){
        lx <- hx
      }
      else{
        ux <- hx
      }
    }
  }
  
  return (l+u)/2
}

pmway.table.botev.sig <- function(stat,sig,lower.tail,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  
  
  u <- rep(stat,nrow(sig))
  pr <- 1 - mvNcdf(-u,u,sig,100)$prob
  
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr)/2
    }else{
      ret <- pr
    }
  }else{
    if(one.side){
      ret <- (1 - pr)/2
    }else{
      ret <- 1 - pr
    }
  }
  ret
}

qmway.table.botev.lite <- function(p,table,tests,lower.tail,x,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(x))x <- c(1e-4,100)
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  
  d <- dim(table)
  df <- prod(d-1)
  rotate <- calc.rotate(table)$Pinv
  test.vecs <- make.test.vecs(rotate,tests)
  
  sig <- test.vecs %*% t(test.vecs)
  
  l <- x[1]
  u <- x[2]
  
  if(lower.tail){
    lower <- pmway.table.botev.sig(l,sig,lower.tail,one.side,n)
    while(lower < p){
      l <- l*0.5
      lower <- pmway.table.botev.sig(l,sig,lower.tail,one.side,n)
    }
    upper <- pmway.table.botev.sig(u,sig,lower.tail,one.side,n)
    while(upper > p){
      u <- u*2
      upper <- pmway.table.botev.sig(u,sig,lower.tail,one.side,n)
    }
    for(it in 1:10){
      h <- (l+u)/2
      half <- pmway.table.botev.sig(h,sig,lower.tail,one.side,n)
      if(half < p){
        u <- h
      }
      else{
        l <- h
      }
    }
  }
  else{
    lower <- pmway.table.botev.sig(l,sig,lower.tail,one.side,n)
    while(lower > p){
      l <- l*0.5
      lower <- pmway.table.botev.sig(l,sig,lower.tail,one.side,n)
    }
    upper <- pmway.table.botev.sig(u,sig,lower.tail,one.side,n)
    while(upper < p){
      u <- u*2
      upper <- pmway.table.botev.sig(u,sig,lower.tail,one.side,n)
    }
    
    for(it in 1:10){
      h <- (l+u)/2
      half <- mway.table.botev.sig(h,sig,lower.tail,one.side,n)
      if(half < p){
        l <- h
      }
      else{
        u <- h
      }
    }
  }
  
  for(it in 1:100){
    h <- (l+u)/2
    half <- pmway.table.null.botev(h,table,tests,lower.tail,one.side,n)
    if(half < p){
      l <- h
    }
    else{
      u <- h
    }
  }
  return (l+u)/2
}

power.botev <- function(table, tests, alpha, n){
  n.tests <- nrow(tests)
  stat <- qmway.table.botev.lite(1-alpha,table,tests)
  #alpha.est <- 1 - pmway.table.null.botev(stat, table,tests)
  pval <- pmway.table.botev(stat, table, tests,n=n)
  pval
}
