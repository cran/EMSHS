
#' Interal Use Only
#' @keywords internal
DWLENV <- new.env(parent=emptyenv())

#' Interal Use Only
#' @keywords internal
DWL.Init <- function()
{
  DWLENV$n <- nrow(DWLENV$X)
  DWLENV$p <- ncol(DWLENV$X)
  DWLENV$m <- min(DWLENV$n,DWLENV$p)

  DWLENV$XxMax <- min(3*DWLENV$n,DWLENV$p)
  DWLENV$Xx <- matrix(0,DWLENV$p,DWLENV$XxMax)
  DWLENV$XxIdx <- rep(0,DWLENV$p)
  DWLENV$XxCnt <- 0

  DWLENV$Xy <- as.vector(t(DWLENV$X) %*% DWLENV$y)
  DWLENV$yy <- sum(DWLENV$y^2)

  DWLENV$lam <- rep(max(abs(DWLENV$Xy))*1.2,DWLENV$p)
  DWLENV$A <- c()
  DWLENV$nA <- 0
  DWLENV$B <- c()
  DWLENV$S <- c()

  DWLENV$C <- DWLENV$Xy
  DWLENV$iXXa <- matrix(0,DWLENV$m,DWLENV$m)
  DWLENV$Idx <- rep(0,DWLENV$p)
}



#' Interal Use Only
#' @keywords internal
DWL <- function(lam)
{
  # Eat free lunch
  for ( i in 1:DWLENV$p )
  {
    if ( DWLENV$Idx[i] == 0 & DWLENV$lam[i] < lam[i] )
      DWLENV$lam[i] <- lam[i]
  }

  niter = 0
  repeat
  {
    niter = niter + 1
    dlam = lam - DWLENV$lam

    # calculate dB/dalpha and dC/dalpha
    if ( DWLENV$nA > 0 )
    {
      dB = -DWLENV$iXXa[1:DWLENV$nA,1:DWLENV$nA] %*% (DWLENV$S * dlam[DWLENV$A])
      dC = -DWL.getXXa(DWLENV$A) %*% dB
    }
    else
    {
      dB = numeric()
      dC = rep(0,DWLENV$p)
    }

    # find breakpoint
    alpha = 1

    if ( DWLENV$nA > 0 )
    {
      pbp0 = -DWLENV$B/dB
      for ( l in 1:DWLENV$nA )
        if ( (DWLENV$B[l]+dB[l])*DWLENV$S[l] < 0 & pbp0[l] < alpha )
        {
          alpha = pbp0[l]
          type = 0
          idx = DWLENV$A[l]
        }
    }

    pbp1 = (DWLENV$lam-DWLENV$C)/(dC-dlam)
    pbp2 = -(DWLENV$lam+DWLENV$C)/(dC+dlam)
    for ( k in 1:DWLENV$p )
      if ( DWLENV$Idx[k] == 0 )
      {
        if ( DWLENV$C[k]+dC[k] > DWLENV$lam[k]+dlam[k] & pbp1[k] < alpha )
        {
          alpha = pbp1[k]
          type = 1
          idx = k
        }
        if ( DWLENV$C[k]+dC[k] < -DWLENV$lam[k]-dlam[k] & pbp2[k] < alpha )
        {
          alpha = pbp2[k]
          type = -1
          idx = k
        }
      }


    # Add or remove var
    if ( alpha < 1 )
      if ( type == 0 )
        DWL.remove(idx)
      else
      {
        DWLENV$B <- DWLENV$B + alpha*dB
        DWL.add(idx,type)
      }


    # compute B and C at alpha with new A and S
    DWLENV$lam <- DWLENV$lam + dlam*alpha
    if ( DWLENV$nA > 0 )
    {
      DWLENV$B <- DWLENV$iXXa[1:DWLENV$nA,1:DWLENV$nA] %*% ( DWLENV$Xy[DWLENV$A] - DWLENV$S*DWLENV$lam[DWLENV$A] )
      DWLENV$C <- DWLENV$Xy - DWL.getXXa(DWLENV$A) %*% DWLENV$B
    }
    else
    {
      DWLENV$B <- c()
      DWLENV$C <- DWLENV$Xy
    }

    if ( alpha ==  1 )
      break
  }

  coef = rep(0,DWLENV$p)
  coef[DWLENV$A] = DWLENV$B

  list(coef=coef,niter=niter)
}

#' Interal Use Only
#' @keywords internal
DWL.add <- function(k,sgn)
{
  b = DWL.getXXa(k)

  if ( DWLENV$nA > 0 )
  {
    a = DWLENV$iXXa[1:DWLENV$nA,1:DWLENV$nA]
    del = drop(a %*% b[DWLENV$A])
    d = drop(b[k] - crossprod(del,b[DWLENV$A]))

    if ( d < 1e-8 )
    {
#      message("Warning: numerical instability")

      pos = which.max(del*sgn/DWLENV$B)
      DWL.remove(DWLENV$A[pos])

      if ( DWLENV$nA > 0 )
      {
        a = DWLENV$iXXa[1:DWLENV$nA,1:DWLENV$nA]
        del = drop(a %*% b[DWLENV$A])
        d = drop(b[k] - crossprod(del,b[DWLENV$A]))
      }
    }
  }

  # Now add k
  if ( DWLENV$nA > 0 )
  {
    DWLENV$iXXa[1:DWLENV$nA,1:DWLENV$nA] <- a + del %*% t(del) / d
    DWLENV$iXXa[1:DWLENV$nA,DWLENV$nA+1] <- -del / d
    DWLENV$iXXa[DWLENV$nA+1,1:DWLENV$nA] <- -del / d
    DWLENV$iXXa[DWLENV$nA+1,DWLENV$nA+1] <- 1/d
  }
  else
  {
    DWLENV$iXXa[1] <- 1/b[k]
  }
  DWLENV$Idx[k] <- DWLENV$nA+1
  DWLENV$nA <- DWLENV$nA+1
  DWLENV$A <- c(DWLENV$A,k)
  DWLENV$S <- c(DWLENV$S,sgn)
}


#' Interal Use Only
#' @keywords internal
DWL.remove <- function(k)
{
  l = DWLENV$Idx[k]
  m = DWLENV$nA
  DWLENV$Idx[k] <- 0
  if ( l<m )
    DWLENV$Idx[DWLENV$A[(l+1):m]] <- DWLENV$Idx[DWLENV$A[(l+1):m]] - 1
  DWLENV$nA <- m-1
  DWLENV$A <- DWLENV$A[-l]
  DWLENV$S <- DWLENV$S[-l]

  if ( m>1 )
  {
    a = DWLENV$iXXa[1:m,1:m]
    b = a[,l]
    DWLENV$iXXa[1:(m-1),1:(m-1)] <- a[-l,-l] - b[-l] %*% t(b[-l]) / b[l]
  }

  DWLENV$iXXa[,m] <- 0
  DWLENV$iXXa[m,] <- 0
}


#' Interal Use Only
#' @keywords internal
DWL.getXXa <- function(A)
{
  for ( k in A )
    if ( DWLENV$XxIdx[k] == 0 )
    {
      DWLENV$XxCnt <- DWLENV$XxCnt + 1
      if ( DWLENV$XxCnt > DWLENV$XxMax )
      {
        oldmax = DWLENV$XxMax
        oldXx = DWLENV$Xx
        DWLENV$XxMax <- min(oldmax*2,DWLENV$p)
        DWLENV$Xx <- matrix(0,DWLENV$p,DWLENV$XxMax)
        DWLENV$Xx[,1:oldmax] <- oldXx
      }
      DWLENV$XxIdx[k] <- DWLENV$XxCnt
      DWLENV$Xx[,DWLENV$XxCnt] <- t(DWLENV$X) %*% DWLENV$X[,k]
    }

  DWLENV$Xx[,DWLENV$XxIdx[A]]
}


#' Interal Use Only
#' @keywords internal
DWL.version <- function()
{
  print("0.2")
}
