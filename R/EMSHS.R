#' EM Estimator for Bayesian Shrinkage approach with Structural Information incorporated
#'
#' \code{EMSHS} implements the EM algorithm for Bayesian shrinkage approach that incorporates
#' structural information. Users are referred to Chang et al (2018).
#'
#' @param y       An n by 1 response vector
#' @param X       An n by p design matrix
#' @param mus     A vector of shrinkage parameters
#' @param nu      The adaptivity parameter
#' @param E       An e by 2 matrix with edges. Edges must be sorted by the first column
#'                and then the second column.
#'                A single edge (j,k) must be duplicated with (k,j) in E.
#'                NULL if no edge.
#' @param a_sigma The shape parameter of the prior for residual variance.
#' @param b_sigma The rate parameter of the prior for residual variance.
#' @param a_omega The shape parameter of the prior for nonzero omega values.
#' @param b_omega The rate parameter of the prior for nonzero omega values.
#' @param w       A weight vector for samples.
#' @param eps     The algorithm stops if relative improvement goes below eps.
#'
#' @return A list that contains the number of EM iterations (niter),
#'         the estimated coefficients (beta), the estimated residual variance (sigma),
#'         the estimated shrinkage parameter (lambda),
#'         and the imputed correlations for the shrinkage parameter (omega) is returned.
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats diffinv
#' @references
#' \insertRef{chang2018}{EMSHS}
#'
#' @examples
#'
#'  # Example with no edges for a high-dimensional data
#'  # with n = 25 observations and p = 50 predictors
#'
#'  set.seed(100)
#'
#'  X <- matrix(rnorm(25*50), ncol = 50)
#'  B <- matrix(c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0), ncol = 1)
#'  e <- matrix(rnorm(25*1), ncol = 1)
#'  y <- matrix(X %*% B + e, ncol = 1)
#'  mus <- 2.3
#'  nu <- 0.3
#'
#'
#'  em_no_edge <- EMSHS(y, X, mus, nu, E = NULL,
#'                      a_sigma = 1, b_sigma = 1, a_omega = 2, b_omega = 1,
#'                      w = 1, eps = 1e-5)
#'
#'
#'  # Example with user-defined set of sorted, undirected
#'  # edges (E) for a high-dimensional data with n = 25
#'  # and p = 50
#'
#'  X <- matrix(rnorm(25*50), ncol = 50)
#'  B <- matrix(c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'                0,0,0,0,0,0,0,0,0,0), ncol = 1)
#'  e <- matrix(rnorm(25*1), ncol = 1)
#'  y <- matrix(X %*% B + e, ncol = 1)
#'  mus <- 2.3
#'  nu <- 0.3
#' EE <- matrix(c(1,4,
#'               4,1,
#'               1,2,
#'               2,1,
#'               1,5,
#'               5,1,
#'               2,3,
#'               3,2,
#'               3,5,
#'               5,3,
#'               10,11,
#'               11,10,
#'               19,11,
#'               11,19,
#'               36,35,
#'               35,36,
#'               31,35,
#'               35,31,
#'               31,22,
#'               22,31,
#'               22,45,
#'               45,22,
#'               45,32,
#'               32,45,
#'               22,21,
#'               21,22,
#'               31,21,
#'               21,31,
#'               21,25,
#'               25,21,
#'               21,18,
#'               18,21,
#'               18,49,
#'               49,18,
#'               49,47,
#'               47,49,
#'               47,37,
#'               37,47,
#'               37,21,
#'               21,37,
#'               18,25,
#'               25,18), nrow = 42, ncol = 2, byrow = TRUE)
#'
#'  # Sort edges by first column then second column
#'
#'  E <- EE[do.call(order, lapply(1:ncol(EE), function(i) EE[,i])),]
#'
#'  em_edge <- EMSHS(y, X, mus, nu, E,
#'                   a_sigma = 1, b_sigma = 1, a_omega = 2, b_omega = 1,
#'                   w = 1, eps = 1e-5)

#'@export
#'

EMSHS <- function(y,X,mus,nu,E=NULL,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1,w=1,eps=1e-5)
{
  n = nrow(X)
  p = ncol(X)

  # Initialize graph information
  # nadj[i] has the number of neighboring variables of variable i
  # idx[i] keeps the beginning and ending index in E[,2] for neighboring variables of variable i
  nadj = rep(0,p)
  if ( is.null(E) )
  {
    e = 0
  }
  else
  {
    e = nrow(E)
    for ( i in 1:e )
      nadj[E[i,1]] = nadj[E[i,1]] + 1
  }
  idx = diffinv(nadj)


  # initialize for DWL
  DWLENV$X <- as.matrix(X)*sqrt(w)
  DWLENV$y <- as.vector(y)*sqrt(w)

  DWL.Init()



  # Initialize beta
  beta = rep(0,p)



  # Initialize sigma
  c1 = DWLENV$yy/2 + b_sigma
  c2 = 0
  c3 = n+p+2*a_sigma+2
  sigma = (c2+sqrt(c2^2+8*c1*c3))/(2*c3)


  # Initialize storages
  M = length(mus)

  niters = rep(0,M)
  Qs = rep(0,M)
  betas = matrix(0,p,M)
  sigmas = rep(0,M)
  alphas = matrix(0,p,M)
  omegas = matrix(0,e,M)


  for ( mm in 1:M )
  {
    mu = mus[mm]

    # Initialize alpha
    alpha = rep(mus[mm],p)

    niter = 0

    repeat
    {
      niter = niter + 1

      ## E-Steps

      # E-Step: omega
      Eomega = 2*nu*a_omega / (2*nu*b_omega + (alpha[E[,1]]-alpha[E[,2]])^2)


      pQ = c3*log(sigma) + c1/sigma^2 + c2/sigma - sum(alpha) + sum((alpha-mu)^2)/2/nu + sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4/nu


      ## M-Step

      # M-Step: beta
      res = DWL(sigma*exp(alpha))
      beta = res$coef
      c1 = (DWLENV$yy - sum((DWLENV$Xy+DWLENV$C)*beta))/2 + b_sigma
      c2 = sum(exp(alpha)*abs(beta))


      # M-Step: sigma
      sigma = (c2+sqrt(c2^2+8*c1*c3))/(2*c3)



      # M-Step: alpha
      h = rep(0,p)
      g = rep(0,p)

      for ( i in 1:p )
      {
        if ( idx[i] < idx[i+1] )
        {
          adjidx = (idx[i]+1):idx[i+1]
          adjvar = E[adjidx,2]
          sEomega = sum(Eomega[adjidx])
          sEomegaalpha = sum(Eomega[adjidx]*alpha[adjvar])
        }
        else
        {
          sEomega = 0
          sEomegaalpha = 0
        }

        h[i] = sigma*(1+sEomega) + nu*exp(alpha[i])*abs(beta[i])
        g[i] = sigma*((1+sEomega)*alpha[i] - mu - sEomegaalpha) - nu*sigma + nu*exp(alpha[i])*abs(beta[i])
      }
      f = -sigma*(mu+nu)*sum(alpha) + nu*c2 + sigma*sum(alpha^2)/2 + sigma*sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4

      dir = g/h
      maxdir = max(abs(dir))
      m = sum(g*dir)

      ss = 1
      repeat
      {
        nalpha = alpha - ss*dir
        nc2 = sum(exp(nalpha)*abs(beta))
        nf = -sigma*(mu+nu)*sum(nalpha) + nu*nc2 + sigma*sum(nalpha^2)/2 + sigma*sum(Eomega*(nalpha[E[,1]]-nalpha[E[,2]])^2)/4

        if ( ss*maxdir < eps | f-nf > ss*m*0.05 )
        {
          alpha = nalpha
          c2 = nc2
          break
        }
        ss = ss/2
      }


      Q = c3*log(sigma) + c1/sigma^2 + c2/sigma - sum(alpha) + sum((alpha-mu)^2)/2/nu + sum(Eomega*(alpha[E[,1]]-alpha[E[,2]])^2)/4/nu

      if ( (pQ-Q)/abs(Q) < eps )
        break

    }

    niters[mm] = niter
    betas[,mm] = beta
    sigmas[mm] = sigma
    alphas[,mm] = alpha
    omegas[,mm] = Eomega
  }

  list(niter=niters,beta=betas,sigma=sigmas,lambda=exp(alphas),omega=omegas)
}



