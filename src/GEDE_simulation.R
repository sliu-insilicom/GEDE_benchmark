# store all helper functions in the simulation folder: https://www.dropbox.com/home/new_GEDE/simulations
print("GEDE_simulation sourced!")



## apply observation-level Hampel filter to subgroups of data. When
## Grp is not given, it takes the entire X as one group. When
## winsorize=TRUE, it replaces the outliers by upper/lower bounds;
## otherwise, with NA. If with.out.ids=TRUE, it returns a list of two
## objects, X is the vector of results after outlier
## removal/winsorization, out.ids are the indices of outliers.
GrpHampel <- function(X, Grp=NULL, nMAD=3, type=c("two.sided", "less", "greater"), winsorize=FALSE, with.out.ids=FALSE) {
  type <- match.arg(type)
  if (is.null(Grp)) Grp <- rep(0, length(X))
  out.ids <- c()
  for (gp in unique(Grp)){
    ids <- which(Grp==gp); X.gp <- X[ids]
    Xmed <- median(X.gp, na.rm=TRUE); Xmad <- mad(X.gp, na.rm=TRUE)
    U <- Xmed+nMAD*Xmad; L <- Xmed-nMAD*Xmad
    out.lower.ids <- which(X.gp<L); out.upper.ids <- which(X.gp>U)
    if (type=="two.sided"){
      if (winsorize) {
        X.gp[out.lower.ids] <- L; X.gp[out.upper.ids] <- U
      } else {
        X.gp[out.lower.ids] <- NA; X.gp[out.upper.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.lower.ids], ids[out.upper.ids])
    } else if (type=="less") {
      if (winsorize) {
        X.gp[out.lower.ids] <- L
      } else {
        X.gp[out.lower.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.lower.ids])
    } else if (type=="greater") {
      if (winsorize) {
        X.gp[out.upper.ids] <- U
      } else {
        X.gp[out.upper.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.upper.ids])
    } else {
      stop("Valid choice of type are: two.sided, less, and greater.")
    }
    X[ids] <- X.gp
  }
  if (with.out.ids) {
    return(list(X=X, out.ids=sort(out.ids)))
  } else {
    return(X)
  }
}



## a function to generate simulated data (sim1 and sim2)
# the "Est" is previously estimated eigen and cov matrix stored in object "oracle"
PPCA.datagen <- function(n, Est, X=NULL, miss.prop=0, out.prop=0, out.L=5, out.U=10) {
  K <- Est$K; Lk <- Est$Lk; sigma2 <- Est$sigma2; Tk <- Est$Tk;
  m <- nrow(Est$Tk); N <- m*n
  ## generate the mean values
  betahat <- Est$betahat
  if (is.null(X)) { #only use the intercept
    mumat <- rep(1,n)%*%t(betahat[1,])
  } else {
    mumat <- cbind(1, X)%*%betahat
  }
  ## noise-free data
  T1 <- sweep(Tk, 2, sqrt(Lk), "*")
  Ytrue <- matrix(rnorm(K*n),n)%*%t(T1) +mumat
  ## add noise
  epsilon <- sqrt(sigma2)*matrix(rnorm(N),n)
  Y <- Ytrue +epsilon
  ## create artificial missingness and outliers
  n.miss <- N*miss.prop; n.out <- N*out.prop
  miss.out.idx <- sample(1:N, size=n.miss+n.out)
  miss.idx <- miss.out.idx[1:n.miss]
  out.idx <- setdiff(miss.out.idx, miss.idx); N.out <- length(out.idx)
  ## artificial missingness
  Ym <- Y; Ym[miss.idx] <- NA
  ## artificial outliers
  u <- runif(N.out, min=out.L, max=out.U)
  Ymo <- Ym; Ymo[out.idx] <- Ymo[out.idx] +sign(epsilon[out.idx])*u
  ##
  colnames(Ytrue) <- colnames(Y) <- colnames(Ym) <- colnames(Ymo) <- paste0("x", 1:m)
  return(list(Ytrue=Ytrue, Y=Y, Ym=Ym, Ymo=Ymo, miss.idx=miss.idx, out.idx=out.idx))
}
