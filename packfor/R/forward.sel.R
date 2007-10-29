"forward.sel" <-
  function(Y,X,K=nrow(X)-1, R2thresh=.99,adjR2thresh=.99,nperm=999,R2more=0.001,alpha=0.05,Xscale=TRUE,Ycenter=TRUE,Yscale=FALSE){
    X <- as.data.frame(X)
    Y <- as.data.frame(Y)
    if (any(is.na(X))|any(is.na(X))) stop("na entries in table")
    if(nrow(X)!=nrow(Y)) stop("different number of rows")
    if (any(apply(X,2,is.factor))|any(apply(Y,2,is.factor))) stop("not yet implemented for factors")
    X <- apply(X,2,scale,scale=Xscale)
    Y <- apply(Y,2,scale,scale=Yscale,center=Ycenter)
    nbcovar <- 0
    W <- NULL
    if(!is.null(W)){
      W <- as.data.frame(W)
      W <- apply(W,2,scale,scale=Wscale)
      if(nrow(W)!=nrow(Y)) stop("different number of rows")
      Yori <- Y
      Xori <- X
      X <- as.data.frame(residuals(lm(as.matrix(X)~as.matrix(W))))
      Y <- as.data.frame(residuals(lm(as.matrix(Y)~as.matrix(W))))
      nbcovar=ncol(W)
    }
    pval <- rep(1,ncol(X))
    ordre <- rep(0,ncol(X))
    R2 <- rep(0,ncol(X))
    adjR2 <- rep(0,ncol(X))
    Fvalue <- rep(0,ncol(X))
    res <- list()
    res <- .C("forwardsel",as.double(t(X)),as.double(t(Y)),as.integer(nrow(X)),as.integer(ncol(X)),as.integer(ncol(Y)),pval=as.double(pval), ord=as.integer(ordre),Fval=as.double(Fvalue),as.integer(nperm),R2=as.double(R2),adjR2=as.double(adjR2),as.integer(K),as.double(R2thresh),as.double(adjR2thresh),as.double(R2more),as.integer(nbcovar),as.double(alpha),PACKAGE="packfor")[c("ord","Fval","pval","R2","adjR2")]
    lambdA <- c(res$R2[1],diff(res$R2))
    resmat <- data.frame(res$ord,lambdA,res$R2,res$adjR2,res$Fval,res$pval)
    if(sum(res$ord>0)==0) stop("No variables selected. Please change your parameters.")
    resmat <- resmat[res$ord>0,]
    resmat <- cbind(I(colnames(X)[resmat[,1]]),resmat)
    if(!is.null(W)){
      Yori <- as.matrix(Yori)
      Y <- as.matrix(Y)
      trori <- sum(diag(crossprod(Yori)))
      trdt <- sum(diag(crossprod(Y)))
      resmat[,3] <- resmat[,3]*trdt/trori
      resmat[,4] <- resmat[,4]*trdt/trori
    }
    names(resmat) <- c("variables","order","R2","R2Cum","AdjR2Cum","F","pval")
    return(resmat)
  }
