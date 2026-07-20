BVAR <- function(nr,gW=gW,bigx=xglobal,Daten,cN,nsave=1500,nburn=1500,thin_chain=thin,ext.inst=NULL,parms=NULL){
#nr <- 1;gW=gW;bigx=xglobal;Daten;cN;nsave=saves;nburn=burns;thin_chain=0.1; ext.inst=ext.inst
  require(coda)
  require(threshtvp)
  library(MASS) 
  require(bayesm)
  require(psych)
  require(MCMCpack)
  require(mvtnorm)
  require(mnormt)
  library(compiler)
  #----------------------user input stuff---------------- 
  xglobal <- bigx
  gW <- gW
  nr <- nr
  p <<- 1#number of lags of the dependent variable
  pwex <- 1 #number of lags of weakly exogenous variables (will be included later) 
  cons <- 1 #1 includes constant, 0 exclude
  trend <- 1 #1 includes trend component i.e. alpha*t
  #Starting values for the hyperparameters used later on
  #--------------Gibbs-related prelims---------------------
  ntot <- nsave+nburn # total draws
  #--------------------Reading in the Data-----------------
  Names <- cN
  End <- xglobal[,substr(colnames(xglobal),1,2)==Names[[nr]]]
  class(End) <- "numeric"
  PIP <- NULL
  class(xglobal) <- "numeric"
  W <- gW[[nr]]
  
  #Construct wWeakly exogenous variables
  all <- W%*%t(xglobal)
  Wex <- all[(ncol(End)+1):nrow(all),]
  Wex <- t(Wex)
  class(Wex) <- "numeric"  
  
  X <- End
  Yraw <- X
  Yraw <- as.matrix(Yraw)
  class(Yraw) <- "numeric"
  #--------------------Data Preparation-------------------------
  dimYraw <- dim(Yraw) #get dimensions of the dependent variable
  Traw <- dimYraw[1]
  M <- N <- dimYraw[2]
  #-------------------------------------------------------------
  #uses the function above
  Ylag <- mlag(Yraw,p)
  Wexlag <- mlag(Wex,pwex)
  
  M_ <- ncol(Wex)
  
  colnames(Wex) <- rep("Wex",ncol(Wex))
  
  #creates nametags for the lags, makes it easier to select the corresponding coefficient matrices later on
  nameslags <- NULL
  wexnameslags <- NULL
  for (ii in 1:p){
    nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),ncol(Yraw)))
    wexnameslags <- c(wexnameslags,rep(paste("Wexlag",ii,sep=""),ncol(Wex)))
  }
  colnames(Ylag) <- nameslags
  colnames(Wexlag) <- wexnameslags
  
  X1 <- cbind(Wex,Wexlag,Ylag)
  X1 <- X1[(p+1):nrow(X1),]
  #---------we need to include this for the trend component and constant------------#
  #  colnames(X1)[[1]] <- "trend"
  X1 <- cbind(1,X1)
  colnames(X1)[[1]] <- "constant" 
  if (length(ext.inst)>1){
    X1 <- cbind(ext.inst[(p+1):length(ext.inst)],X1)
    colnames(X1)[[1]] <- "instr" 
  }
  
  #-----------------get size of the final matrix X---------------------#
  dimX <- dim(X1)
  Traw3 <- dimX[1]
  K <- dimX[2]
  T <- Traw-p
  
  #------------------create block diagonal matrix z-------------------#
  Z1 <- kronecker(diag(N),X1)
  Y <- Yraw[(p+1):nrow(Yraw),]
  X <- X1
  Z <- kronecker(diag(N),X)
  #----------------------PRIOR FUN------------------------------
  #----------------------getting OLS estimates------------------ 
  A_OLS  <-  ginv(t(X)%*%X)%*%(t(X)%*%Y)
  a_OLS <- as.vector(A_OLS) #vectorizes A_OLS, i.e. a_OLS=vec(A_OLS)
  SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <-  SSE/(T-K+1)
  
  #-------------Initialize Bayesian Posteriors using OLS values-
  alpha <- a_OLS # single draw from the posterior of alpha
  ALPHA <- A_OLS # -,,- ALPHA
  SSE_Gibbs <- SSE # of SSE
  SIGMA <- SIGMA_OLS #of SIGMA
  
  #store the draws:
  alpha_draws <- matrix(0,nsave,K*N)
  ALPHA_draws <- array(0,dim=c(nsave,K,N))
  SIGMA_draws <- array(0,dim=c(nsave,N,N))
  lambda_draws <- matrix(0,nsave,4)
  PL <- matrix(0,nsave,1)
  
  #-------------------PRIOR Hyperparameters--------------------
  M <- N
  n <- K*N
  post_draws <-D.list<-thresh_draws<- list()
  if (is.null(parms)){
      parms <- list(B_1=3, B_2=0.03,kappa0=-1e-05)
  }
  
  
  
  for (ii in 1:M){
  if (ii==1) slct <- NULL else slct <- 1:(ii-1)
    Y__ <- Y[,ii]
    X__ <- cbind(Y[,slct],X)
    
    post.iii  <- try(estimate_tvp(Y__,X__,save=nsave,burn=nburn,p=p,sv_on = TRUE,thin = thin_chain,priorbtheta = parms,priormu=c(0,10),h0prior="stationary", grid.length = 150, thrsh.pct = 0.1,thrsh.pct.high = 1.5,TVS=TRUE,CPU=1),silent=TRUE)
    
    if (is(post.iii,"try-error")){
         print(post.iii)
            post.ii  <- estimate_tvp(Y__,X__,save=nsave,burn=nburn,p=p,sv_on = TRUE,thin = thin_chain,priorbtheta = parms,priormu=c(0,10),h0prior="stationary", grid.length = 150, thrsh.pct = 0.1,thrsh.pct.high = 1.5,TVS=TRUE,CPU=1)
    }else{
        post.ii <- post.iii
    }
    
    post.ii$posterior$slct <- slct
    post.ii$posterior$A <- aperm(post.ii$posterior$A,c(2,3,1))
    post_draws[[ii]] <- post.ii
    thresh_draws[[ii]]<-post.ii$posterior$thresholds
    D.list[[ii]]<-post.ii$posterior$D_dyn
  }
  nthin <- round(thin_chain*nsave)
  #Storage matrices for the full system
  H_store <- array(0,c(T,M,nthin))
  ALPHA_draws <- array(0,c(T,K,M,nthin))
  A0_store <- array(0,c(T,M,M,nthin))  
  S_draws <- array(0,c(T,M,M,nthin))
  res_store<-matrix(0,nrow(X),M)
  time.var.process <- array(NA,c(nthin,K,M,T))
  
  for (ii in 1:nthin){
    for (jj in 1:M){
      slct <- post_draws[[jj]]$posterior$slct
      #split and create structural matrix
      A0_store[,jj,slct,ii] <- (-1*post_draws[[jj]]$posterior$A[,slct,ii])
      if (jj==1){
        ALPHA_draws[,,jj,ii] <- (post_draws[[jj]]$posterior$A[,,ii])
        time.var.process[ii,,jj,] <- post_draws[[jj]]$posterior$Omega[ii,,]
      }else{
        ALPHA_draws[,,jj,ii] <- (post_draws[[jj]]$posterior$A[,-slct,ii])
        time.var.process[ii,,jj,] <- post_draws[[jj]]$posterior$Omega[ii,-slct,]
        
      }
      H_store[,jj,ii] <- post_draws[[jj]]$posterior$H[ii,]
    }
    res<-NULL
    for (nn in 1:T){
      A0 <- A0_store[nn,,,ii]
      diag(A0) <- 1
      A0inv <- try(solve(A0),silent=TRUE)
      if (is(A0inv,"try-error")) A0inv <- ginv(A0inv)
      ALPHA_draws[nn,,,ii] <- t(A0inv%*%t(ALPHA_draws[nn,,,ii]))
      S_draws[nn,,,ii] <- A0inv%*%diag(exp(H_store[nn,,ii]))%*%t(A0inv)
      res<-rbind(res,Y[nn,]-X[nn,,drop=FALSE]%*%ALPHA_draws[nn,,,ii])
    }
    res_store<-res_store+res
  }
 
  
  ALPHA <- apply(ALPHA_draws,c(1,2,3),mean)
  D.list<-lapply(D.list,function(x) apply(x,c(2,3),mean,na.rm=TRUE))
  fit <- NULL
  for (nn in 1:T){
    fit <- rbind(fit,Y[nn,]-X[nn,,drop=FALSE]%*%ALPHA[nn,,])
  }
  
  dimnames(ALPHA_draws)=list(NULL,colnames(X),colnames(A_OLS))
  
  A_post <- apply(ALPHA_draws, c(1,2,3), mean)
  S_post <- apply(S_draws,c(1,2,3),mean)
  H_store<-exp(0.5*H_store)
  H_post<-apply(H_store,c(1,2),median)
  Res <- fit
  
  slct <- round(seq(1,nsave,length.out = thin_chain))#sample(1:nsave,0.3*nsave)
  dims <- dimnames(ALPHA_draws)[[2]] 
  #splitting the output of the gibbs sampler for predictive density later on
  a0post <- ALPHA_draws[,which(dims=="constant"),,]
  a1post <- ALPHA_draws[,which(dims=="trend"),,] #coefficients on the trend
  a2post <- NULL
  alpha_var <- NULL
  dummiespost <- ALPHA_draws[,which(dims=="Dummies"),,] #coefficients for the dummies
  postExpost <- ALPHA_draws[,which(dims=="Exogenous"),,] #coefficients for the cont.exogenous
  postExlpost <- ALPHA_draws[,which(dims=="ExogenousLag"),,] #coefficients for the lagged. ex
  Lambda0post <- ALPHA_draws[,which(dims=="Wex"),,]#coefficients on WEX
  SIGMApost <- S_draws
  INST_post <- ALPHA_draws[,which(dims=="instr"),,]
  cc.res<-res_store/nthin
  
  Lambdapost <- NULL
  Thetapost <- NULL
  
  for (jj in 1:p) {
    Lambdapost[[jj]] <- ALPHA_draws[,which(dims==paste("Wexlag",jj,sep="")),,]
    Thetapost[[jj]]<- ALPHA_draws[,which(dims==paste("Ylag",jj,sep="")),,]
  }
  #--------------Splits the posterior matrix & creates residuals---------------------
  dims <- dimnames(A_post)[[2]] 
  ML <- 0#-0.5*(-2*Lik+ncol(ALPHA)*nrow(ALPHA)*log(T))
  
  
  a0 <- A_post[,which(dims=="constant"),] #coefficients on the constant
  a1 <- A_post[,which(dims=="trend"),] #coefficients on the trend
  a2 <- NULL
  dummies <- A_post[,which(dims=="Dummies"),] #coefficients for the dummies
  postEx <- A_post[,which(dims=="Exogenous"),] #coefficients for the cont.exogenous
  postExl <- A_post[,which(dims=="ExogenousLag"),] #coefficients for the lagged. ex
  Lambda0 <- A_post[,which(dims=="Wex"),]#coefficients on WEX
  Lambda <- A_post[,which(dims=="Wexlag1"),] #coefficients on Wexlag
  Theta <- A_post[,which(dims=="Ylag1"),]#coefficients on endogenous lagged
  
  Lambda <- NULL
  Theta <- NULL
  
  for (jj in 1:p) {
    Lambda[[jj]] <- A_post[,which(dims==paste("Wexlag",jj,sep="")),]
    Theta[[jj]]<- A_post[,which(dims==paste("Ylag",jj,sep="")),]
  }
  sfCat(paste("Iteration ", nr), sep="\n")

  
  #Data Splitting cbind(Dummies,Exogenous,Exlag,Wex,Wexlag,Ylag)
  DDummies <- X[,which(dims=="Dummies")]
  DExogenous <- X[,which(dims=="Exogenous")]
  DExogenousL <- X[,which(dims=="ExogenousLag")]
  DWex <- X[,which(dims=="Wex")]
  DWexlag <- X[,which(dims=="Wexlag")]
  DEndlag <- Ylag[(p+1):nrow(Ylag),]
  # return(list(Dummies=DDummies,End=Y,Exogenous=DExogenous,Exlag=DExogenousL,Wex=DWex,Wexlag=DWexlag,Ylag=DEndlag,A_post=A_post,
  #             Res=Res,a0=a0,a1=a1,a2=a2,dummies=dummies,postEx=postEx,postExl=postExl,Lambda0=Lambda0,Lambda=Lambda,Theta=Theta,p=p,bigX=X,X=X,W=W,a0post=a0post,a1post=a1post,a2post=a2post,
  #             dummiespost=dummiespost,postExpost=postExpost,postExlpost=postExlpost,Lambda0post=Lambda0post,Lambdapost=Lambdapost,Thetapost=Thetapost,SIGMApost=SIGMApost,alpha_Var=alpha_var,PL=PL,ALPHA_draws=NULL,
  #             ML=ML,SIGMA=S_post,PIP=PIP,lambda=lambda_draws,inst=INST_post,time.var.process=time.var.process))
   return(list(ALPHA=ALPHA_draws,SIGMApost=SIGMApost,SIGMA=S_post,PIP=PIP,lambda=lambda_draws,inst=INST_post,time.var.process=time.var.process,W=W,thresholds=thresh_draws,D.list=D.list,H_post=H_post,cc.res=cc.res))

}
