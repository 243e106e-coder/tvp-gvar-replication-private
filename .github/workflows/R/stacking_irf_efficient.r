
irf.mcmc <- function(draw_i, Sig_draw_i, x, globalG){
  rm(A.list)
  rm(Sigma.posterior)
  
  bigT<- ncol(x)  
  pts <- seq(1,bigT)
  IRF_post <- array(NA,dim=c(length(pts),ncol(xglobal),nhor+1))
  F.eigen<-rep(NA,bigT)
  S_post <- list()
  
  
  #start the loop over time
  tt <- 0
  for (t in 1:bigT){
    if (t %in% pts){
      tt <- tt+1
      
      VAR0 <- split.function(draw_i[[1]])
      #stacking/solving like in pesaran
      #parameter stuff
      k_i <- dim(VAR0$Lambda0post)[[3]]
      A0 <- cbind(diag(k_i),-t(VAR0$Lambda0post[t,,]))
      W0 <- globalG[[1]]
      
      #creates matrices according to PSW 2004 for every lag
      for (kk in 1:length(VAR0$Theta)){
        assign(paste("B0",kk,sep=""),cbind(t(VAR0$Thetapost[[kk]][t,,]),t(VAR0$Lambdapost[[kk]][t,,])))
        assign(paste("H0",kk,sep=""),get(paste("B0",kk,sep=""))%*%W0)
      }
      
      if (ext.inst) kappa0 <- matrix(VAR0$inst[t,])
      G <- A0%*%W0
      S_post[[1]] <- Sig_draw_i[[1]][t,,]
      
      for (i in 2:length(Sig_draw_i)){
        S_post[[i]] <- Sig_draw_i[[i]][t,,]
        #N might create a list containing all country specific W matrices and use that in here (should be much faster)
        VAR1 <- split.function(draw_i[[i]])
        W1 <- globalG[[i]]
        k_i <- dim(VAR1$Lambda0post)[[3]]
        
        A1 <- cbind(diag(k_i),-1*t(VAR1$Lambda0post[t,,]))
        
        for (kk in 1:length(VAR0$Theta)){
          assign(paste("B1",kk,sep=""),cbind(t(VAR1$Thetapost[[kk]][t,,]),t(VAR1$Lambdapost[[kk]][t,,])))
          assign(paste("H0",kk,sep=""),rbind(get(paste("H0",kk,sep="")),get(paste("B1",kk,sep=""))%*%W1))
        }
        G <- rbind(G,A1%*%W1)
        if (ext.inst)  kappa0 <- rbind(kappa0,matrix(VAR1$inst[t,,irep]))
      }
      
      G.inv<-solve(G)

      if (ext.inst)  zeta <- G.inv%*%kappa0
      for (kk in 1:length(VAR0$Theta)){
        assign(paste("F",kk,sep=""),G.inv%*%get(paste("H0",kk,sep="")))
      }
      
      #err <- solve(G)%*%resids
      F_1 <- list()
      for (kk in 1:length(VAR0$Theta)){
        F_1[[kk]] <- get(paste("F",kk,sep=""))
      }
      
      #This code is used to construct the responses to a shock to the external instrument
      compMat <- F_1[[1]] #NEED COMPANION FORM HERE IF WE INCLUDE MORE LAGS
      # aux<-suppressWarnings(eigs(compMat,k=1,silent=TRUE))
      # if(length(aux$values)>0){
      #   F.eigen[t]<-Mod(aux$values)
      # }
      # 
      # if (Re(aux$values)>1.05) next
      
      plag <- length(F_1)
      M <- nrow(x)
      if (ext.inst){
        irf.1 <- array(0,c(M*plag,nhor+1))
        for (jj in 1:(nhor+1)){
          irf.1[,jj] <- c(zeta,matrix(0,M*((plag-1)),1))
          if (jj>1){
            irf.1[,jj] <- compMat%*%irf.1[,jj-1]
          }
        }
        
        rownames(irf.1) <- rownames(x)
        scal<-irf.1["US.stir",1];irf.1<-(irf.1/scal)*0.01
        IRF_post[tt,,] <- irf.1
      }else{
        irf_draw <- sign.irf(maxlag=length(F_1),slct_var="US.stir",G=G,F=F_1,sig=S_post,x=x,horizon=20)
        
        rownames(irf_draw) <- rownames(x)
        IRF_post[tt,,] <- irf_draw#/irf_draw["US.stir",1]*0.0025
        print(tt)
      }
    }
  }
   return(list(IRF_post=IRF_post,F.eigen=F.eigen))
 }  


split.function <- function(ALPHA.draw){
  dims <- dimnames(ALPHA.draw)[[2]]
  
  a0post <- ALPHA.draw[,which(dims=="constant"),]
  a1post <- ALPHA.draw[,which(dims=="trend"),] #coefficients on the trend
  a2post <- NULL
  alpha_var <- NULL
  dummiespost <- ALPHA.draw[,which(dims=="Dummies"),] #coefficients for the dummies
  postExpost <- ALPHA.draw[,which(dims=="Exogenous"),] #coefficients for the cont.exogenous
  postExlpost <- ALPHA.draw[,which(dims=="ExogenousLag"),] #coefficients for the lagged. ex
  Lambda0post <- ALPHA.draw[,which(dims=="Wex"),]#coefficients on WEX
  #  SIGMApost <- S_draws
  INST_post <- ALPHA.draw[,which(dims=="instr"),]
  
  Lambdapost <- NULL
  Thetapost <- NULL
  
  for (jj in 1:1) {
    Lambdapost[[jj]] <- ALPHA.draw[,which(dims==paste("Wexlag",jj,sep="")),]
    Thetapost[[jj]]<- ALPHA.draw[,which(dims==paste("Ylag",jj,sep="")),]
  }
  return(list(a0post=a0post,a1post=a1post,a2post=a2post,dummiespost=dummiespost,Lambda0post=Lambda0post,Lambdapost=Lambdapost,Thetapost=Thetapost))
}

get.irfa.t <-   function(tt,draw_i, Sig_draw_i, x, globalG,horz=nhor){
      t <- tt
      bigT<- ncol(x)  
      IRF_post <- array(NA,dim=c(ncol(xglobal),horz+1))
      S_post <- list()
      #start the loop over time
      VAR0 <- split.function(draw_i[[1]])
      #stacking/solving like in pesaran
      #parameter stuff
      k_i <- dim(VAR0$Lambda0post)[[3]]
      A0 <- cbind(diag(k_i),-t(VAR0$Lambda0post[t,,]))
      W0 <- globalG[[1]]
      
      #creates matrices according to PSW 2004 for every lag
      for (kk in 1:length(VAR0$Theta)){
        assign(paste("B0",kk,sep=""),cbind(t(VAR0$Thetapost[[kk]][t,,]),t(VAR0$Lambdapost[[kk]][t,,])))
        assign(paste("H0",kk,sep=""),get(paste("B0",kk,sep=""))%*%W0)
      }
      
      if (ext.inst) kappa0 <- matrix(VAR0$inst[t,])
      G <- A0%*%W0
      S_post[[1]] <- Sig_draw_i[[1]][t,,]
      
      for (i in 2:length(Sig_draw_i)){
        S_post[[i]] <- Sig_draw_i[[i]][t,,]
        #N might create a list containing all country specific W matrices and use that in here (should be much faster)
        VAR1 <- split.function(draw_i[[i]])
        W1 <- globalG[[i]]
        k_i <- dim(VAR1$Lambda0post)[[3]]
        
        A1 <- cbind(diag(k_i),-1*t(VAR1$Lambda0post[t,,]))
        
        for (kk in 1:length(VAR0$Theta)){
          assign(paste("B1",kk,sep=""),cbind(t(VAR1$Thetapost[[kk]][t,,]),t(VAR1$Lambdapost[[kk]][t,,])))
          assign(paste("H0",kk,sep=""),rbind(get(paste("H0",kk,sep="")),get(paste("B1",kk,sep=""))%*%W1))
        }
        G <- rbind(G,A1%*%W1)
        if (ext.inst)  kappa0 <- rbind(kappa0,matrix(VAR1$inst[t,,irep]))
      }
      
      G.inv<-solve(G)
      
      if (ext.inst)  zeta <- G.inv%*%kappa0
      for (kk in 1:length(VAR0$Theta)){
        assign(paste("F",kk,sep=""),G.inv%*%get(paste("H0",kk,sep="")))
      }
      
      #err <- solve(G)%*%resids
      F_1 <- list()
      for (kk in 1:length(VAR0$Theta)){
        F_1[[kk]] <- get(paste("F",kk,sep=""))
      }
      
      #This code is used to construct the responses to a shock to the external instrument
      compMat <- F_1[[1]] #NEED COMPANION FORM HERE IF WE INCLUDE MORE LAGS

      plag <- length(F_1)
      M <- nrow(x)
      if (ext.inst){
        irf.1 <- array(0,c(M*plag,horz+1))
        for (jj in 1:(horz+1)){
          irf.1[,jj] <- c(zeta,matrix(0,M*((plag-1)),1))
          if (jj>1){
            irf.1[,jj] <- compMat%*%irf.1[,jj-1]
          }
        }
        
        rownames(irf.1) <- rownames(x)
        scal<-irf.1["US.stir",1];irf.1<-(irf.1/scal)*0.01
        IRF_post[tt,,] <- irf.1
      }else{
        irf_draw <- sign.irf(maxlag=length(F_1),slct_var="US.stir",G=G,F=F_1,sig=S_post,x=x,horizon=horz)
        
        rownames(irf_draw) <- rownames(x)
        IRF_post <- irf_draw#/irf_draw["US.stir",1]*0.0025 # this is the normalization
       
      }
  
  return(list(IRF_post=IRF_post))
}

get.gpr.irfa.t <- function(tt, draw_i, Sig_draw_i, x, globalG, horz=12) {
  t <- tt
  S_post <- list()
  VAR0 <- split.function(draw_i[[1]])
  k_i <- dim(VAR0$Lambda0post)[[3]]
  A0 <- cbind(diag(k_i), -t(VAR0$Lambda0post[t,,]))
  W0 <- globalG[[1]]
  H <- list()
  for (kk in seq_along(VAR0$Theta)) {
    B0 <- cbind(t(VAR0$Thetapost[[kk]][t,,]), t(VAR0$Lambdapost[[kk]][t,,]))
    H[[kk]] <- B0 %*% W0
  }
  G <- A0 %*% W0
  S_post[[1]] <- Sig_draw_i[[1]][t,,]
  if (length(Sig_draw_i) >= 2) for (i in 2:length(Sig_draw_i)) {
    VAR1 <- split.function(draw_i[[i]])
    k_i <- dim(VAR1$Lambda0post)[[3]]
    A1 <- cbind(diag(k_i), -t(VAR1$Lambda0post[t,,]))
    W1 <- globalG[[i]]
    for (kk in seq_along(VAR1$Theta)) {
      B1 <- cbind(t(VAR1$Thetapost[[kk]][t,,]), t(VAR1$Lambdapost[[kk]][t,,]))
      H[[kk]] <- rbind(H[[kk]], B1 %*% W1)
    }
    G <- rbind(G, A1 %*% W1)
    S_post[[i]] <- Sig_draw_i[[i]][t,,]
  }
  invG <- solve(G)
  F <- lapply(H, function(h) invG %*% h)
  list(IRF_post = gpr_girf(length(F), G, F, S_post, x,
                           horizon=horz, shock_var="US_gpr"))
}
