sign.irf <- function(maxlag=length(gvarobj$F),slct_var="US.stir",G=gvarobj$G,F=gvarobj$F,sig=sigs,x=gvarobj$x,horizon=40){  
 # maxlag=length(F_1);slct_var="US.stir";G=G;F=F_1;sig=S_post;x=x;horizon=20
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  #resids <- G%*%err
  #rownames(resids) <- rownames(err) <- rownames(x)
  shock_names <- substr(slct_var,1,2)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
    if (names[[i]] %in% shock_names){ 
      k_us <- ncol(sig[[i]])
      covlist[[i]] <- t(chol(sig[[i]]))
    }else{
      covlist[[i]] <- sig[[i]]  
    }   
  }
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  P0G <- diag(nrow(x))
  for (jj in 1:length(slct_var)){
  P0 <- t(chol(sig[[which(names==substr(slct_var,1,2)[[jj]])]]))
  rownames(P0G) <- colnames(P0G) <- rownames(x)
  #slct US
  P0G[which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]]),which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]])] <- solve(P0)
  }
  K <- nrow(x) #gvarobj$x
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
  
  for (i in 1:maxlag){
    PHIx[,,i]  <-  matrix(0,K,K)
  }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  
  invG <- solve(G)
  # preallocate matrix containing GIRFs
  eslct <- matrix(0,K,1)
  
  irfa  <- array(0,c(K,K,horizon+1))
  
  colnames(irfa) <- rownames(irfa) <- rownames(x)
  sirf <- NULL
  rotstore <- matrix(0,c(nrow(x),nrow(x)))
  names <- substr(rownames(x),1,2)[!duplicated(substr(rownames(x),1,2))]
  h_sign <- 1:1
  #sirf <- array(0,c(nrow(x),length(shockc),horizon+1,ssave))
  rotS <- diag(nrow(x));colnames(rotS) <- rownames(rotS) <- substr(rownames(x),1,2)
  
  icounter <- 0
  icheck <- 0
  max.counter <- 5000
  a <- 0;b <- 0;c <- 0;d <- 0
  while(icheck<1 && icounter < max.counter){
    icounter <- icounter+1
    #print(c(a,b,c,d,icounter))
    #icounter <- icounter+1
    #Step 1: Draw a rotation matrix
    k <- length(which(substr(rownames(x),1,2)=="US")) #where to impose the sign restrictions?
    A <- matrix(rnorm(k*k,0,1),k,k)
    qA <- qr(A)
    rotA <- qr.Q(qA)
    rotA <- rotA%*%diag(((diag(rotA)>0)-(diag(rotA)<0)))
    rotS[rownames(rotS)=="US",rownames(rotS)=="US"] <- rotA
    temp <- (rotS)%*%Cm
    
    temp  <-  invG%*%Cm%*%t(rotS)#%*%Cm
    slct <- diag(nrow(x))
    impact <- as.matrix(PHI[grep("US",rownames(x)),,1]%*%(temp))[,grep("US",rownames(x))]
    colnames(impact) <- rownames(impact) #SHITTY 
    # R>0, M<0, y<0, and P<0 for the irs periods.
    #a = (imf3hat(1:irs,3,3) > 0) .* (imf3hat(1:irs,4,3) < 0) .* (imf3hat(1:irs,1,3) < 0) .* (imf3hat(1:irs,2,3) < 0);
    #Impose AD shock
    a <- (impact["US.y.qoq","US.y.qoq"]>0)*(impact["US.Dp.qoq","US.y.qoq"]>0)*(impact["US.stir","US.y.qoq"]>0)
    if (!(all(a)==1)){
      next
    }
    #impose restrictions on AS shock
    b <- (impact["US.y.qoq","US.Dp.qoq"]<0)*(impact["US.Dp.qoq","US.Dp.qoq"]>0)*(impact["US.stir","US.Dp.qoq"]>0)
    if (!(all(b)==1)){
      next
    }
    #impose restrictions on MP shock
    c <- (impact["US.y.qoq","US.stir"]<0)*(impact["US.Dp.qoq","US.stir"]<0)*(impact["US.stir","US.stir"]>0)*(impact["US.sp","US.stir"]<0)*(impact["US.eq.qoq","US.stir"]<0)
    if (!(all(c)==1)){
      next
    }

    icheck <- 1
  }
  irfa  <- array(0,c(K,K,horizon+1))
  
  temp  <- invG%*%Cm%*%t(rotS)
  for (ii in 1:(horizon+1)){
    irfa[,,ii]  <-  PHI[,,ii]%*%as.matrix(temp) #checks for the first three periods
  }

  dimnames(irfa)[[1]] <- dimnames(irfa)[[2]]  <- rownames(x)
  rot.max <- ifelse(icounter==max.counter,NA,1)
  irfa <- rot.max*irfa[,"US.stir",] #CAN SELECT AS AD AND MP SHOCK HERE; if wanted
  return(irfa)
}


irf <- function(maxlag=length(gvarobj$F),slct_var="US.stir",G=gvarobj$G,F=gvarobj$F,sig=sigs,x=gvarobj$x,horizon=40){  
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  #resids <- G%*%err
  #rownames(resids) <- rownames(err) <- rownames(x)
  shock_names <- substr(slct_var,1,2)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
    if (names[[i]] %in% shock_names){ 
      k_us <- ncol(sig[[i]])
      covlist[[i]] <- diag(k_us) 
    }else{
      covlist[[i]] <- sig[[i]]  
    }   
  }
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  P0G <- diag(nrow(x))
  for (jj in 1:length(slct_var)){
    P0 <- t(chol(sig[[which(names==substr(slct_var,1,2)[[jj]])]]))
    rownames(P0G) <- colnames(P0G) <- rownames(x)
    #slct US
    P0G[which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]]),which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]])] <- solve(P0)
  }
  K <- nrow(x) #gvarobj$x
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
  
  for (i in 1:maxlag){
    PHIx[,,i]  <-  matrix(0,K,K)
  }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  
  irfa  <- array(0,c(K,horizon+1))
  
  invGSigma_u  <-  solve(G)%*%solve(P0G)%*%Cm
  eslct <- matrix(0,nrow(x),1)
  eslct[rownames(x) %in% slct_var,1] <- 1 
  
  for (i in 1:(horizon+1)){
    irfa[,i]  <-  (PHI[,,i]%*%(invGSigma_u)%*%eslct)*as.numeric(1/sqrt(t(eslct)%*%Cm%*%eslct))
  }
  dimnames(irfa)[[1]]  <- rownames(x)
  
  
  return(irfa)
}

irf_all <- function(maxlag=length(gvarobj$F),slct_var="US.stir",G=gvarobj$G,F=gvarobj$F,sig=sigs,x=gvarobj$x,horizon=40,group_select=Grps){  
 #maxlag=1;G=G;F=F_1;sig=S_post;x=x;horizon=20;group_select=Grps
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  #resids <- G%*%err
  #rownames(resids) <- rownames(err) <- rownames(x)
  #shock_names <- substr(slct_var,1,2)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
        k_us <- ncol(sig[[i]])
        covlist[[i]] <- diag(k_us) 
 
  }
  slct_var <- rownames(x)
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  P0G <- diag(nrow(x))
  for (jj in 1:length(slct_var)){
    P0 <- t(chol(sig[[which(names==substr(slct_var,1,2)[[jj]])]]))
    rownames(P0G) <- colnames(P0G) <- rownames(x)
    #slct US
    P0G[which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]]),which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]])] <- solve(P0)
  }
  K <- nrow(x) #gvarobj$x
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
#   
#   for (i in 1:maxlag){
#     PHIx[,,i]  <-  matrix(0,K,K)
#   }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  
  irfa  <- array(0,c(K,horizon+1))
  stir_irf  <- array(0,c(horizon+1,length(group_select)))

  invGSigma_u  <-  solve(G)%*%solve(P0G)%*%Cm
  var_slct <- c("stir")
  impacts <- array(NA,c(nrow(x),length(group_select)),dimnames=list(rownames(x),names(group_select)))
  for (ik in 1:length(group_select)){
  for (jj in 1:length(var_slct)){
  sl_var <- paste(group_select[[ik]],var_slct[[jj]],sep=".")
  sl_var <- intersect(sl_var,rownames(x))
  eslct <- matrix(0,nrow(x),1,dimnames=list(rownames(x),NULL))
  eslct[sl_var,1] <- 1 
  
  for (i in 1:(horizon+1)){
    temp <- (PHI[,,i]%*%(invGSigma_u)%*%eslct)*as.numeric(1/sqrt(t(eslct)%*%Cm%*%eslct))
    irfa[,i]  <-  temp
    if (i==1){
      impacts[sl_var,ik] <- temp[sl_var,]
    }
  }
  dimnames(irfa)[[1]]  <- rownames(x)
  }
  
  stir_resp <- irfa["US.stir",]
  stir_irf[,ik] <- stir_resp
  }
  stir_all <- irfa
  dimnames(stir_irf)[[2]] <- names(group_select)
  
  return(list(stir_irf,stir_all,impacts))
}

irf_global <- function(maxlag=length(gvarobj$F),slct_var="US.stir",G=gvarobj$G,F=gvarobj$F,sig=sigs,x=gvarobj$x,horizon=40){  
  # maxlag=1;slct_var=paste(EA,"stir",sep=".");G=G;F=F_1;sig=S_post;x=x;horizon=20
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  #resids <- G%*%err
  #rownames(resids) <- rownames(err) <- rownames(x)
  shock_names <- substr(slct_var,1,2)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
    # if (names[[i]] %in% shock_names){ 
    k_us <- ncol(sig[[i]])
    covlist[[i]] <- diag(k_us) 
    #     }else{
    #       covlist[[i]] <- sig[[i]]  
    #     }   
  }
  slct_var <- cN
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  P0G <- diag(nrow(x))
  for (jj in 1:length(slct_var)){
    P0 <- t(chol(sig[[which(names==substr(slct_var,1,2)[[jj]])]]))
    rownames(P0G) <- colnames(P0G) <- rownames(x)
    #slct US
    P0G[which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]]),which(substr(rownames(x),1,2)==substr(slct_var,1,2)[[jj]])] <- solve(P0)
  }
  K <- nrow(x) #gvarobj$x
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
  
  for (i in 1:maxlag){
    PHIx[,,i]  <-  matrix(0,K,K)
  }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  
  irfa  <- array(0,c(K,horizon+1,4))
  
  invGSigma_u  <-  solve(G)%*%solve(P0G)%*%Cm
  var_slct <- c("stir","Dp","y","rer")
  
  for (jj in 1:4){
    sl_var <- rownames(x)[grepl(var_slct[[jj]],rownames(x))]
    sl_var <- sl_var[!grepl("US",sl_var)]   #remove US
    eslct <- matrix(0,nrow(x),1,dimnames=list(rownames(x),NULL))
    eslct[sl_var,1] <- 1 
    
    for (i in 1:(horizon+1)){
      irfa[,i,jj]  <-  (PHI[,,i]%*%(invGSigma_u)%*%eslct)*as.numeric(1/sqrt(t(eslct)%*%Cm%*%eslct))
    }
    dimnames(irfa)[[1]]  <- rownames(x)
  }
  stir_resp <- irfa["US.stir",,]
  stir_all <- irfa
  dimnames(stir_resp)[[2]] <- var_slct
  
  return(list(stir_resp,stir_all))
}

irf_girf <- function(maxlag=length(gvarobj$F),slct_var="US.stir",G=gvarobj$G,F=gvarobj$F,sig=sigs,x=gvarobj$x,horizon=40,group_select=Grps){  
 # maxlag=1;slct_var="US.stir";G=G;F=F_1;sig=S_post;x=x;horizon=20;group_select=Grps
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  #resids <- G%*%err
  #rownames(resids) <- rownames(err) <- rownames(x)
  shock_names <- substr(slct_var,1,2)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
    k_us <- ncol(sig[[i]])
    covlist[[i]] <- sig[[i]]
  }
  slct_var <- cN
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  P0G <- diag(nrow(x))
  K <- nrow(x) #gvarobj$x
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
  
  #   
  #   for (i in 1:maxlag){
  #     PHIx[,,i]  <-  matrix(0,K,K)
  #   }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  
  irfa  <- array(0,c(K,horizon+1,4))
  stir_irf  <- array(0,c(horizon+1,4,length(group_select)))
  
  invGSigma_u  <-  solve(G)%*%Cm
  var_slct <- c("stir","Dp","y","rer")
  impacts <- array(NA,c(nrow(x),4,length(group_select)),dimnames=list(rownames(x),NULL,names(group_select)))
  for (ik in 1:length(group_select)){
    for (jj in 1:4){
      sl_var <- paste(group_select[[ik]],var_slct[[jj]],sep=".")
      sl_var <- intersect(sl_var,rownames(x))
      eslct <- matrix(0,nrow(x),1,dimnames=list(rownames(x),NULL))
      eslct[sl_var,1] <- 1 
      
      for (i in 1:(horizon+1)){
        temp <- (PHI[,,i]%*%(invGSigma_u)%*%eslct)*as.numeric(1/sqrt(t(eslct)%*%Cm%*%eslct))
        irfa[,i,jj]  <-  temp
        if (i==1){
          impacts[sl_var,jj,ik] <- temp[sl_var,]
        }
      }
      dimnames(irfa)[[1]]  <- rownames(x)
    }
    
    stir_resp <- irfa["US.stir",,]
    stir_irf[,,ik] <- stir_resp
  }
  stir_all <- irfa
  dimnames(stir_irf)[[2]] <- var_slct
  dimnames(stir_irf)[[3]] <- names(group_select)
  
  return(list(stir_irf,stir_all,impacts))
}



