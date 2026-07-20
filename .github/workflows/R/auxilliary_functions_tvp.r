bvartvpm <- function(nr,Y,X,nburn,nsave,B0inv,b0,Qprmean=0.0001,Q_prvar=100){
  require(bayesm)
  require(stochvol)
  require(snowfall)
  require(rcpp01)
  require(GIGrvg)
  
  #Get data/model dimensions (full TVP-VAR)
  T <- nrow(Y)
  K <- ncol(X)
  M <- ncol(Y)
  #Prior stuff for SV
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,T))
  hv <- svdraw$latent
  para <- list(mu=-10,phi=.9,sigma=.2)
  H <- matrix(-10,T,M)
  #Create full-data matrices
  if (nr==1) slct <- NULL else slct <- 1:(nr-1)
  Y__ <- Y[,nr]
  X__ <- cbind(Y[,slct],X)
  K_ <- ncol(X__)
  M_ <- M-length(slct)
  b00 <- matrix(0,K_,1)#b0[(M_):(K+M-1),nr,drop=FALSE]
  
  
  B_gamma <- .1  #standard value 0.1
  gamma_alpha <- 1/2
  gamma_beta <- 1/(2*B_gamma)
  omega_tau <- 10
  
  a_tau <- .1 #standard value 0.1
  c_tau <- 0.01
  d_tau <- 0.01
  V_prior <- diag(K_)*0.01
  tau2 <- matrix(0.01,K_,1)
  Qdraw <- diag(K_)*0.1#diag(as.numeric(tau2))
  xistore <- matrix(0.1,K_,1)
  #storage matrices
  H_store <- matrix(0,T,nsave)
  ALPHA_store <- array(0,c(T,K_,nsave))
  svparms_store <- matrix(0,3,nsave)
  Q_store <- matrix(NA,nsave,K_)#array(0,c(nsave,K_,K_))
  u_ <- matrix(0,T,1)
  ntot <- nburn+nsave
  for (irep in 1:ntot){
    #------------------ normalize data-------------------#
    #rescaled data, only needed in the constant VAR case
    X_ <- X__*as.numeric(exp(-hv/2))
    Y_ <- Y__*as.numeric(exp(-hv/2))

     ALPHA0 <- try(KF(t(as.matrix(Y__)),X__,as.matrix(exp(hv)),Qdraw,K_,1,T,b00,V_prior),silent=TRUE)
     if (is(ALPHA0,"try-error")){ 
       ALPHA <- try(KF_alternative(t(as.matrix(Y__)),X__,as.matrix(exp(hv)),Qdraw,K_,1,T,b00,V_prior),silent=TRUE)
       }else{
         ALPHA <- ALPHA0
       } 
    #sample variances
    omega_tau <-  rgamma(1,c_tau+a_tau*K_,d_tau+a_tau/2*sum(xistore))
     
    Em <- diff(t(ALPHA))
    for (jj in 1:K_){
      xi2 <- rgig(1,gamma_alpha-1/2,Qdraw[jj,jj],gamma_alpha*omega_tau)
      sig_q <- try(rgig(1,gamma_alpha-T/2,sum(Em[,jj]^2),1/(2*xi2)),silent=TRUE)
      if (!is(sig_q,"try-error")){
        if (sig_q<1e-07) sig_q <- 1e-07
        if (sig_q>1e+5) sig_q <- 1e+5
        Qdraw[jj,jj] <- sig_q
        xistore[jj,1] <- xi2
      } 
    }
    
    #Sample hyperparameter lambda from G(a,b)
    lambda2_tau <- rgamma(1,c_tau+a_tau*K_,d_tau+a_tau/2*sum(tau2))

    #Sample variances of the time invariant part first
    for (nn in 1:K_){
       scale1 <- ALPHA[nn,1]^2
       tautemp <- try(rgig(n=1,lambda=a_tau-0.5,ALPHA[nn,1]^2,a_tau*lambda2_tau),silent=TRUE)
       if (!is(tautemp,"try-error")){
         if (tautemp<1e-7) tautemp <- 1e-7
         if (tautemp>1e+5) tautemp <- 1e+5

         tau2[nn] <-tautemp
       }
    }
    V_prior <- diag(as.numeric(tau2))

    for (i in 1:T){
      u_[i,1] <- Y__[i]-X__[i,]%*%ALPHA[,i]
    }
    u_[abs(u_)<1e-08] <- 1e-08
    
    #-----------------sample log volas-------------------#
    svdraw <- svsample2(u_,startpara=para(svdraw),startlatent=latent(svdraw),priorphi=c(25,5))
    hv <- latent(svdraw)
    if (irep>nburn){
      H_store[,irep-nburn] <- hv
      ALPHA_store[,,irep-nburn] <- t(ALPHA)
      svparms_store[,irep-nburn] <- para(svdraw)
      Q_store[irep-nburn,] <- ifelse(runif(1,0,1)>0.5,-1,1)*sqrt(diag(Qdraw))
    }
    print(irep)
    sfCat(paste("Iteration ", irep), sep="\n")
    
  }
  return(list(A=ALPHA_store,H=H_store,svparms=svparms_store,Q=Q_store,slct=slct))
}
mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag) 
}

dmean <- function(x){
  
  xnew <- (x-mean(x))/sd(x)
  return(xnew)
}

get_companion <- function(Beta_,varndxv){
 # Beta_ <- Atilda
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
  
  return(list(MM=MM,Jm=Jm))
}




datahandling <- function(Data,nr){
  #---------------------------------------------------------------------------------------------------------------#  
  #Input:Data and the number of the country in the dataset which is the current home country
  #Output: Data_Country , a list, where Data_Country[[2]] refers to the true exogenous variables, Data_Country[[3]] 
  #are the weak exogenous variables and Data_Country[[4]] are the endogenous variables and Country[[1]] are the mn
  #This function feeds in DATA and creates a 1 x K vector which indicates which variables are exogenous, own lags 
  #lags of other variables
  #indicator 0 means home country, 1 denotes the rest of the world and 2 are truly exogenous variables
  #WARNING: original data needs labels for the countries in order to create the indicator vector
  #---------------------------------------------------------------------------------------------------------------# 
  Daten <- as.data.frame(Data) #reads in  data and converts it to a data.frame
  names  <-  colnames(Daten)
  cols <- ncol(Daten)
  indicator  <-  matrix("NA",1,ncol(Daten))
  ind  <- substr(names,1,2)
  
  number  <-nr #indicates which country is selected
  
  IND <- ind[!duplicated(ind)]
  
  for (i in 1:cols){
    
    if (grepl(IND[number],names[i]))  {
      
      indicator[i]=0 
      
    }else indicator[i]=1
    if(grepl(IND[length(IND)],names[i])) indicator[i]=2
    
  }
  
  sortingM <- rbind(indicator,as.matrix(Daten))
  
  #optional: would put the choosen country's columns at the start of the matrix, using indicator would be more efficient i guess
  sort1.data <- sortingM[,order(sortingM[1,],decreasing="TRUE")]
  
  #splitting the data into (true) exogenous, weak exogenous and system variables
  #could be more elegant, but it works
  ex <- as.matrix(which(sort1.data[1,]==2, arr.ind=T))
  wex <- as.matrix(which(sort1.data[1,]==1, arr.ind=T))
  end <- as.matrix(which(sort1.data[1,]==0, arr.ind=T))
  #Exogenous variables need extra treatment to extract dummys
  Exogenous <- sort1.data[2:nrow(sort1.data),ex]
  Exnames <- colnames(Exogenous)
  inddummy <- substr(Exnames,7,10)
  dummies <- matrix("NA",1,ncol(Exogenous))
  
  for (j in 1:ncol(Exogenous)){
    
    if (grepl("Dumm",Exnames[j])) dummies[j]=1 else dummies[j]=0
  }
  class(dummies) <- "numeric"
  
  Exogenous <- rbind(dummies,Exogenous)
  sort.Ex <- Exogenous[,order(Exogenous[1,],decreasing="TRUE")]
  
  tex <- as.matrix(which(sort.Ex[1,]==0),arr.ind=T)
  dum <- as.matrix(which(sort.Ex[1,]==1),arr.ind=T)
  
  Dummies <- sort.Ex[2:nrow(sort.Ex),dum]
  Exogenous <- sort.Ex[2:nrow(sort.Ex),tex]
  
  
  WExogenous <- sort1.data[2:nrow(sort1.data),wex]
  Endogenous <- sort1.data[2:nrow(sort1.data),end]
  
  assign(paste("WEx_",IND[number],sep=""),WExogenous) #creates a matrix WEx_Countryname
  assign(paste("END_",IND[number],sep=""),Endogenous) # same, we dont need this for Exogenous because it stays the same for all countries
  assign(paste("DataFULL_",IND[number],sep=""),rbind(indicator,as.matrix(Daten))) #returns the full matrix with a indicator row
  
  assign(paste("",IND[number],sep=""),list(IND,Dummies,Exogenous,assign(paste("WEx_",IND[number],sep=""),WExogenous),assign(paste("END_",IND[number],sep=""),Endogenous)))
  
  return(get(IND[nr]))
}
