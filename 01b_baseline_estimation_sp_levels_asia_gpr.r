  rm(list=ls())
  #-------------------------------------Load necessary packages ---------------------------------------------------------#

  library(compiler)
  library(rARPACK)
  require(snowfall)
  require(zoo)
  require(Matrix)
  require(mvtnorm)
  require(threshtvp)
  
#-------------------------------------Source helper scripts ---------------------------------------------------------#
  source("BVAR_ttvp.r")
  source("Datahandling.r")
  source("stacking_irf_efficient.r")
  source("auxilliary_functions_tvp.r")
  source("irf_cholesky.r")
  
  #------------------------------------Some preliminaries------------------------------------------------------------#
  
  foldername <- paste0("Results",sep="")
  dir.create(foldername, showWarnings = FALSE)
  
  glob.shrink.list <- list(one=list(B_1=3, B_2=0.03,kappa0=-0.1/20)) #selects key hyperparameters for the TTVP specification, in the paper the off-setting constant kappa0 is defined as kappa0^2
  
  comb.grid <- expand.grid(1:length(glob.shrink.list),1)
  zz <- 1  #in case length(glob.shrink.list)>1

  CPU=1 #set equal to the total number of CPU cores available on your desktop/server/cluster ;
  saves <- 5 #number of retained draws (should be really large, >10,000)
  burns <- 5 #number of burned draws (should be even larger, >30,000)
  thin <- 1 #thinning (i.e. retain every 10th draw, saves memory)

  #-------------------------------------Loads the dataset ---------------------------------------------------------#
  long.sl <- "long"
  if (long.sl=="long"){
    load("asia_gpr_3var_dataset_for_cluster.rda")
  }else{
    load(paste0("Data_sets.rda"))
  }
  enableJIT(2)
  
  West <- c("DE","ES","FR","NO","GB")
  Rest <- c("US","CA","JP","AU","NZ")
  Asia <- c("CN","IN","ID","KR","TH")
  LA <- c("BR","CL","PE","MX","AR")
  EA<-c("AT","BE","DE","ES","FI","FR","GR", "IT","NL","PT")
  Grps <- list(West=West,Rest=Rest,Asia=Asia,LA=LA)
  
  Data.setup <- data.sets$asia_gpr_3var
  
  Daten <- NULL
  
  xglobal <- ts(
    Data.setup$bigx / 100,
    start = c(1998, 2),
    frequency = 4
  )
  

  #-------------------------------------Estimating country-specific TTVP VARs ---------------------------------------------------------#
  
    shrink.parm <- glob.shrink.list[[as.numeric(comb.grid[zz,1])]]
    #Data.setup$bigx
    gW<-Data.setup$gW
    nhor <- 14
    plag<-1
    cN <- names(Data.setup$gW)
    BVAR <- cmpfun(BVAR)
    

    ext.inst <- FALSE
    
    x <- t(xglobal)
#---------------------------------------------------------------------------------#
    
    #sfInit(parallel=TRUE,cpus=CPU)
    #sfExport(list=list("mlag","BVAR","datahandling","xglobal","gW","Daten","cN","bvartvpm","saves","burns","thin","ext.inst","shrink.parm"))
    predDens <- lapply(1:length(cN),function(i) BVAR(i,gW=gW,bigx=xglobal,Daten,cN,nsave=saves,nburn=burns,thin_chain=thin,ext.inst=ext.inst,parms=shrink.parm))#,c_tau=shrink.parm[[1]],d_tau=shrink.parm[[2]]
    #sfStop()
    
    save(predDens,file="ttvp_gvar_ssr_sp_level.rda")
    q(save="no")
    
    #Creates a set of lists that only include coefficients and VC matrices + Weights
    #Stuff we have to compute once for the full system
    Sigma.posterior <- A.list <- globalG <- list()
    for (i in 1:length(predDens)){ 
      globalG[[i]] <- predDens[[i]]$W
      Sigma.posterior[[i]] <- predDens[[i]]$SIGMApost
      A.list[[i]] <- predDens[[i]]$ALPHA
    }
    
  thin.fac <- round(thin*saves)
  #-------------------------------------Compute IRFs based on sign-restrictions using parallel computing---------------------------------------------#
  IRF_post <-  array(NA,dim=c(nrow(xglobal)-1,ncol(xglobal),nhor+1,thin.fac))
  t.names <- time(xglobal)[-1]
  dimnames(IRF_post) <- list(t.names,colnames(xglobal),c(0:nhor),NULL)

  pb <- txtProgressBar(min = 0, max = thin.fac, style = 3) #start progress bar
  start <- Sys.time()
  for (irep in 1:thin.fac){  
    A.i <- rapply(A.list, classes = 'array', how = 'list', f = function(x) x[,,,irep])
    S.i <- rapply(Sigma.posterior, classes = 'array', how = 'list', f = function(x) x[,,,irep])
    #sfInit(parallel=TRUE,cpus=CPU) #can put whatever here
    #sfExport(list=list("irf.mcmc","get.irfa.t","xglobal","nhor","irf","ext.inst","sign.irf","split.function",'A.i','globalG','S.i'))
    IRF.big<- sfLapply(1:(nrow(xglobal)-1),function(i) get.irfa.t(i, A.i, S.i  , t(xglobal), globalG, horz=nhor))
    #sfStop()
    
    for (ss in 1:(nrow(xglobal)-1)){
      IRF_post[ss,,,irep] <- IRF.big[[ss]]$IRF_post
    }
    setTxtProgressBar(pb, irep)
  }  
  end <- Sys.time()
  print(end-start)   
  
  save(IRF_post,file="irf_ttvp_gvar_ssr.rda")
  q(save="no")
  
  #-------------------------------------Do a lot of plots ---------------------------------------------------------#
  
  foldername <- paste0("Results/Results_","shrink",names(glob.shrink.list)[[zz]],sep="")
  dir.create(foldername, showWarnings = T)

  
  IRF.mean <- apply(IRF_post,c(1,2,3),median,na.rm=TRUE)
  IRF.low <- apply(IRF_post,c(1,2,3),quantile,0.16,na.rm=TRUE)
  IRF.high <- apply(IRF_post,c(1,2,3),quantile,0.84,na.rm=TRUE)
  
  
  pdf(paste0(foldername,"/irfs.pdf"),width = 8,height=22)
  for (jj in 1:length(cN)){
    #Slct. plots across countries
    sl.countries <- grep(cN[[jj]],colnames(xglobal))
    names.sl <- colnames(xglobal)[sl.countries]
    par(mfrow=c(length(sl.countries),4))
    for (ss in 1:length(sl.countries)){
      slct <- sl.countries[[ss]]
      for (fhorz in seq(1,nrow(xglobal)-1,length.out=4)){
        matplot(cbind(IRF.mean[fhorz,slct,],IRF.high[fhorz,slct,],IRF.low[fhorz,slct,]),type="l",col=c("black","black","black"),lty=c(1,2,2),ylab=names.sl[[ss]])
        abline(h=0,col="red")
      }
    }
  }
  dev.off()
  
  
  pdf(paste0(foldername,"/irfs_sliced.pdf"),width = 8,height=22)
  for (jj in 1:length(cN)){
    #Slct. plots across countries
    sl.countries <- grep(cN[[jj]],colnames(xglobal))
    names.sl <- colnames(xglobal)[sl.countries]
    par(mfrow=c(length(sl.countries),4))
    for (ss in 1:length(sl.countries)){
      slct <- sl.countries[[ss]]
      for (fhorz in c(1,3,5,9)){
        matplot(cbind(IRF.mean[,slct,fhorz],IRF.high[,slct,fhorz],IRF.low[,slct,fhorz]),type="l",col=c("black","black","black"),lty=c(1,2,2),ylab=names.sl[[ss]])
        abline(h=0,col="red")
      }
    }
  }
  dev.off()

  
  irf.list <- list(median=IRF.mean, low=IRF.low, high=IRF.high)
  save(irf.list,file=paste0(foldername,"/irfs.RData",sep=""))
  
