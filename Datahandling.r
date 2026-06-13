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
