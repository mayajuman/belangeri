library(mice)
library(shapes)

#set up data
setwd("~/Desktop/belangeri")
bel <- read.csv("belangeri.csv")

#dealing with NAs outside of measurements
bel$Island.Area[is.na(bel$Island.Area)] <- ""
bel$Mainland.Distance[is.na(bel$Mainland.Distance)] <- ""
bel$Sea.Depth[is.na(bel$Sea.Depth)] <- ""
bel$Island.Area <- as.factor(bel$Island.Area)
bel$Mainland.Distance <- as.factor(bel$Mainland.Distance)
bel$Sea.Depth <- as.factor(bel$Sea.Depth)
bel[bel == 0.00] <- NA

bel <- subset(bel, select = c(Museum, MuseumCat., Genus, species, 
                              Subspecies, CPL, CIL, UTL, MTL, EPL, PPL, 
                              EB, MB, LB, LIB, ZB, BB, LPL, CNL, PBPL, 
                              LTPL, LCH, MH, MCH, MCW, MCIL, LTL, Sex, 
                              Lat, Long, Locality, General.Locality, Island, 
                              Island.Area, Mainland.Distance, Sea.Depth, 
                              Collection.Year, Collection.Date, TL, TV, HF, 
                              E, WT, Type.information, elev))

table(bel$Island) #111 island, 729 mainland individuals
(colMeans(is.na(bel[,6:27])))*100

bel$Sex[bel$Sex == "sex unknown"] <- NA
bel$Sex <- droplevels(bel$Sex)
summary(bel$Sex)    #363 F, 440 M, 37 unknown

#potentially problematic individuals (see outlier plot below):
bel[817,] #no obvious issues
bel[199,] #specimen is on the larger side
bel[709,] #specimen is on the larger side but measurements are almost complete
bel <- bel[-185,] #no obviously weird measurements, but missing MANY

#drop 2 specimens without locality
bel <- bel[-c(94,142),]

#dump empty rows
bel <- bel[-which(bel$MuseumCat. == 35822),]
bel <- bel[-which(bel$MuseumCat. == 1990512),]

#cut out LACM:
#bel <- bel[-c(1,836),]

sum(is.na(bel[,6:27]))/prod(dim(bel[,6:27])) #18.63% missing data overall

#<20% missing data cutoff; conservative but we lose more cranial variables
#bel20 <- bel[,6:27][,colMeans(is.na(bel[,6:27]))<0.20]

#<25% missing data cutoff
bel25 <- bel[,6:27][,colMeans(is.na(bel[,6:27]))<0.25]
sum(is.na(bel25))/prod(dim(bel25)) #15% missing data now; better


######IMPUTATION FUNCTIONS (Clavel et al. 2014)

#agglomerate.data: function to store the m imputed datasets and to compute the averaged dataset

agglomerate.data<-function(data,imp,Mimp,Method="mice"){
  
  Moy<-Mimp+1
  redata<-as.matrix(data)
  ximp<-array(redata,dim=c(nrow(redata),ncol(redata),Moy))
  #####################
  
  if(any(is.na(redata))==TRUE){
    if(Method=="mice" || Method=="amelia" || Method=="missmda" || Method=="hmisc" || Method=="norm"){
      #####################MICE
      if(Method=="mice"){
        
        for(i in 1:Mimp){
          ximp[,,i]<-as.matrix(complete(imp,i))
          
        }
        ##Averaged dataset
        ximp[,,Moy]<-apply(ximp[,,1:Mimp],c(1,2),mean)
      }
      
    }else{
      ####################
      ##Warning messages
      cat("Error! You must indicate if you are using Mice, Amelia, missMDA, NORM, or Hmisc package","\n")
    }}else{
      ## Warning messages
      cat("There is no missing value in your dataset","\n")
    }
  #return(ximp)
  tabM<-ximp[,,Moy]
  colnames(tabM)<-colnames(redata)
  list("ImpM"=tabM,"Mi"=ximp[,,1:Mimp],"nbMI"=Mimp, "missing"=as.data.frame(redata))
}

#ELLI: function to draw confidence ellipses
ELLI<-function(x,y,conf=0.95,np)
{centroid<-apply(cbind(x,y),2,mean)
ang <- seq(0,2*pi,length=np)
z<-cbind(cos(ang),sin(ang))
radiuscoef<-qnorm((1-conf)/2, lower.tail=F)
vcvxy<-var(cbind(x,y))
r<-cor(x,y)
M1<-matrix(c(1,1,-1,1),2,2)
M2<-matrix(c(var(x), var(y)),2,2)
M3<-matrix(c(1+r, 1-r),2,2, byrow=T)
ellpar<-M1*sqrt(M2*M3/2)
t(centroid + radiuscoef * ellpar %*% t(z))}

#plot.MI: function to Plot MI confidence ellipses using procruste superimposition

plot.MI<-function(IM,symmetric=FALSE,DIM=c(1,2),scale=FALSE,web=FALSE,ellipses=TRUE,...){
  if(any(is.na(IM$ImpM)==TRUE))
  { cat("There is still missing values in the imputed dataset, please check your imputation")
    break
  }else{
    Mo<-IM$nbMI+1
    pcaM<-princomp(IM$ImpM)
    cpdimM<-as.matrix(pcaM$scores[,DIM])   
    opa<-array(cpdimM,dim=c(nrow(cpdimM),ncol(cpdimM),Mo)) 
    for(i in 1:IM$nbMI){
      pca<-princomp(IM$Mi[,,i])
      opa[,,i]<-as.matrix(pca$scores[,DIM])
    }
    if(symmetric==TRUE){
      for (i in 1:IM$nbMI+1){
        trace<-sum(opa[,,i]^2)
        opa[,,i]<-opa[,,i]/sqrt(trace) 
      }
    }
    #Ordinary Procrustes Analysis (library(shapes))
    for(k in 1:IM$nbMI){
      analyse<-procOPA(opa[,,Mo],opa[,,k], reflect=TRUE)
      opa[,,k]<-analyse$Bhat
    }
    opa[,,Mo]<-analyse$Ahat
    #Principal component explained variance
    pvar<-pcaM$sdev^2
    tot<-sum(pvar)
    valX<-pvar[DIM[1]]
    valY<-pvar[DIM[2]]
    valX<-round(valX*100/tot,digits=2)
    valY<-round(valY*100/tot, digits=2)
    ######################## Plot function
    op <- par(no.readonly=TRUE)
    if(scale==TRUE){
      plot(opa[,1,Mo],opa[,2,Mo], type="p", pch=3, col=c(as.factor(ifelse(complete.cases(IM$missing) ==T, 1, 5))),lwd=1,xlim=range(opa[,1,Mo]),ylim=range(opa[,1,Mo]),xlab=paste("DIM",DIM[1],valX,"%",sep=" "),ylab=paste("DIM",DIM[2],valY,"%",sep=" "))
    }
    if(scale==FALSE){
      plot(opa[,1,Mo],opa[,2,Mo], type="p", pch=3, col=c(as.factor(ifelse(complete.cases(IM$missing) ==T, 1, 5))),lwd=1,xlab=paste("DIM",DIM[1],valX,"%",sep=" "),ylab=paste("DIM",DIM[2],valY,"%",sep=" "))
    }
    title("MI effect on Multivariate Analysis", font.main=3, adj=1)
    ## Store row names
    NR<-IM$missing
    rownames(IM$missing)<-NULL
    ##
    if(ellipses==TRUE){                                      
      coul<-as.numeric(rownames(IM$missing[complete.cases(IM$missing),]))
      for (j in coul){
        lines(ELLI(opa[j,1,],opa[j,2,],np=Mo), col="black", lwd=1)}
      coul<-as.numeric(rownames(IM$missing[!complete.cases(IM$missing),]))
      for (j in coul){
        lines(ELLI(opa[j,1,],opa[j,2,],np=Mo), col="red", lwd=1)}
    }else{ points(opa[,1,],opa[,2,],cex=0.5) }
    if(web==TRUE){
      coul<-as.numeric(rownames(IM$missing[complete.cases(IM$missing),]))
      for (j in coul){ 
        for(f in 1:IM$nbMI){
          segments(opa[j,1,Mo],opa[j,2,Mo], opa[j,1,f],opa[j,2,f], col="black", lwd=1) }
      } 
      coul<-as.numeric(rownames(IM$missing[!complete.cases(IM$missing),]))
      for (j in coul){ 
        for(f in 1:IM$nbMI){
          segments(opa[j,1,Mo],opa[j,2,Mo], opa[j,1,f],opa[j,2,f], col="red", lwd=1)}
      } 
      points(opa[,1,],opa[,2,],cex=0.5) 
    } 
    nom<-rownames(NR)
    text(opa[,1,Mo],opa[,2,Mo],nom, pos=1)
    
    abline(h=0,v=0, lty=3)
    par(xpd=TRUE)  # Do not clip to the drawing area
    lambda <- .025
    legend(par("usr")[1], (1 + lambda) * par("usr")[4] - lambda * par("usr")[3],c("Complete", "Missing"), xjust = 0, yjust = 0,lwd=3, lty=1, col=c(par('fg'), 'red'))
    par(op)      
  }
}

######IMPUTING

#only variables with <25% missing data

m = 100

imp_bel25 <- mice(bel25, m = m, method="pmm")
IM_mice25<-agglomerate.data(data=bel25, imp=imp_bel25, Mimp=100, Method="mice")

#IM_mice25$ImpM   # average dataset
#IM_mice25$Mi #list with the m imputed datasets

write.csv(IM_mice25$ImpM, file="imputed belangeri with 2 LACM.csv")
#imp <- read.csv("imputed belangeri.csv")


######CHECK FOR OUTLIERS
#procrustes superimposition of the m imputed datasets onto the principal components calculated from the average MI-dataset.
par(mar=c(3,3,3,0))
plot.MI(IM_mice25, symmetric=TRUE, DIM=c(1,2), web=FALSE, ellipses=TRUE)


###remove outliers + repeat if necessary


###references

#Clavel, J., Merceron, G., & Escarguel, G. Syst. Biol. 63, 203–218 (2014).

