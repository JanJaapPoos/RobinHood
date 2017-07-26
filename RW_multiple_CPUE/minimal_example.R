## ----echo=FALSE----------------------------------------------------------
library(knitr);  library(TMB); library(corrplot);
##
opts_chunk$set(size="footnotesize")
## set the working directory
##setwd('w:/imares/ijmuiden/afdeling/projecten/data poor mixed fisheries/robin hood/include data rich'); 
setwd('~/wur/N/Projecten/data poor mixed fisheries/robin hood/RW_multiple_CPUE'); 
compile("sprobinhood_model.cpp")
dyn.load(dynlib("sprobinhood_model"))

cuslogit <- function(x){return(log(x/(2-x)))}
invcuslogit <- function(x) {return( 2/(1+exp(-x)))} 

## ------------------------------------------------------------------------
#years
years <- 1976:2016
## indices first is tur, second is brill (one more I obs than C obs per spec )
Im <-   read.csv("~/wur/N/Projecten/data poor mixed fisheries/Robin Hood/data/Itur_bll.csv")
# Catches
Cm <- read.csv("~/wur/N/Projecten/data poor mixed fisheries/Robin Hood/data/Cple_sol_tur_bll_lem.csv")
## F estimates from data rich assessments (ple, sol),same length as I )
Fm <- read.csv("~/wur/N/Projecten/data poor mixed fisheries/Robin Hood/data/Fple_sol_cod.csv")
## print data for checks
merge(merge(Fm,Cm),Im, all=T)

#plot Fs
#plot(x=Fm[,1], y=Fm[,2], ylim=c(0,1.2), xlab="Year", ylab="Fishing mortality (year-1))", type="b", col="black", pch=15, yaxs="i", las=1, mgp=c(3,0.7,0) )
#lines(x=Fm[,1], y=Fm[,3],  type="b", col="red", pch=16)
#lines(x=Fm[,1], y=Fm[,4],  type="b", col="red", pch=16)
#legend("topright",legend=c("European plaice", "sole"), col=c("black","red"), pch=c(15,16))

                                        #plot catches
#options(scipen=8)
#plot(x=Cm[,1],  y=Cm[,2] + Cm[,3] , ylim=c(600,600000 ),  xlab="Year", ylab="Catch (tonnes)", type="b",log="y", col="black", pch=15, yaxs="i", las=1, mgp=c(3,0.7,0) )
#lines(x=Cm[,1], y=Cm[,4] + Cm[,5],  type="b", col="red", pch=16)
#lines(x=Cm[,1], y=Cm[,6],  type="b", col="blue", pch=17, cex=1.2)
#lines(x=Cm[,1], y=Cm[,7],  type="b", col="magenta", pch=18)
#grid()
#legend("topright",legend=c("European plaice", "sole", "turbot", "brill"), col=c("black","red", "blue","magenta"), pch=c(15,16,17,18))
#box()

##remove data from start of TS
#for( startcut in 1:20)
startcut <- 0 #it works with startcut 13

#do actual removal NOTE THAT FOR LEM BUT NOW BOTH INDICES
Im_cut <- as.matrix(Im[(startcut+1):(dim(Im)[1]),c(2,3,4,5)])
#note that we make catches from L + D for Lem here
Cm_cut <- as.matrix(Cm[(startcut+1):(dim(Cm)[1]),c(6,7,8,9)])
Cm_cut[,3] <- Cm_cut[,3] +   Cm_cut[,4]
Cm_cut <- Cm_cut[,1:3]
Fm_cut <- as.matrix(Fm[(startcut+1):(dim(Fm)[1]),2:4])

# number of stocks
mc  <- dim(Cm_cut)[2]
# number of indices
mi <- dim(Im_cut)[2]
# matrix of number of indices for each species  (the sum of this matrix must be equal to the columns of Imat and Ipresent)
nIndices <- matrix(c(1,1,2))

## matrix of ones and zeros for whether index is present for that year
Ipresent <- matrix(as.numeric(!is.na(Im_cut)), ncol = mi)
## calc sigmas for random walks data rich (sdrwdr)
sdrwdr <- matrix(apply(Fm_cut,2, FUN=function(x){sd(diff(x),na.rm=T)}))
# set bounds lm         lK         lq               lsdi           lsdc       lsdrw     logitalpha
lower <- c(7,7,7 ,    7,7,7,     -9.5,-9.5,-9.5,-9.5,  -5,-5,-5,-5,  -10,-10,-7,  -5,-5,-5,  -20,-20,-20 )
upper <- c(13,13,13,  15,15,15,    1,1,1,1,          2,2,2,2,      1,1,1,      1,1,1,     20,20,20)

corr.mode <- F

## make  corr matrix set off diag to 0 (default
corrmatrix <- matrix(0.0, nrow=6,ncol=6)

if(corr.mode == T){
    ## make  corr matrix set according to Figure in into (assumed that first DP, then DR)
    corrmatrix[2,1] <- corrmatrix[1,2] <-  0.99 #TUR-BLL
    corrmatrix[3,1] <- corrmatrix[1,3] <-  0.13 #TUR-LEM
    corrmatrix[4,1] <- corrmatrix[1,4] <-  0.94 #TUR-PLE
    corrmatrix[5,1] <- corrmatrix[1,5] <-  0.99 #TUR-SOL
    corrmatrix[6,1] <- corrmatrix[1,6] <- -0.06 #TUR-COD
    corrmatrix[2,3] <- corrmatrix[3,2] <-  0.07 #BLL-LEM
    corrmatrix[2,4] <- corrmatrix[4,2] <-  0.93 #BLL-PLE
    corrmatrix[2,5] <- corrmatrix[5,2] <-  0.99 #BLL-SOL
    corrmatrix[2,6] <- corrmatrix[6,2] <- -0.12 #BLL-COD
    corrmatrix[3,4] <- corrmatrix[4,3] <-  0.44 #LEM-PLE
    corrmatrix[3,5] <- corrmatrix[5,3] <- -0.02 #LEM-SOL
    corrmatrix[3,6] <- corrmatrix[6,3] <-  0.94 #LEM-COD
    corrmatrix[4,5] <- corrmatrix[5,4] <-  0.89 #PLE-SOL
    corrmatrix[4,6] <- corrmatrix[6,4] <-  0.24 #PLE-COD
    corrmatrix[5,6] <- corrmatrix[6,5] <- -0.18 #SOL-COD
    #  corrmatrix <- abs(corrmatrix)
}

corrmatrix <- corrmatrix *0.81

diag(corrmatrix) <- 1


obj <- MakeADFun(
    data = list(
        Imat = Im_cut,
        Cmat = Cm_cut,
        Fmat = Fm_cut,
        nIndices = nIndices,
        Ipresent = Ipresent,
        corrmatrix = corrmatrix,
        sdrwdr = sdrwdr ),
    parameters = list(
        logm  = c(8.4,7.4,8.3),
        logK = c(10.8,9.05,10.3),
        logq = c(-5.9,-4.6,-7.5,-8),
        logsdi = c(-2.4,-2.2,-1.,-1.),
        logsdc = c(-2.5,-2.5,-2.5),
        logsdrw = matrix(log(0.05),nrow=mc ),
        logitalpha = c(0.5,1.1,-0.6),
        ## CM: changed number of rows
        logitU = cbind(cuslogit(0.5),cuslogit(0.5),cuslogit(0.5), cuslogit(Fm_cut))),
    random = "logitU", 
    map = list(logsdc = rep(as.factor(NA),3 ) ),
    DLL = "sprobinhood_model",
    silent = F)

## fit    
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, lower=lower, upper=upper,  control = list(iter.max = 1e7, eval.max = 1e7))
(rep <- sdreport(obj))

par.corr <- cov2cor(rep$cov.fixed)
corrplot.mixed(par.corr, lower = "ellipse", upper = "number")

rep.summ <- summary(rep)

Bhat   <- matrix( rep.summ[rownames(rep.summ) == "B", "Estimate"], nc = mc)
Blow   <- Bhat - 2 * matrix(rep.summ[rownames(rep.summ) == "B", "Std. Error"],  nc = mc)
Bhigh  <- Bhat + 2 * matrix(rep.summ[rownames(rep.summ) == "B", "Std. Error"], nc = mc)
BoBMSYhat  <- matrix( rep.summ[rownames(rep.summ) == "BoBMSY", "Estimate"], nc = mc)
BoBMSYlow  <- BoBMSYhat - 2 * matrix(rep.summ[rownames(rep.summ) == "BoBMSY", "Std. Error"],  nc = mc)
BoBMSYhigh <- BoBMSYhat + 2 * matrix(rep.summ[rownames(rep.summ) == "BoBMSY", "Std. Error"], nc = mc)
FoFMSYhat  <- matrix( rep.summ[rownames(rep.summ) == "FoFMSY", "Estimate"], nc = mc)
FoFMSYlow  <- FoFMSYhat - 2 * matrix(rep.summ[rownames(rep.summ) == "FoFMSY", "Std. Error"],  nc = mc)
FoFMSYhigh <- FoFMSYhat + 2 * matrix(rep.summ[rownames(rep.summ) == "FoFMSY", "Std. Error"], nc = mc)
Uhat  <- matrix( rep.summ[rownames(rep.summ) == "logitU", "Estimate"], nc = mc+3)
Ulow  <- invcuslogit(Uhat - 2 * matrix(rep.summ[rownames(rep.summ) == "logitU", "Std. Error"], nc = mc+3))
Uhigh <- invcuslogit(Uhat + 2 * matrix(rep.summ[rownames(rep.summ) == "logitU", "Std. Error"], nc = mc+3))
Uhat <- invcuslogit(Uhat)
BMSY <-  0.5* exp( matrix( rep.summ[rownames(rep.summ) == "logK", "Estimate"], nc=mc))

## check fit to the indices
Ihat <- matrix( rep.summ[rownames(rep.summ) == "Ihat", "Estimate"], nc = mi)
## ses for mean CIs
Ihat.se <- matrix( rep.summ[rownames(rep.summ) == "Ihat", "Std. Error"], nc = mi)

##plot data and what is removed
colscheme <- c("black","red", "blue","magenta", "green", "green")
pchscheme <- c(1,1,1,16)
par(mfrow=c(1,3), mar=c(5.1,4.1,1,1))
plot(x=years,y=Im[,2], log="y", type="b", ylim=c(0.4,150), col="black",xlab="Year", ylab="Index values", las=1, yaxs="i",panel.first=grid() )
for(ii in 1:4){
  lines(x=years,y=Im[,ii+1], type="b", col=colscheme[ii+2], pch=pchscheme[ii])
  lines(x=years,y=c(rep(NA,startcut ), Ihat[,ii] ), lwd=2, col=colscheme[ii+2])
  lines(x=years,y=c(rep(NA,startcut ), Ihat[,ii] +2 * Ihat.se[,ii]), lty=2, col=colscheme[ii+2]) 
  lines(x=years,y=c(rep(NA,startcut ), Ihat[,ii] -2 * Ihat.se[,ii]), lty=2, col=colscheme[ii+2])
}
if (startcut >0) abline(v=years[1]+startcut+0.5, lty=2)
#legend("topleft", legend=c("Turbot","Brill"), border=NA, col=c("black","red"), pch=1)


##
plot(x=years,y=c(rep(NA,startcut),FoFMSYhat[,1]), type="l", ylim=c(0,4), col="black", xlab="Year", ylab="F / FMSY",  lwd=2, las=1, yaxs="i", panel.first=grid())
for(ii in 1:3){
    lines(x=years,y=c(rep(NA,startcut ), FoFMSYhat[,ii]),  lwd=2, col=colscheme[ii+2])
    lines(x=years,y=c(rep(NA,startcut ), FoFMSYlow[,ii]),  lty=2, col=colscheme[ii+2]) 
    lines(x=years,y=c(rep(NA,startcut ), FoFMSYhigh[,ii]), lty=2, col=colscheme[ii+2])
}
if (startcut >0) abline(v=years[1]+startcut+0.5, lty=2)
lines(x=years,y=c(rep(NA,startcut), Uhat[,4]/0.21), lwd=2, col=colscheme[1]) #ple
lines(x=years,y=c(rep(NA,startcut), Uhat[,5]/0.20), lwd=2, col=colscheme[2]) # sol
lines(x=years,y=c(rep(NA,startcut), Uhat[,6]/0.36), lwd=2, col="orange") #cod

abline(h=1)
legend("topright", legend=c("European plaice","sole", "cod", "turbot","brill","lemon sole"), border=NA, lty=1, col=c("black","red","orange", "blue", "magenta", "green"),  bty="n", pch=1)

##
plot(x=years,y=c(rep(NA,startcut),BoBMSYhat[,1]), type="l", ylim=c(0,2), col="black", xlab="Year", ylab="B / BMSY",  lwd=2, las=1, yaxs="i", panel.first=grid())
for(ii in 1:3){
    lines(x=years,y=c(rep(NA,startcut ), BoBMSYhat[,ii]), type="l", lwd=2, col=colscheme[ii+2])
    lines(x=years,y=c(rep(NA,startcut ), BoBMSYlow[,ii]),  lty=2,  col=colscheme[ii+2])
    lines(x=years,y=c(rep(NA,startcut ), BoBMSYhigh[,ii]), lty=2,  col=colscheme[ii+2])
}
if (startcut >0) abline(v=years[1]+startcut+0.5, lty=2)
abline(h=1)


## plot U vals for DL stocks in assessment and for DR inputs 
#plot(x=years,y=c(rep(NA,startcut),Uhat[,1]), type="l", ylim=c(0,1.2), col="black", ylab="Harvest ratios", las=1,lwd=2,  yaxs="i", panel.first=grid())
#for(ii in 1:3){
#    lines(x=years,y=c(rep(NA,startcut), Uhat[,ii] ),lwd=2, col=colscheme[ii+2])
#    lines(x=years,y=c(rep(NA,startcut), Ulow[,ii] ), lty=2, col=colscheme[ii+2])
#    lines(x=years,y=c(rep(NA,startcut), Uhigh[,ii]), lty=2, col=colscheme[ii+2])
#}
#abline(v=years[1]+startcut+0.5, lty=2)
#lines(x=years,y=c(rep(NA,startcut), Uhat[,4]),lwd=2, col=colscheme[1])
#lines(x=years,y=c(rep(NA,startcut), Uhat[,5]),lwd=2, col=colscheme[2])
#lines(x=years,y=c(rep(NA,startcut), Uhat[,6]),lwd=2, col="orange")

#legend("topleft", legend=c("Turbot","Brill"), border=NA, col=c("black","red"), pch=1)

## plot Biomass estimates
#plot(x=years,y=c(rep(NA,startcut),Bhat[,1]), type="l", ylim=c(200,100000), log="y", col="black", ylab="Biomass (tonnes)", las=1, lwd=2,  yaxs="i", panel.first=grid())
#for(ii in 1:3){
#    lines(x=years,y=c(rep(NA,startcut), Bhat[,ii]), type="l", lwd=2,  col=colscheme[ii+2])
#    lines(x=years,y=c(rep(NA,startcut), Blow[,ii]),  lty=2, col=colscheme[ii+2])
#    lines(x=years,y=c(rep(NA,startcut), Bhigh[,ii]), lty=2, col=colscheme[ii+2])
#}
#abline(v=years[1]+startcut+0.5, lty=2)
#legend("topleft", legend=c("Turbot","Brill"), border=NA, col=c("black","red"), pch=1)













