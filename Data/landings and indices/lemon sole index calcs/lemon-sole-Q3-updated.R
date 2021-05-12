##install.packages("remotes")
##install.packages("icesDatras")
##remotes::install_github("DTUAqua/DATRAS/DATRAS")
##remotes::install_github("casperwberg/surveyIndex/surveyIndex")

library(DATRAS)
library(mgcv);
library(parallel);
library(maps); library(mapdata);
library(maptools)
library(surveyIndex);
## Species specific parameters:
cmSize=1;
agesQ3=0:10
years=2004:2019 
setwd("M:\\robinhood\\data\\landings and indices\\lemon sole index calcs")
outFolder <- "."
genus="Microstomus"
bfamily="kitt";

datafile<-paste0("lemonSole-Q3-",max(years),".RData")

if(!file.exists(datafile)){
    IBTS <- getDatrasExchange("NS-IBTS",years=years,quarters=3,strict=FALSE)
    BTS <- getDatrasExchange("BTS",years=years,quarters=2:4,strict=FALSE)
    dAll <- c(IBTS,BTS)
    dAll <-addSpatialData(dAll,"~/Documents/shapefiles/ICES_areas.shp")
    dQ3 <- subset(dAll,Species==paste(genus,bfamily),Year %in% years,HaulVal=="V",StdSpecRecCode==1)
    dAll <- IBTS <- BTS <- NULL
    save(dQ3,file=datafile)
} else {
    load(datafile)
}


dQ3=addSpectrum(dQ3,by=cmSize)



## impute missing depths
summary(dQ3$Depth)

dmodel=gam(log(Depth) ~ s(lon,lat,k=200),data=dQ3[[2]])
sel=subset(dQ3,is.na(Depth))
sel$Depth=0; ## Guard against NA-error
dQ3$Depth[is.na(dQ3$Depth)]=exp(predict(dmodel,newdata=sel[[2]]))
dmodel=NULL
sel=NULL
gc()

dQ3=subset(dQ3,ICESAREA %in% as.character(c("IIIa20","IIIa21","IVa","IVb","IVc","VIId")), Month>6, Month<11)

##tmp=merge(dQ3[[1]],dQ3[[2]][,c("haul.id","ICESAREA")],all.x=TRUE)
##tmp=subset(tmp,!is.na(Age))
##xtabs(~Year+ICESAREA,data=tmp)

dQ3=addWeightByHaul(dQ3)

dQ3.BTS = subset(dQ3,Survey=="BTS")
dQ3.IBTS = subset(dQ3,Survey=="NS-IBTS")

mybubblePlot<-function (d, response = "HaulWgt", scale = NULL, col.zero = "red", 
    pch.zero = "+", ...) 
{
    d[[2]]$resp.var <- d[[2]][[response]]
    if (is.null(scale)) 
        scale = mean(d[[2]]$resp.var, na.rm = TRUE)/max(d[[2]]$resp.var, 
            na.rm = TRUE)
    plot(d$lon, d$lat, type = "n", xlab = "Longitude", ylab = "Latitude",...)
    map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
    points(d$lon, d$lat, pch = 16, cex = scale * sqrt(d[[2]]$resp.var), 
        ...)
    zero = subset(d, resp.var == 0)
    points(zero$lon, zero$lat, pch = pch.zero, col = col.zero)
}

pdf("bubbles-lemonSole-Q3.pdf")
par(mfrow=c(1,2))
mybubblePlot(dQ3.BTS,scale=1/20,ylim=c(51,62),main="BTS Q3")
mybubblePlot(dQ3.IBTS,scale=1/20,ylim=c(51,62),main="IBTS Q3")
dev.off()

## Gear subsetting 
xtabs(~Year+Gear,data=dQ3[[2]])

dQ3$Gear[dQ3$Gear=="BT4S"]="BT4A"
mytab=xtabs(~Gear,dQ3[["HH"]])
goodGears=names(mytab[mytab>120])
dQ3<-subset(dQ3, Gear %in% goodGears)



## Check for enough age data
xtabs(NoAtALK~Year+Survey,data=dQ3[[1]])
## Only enough from 2004 and onwards
dQ3<-subset(dQ3, Year %in% as.character(2004:max(years)))

## But no, all ages are NA in 2004 and 2017...
## Let's require at least 100 age samples per year:
noagesamples<-xtabs(!is.na(Age)~Year,data=dQ3[[1]])
goodYears = names( noagesamples[ noagesamples>100 ] ) 
dQ3 <- subset( dQ3, Year %in% goodYears)

removeAgeNAs<-function(x) {
    x[[1]]=subset(x[[1]],!is.na(x[[1]]$Age))
    x[[1]]=subset(x[[1]],!is.na(x[[1]]$NoAtALK))
    x
}

dQ3=removeAgeNAs(dQ3)

sink("Q3noages.txt")
xtabs(NoAtALK~Year+Age,data=dQ3[[1]])
sink()
####################
## Age-length key 
####################
## Declare settings for ALK model
mf = "" 
ack=TRUE;
useBICs=TRUE;
varCofs=FALSE;
maxKs=50;
mc.cores=1 ## Windows users should use mc.cores=1

add.ALK<-function(d){

    ages=agesQ3
    
    if(d$Quarter[1]=="1"){
        ages=agesQ1
    }
    d[[1]]=subset(d[[1]],Age>=min(ages))
    for(aa in ages){
        d=fixAgeGroup(d,age=aa,n=1,fun=ifelse(aa==min(ages),min,mean))
    }
    
        
    
    d=addSpectrum(d,by=cmSize)
    
    d.ysplit = split(d,d$Year)
    
    d.ALK= mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,mc.cores=mc.cores)

    d.Nage=mclapply(d.ALK,predict,mc.cores=mc.cores)
    for(i in 1:length(d.ALK)) d.ysplit[[i]]$Nage=d.Nage[[i]];
    dd <- do.call("c",d.ysplit)
    dd    
}

dQ3=add.ALK(dQ3)

#######################
## Model
#######################

## Stationary model
modelsStatZ=rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecZ[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(agesQ3))

modelsStatP=rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecP[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(agesQ3))

grid=getGrid(dQ3,nLon=40)

SI = getSurveyIdx(dQ3,ages=agesQ3,myids=grid[[3]],cutOff=0.1,fam="LogNormal",mc.cores=mc.cores,modelZ=modelsStatZ,modelP=modelsStatP)

# Removed due to a change in getSurveyIdx (email from Casper 15/04/20)
# SI$idx[ SI$idx==0] = NA

SI.alt = getSurveyIdxStratMean(dQ3,agesQ3+1)

surveyIdxPlots(SI,dQ3,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,4),mar=c(4,1,1,1)),select=c("index"),plotByAge=FALSE)
surveyIdxPlots(SI,dQ3,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,4),mar=c(4,1,1,1)),select=c("2"),plotByAge=FALSE)
surveyIdxPlots(SI,dQ3,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,4),mar=c(4,1,1,1)),select=c("residuals"),plotByAge=FALSE)
surveyIdxPlots(SI,dQ3,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,4),mar=c(4,1,1,1)),select=c("map"),plotByAge=FALSE,colors=rev(heat.colors(5)),legend=FALSE)

gearEffects = getEffect(SI,dQ3,parName="Gear",cutOff=0.1) ## compared to GOV=1.0

internalCons(SI$idx)

exportSI<-function(x,ages,years,toy,file,nam="",exclude=c()){
  cat(nam,"\n",file=file)
  cat(range(as.numeric(as.character(years))),"\n",file=file,append=TRUE)
  cat("1 1 ",rep(toy,2),"\n",file=file,append=TRUE)
  cat(min(ages),max(ages),"\n",file=file,append=TRUE)
  write.table(round(cbind(1,x[,]),4),file=file,row.names=FALSE,col.names=FALSE,append=TRUE)
}

ynum = as.numeric(as.character(dQ3$Year))

dropAgeGroups=c(1)

exportSI(SI$idx[,-dropAgeGroups],agesQ3[-dropAgeGroups],min(ynum):max(ynum),toy=mean(dQ3[[2]]$timeOfYear,na.rm=TRUE),file="survey-lemonSole-IBTS-BTS-Q3.dat",nam=paste("NS Lemon Sole ; Last age is plus group, calculated",Sys.time())); 

