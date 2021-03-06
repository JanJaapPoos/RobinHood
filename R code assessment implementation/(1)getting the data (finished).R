#getting the data
##---data input---------------------------------------------------------------------
#define start year
startyear <- 1975
endyear <- 2018

#####################################################################################
## Read Catches and indices data for data poor
#####################################################################################
#rootPath <- "m:/robinhood/"
rootPath <- "C:/Users/jappi/Bureaublad/robinhood"


#####################################################################################
#depending on wether discard estimates should be mean proportion landings and discards or also from linear regression:
#Two choices: either dab and witch flounder with lm and rest proportions (recalculated catches data and indices)
# or all discards reconstructions as proportions (recalculated catches data and indices_only_mean_prop)
#####################################################################################
filename <- "recalculated catches data and indices_only_mean_prop"

alldat <-   read.csv(file.path(rootPath,"Data/landings and indices",paste0(filename,".csv")), stringsAsFactors = F)

names(alldat)[1] <- "years"

# start all analyses in 1975, that's the first year with Bll landings
alldat <- alldat[!is.na(alldat$years),]
alldat <- alldat[alldat$years>(startyear-1),]
alldat <- alldat[alldat$years<(endyear+1),]
alldat[alldat==""] <- NA

# Extract indices
Im <- alldat[,names(alldat) %in% c("years","Ibll","Itur", "Iwit","Iwit.1","Idab","llem.3","Ilem")]
# some indices are in kg /day (Witch.1 (IBTSQ3), Ibll (commercial), and Turbot. Note that the units for the turbot index are kg per kwday.
#So there we need an additional conversion. 
Im$Iwit.1 <- as.numeric(Im$Iwit.1)/24
Im$Ibll  <- as.numeric(Im$Ibll)/24
Im$Itur <- (as.numeric(Im$Itur) *1491) /24

Im_lemonsole<-read.csv(file.path(rootPath,"Data/landings and indices/lemon sole index calcs/Resulting lemon sole indices from IBTS.csv"), stringsAsFactors = F)
Im_lemonsole<-Im_lemonsole[,2:4]
Im_lemonsole<-tidyr::spread(Im_lemonsole, Quarter, kgperhour.x)
top_Im_lem<-cbind(c(1975:1979),rep(NA,5),rep(NA,5))
colnames(top_Im_lem)<-names(Im_lemonsole)
Im_lemonsole<-rbind(top_Im_lem,Im_lemonsole)
Im$Ilem  <-Im_lemonsole[1:44,2]
Im$llem.3<-Im_lemonsole[1:44,3]
Im
colnames(Im)<-c("years","Itur","Ibll","Iwit","Iwit.1","Idab","Ilem","Ilem.1")
Im
# Extract catches
Lm <- alldat[,names(alldat) %in% c("years","Ltur","Lbll", "Lwit","Ldab","Llem") ]
names(Lm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")
Dm <- alldat[,names(alldat) %in% c("years","Dtur","Dbll", "Dwit","Ddab","Dlem") ]
names(Dm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")
Cm<- alldat[,names(alldat) %in% c("years","Ctur","Cbll", "Cwit","Cdab","Clem") ]
names(Cm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")

#######################################################################################
## Read F estimates from data rich assessments (ple, sol),same length as I )
#######################################################################################
allF <- read.csv(file.path(rootPath,"Data/F data rich/ple_sol_lem_wit_tur_cod_fishing mortalities.csv"), sep=";",stringsAsFactors = F)

names(allF)[1] <- "years" 
allF <- allF[allF$years>(startyear-1),]
allF <- allF[allF$years<(endyear+1),]
allF[allF==""] <- NA
allF$Plaice<-as.numeric(allF$Plaice)
allF$Sole<-as.numeric(allF$Sole)
allF$Cod<-as.numeric(allF$Cod)
allF$Whiting<-as.numeric(allF$Whiting)
allF$Haddock<-as.numeric(allF$Haddock)
## extract Fs
Fm <- allF[,names(allF) %in% c("years","Plaice","Sole", "Cod","Whiting","Haddock") ]


#######################################################################################
# Read landings data by gear rect spec from STECF website
# this data is downloaded from STECF website
##############################################################################

landings <- read.csv(file.path(rootPath,"data/stecf data/map__landings_by_rectangle_data.csv"), stringsAsFactors = F)

save(Cm, Fm, Im, Lm, Dm, file=file.path(rootPath,"Data/intermediate results/original data.Rdata"))



