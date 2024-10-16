
#getting the data
library(readxl)
library(dplyr)
##---data input---------------------------------------------------------------------
#define start year
startyear <- 1975
endyear <- 2023

#####################################################################################
## Read Catches and indices data for data poor
#####################################################################################
#rootPath <- "m:/robinhood/"
rootPath <- "C:/Users/bleij009/OneDrive - Wageningen University & Research/new_git_clones/RobinHood/"


#####################################################################################
#depending on wether discard estimates should be mean proportion landings and discards or also from linear regression:
#Two choices: either dab and witch flounder with lm and rest proportions (recalculated catches data and indices)
# or all discards reconstructions as proportions (recalculated catches data and indices_only_mean_prop)
#####################################################################################
filename <- "recalculated catches data and indices_only_mean_prop"


extract_yrs <- function(x, startyear, endyear){
    sel<-as.vector(x[,1])$Species %in% as.character((seq(startyear,endyear)))
    
    df_extract <- x[sel,]
    df_extract <- as.data.frame(df_extract)
    df_extract <- data.frame(lapply(df_extract, function(x) as.numeric(as.character(x))))
    return(df_extract)
}

#turbot
Tur <- read_excel(
           path = paste0(rootPath,"/Data/landings and indices/",filename,".xlsx"),
           range = "B2:F56",
           col_names =T,
           sheet = 'Turbot'
           )
tur_extract <- extract_yrs(Tur, startyear = startyear , endyear = endyear)

head(Tur)
names(tur_extract) <- c("Year","Landings", "Discards", "Catches", "Survey1")

#lemon sole
Lem <- read_excel(
           path = paste0(rootPath,"/Data/landings and indices/",filename,".xlsx"),
           range = "S3:AD82",
           col_names =T,
           sheet = 'lemon sole'
           )

lem_extract <- extract_yrs(Lem, startyear = startyear , endyear = endyear)
lem_extract <- lem_extract[,c(1,3,5,6:12)]

head(Lem)
names(lem_extract) <- c("Year","Landings", "Discards", "Catches", "Survey1","Survey2","Survey3","Survey4","Survey5","Survey6")

#witch
Wit <- read_excel(
           path = paste0(rootPath,"/Data/landings and indices/",filename,".xlsx"),
           range = "B3:L82",
           col_names =T,
           sheet = 'Witch'
           )

wit_extract <- extract_yrs(Wit, startyear = startyear , endyear = endyear)
wit_extract <- wit_extract[,c(1,3,6, 7:11)]

head(Wit)
names(wit_extract) <- c("Year","Landings", "Discards", "Catches", "Survey1","Survey2","Survey3","Survey4")

#dab
Dab <- read_excel(
           path = paste0(rootPath,"/Data/landings and indices/",filename,".xlsx"),
           range = "B2:H81",
           col_names =T,
           sheet = 'Dab'
           )

dab_extract <- extract_yrs(Dab, startyear = startyear , endyear = endyear)
dab_extract <- dab_extract[,c(1,3:7)]

head(Dab)
names(dab_extract) <- c("Year","Landings", "Discards", "Catches", "Survey1","Survey2")

#brill
Bll <- read_excel(
           path = paste0(rootPath,"/Data/landings and indices/",filename,".xlsx"),
           range = "B3:J82",
           col_names =T,
           sheet = 'brill'
           )

bll_extract <- extract_yrs(Bll, startyear = startyear , endyear = endyear)

bll_extract <- bll_extract[,c(1,3:9)]

head(Bll)
names(bll_extract) <- c("Year","Landings", "Discards", "Catches", "Survey1","Survey2", "Survey3","Survey4" )


######
extract_cols <- function(specs, colname) {
    # Dynamically construct the data frame name
    df_name <- paste0(specs, "_extract")
    
    # Retrieve the data frame
    int_df <- get(df_name)
    
    # Extract the specified column
    int_vec <- int_df  %>%  select(starts_with(colname))
    
    return(int_vec)
}

#create Lm
Lm <- unname(as.data.frame(lapply(c("tur","bll", "wit","dab","lem"), function(x) extract_cols(x, "Landings"))))
Lm$years <- startyear:endyear
names(Lm) <- c("Ltur","Lbll", "Lwit","Ldab","Llem","years")
Lm <- Lm[,c(6,1:5)]
names(Lm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")

#create Dm
Dm <- unname(as.data.frame(lapply(c("tur","bll", "wit","dab","lem"), function(x) extract_cols(x, "Discards"))))
Dm$years <- startyear:endyear
names(Dm) <- c("Dtur","Dbll", "Dwit","Ddab","Dlem","years")
Dm <- Dm[,c(6,1:5)]
names(Dm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")

#create Cm
Cm <- unname(as.data.frame(lapply(c("tur","bll", "wit","dab","lem"), function(x) extract_cols(x, "Catches"))))
Cm$years <- startyear:endyear
names(Cm) <- c("Ctur","Cbll", "Cwit","Cdab","Clem","years")
Cm <- Cm[,c(6,1:5)]
names(Cm)<-c("years","Turbot","Brill","Witch","Dab","Lemon sole")

#create Index df
Im <- unname(as.data.frame(lapply(c("tur","bll", "wit","dab","lem"), function(x) extract_cols(x, "Survey"))))
Im$years <- startyear:endyear
names(Im) <- c("Itur1","Ibll1","Ibll2","Ibll3","Ibll4","Iwit1","Iwit2","Iwit3","Iwit4","Idab1","Idab2","Ilem1","Ilem2","Ilem3","Ilem4","Ilem5","Ilem6", "years")

#######################################################################################
## Read F estimates from data rich assessments (ple, sol),same length as I )
#######################################################################################
allF <- read.csv(file.path(rootPath,"Data/F data rich/ple_sol_lem_wit_tur_cod_fishing mortalities.csv"),stringsAsFactors = F)

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

save(Cm, Fm, Im, Lm, Dm, file=file.path(rootPath,"Data/intermediate results/original data.Rdata"))