#load proper libraries-----
library(TMB); library(corrplot);library(knitr); library(ggplot2);library(wesanderson); library(cowplot)
library(adnuts); library(tidyr) 

# set the working directory & compile C++ model----
#rootPath <- "m:/robinhood/"
rootPath <- "C:/Users/jappi/Bureaublad/robinhood"

setwd(file.path(rootPath,'code/Robinhood')); 

#compile("sprobinhood_model.cpp")
dyn.load(dynlib("sprobinhood_model"))

load(file.path(rootPath,"Data/intermediate results/original data.Rdata"))

#masterfunctions to keep results within biological limits--------------
source(file=file.path(rootPath,"code/Jasper final/assessment implementation/(0)helper functions.R"))

##---------make objects that correspond with C++
##----------make tamplate of corrmatrix-----------
##remove data from start of TS
# startcut at 3 because whiting Fishing mortality estimate from data rich whiting starts at 1978 
startcut <- 3 #it works with star tcut 3

#do actual removal NOTE THAT FOR LEM BUT NOW BOTH INDICES
Im_cut <- as.matrix(Im[(startcut+1):(dim(Im)[1]),-1])
Cm_cut <- as.matrix(Cm[(startcut+1):(dim(Cm)[1]),-1])

#transpose F values datarich so that they are expressed as fractions of biomass
Fm_cut<- 1-exp(-Fm[2:6])
#here we remove the start of the timeseries
Fm_cut <- as.matrix(Fm_cut[(startcut+1):(dim(Fm_cut)[1]),])

years <- 1975:2023

# number of data limited stocks
print(mc  <- dim(Cm_cut)[2])
# number of data rich stocks
print(mcdr  <- dim(Fm_cut)[2])

# number of indices for our DLS stocks (note that we can have more than 1 index for a stock)
print(mi <- dim(Im_cut)[2])
# matrix of number of indices for each species  (the sum of this matrix must be equal to the columns of Im and Ipresent)
nIndices <- matrix(c(1,1,2,1,2))

## matrix of ones and zeros for whether index is present for that year
Ipresent <- matrix(as.numeric(!is.na(Im_cut)), ncol = mi)
## calc sigmas for random walks data rich (sdrwdr)
sdrwdr <- matrix(apply(Fm_cut,2, FUN=function(x){sd(diff(as.numeric(x)),na.rm=T)}))

## make  corr matrix set off diag to 0 (default)
## number of rows and columns are equal to number of data poor and data rich
corrmatrix <- matrix(0.0, nrow=(mc+mcdr),ncol=(mc+mcdr))
diag(corrmatrix) <- 1

#make sure Im_cut and Cm_cut are numeric
Im_cut <- matrix(as.numeric(Im_cut), nrow=dim(Im_cut)[1])
Cm_cut <- matrix(as.numeric(Cm_cut), nrow=dim(Cm_cut)[1])
#Fm_cut <- matrix(as.numeric(Fm_cut), nrow=dim(Fm_cut)[1])

# Set boundaries, these come in sets that are equal to the number data poor stocks, the optimizer can use this info
#            log(m)      log(K)           log(q)         log(sd_i)     log(sd_c)    log(sd_rw)  logit(alpha)
lower <- c(rep(7,mc),  rep(log(1e3),mc), rep(-18,mi),  rep(-3,mi), rep(-13,mc), rep(-6,mc), rep(-20,mc) )
upper <- c(rep(20,mc), rep(log(9.0e6),mc), rep(0,mi),    rep(1,mi),  rep(1,mc),   rep(1,mc),  rep(20,mc))

############################
############################
############################
#step #1:----
###---------part for probing for best starting parameter values to find lowest NLL----
#loopresult for m and k of species (current one is for witch and lemon sole)----------
#  loopres_d_l <- data.frame("mm_d"=0, "KK_d"=0,"mm_w"=0,"KK_w"=0, "nll" = NA )
#  bestnll <- 1e20
#  bestpar <- NULL
#  
#   for(mm_d in seq(7.4,12.4,1)){ 
#    for(KK_d in seq(mm_d+1,14,1)){
#      for (mm_w in seq(7.4,12.4,1)){ 
#        for(KK_w in seq(mm_w+1,14,1)){
#         # for(mm_w in seq(7.4,12.4,1)){
#          #  for(KK_w in seq(12.4,14,1)){
#            
#            obj <- MakeADFun(
#              data = list(
#                Imat = Im_cut,
#                Cmat = Cm_cut,
#                Fmat = Fm_cut,
#                nIndices = nIndices,
#                Ipresent = Ipresent,
#                corrmatrix = corrmatrix,
#                sdrwdr = sdrwdr ),
#              parameters = list(
#                ##order is tur bll wit dab lem
#                logm   = c( 8.54,  7.95,  mm_w,  mm_d, 9.01),
#                logK   = c(10.99, 9.11, KK_w, KK_d, 11.5),
#                logq   = c(-8.55, -4.79, -11.17, -8.15, -10.38, -8,-8),       #we have as many log(q) as we have columns in Im_cut (now 7)
#                logsdi = c(-2.4, -2.2, -1.0, -1.0, -2.0, -2.0,-2.0),   #we have as many sd(I) as we have columns in Im_cut (now 7)
#                logsdc = c(-2.5, -2.5, -2.5, -2.5, -2.5),     #we have as many sd(C) as we have columns in Cm_cut (now 5)
#                logsdrw = matrix(-2.2,nrow=mc ),
#                logitalpha = c(0.07,-1,-0.6,0.5,-1.5),
#              ## CM: changed number of rows
#                logitU = cbind(cuslogit(0.1),cuslogit(0.1),cuslogit(0.1), cuslogit(0.1),cuslogit(0.1), cuslogit(Fm_cut))),
#              random = "logitU",  #this is the part that follows the random walk
#              map = list(logsdc = rep(as.factor(NA),5 ) ),
#              DLL = "sprobinhood_model",
#              silent = T)
#            
#            ## fit    
#            opt <- try( nlminb(start     = obj$par,
#                               objective = obj$fn,
#                               gradient  = obj$gr,
#                               lower     = lower,  upper     = upper,
#                               control   = list(iter.max = 1e7, eval.max = 1e7)))
#            if (class(opt) == "try-error"){
#              loopres_d_l <- rbind(loopres_d_l, data.frame( "mm_d" = mm_d, "KK_d"=KK_d,"mm_w" = mm_w,"KK_w"=KK_w, "nll" = NA ))#"mm_w" = mm_w, "KK_w"= KK_w,
#            } else { 
#              loopres_d_l <- rbind(loopres_d_l, data.frame( "mm_d" = mm_d, "KK_d"=KK_d,"mm_w" = mm_w,"KK_w"=KK_w, "nll" = opt$objective ))#"mm_w" = mm_w, "KK_w"= KK_w,
#              if (opt$objective < bestnll) bestpar <- opt$par
#            }
#            print(loopres_d_l)
#          } 
#        } 
#      }
#    }
# 
#   #}
#  #}
#  
# #step #2:----- ASSESSMENT START WITH MIX LM AND PROPORTION DISCARDS
# ###---------create the object for C++ for decided parameters--------------        
#   obj <- MakeADFun(
#     data = list(
#       Imat = Im_cut,
#       Cmat = Cm_cut,
#       Fmat = Fm_cut,
#       nIndices = nIndices,
#       Ipresent = Ipresent,
#       corrmatrix = corrmatrix,
#       sdrwdr = sdrwdr ),
#     parameters = list(
#       ##order is tur bll wit dab lem
#       logm   = c( 8.54,  7.95,  8.17,  11.1,  9.01),
#       logK   = c( 10.99, 9.11, 10.41, 14, 11.5),
#       logq   = c(-8.55, -4.79, -11.17, -8.15, -10.38, -8,-8), #we have as many log(q) as we have columns in Im_cut (now 7)
#       logsdi = c(-2.4, -2.2, -1.0, -1.0, -2.0, -2.0,-2.0), #we have as many sd(I) as we have columns in Im_cut (now 7)
#       logsdc = c(-2.5, -2.5, -2.5, -2.5, -2.5),       #we have as many sd(C) as we have columns in Cm_cut (now 5)
#       logsdrw = matrix(-2.2,nrow=mc ),
#       logitalpha = c(0.07,-1,-0.6,0.5,-1.5),
#       ## CM: changed number of rows
#       logitU = cbind(cuslogit(0.1),cuslogit(0.1),cuslogit(0.1), cuslogit(0.1),cuslogit(0.1), cuslogit(Fm_cut))),
#     random = "logitU", 
#     map = list(logsdc = rep(as.factor(NA),5 ) ),
#     DLL = "sprobinhood_model",
#     silent = F)
#  
 #step #2:----- ASSESSMENT START WITH ONLY PROPORTION DISCARDS. ALL STARTING VALUES BASED ON ASSESSMENT WITH MIX
 ###---------create the object for C++ for decided parameters--------------        
 
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
     ##order is tur bll wit dab lem
     logm   = c( 8.5424, 7.9718,  8.1748, 11.0313,  9.529),
     logK   = c( 10.982, 9.4469, 10.4173, 13.8124, 12.8577),
     logq   = c(-9.0152, -8.3071, -11.1704, -11.3254, -10.1460, -11.9822,-10.9759), #we have as many log(q) as we have columns in Im_cut (now 7)
     logsdi = c(-2.35, -2.11, -1.09, -1.31, -1.41, -1.15,-1.81), #we have as many sd(I) as we have columns in Im_cut (now 7)
     logsdc = c(-2.5, -2.5, -2.5, -2.5, -2.5),       #we have as many sd(C) as we have columns in Cm_cut (now 5)
     logsdrw = matrix(c(-2.07,-1.38,-2.06,-2.21,-1.88),nrow=mc ),
     logitalpha = c(0.077,-1.15,-0.59,-0.56,-0.2002),
     ## CM: changed number of rows
     logitU = cbind(cuslogit(0.1),cuslogit(0.1),cuslogit(0.1), cuslogit(0.1),cuslogit(0.1), cuslogit(Fm_cut))),
   random = "logitU", 
   map = list(logsdc = rep(as.factor(NA),5 ) ),
   DLL = "sprobinhood_model",
   silent = F)
  
 
############################
############################
############################
#step #3:-----  
###---------optimize parameters by fitting model to available data-----    
opt <- try( nlminb(start     = obj$par,
                     objective = obj$fn,
                     gradient  = obj$gr,
                     lower     = lower,
                     upper     = upper,
                     control   = list(iter.max = 1e7, eval.max = 1e7)))
 
 
  
#create solution of optimization---------
(rep <- sdreport(obj))
rep$value
opt$objective

latests_opt<-summary(rep)[1:34,1]  
  
#put solution in matrix and table  ----------------
paramatrix<-matrix(c(rep[[4]],
            sqrt(diag(rep[[5]])),
            signif(lower[c(1:24, 30:39)],4),
            signif(upper[c(1:24, 30:39)],4)),
            ncol =4)
  
paraframe<-as.data.frame(list( c(rep(c("tur", "bll", "wit", "dab", "lem"),2),
                                   rep(c("tur", "bll", "wit", "wit.1", "dab", "lem", "lem.1"),2),
                                   rep(c("tur", "bll", "wit", "dab", "lem"),2)),
                                 names(rep[[4]]),
                                 paramatrix))
  
colnames(paraframe)<-c("species","parameter names", "mean","sd"," lower","upper" )
  
paraframe<-paraframe[order(paraframe$species ),]
paraframe

######################################################################
### Estimates +sd for parameters:
### logm","logK","logsdrw","logitalpha"
### Final format can be copied to excel to make table
######################################################################
rep.summ <-summary(rep)
est_res  <- cbind(rep.summ[1:5,1],rep.summ[6:10,1],rep.summ[25:29,1],rep.summ[30:34,1])
sd_res   <- cbind(rep.summ[1:5,2],rep.summ[6:10,2],rep.summ[25:29,2],rep.summ[30:34,2])

colnames(est_res) <- colnames(sd_res) <- c("logm","logK","logsdrw","logitalpha")
rownames(est_res) <- rownames(sd_res) <- c("tur","bll", "wit","dab","lem")
est_res; sd_res

est_res_sig2<- apply(est_res,c(1,2),function(x){ round(x,2)})
sd_res_sig2<-apply(sd_res,c(1,2),function(x){ round(x,2)})

collap_est_sd_unconstrained_pop <- paste0(est_res_sig2," (", sd_res_sig2,")")
dim(collap_est_sd_unconstrained_pop) <- c(5,4)

colnames(collap_est_sd_unconstrained_pop) <- c("logm","logK","logsdrw","logitalpha")
rownames(collap_est_sd_unconstrained_pop) <- c("tur","bll", "wit","dab","lem")

# witch order of spec in table to match orderdab bll lem tur witch and transpose
collap_est_sd_unconstrained_pop <- t(collap_est_sd_unconstrained_pop[c("dab","bll","lem","tur","wit"),])

collap_est_sd_unconstrained_pop
#to make a nice table, use excel to separate text to columns

###########################################################################
### same (output of estimates for excel) for estimates for parameters "logq","logsdi"
###
############################################################################
est_res <- cbind(rep.summ[11:17,1],rep.summ[18:24,1])
sd_res <- cbind(rep.summ[11:17,2],rep.summ[18:24,2])

colnames(est_res) <- colnames(sd_res) <- c("logq","logsdi")
rownames(est_res) <- rownames(sd_res) <- c("tur","bll", "wit","wit.1","dab","lem","lem.1")
est_res; sd_res

est_res_sig2<- apply(est_res,c(1,2),function(x){ round(x,2)})
sd_res_sig2<-apply(sd_res,c(1,2),function(x){ round(x,2)})

collap_est_sd_unconstrained_surveys <- paste0(est_res_sig2," (",sd_res_sig2,")")
dim(collap_est_sd_unconstrained_surveys)<-c(7,2)

colnames(collap_est_sd_unconstrained_surveys) <- c("logq","logsdi")
rownames(collap_est_sd_unconstrained_surveys) <- c("tur","bll", "wit","wit.1","dab","lem","lem.1")
collap_est_sd_unconstrained_surveys <- t(collap_est_sd_unconstrained_surveys[c("dab","bll","lem","lem.1","tur","wit","wit.1"),])

collap_est_sd_unconstrained_surveys
#to make a nice table, use excel to separate text to columns

##########################################################################
#Correlation plot of the fixed effects parameters
#For now these are not used in the Manuscript, but this
#should be considered, unfortunatel because of the large 
# number of pars, the matrix is too large to be printed easily
###########################################################################
par.corr <- cov2cor(rep$cov.fixed)
corrplot.mixed(par.corr, lower = "ellipse", upper = "number") 
# colors can be changes using upper.col and lower.col

############################
############################
############################
#step #4:-----  
###---------extract the calculated estimates--------------
#B,BoBMSY,FoFMSY,Ihat,Sigma

Bhat   <- matrix( rep.summ[rownames(rep.summ) == "B", "Estimate"], nc = mc)
Blow   <- Bhat - 2 * matrix(rep.summ[rownames(rep.summ) == "B", "Std. Error"],  nc = mc)
Blow   <- apply(Blow,c(1,2),function(x)if(x<0){0}else{x})

Bhigh  <- Bhat + 2 * matrix(rep.summ[rownames(rep.summ) == "B", "Std. Error"], nc = mc)

BoBMSYhat  <- matrix( rep.summ[rownames(rep.summ) == "BoBMSY", "Estimate"], nc = mc)
BoBMSYlow  <- BoBMSYhat - 2 * matrix(rep.summ[rownames(rep.summ) == "BoBMSY", "Std. Error"],  nc = mc)
#BoBMSYlow <- apply(BoBMSYlow,c(1,2),function(x)if(x<0){0}else{x})

BoBMSYhigh <- BoBMSYhat + 2 * matrix(rep.summ[rownames(rep.summ) == "BoBMSY", "Std. Error"], nc = mc)

FoFMSYhat  <- matrix( rep.summ[rownames(rep.summ) == "FoFMSY", "Estimate"], nc = mc)
FoFMSYlow  <- FoFMSYhat - 2 * matrix(rep.summ[rownames(rep.summ) == "FoFMSY", "Std. Error"],  nc = mc)
FoFMSYlow <- apply(FoFMSYlow,c(1,2),function(x)if(x<0){0}else{x})

FoFMSYhigh <- FoFMSYhat + 2 * matrix(rep.summ[rownames(rep.summ) == "FoFMSY", "Std. Error"], nc = mc)
  
logitUhat  <- matrix( rep.summ[rownames(rep.summ) == "logitU", "Estimate"], nc = mc+mcdr)
Uhat <- invcuslogit(logitUhat)
Ulow  <- invcuslogit(logitUhat - 2 * matrix(rep.summ[rownames(rep.summ) == "logitU", "Std. Error"], nc = mc+mcdr))
Ulow <- apply(Ulow,c(1,2),function(x)if(x<0){0}else{x})
Uhigh <- invcuslogit(logitUhat + 2 * matrix(rep.summ[rownames(rep.summ) == "logitU", "Std. Error"], nc = mc+mcdr))

BMSY <-  0.5* exp( matrix( rep.summ[rownames(rep.summ) == "logK", "Estimate"], nc=mc))

## check fit to the catches
Chat <- matrix( rep.summ[rownames(rep.summ) == "Chat", "Estimate"], nc = mc)

## ses for mean Cms
Chat.se <- matrix( rep.summ[rownames(rep.summ) == "Chat", "Std. Error"], nc = mc)

## check fit to the indices
Ihat<-matrix( rep.summ[rownames(rep.summ) == "Ihat", "Estimate"], nc = mi)
## ses for mean CIs
Ihat.se <- matrix( rep.summ[rownames(rep.summ) == "Ihat", "Std. Error"], nc = mi)


###########################################################
###plot the new findings
###################################################################

#####to plot Chat results of optimization vs Cm
#-----first: data manipulation of catches to long format----
colnames(Cm) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
Cm_long<-gather(Cm,key=species, value="Catches", -years)
mis_data_Chat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA)))
colnames(mis_data_Chat) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")

years_hat<-1978:2018

Chat_incl<-cbind(years_hat,Chat)
colnames(Chat_incl) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
Chat_incl_dataframe<-as.data.frame(rbind(mis_data_Chat,Chat_incl))
C_hat_long<-gather(Chat_incl_dataframe, key=species,value="Est_Catches", -years)
Chat.se_min<- Chat - (2*Chat.se)
Chat.se_min <- apply(Chat.se_min,c(1,2),function(x)if(x<0){0}else{x})

Chat.se_min_incl<-cbind(years_hat,Chat.se_min)
colnames(Chat.se_min_incl) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
Chat.se_min_dataframe<-as.data.frame(rbind(mis_data_Chat,Chat.se_min_incl))
Chat.se_min_long<-gather(Chat.se_min_dataframe, key=species,value="SE.min", -years)

Chat.se_max<- Chat + (2*Chat.se)
Chat.se_max_incl<-cbind(years_hat,Chat.se_max)
colnames(Chat.se_max_incl) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
Chat.se_max_dataframe<-as.data.frame(rbind(mis_data_Chat,Chat.se_max_incl))
Chat.se_max_long<-gather(Chat.se_max_dataframe, key=species,value="SE.max", -years)

#-----then merge all catches into single dataframe for ggplot----
merge_catches<-as.data.frame(cbind(Cm_long,C_hat_long$Est_Catches,Chat.se_min_long$SE.min, Chat.se_max_long$SE.max ))
merge_catches[merge_catches$years==2018,4]<- NA #c(4,5,6,9,10,11)
merge_catches[merge_catches$years==2018,5]<- NA 
merge_catches[merge_catches$years==2018,6]<- NA 
colnames(merge_catches) <- c("years", "species", "Catches", "est_Catches","Min","Max")
merge_catches$species.1<-merge_catches$species
merge_catches$species.1<- factor(merge_catches$species.1, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))
merge_catches$Catches.1<-as.numeric(merge_catches$Catches)/1000
merge_catches$est_catches.1<-as.numeric(merge_catches$est_Catches)/1000
merge_catches$Min.1<-as.numeric(merge_catches$Min)/1000
merge_catches$Max.1<-as.numeric(merge_catches$Max)/1000

#-----then plot the catch data ----------

#with facets
ggplot(data = merge_catches) +
  ylab("Catch (x 1000 tonnes)")+
  xlab("Years")+
  geom_point(aes(x=years, y=as.numeric(merge_catches$Catches.1), shape="16", colour=species.1),show.legend = F,size=1.2,alpha=1)+
  facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
   theme_bw()+
  scale_colour_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  geom_line(aes(x=years,y=as.numeric(merge_catches$est_catches.1),linetype="dashed",colour=species.1),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
  geom_ribbon(aes(x=years,ymin=as.numeric(merge_catches$Min.1), ymax=as.numeric(merge_catches$Max.1),fill=species.1), linetype="22",alpha=0.4) +
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  ggtitle("Estimates of Catches")

#########################################################################
#calc residuals for Indices & Catches-----------
#########################################################################
Iresid <- log(Im_cut)-log(Ihat)##logtransformed because that is what we did the nll with
hist(Iresid,20)

(Cresid<- log(Cm_cut)-log(Chat))
hist(Cresid,20)

################################################################################
#####to plot Ihat results of optimization vs Im
###################################################################################
#-----first: data manipulation of Indices to long format----
##values are not in the same units but B is correctly calculated from Ihat trough q (catchability)
colnames(Im)  <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
Im_long       <- gather(Im,key=species, value="Indices", -years)
mis_data_Ihat <- data.frame(rbind(c(1975,rep(NA,mi)),
                                c(1976,rep(NA,mi)),
                                c(1977,rep(NA,mi))))
colnames(mis_data_Ihat) <- colnames(Im)

#Ihat_incl combines years and estimated indices, the colnames should be equal to Im colnames
Ihat_incl           <- cbind(years_hat,Ihat)
colnames(Ihat_incl) <- colnames(Im)

#Ihat_incl_dataframe includes NAs for missing years (1975-1977)
Ihat_incl_dataframe <- as.data.frame(rbind(mis_data_Ihat,Ihat_incl))
I_hat_long          <- gather(Ihat_incl_dataframe, key=species,value="Est_indices", -years)

Ihat.se_min <- Ihat - (2*Ihat.se)
Ihat.se_min <- apply(Ihat.se_min,c(1,2),function(x)if(x<0){0}else{x})
Ihat.se_min_incl<-cbind(years_hat,Ihat.se_min)
colnames(Ihat.se_min_incl) <- colnames(Im)
Ihat.se_min_dataframe<-as.data.frame(rbind(mis_data_Ihat,Ihat.se_min_incl))
Ihat.se_min_long<-gather(Ihat.se_min_dataframe, key=species,value="SE.min", -years)

Ihat.se_max<- Ihat + (2*Ihat.se)
Ihat.se_max_incl<-cbind(years_hat,Ihat.se_max)
colnames(Ihat.se_max_incl) <- colnames(Im)
Ihat.se_max_dataframe<-as.data.frame(rbind(mis_data_Ihat,Ihat.se_max_incl))
Ihat.se_max_long<-gather(Ihat.se_max_dataframe, key=species,value="SE.max", -years)

#-----then merge all indices into single dataframe for ggplot----
merge_indices <- as.data.frame(cbind(Im_long,I_hat_long$Est_indices,Ihat.se_min_long$SE.min, Ihat.se_max_long$SE.max,Ihat.se_max_long$species ))
colnames(merge_indices) <- c("years", "species", "indices", "est_indices","Min","Max","species.1")
merge_indices$species.1 <- factor(merge_indices$species.1, levels=c("Dab","Brill","Lemon sole","Lemon sole.1", "Turbot","Witch","Witch.1"))

#-----then plot the data (facet not on log scale)----------

ggplot(data = merge_indices) +
  ylab("Indices (units differ)")+
  xlab("Years")+
  geom_point(aes(x = years, y = as.numeric(merge_indices$indices),shape="16",colour=species.1),show.legend = F,size=1.2, alpha=1)+
  facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
  #scale_y_log10(limits = c(0.01,1e2))+
  #maybe change the scales manually
  theme_bw()+
  scale_colour_manual(values = wes_palette("Zissou1",7, type = "continuous"),guide=F)+
  geom_line(aes(x=years,y=as.numeric(merge_indices$est_indices),linetype="dashed", colour=species.1),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
  geom_ribbon(aes(x=years,ymin=as.numeric(merge_indices$Min), ymax=as.numeric(merge_indices$Max),fill=species.1),alpha=0.4 ,show.legend = F)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = wes_palette("Zissou1",7, type = "continuous"),guide=F)+
  #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
  #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
  ggtitle("Estimates of Indices")

###########################################################################
#####to plot Uhat results (fishing mortality)
###########################################################################
#-----first: data manipulation to long format-------    
Uhat                    <- cbind(years_hat,Uhat)
colnames(Uhat)          <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
Uhat                    <- as.data.frame(Uhat)
mis_data_Uhat           <- data.frame(rbind(c(1975, rep(NA,mc)),
                                c(1976, rep(NA,mc)),
                                c(1977, rep(NA,mc))))
mis_data_Uhat           <- data.frame(cbind(mis_data_Uhat,Fm[1:3,c("Plaice","Sole","Cod","Whiting","Haddock")]))
colnames(mis_data_Uhat) <- colnames(Uhat)
  
Uhat_incl_dataframe <- as.data.frame(rbind(mis_data_Uhat,Uhat))
Uhat_long           <- gather(Uhat_incl_dataframe, key=species,value=U, -years)

Uhat.se_min <- Ulow
Uhat.se_min_incl <- cbind(years_hat,Uhat.se_min)

colnames(Uhat.se_min_incl) <- colnames(Uhat)
mis_data_Uhat              <- data.frame(rbind(c(1975,rep(NA,mc+mcdr)),
                                               c(1976,rep(NA,mc+mcdr)),
                                               c(1977,rep(NA,mc+mcdr))))
colnames(mis_data_Uhat) <- colnames(Uhat)
Uhat.se_min_dataframe<-as.data.frame(rbind(mis_data_Uhat,Uhat.se_min_incl))
Uhat.se_min_long<-gather(Uhat.se_min_dataframe, key=species,value="SE.min", -years)

Uhat.se_max                <- Uhigh 
Uhat.se_max_incl           <- cbind(years_hat,Uhat.se_max)
colnames(Uhat.se_max_incl) <- colnames(Uhat)
Uhat.se_max_dataframe      <- as.data.frame(rbind(mis_data_Uhat,Uhat.se_max_incl))
Uhat.se_max_long           <- gather(Uhat.se_max_dataframe, key=species,value="SE.max", -years)
  
#-----then merge all indices into single dataframe for ggplot----
merge_Uhat<-as.data.frame(cbind(Uhat_long,Uhat.se_min_long$SE.min, Uhat.se_max_long$SE.max,Uhat.se_max_long$species ))
colnames(merge_Uhat) <- c("years", "species", "U","Min","Max","species.1")
merge_Uhat$species.1<- factor(merge_Uhat$species.1, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch","Cod","Plaice","Sole","Whiting","Haddock"))

merge_Uhat_data_poor <- merge_Uhat[merge_Uhat$species %in% names(Cm)[-1],]

###################################################################################
#-----then plot the data both complete and data poor (facet not on log scale)----------
#########################################################################
ggplot(data=merge_Uhat)+
  ylab("U")+
  xlab("Years")+
  geom_line(aes(x = years, y = as.numeric(merge_Uhat$U),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.6,show.legend=F)+
  facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+#ncol =1
  theme_bw()+
  scale_colour_manual(values = wes_palette("Zissou1",mc+mcdr, type = "continuous"),guide=F)+
  geom_ribbon(aes(x=years, ymin=as.numeric(merge_Uhat$Min), ymax=as.numeric(merge_Uhat$Max),fill=species.1),linetype="solid",alpha=0.2)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = wes_palette("Zissou1", mc+mcdr, type = "continuous"),guide=F)+
  #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
  ggtitle("Estimates of Fishing mortality")

ggplot(data=merge_Uhat_data_poor)+
  ylab("U")+
  xlab("Years")+
  geom_line(aes(x = years, y = as.numeric(merge_Uhat_data_poor$U),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
  facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+#,ncol =1
  theme_bw()+
  scale_colour_manual(values = c(wes_palette("Zissou1",mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  geom_ribbon(aes(x=years, ymin=as.numeric(merge_Uhat_data_poor$Min), ymax=as.numeric(merge_Uhat_data_poor$Max),fill=species.1),linetype="solid",alpha=0.4)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = c(wes_palette("Zissou1",mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
  ggtitle("Estimates of Fishing mortality")

#############################################################################################
#####to plot Bhat results
##############################################################################################
#-----first: data manipulation to long format----

Bhat <- cbind(years_hat,Bhat)
colnames(Bhat) <- names(Cm)
Bhat <- as.data.frame(Bhat)
Bhat_long <- gather(Bhat,key=species, value=B, -years)

Bhat.se_min <- Blow
Bhat.se_min_incl <- cbind(years_hat,Bhat.se_min)
colnames(Bhat.se_min_incl) <- names(Cm)
Bhat.se_min_incl <- as.data.frame(Bhat.se_min_incl)
Bhat.se_min_long <- gather(Bhat.se_min_incl, key=species,value="SE.min", -years)

Bhat.se_max <- Bhigh 
Bhat.se_max_incl <- cbind(years_hat,Bhat.se_max)
colnames(Bhat.se_max_incl) <- names(Cm)
Bhat.se_max_incl <- as.data.frame(Bhat.se_max_incl)
Bhat.se_max_long <- gather(Bhat.se_max_incl, key=species,value="SE.max", -years)


#-----then merge all Biomass estimates into single dataframe for ggplot----
merge_Bhat<-as.data.frame(cbind(Bhat_long,Bhat.se_min_long$SE.min, Bhat.se_max_long$SE.max,Bhat.se_max_long$species ))
colnames(merge_Bhat) <- c("years", "species", "B","Min","Max","species.1")
merge_Bhat$species.1<- factor(merge_Bhat$species.1, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))
merge_Bhat$B.1<-as.numeric(merge_Bhat$B)/1000
merge_Bhat$Min.1<-as.numeric(merge_Bhat$Min)/1000
merge_Bhat$Max.1<-as.numeric(merge_Bhat$Max)/1000

#-----then plot the data------
ggplot(data=merge_Bhat)+
  ylab("Biomass (x 1000 tonnes)")+
  xlab("Years")+
  geom_line(aes(x = years, y = as.numeric(merge_Bhat$B.1),linetype="dashed",colour=species.1),linetype="22",size=0.75, alpha=0.8,show.legend=F)+
  facet_wrap(facets = vars(species.1),scales="free_y", ncol=1)+#ncol=1
  theme_bw()+
  scale_colour_manual(values = c(wes_palette("Zissou1",mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  geom_ribbon(aes(x=years, ymin=as.numeric(merge_Bhat$Min.1), ymax=as.numeric(merge_Bhat$Max.1),fill=species.1),linetype="22",alpha=0.4)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = c(wes_palette("Zissou1",mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
  ggtitle("Estimates of Biomasses")

##############################################################################
#####to plot BoBMSYhat results-----
############################################################################
#-----first: data manipulation to long format----
  
BoBMSYhat<-cbind(years_hat,BoBMSYhat)
colnames(BoBMSYhat) <- names(Cm)
BoBMSYhat<-as.data.frame(BoBMSYhat)
BoBMSYhat_long<-gather(BoBMSYhat,key=species, value=BoBMSY, -years)

BoBMSYhat.se_min<-BoBMSYlow
BoBMSYhat.se_min_incl<-cbind(years_hat,BoBMSYhat.se_min)
colnames(BoBMSYhat.se_min_incl) <- names(Cm)
BoBMSYhat.se_min_incl<-as.data.frame(BoBMSYhat.se_min_incl)
BoBMSYhat.se_min_long<-gather(BoBMSYhat.se_min_incl, key=species,value="SE.min", -years)

BoBMSYhat.se_max<-BoBMSYhigh 
BoBMSYhat.se_max_incl<-cbind(years_hat,BoBMSYhat.se_max)
colnames(BoBMSYhat.se_max_incl) <- names(Cm)
BoBMSYhat.se_max_incl<-as.data.frame(BoBMSYhat.se_max_incl)
BoBMSYhat.se_max_long<-gather(BoBMSYhat.se_max_incl, key=species,value="SE.max", -years)

#------then merge all BoBMSY estimates into single dataframe for ggplot----
merge_BoBMSYhat <- as.data.frame(cbind(BoBMSYhat_long,BoBMSYhat.se_min_long$SE.min, BoBMSYhat.se_max_long$SE.max,BoBMSYhat.se_max_long$species ))
colnames(merge_BoBMSYhat) <- c("years", "species", "BoBMSY","Min","Max","species.1")
merge_BoBMSYhat$species.1 <- factor(merge_BoBMSYhat$species.1, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))

#-----then plot of BoBMSYhat, plot to object because BoBMSY and FoFMSY will be combined later-------
BoBMSY_plot_polite <-  ggplot(data=merge_BoBMSYhat)+
  ylab("B/BMSY")+
  xlab("Years")+
  geom_line(aes(x = years, y = as.numeric(merge_BoBMSYhat$B),linetype="dashed",colour=species.1),linetype="22",size=0.75, alpha=1,show.legend=F)+
  facet_wrap(facets = vars(species.1),nrow=1)+#,scales="free_y"
  theme_bw()+
  geom_hline(aes(yintercept=1,colour=species.1),linetype="twodash", size=0.8, alpha=0.4 )+
  scale_colour_manual(values = c(wes_palette("Zissou1", mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  geom_ribbon(aes(x=years, ymin=as.numeric(merge_BoBMSYhat$Min), ymax=as.numeric(merge_BoBMSYhat$Max),fill=species.1),linetype="22",alpha=0.4)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = c(wes_palette("Zissou1", mc + 1, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
  #ggtitle("Estimates of B/BMSY")+
  theme(plot.margin = margin(10,15 ,0, 10))+
  theme(panel.spacing = unit(1, "lines"))+
  theme(strip.background = element_blank(),strip.text.x = element_blank() )+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)))
 
##################################################################
#####to plot FoFMSYhat results
################################################################
#-----first: data manipulation to long format----

FoFMSYhat  <- matrix( rep.summ[rownames(rep.summ) == "FoFMSY", "Estimate"], nc = mc)
FoFMSYlow  <- FoFMSYhat - 2 * matrix(rep.summ[rownames(rep.summ) == "FoFMSY", "Std. Error"],  nc = mc)
#FoFMSYlow <- apply(FoFMSYlow,c(1,2),function(x)if(x<0){0}else{x})

FoFMSYhat<-cbind(years_hat,FoFMSYhat)
colnames(FoFMSYhat) <- names(Cm)
FoFMSYhat<-as.data.frame(FoFMSYhat)
FoFMSYhat_long<-gather(FoFMSYhat,key=species, value=FoFMSY, -years)

FoFMSYhat.se_min<-FoFMSYlow
FoFMSYhat.se_min_incl<-cbind(years_hat,FoFMSYhat.se_min)
colnames(FoFMSYhat.se_min_incl) <- names(Cm)
FoFMSYhat.se_min_incl<-as.data.frame(FoFMSYhat.se_min_incl)
FoFMSYhat.se_min_long<-gather(FoFMSYhat.se_min_incl, key=species,value="SE.min", -years)

FoFMSYhat.se_max<-FoFMSYhigh 
FoFMSYhat.se_max_incl<-cbind(years_hat,FoFMSYhat.se_max)
colnames(FoFMSYhat.se_max_incl) <- names(Cm)
FoFMSYhat.se_max_incl<-as.data.frame(FoFMSYhat.se_max_incl)
FoFMSYhat.se_max_long<-gather(FoFMSYhat.se_max_incl, key=species,value="SE.max", -years)

#------merge all FoFMSY estimates into single dataframe for ggplot----
merge_FoFMSYhat<-as.data.frame(cbind(FoFMSYhat_long,FoFMSYhat.se_min_long$SE.min, FoFMSYhat.se_max_long$SE.max,FoFMSYhat.se_max_long$species ))
colnames(merge_FoFMSYhat) <- c("years", "species", "FoFMSY","Min","Max","species.1")
merge_FoFMSYhat$species.1<- factor(merge_FoFMSYhat$species.1, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))

#-----then plot of FoFMSYhat, plot to object because combined with BoBMSY plot later ----
FoFMSY_plot_polite <- ggplot(data=merge_FoFMSYhat)+
  ylab("F/FMSY")+
  #xlab("Years")+
  geom_line(aes(x = years, y = as.numeric(merge_FoFMSYhat$FoFMSY),linetype="dashed",colour=species.1),linetype="22",size=0.75, alpha=1,show.legend=F)+
  facet_wrap(facets = vars(species.1),nrow=1)+#,scales="free_y"
  theme_bw()+
  geom_hline(aes(yintercept=1,colour=species.1),linetype="twodash", size=0.8, alpha=0.4 )+
  scale_colour_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  geom_ribbon(aes(x=years, ymin=as.numeric(merge_FoFMSYhat$Min), ymax=as.numeric(merge_FoFMSYhat$Max),fill=species.1),linetype="22",alpha=0.4)+
  scale_alpha_manual(guide=F)+
  scale_fill_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6,7)]),guide=F)+
  #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
  #ggtitle("Estimates of F/FMSY")+
  theme(plot.margin = margin(10,15 ,0, 10))+
  theme(panel.spacing = unit(1, "lines"))+
  theme(axis.text.x =element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank())
  
############################
############################
############################
#step #6:-----
#####
###create first visual representation of the data (in new script to be optimized): make combined plots of the estimates 
#####
#alles bij elkaar per rij

ggdraw() +
  draw_plot(FoFMSY_plot_polite, x = 0, y = 0.5, width = 1, height = .5) +
  draw_plot(BoBMSY_plot_polite, x = 0, y = 0, width = 1, height = .5)

##############################################################################
# Do we need to save objects to file for comparison with constrained version?
################################################################################
save(latests_opt, Im_cut, Cm_cut, Fm_cut, nIndices, Ipresent, sdrwdr,
     rep.summ, years, startcut, mc, mcdr, mi, lower, upper,  
     years_hat,merge_catches, merge_indices, merge_BoBMSYhat, merge_FoFMSYhat,
     collap_est_sd_unconstrained_pop,collap_est_sd_unconstrained_surveys,
     file=file.path(rootPath, "data/intermediate results/unconstrained results.Rdata"))

#################################
# Can we get ADNUTS going?
##################################

# TMB::runExample("simple")
# #init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)
# init <- function() list((rep.summ[1:444,1]) )
# fit <- sample_tmb(obj=obj, init=init, iter=100, chains=2)
# fit <- sample_tmb(obj=obj, init=init, lower=lower, upper=upper)
# fit <- sample_tmb(obj=obj, init=NULL)
# fit <- sample_tmb(obj=obj, init=init, iter=3000, chains = 1, control = list(adapt_delta = 0.99))
#
# fit <- sample_tmb(obj=obj, init=init, iter=3000, control = list(adapt_delta = 0.99))

# post <- extract_samples(fit)
# sp <- extract_sampler_params(fit)
# 
# pairs_admb(fit)
# launch_shinytmb(fit)
