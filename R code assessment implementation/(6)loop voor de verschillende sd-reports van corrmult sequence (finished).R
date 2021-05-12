library(TMB); library(corrplot);library(knitr); library(ggplot2);library(wesanderson); library(cowplot)
library(adnuts) 

# set the working directory & compile C++ model----
#rootPath <- "m:/robinhood/"
rootPath <- "C:/Users/jappi/Bureaublad/robinhood"

setwd(file.path(rootPath,'code/Robinhood')); 

#get dll to do TMB
#compile("sprobinhood_model.cpp")
dyn.load(dynlib("sprobinhood_model"))

# load correlation matrix and intermediate results from robinhood assessment
load(file.path(rootPath,"Data/intermediate results/corrmatrix.Rdata"))
load(file.path(rootPath,"Data/intermediate results/unconstrained results.Rdata"))

#masterfunctions to keep results within biological limits--------------
source(file=file.path(rootPath,"code/Jasper final/assessment implementation/(0)helper functions.R"))

corrseq <- seq(1,1,0.05)

#define names for parameters to be optimized
names_opt_par<-c(rep("logm", mc), rep("logK",mc),
                  rep("logq",mi),
                  rep("logsdi",mi),rep("logsdc",mc),    
                  rep("logsdrw",mc), rep("logitalpha", mc))

parmatrix <- matrix(NA,  nrow=length(corrseq), ncol=39,dimnames=list("corr" = corrseq,"par"=names_opt_par))


#--------effect of different corrmults----

for(corrmult in seq(1,1,0.05)){
  
  ## make  corr matrix set according to Figure in into (assumed that first Data Poor, then Data Rich)
  ## order was TUR, Brill, Lem, ple, sol, cod 
  ## new order 1:tur, 2:brill, 3:witch, 4:dab, 5:lem, 6:plaice, 7:sole, 8:cod, 9:whiting 10:haddock
  
  
  corrmat <- corrmatrix * corrmult 
  
  diag(corrmat) <- 1
  
  obj_rh <- MakeADFun(
    data = list(
      Imat = Im_cut,
      Cmat = Cm_cut,
      Fmat = Fm_cut,
      nIndices = nIndices,
      Ipresent = Ipresent,
      corrmatrix = corrmat,
      sdrwdr = sdrwdr ),
    parameters = list(
      logm   = latests_opt[names(latests_opt) %in% "logm"],
      logK   = latests_opt[names(latests_opt) %in% "logK"],
      logq   = latests_opt[names(latests_opt) %in% "logq"],       #we have as many log(q) as we have columns in Im_cut (now 6)
      logsdi = latests_opt[names(latests_opt) %in% "logsdi"],   #we have as many sd(I) as we have columns in Im_cut (now 6)
      logsdc = c(-2.5, -2.5, -2.5, -2.5, -2.5),     #we have as many sd(C) as we have columns in Cm_cut (now 5)
      logsdrw = matrix(-2.85,nrow=mc ),
      logitalpha = latests_opt[names(latests_opt) %in% "logitalpha"],
      ## CM: changed number of rows
      ##logitU = cbind(cuslogit(0.2),cuslogit(0.2),cuslogit(0.2), cuslogit(0.2),cuslogit(0.2), cuslogit(Fm_cut))),
      logitU = cbind(rep.summ[rownames(rep.summ) == "logitU", "Estimate"][1:(length(years)-startcut)],
                     rep.summ[rownames(rep.summ) == "logitU", "Estimate"][((1*(length(years)-startcut))+1):(2*(length(years)-startcut))],
                     rep.summ[rownames(rep.summ) == "logitU", "Estimate"][((2*(length(years)-startcut))+1):(3*(length(years)-startcut))],
                     rep.summ[rownames(rep.summ) == "logitU", "Estimate"][((3*(length(years)-startcut))+1):(4*(length(years)-startcut))],
                     rep.summ[rownames(rep.summ) == "logitU", "Estimate"][((4*(length(years)-startcut))+1):(5*(length(years)-startcut))],
                     cuslogit(Fm_cut))),
    random = "logitU", 
    #map = list(logsdc = rep(as.factor(NA),5 ) ),
    DLL = "sprobinhood_model",
    silent = F)
  
  ## fit    
  opt_rh <- try(nlminb(start  = obj_rh$par,
                       objective = obj_rh$fn,
                       gradient  = obj_rh$gr,
                       lower     = lower,
                       upper     = upper,
                       control   = list(iter.max = 1e7, eval.max = 1e7)))
  
  parmatrix[as.character(corrmult),] <- opt_rh$par
  sdrep_obj_rh <- sdreport(obj_rh)
  
    if (class(sdrep_obj_rh) == "try-error"){
    print("i is",corrmult,"gives try-error")
  } else { 
    assign(paste0("rep_rh", corrmult),sdrep_obj_rh)
  }
}



######################################################################
### Estimates +sd for parameters:
### logm","logK","logsdrw","logitalpha"
### Final format can be copied to excel to make table
######################################################################
rep.summ <-summary(rep_rh1)
est_res  <- cbind(rep.summ[1:5,1],rep.summ[6:10,1],rep.summ[25:29,1],rep.summ[30:34,1],rep.summ[35:39,1])
sd_res   <- cbind(rep.summ[1:5,2],rep.summ[6:10,2],rep.summ[25:29,2],rep.summ[30:34,2],rep.summ[35:39,2])

colnames(est_res) <- colnames(sd_res) <- c("logm","logK","logsc","logsdrw","logitalpha")
rownames(est_res) <- rownames(sd_res) <- c("tur","bll", "wit","dab","lem")
est_res; sd_res

est_res_sig2<- apply(est_res,c(1,2),function(x){ round(x,2)})
sd_res_sig2<-apply(sd_res,c(1,2),function(x){ round(x,2)})

collap_est_sd_robinhood_pop <- paste(est_res_sig2," (", sd_res_sig2,")")
dim(collap_est_sd_robinhood_pop) <- c(5,5)

colnames(collap_est_sd_robinhood_pop) <- c("logm","logK","logsc","logsdrw","logitalpha")
rownames(collap_est_sd_robinhood_pop) <- c("tur","bll", "wit","dab","lem")
collap_est_sd_robinhood_pop <- t(collap_est_sd_robinhood_pop[c("dab","bll","lem","tur","wit"),])

collap_est_sd_robinhood_pop

#t(collap_est_sd_robinhood_pop)
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

collap_est_sd_robinhood_surveys <- paste0(est_res_sig2," (",sd_res_sig2,")")
dim(collap_est_sd_robinhood_surveys)<-c(7,2)

colnames(collap_est_sd_robinhood_surveys) <- c("logq","logsdi")
rownames(collap_est_sd_robinhood_surveys) <- c("tur","bll", "wit","wit.1","dab","lem","lem.1")
collap_est_sd_robinhood_surveys

collap_est_sd_robinhood_surveys <- t(collap_est_sd_robinhood_surveys[c("dab","bll","lem","lem.1","tur","wit","wit.1"),])

collap_est_sd_robinhood_surveys
#to make a nice table, use excel to separate text to columns


save(rep_rh1,collap_est_sd_robinhood_pop, collap_est_sd_robinhood_surveys,file=file.path(rootPath,"Data/intermediate results/Robin Hood results.Rdata"))

#-------------------------------------------------
