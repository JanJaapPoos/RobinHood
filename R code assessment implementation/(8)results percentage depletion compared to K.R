#first run
rep.summ


biomass_t0<-Bhat[1,]
biomass_inf<-as.vector(exp(as.vector(rep.summ[rep.summ%>%rownames(rep.summ)=="logK",][,1])))
names(biomass_inf)<-c("Turbot","Brill","Witch","Dab","Lemon sole")


biomass_t0
biomass_inf

percent_K<-(biomass_t0/biomass_inf)*100
percent_K


#robin hood


#------------------------------------------------


#order = turbllwitdablem
rep 
#zie sdrw van witch bijvoorbeeld
#dan:
plot(Uhat[colnames(Uhat)=="Witch"][,])
log(sd(Uhat[colnames(Uhat)=="Witch"][,]))

log(mean(diff(Uhat[colnames(Uhat)=="Witch"][,])))
log(sd(diff(Uhat[colnames(Uhat)=="Witch"][,])))

# sdrwdr is: matrix(apply(Fm_cut,2, FUN=function(x){sd(diff(as.numeric(x)),na.rm=T)}))

log(sdrwdr)

