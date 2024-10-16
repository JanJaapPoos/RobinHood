library(lattice); library(rasterVis); library(latticeExtra);library(corrplot)
library(gridExtra); library(dplyr); library(tidyr);library("RColorBrewer"); library(ggplot2)

if(length(grep("coilin", getwd())) == 1){
 rootPath <- ".."
}else{
    rootPath <- "C://Users/poos001/OneDrive - Wageningen University & Research/projects/RobinHood"
}

#Read data (effort and landings)
landings <- read.csv(file.path(rootPath,"Data/STECF data/map__landings_by_rectangle_data.csv"), stringsAsFactors = F)
effort <- read.csv(file.path(rootPath,"Data/STECF data/effort_(hours_fished)_Full_Data_data.csv"), stringsAsFactors = F)

##########################################
# Which fleets-years effort are constant?
############################################

#make all gears with small contributions into "OTHER"
effort[effort$regulated.gear %in% c("PEL_SEINE","PEL_TRAWL", "NONE","DREDGE","POTS", "TR3","OTTER","LL1","DEM_SEINE","BEAM"),]$regulated.gear <- "OTHER"

#make aggregations (summing) for effort in area 4 and for years 
effortbygearyear <- aggregate(Effective.Effort~regulated.gear  + year, data=subset(effort, area=="4" & year> 2003 & year <2016), FUN="sum")

#for one reason or the other, the amount of effort by pots is enormous (and cannot be correct), so let's remove these
effortbygearyear_noother <- subset(effortbygearyear, !regulated.gear == "OTHER")

ggplot(data =effortbygearyear_noother) +
  ylab("Fishing effort(x hours fished)")+
  xlab("Years")+
  geom_point(aes(x=year, y=Effective.Effort, colour=regulated.gear),size=1.2,alpha=1)+
  theme_bw()+
  geom_line(aes(x=year,y=Effective.Effort,linetype="dashed",colour=regulated.gear),linetype="22",alpha=0.8)+
  scale_alpha_manual(guide=F)+
  ggtitle("Estimates of Fishing effort")+
  geom_text(data = effortbygearyear_noother %>% filter(year == last(year)), aes(label = regulated.gear, 
                                                             x = year + 0.5, 
                                                             y = Effective.Effort, 
                                                             color = regulated.gear))  

#############################################
# Make association matrix from landings
############################################

#select only relevant species from STECF landings data 
relevantspecies <-  c("PLE","SOL", "WHG","COD","HAD","TUR","BLL", "LEM", "WIT", "DAB")
landings <- landings[landings$species %in% relevantspecies,]

#make all gears with small contributions into "OTHER"
landings[landings$regulated.gear %in% c("PEL_SEINE","PEL_TRAWL", "NONE","DREDGE","POTS", "TR3","OTTER","LL1","DEM_SEINE","BEAM"),]$regulated.gear <- "OTHER"

#make aggregations (summing)
lanbygearrectyear <- aggregate(Landings~regulated.gear + rectangle + species + year, data=landings, FUN="sum")
lanbygearveslenyear <- aggregate(Landings~regulated.gear + vessel.length + species + year, data=landings, FUN="sum")
lanbygearyear <- aggregate(Landings~regulated.gear + species +year, data=landings, FUN="sum")
lanbygearveslen <- aggregate(Landings~regulated.gear + vessel.length + species , data=landings, FUN="sum")
lanbygear     <- aggregate(Landings~regulated.gear + species , data=landings, FUN="sum")
lanyear       <- aggregate(Landings~ species + year, data=landings, FUN="sum")
lan           <- aggregate(Landings~ species, data=landings, FUN="sum")

final <- merge(lan, lanbygear, by="species")
final$frac <- final$Landings.y/final$Landings.x

finalyear <- merge(lanyear, lanbygearyear, by=c("species","year"))
finalyear$frac <- finalyear$Landings.y/finalyear$Landings.x

finalyear_incl_veslen <- merge(lanyear, lanbygearveslenyear, by=c("species","year"))
finalyear_incl_veslen$frac <- finalyear_incl_veslen$Landings.y/finalyear_incl_veslen$Landings.x
finalyear_incl_veslen$gear_vl <- paste0(finalyear_incl_veslen$regulated.gear,finalyear_incl_veslen$vessel.length)

finalyear_incl_rect <- merge(lanyear, lanbygearrectyear, by=c("species","year"))
finalyear_incl_rect$frac <- finalyear_incl_rect$Landings.y/finalyear_incl_rect$Landings.x
finalyear_incl_rect$gear_rect <- paste0(finalyear_incl_rect$regulated.gear,finalyear_incl_rect$rectangle)

cortest                   <- tidyr::spread(final[,c("regulated.gear","species","frac")],key= species, value=frac)
cortestyear               <- tidyr::spread(finalyear[,c("regulated.gear","species","year","frac")],key= species, value=frac)
cortestyear_raw           <- tidyr::spread(lanbygearyear,key= species, value=Landings) 

cortestyear_incl_veslen   <- tidyr::spread(finalyear_incl_veslen[,c("gear_vl","species","year","frac")],key= species, value=frac)
cortestyear_incl_rect     <- tidyr::spread(finalyear_incl_rect[,c("gear_rect","species","year","frac")],key= species, value=frac)
cortestyear_incl_rect_raw <- tidyr::spread(lanbygearrectyear,key= species, value=Landings) 

#MAKE COLORS (NOTE THAT THE COL DEFINITION IS different FROM (3)PREL PLOTS)
#NOTE THAT IF WE PUT THIS FIG IN THE MS, THEN CAPTION SHOULD EXPLAIN THAT ORDER IS BASED ON CLUSTERING (SEE ALSO ?CORRPLOT())

display.brewer.all()
col<- colorRampPalette(brewer.pal(n=9,name="RdYlBu"))
col(200)[1:100]
adj_col2 <- c(col(200)[1:100],col(100)[50:100],rep("#FFFFFF",50))
adj_col3<-c(col(200)[1:100]) #,col(100)[50:100])

######################################################################
#FINAL FIGURES FOR PAPER
##########################################################################

###convert data into dataframe for corrplot
######here comes the second graph ##############-------

dat_frame_final<-final[,c(1,3,5)]%>%spread(regulated.gear,frac)
dat_frame_final<-cbind(dat_frame_final[,1], round(dat_frame_final[,-1],2))
names_species<-unique(as.matrix(t(dat_frame_final))[1,])

new_matrix_final<-as.matrix(t(dat_frame_final))[-1,]
colnames(new_matrix_final)<-names_species
names_methods<-rownames(new_matrix_final)

new_matrix_final<-matrix(data = sapply(new_matrix_final, as.numeric), ncol = 10, nrow = 7)
colnames(new_matrix_final)<-names_species
rownames(new_matrix_final)<-names_methods

par(mfrow=c(1,1))
####bij 800*550 -------------
p1_j<-corrplot(new_matrix_final,col=rev(adj_col2),
               diag=T,
               method="color", type = "full",number.cex = 0.8,
               cl.lim= c(0,1),is.corr = F, addCoef.col = "black",  # Add coefficient of correlation
               tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.1) # Text label color and rotation

mtext(text = "Geartype", side = 2, line = 1,  cex=1.3, at = 4 )
mtext(text = "Species", side = 1, line = 3,  cex=1.3, at = 5.5 )

######################################################################
#Select dataset to use for correlation
##########################################################################
#Select dataset to use
rawcatch <- cortestyear_raw #alternative cortestyear_incl_rect_raw

#####################################################################################
#Make corr matrix (of either raw catches or corrected catches) 
#####################################################################################
par(mfrow=c(2,2))
for (source in c("ref_corrected","raw")){
  for (corrtype in c("pearson","spearman")){
    if (source == "ref_corrected"){
      
      #if cortestyear_raw then option to divide by fleet/years with constant effort (defined by refgear): ugly loop
      #that is what is done below: catches are divided by the cacthes of the ref gear. That one ends up as having only 1s
      refgear <- "GN1"
      for(ii in unique(rawcatch$regulated.gear)){
        tmp <- rawcatch[rawcatch$regulated.gear==ii,names(rawcatch) %in% relevantspecies]/rawcatch[rawcatch$regulated.gear==refgear,names(rawcatch) %in% relevantspecies]
        tmp <- cbind(rawcatch[rawcatch$regulated.gear==ii,!names(rawcatch) %in% relevantspecies],tmp)
        if (ii ==unique(rawcatch$regulated.gear)[1]){
          corrcatch <- tmp
        }else{
          corrcatch <- rbind(corrcatch,tmp)
        }
      }
      
      #sum multiplication factors of catches by fleet so that we have a sum by year
      corrcatch2 <- aggregate(cbind(BLL,COD,DAB,HAD,LEM,PLE,SOL,TUR,WHG,WIT) ~ year,data=corrcatch, FUN="sum")
      
      #calc differences, subset years up to 2011 where effort was constant for GN1
      corrcatch3 <- apply(subset(corrcatch2, year <= 2011),2,  FUN=diff)
      ## correlation is from here
      #plot(as.data.frame(corrcatch3[,-1]))
      ## using Spearman more robust for time series
      corrmatrix <- cor(corrcatch3[,-1],use="complete.obs", method = corrtype)
      
    } else {
      
      #use either usedat(for raw catches, or effort corrected catches (tmpc) (in this case for ))
      usedat <-   rawcatch[,names(rawcatch) %in% relevantspecies] # alternative corrcatch[,names(corrcatch) %in% relevantspecies]#
      #usedat <- corrcatch[,names(corrcatch) %in% paste0(relevantspecies,"dif")]# alternative rawcatch[,names(rawcatch) %in% relevantspecies]
      corrmatrix <- cor(usedat,use="complete.obs", method = corrtype)
    
      }
    p2_j<-corrplot(corrmatrix,  col=rev(adj_col2[1:150]),
               diag=FALSE,
               method="color", type = "upper", order = "hclust",number.cex = 0.8,
               cl.lim= c(min(corrmatrix),1), is.corr = FALSE, addCoef.col = "black",  # Add coefficient of correlation
               tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.1, tl.offset = 0.4, # Text label color and rotation
               main=paste(source, corrtype))
  
    mtext(text = "Species", side = 2, line = -4,  cex=1.3, at = 5.75 )
    mtext(text = "Species", side = 1, cex=1.3,line = 3.5, at = 6 )
  }
}

## take a look at a realisation of such a correlated process
library(mvtnorm)
## standard deviations
nspp <- nrow(corrmatrix)
sd_diag <- diag(rep(0.05, nspp))
S <- sd_diag %*% corrmatrix %*% sd_diag
H <- 0.5 + apply(rmvnorm(20, rep(0, nspp), S), 2, cumsum)
colnames(H) <- colnames(corrmatrix)
matplot(H, type = "l")

## correlation of the harvest rates - not sure what's up with colours here
corrplot(cor(H),  col  =rev(adj_col2),
         diag=FALSE,
         method="color", type = "upper", order = "hclust",number.cex = 0.8,
         cl.lim= c(-1,1),is.corr = FALSE, addCoef.col = "black",  # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.1, tl.offset = 0.4) # Text label color and rotation


#####################################################################################
# Save corrmatrix for use in Robin Hood code (6)  
###################################################################################
save(corrmatrix, file=file.path(rootPath,"Data/intermediate results/corrmatrix.Rdata"))

##############################################################################
#compare with web-results of f trends------
# THiS is now not in paper yet, but does correlations in trends DR Fishing mort 
###################################################################################

corr_web<-cor(Fm[,2:6],use="complete.obs")
colnames(corr_web)<-c("PLE","SOL","COD","WHG","HAD")

#Note that the figure below takes the corrmatrix, not the observed corrs between DR Fs
corrplot(corrmatrix[6:10,6:10],  col  =rev(adj_col2),
         diag=FALSE,
         method="color", type = "upper", order = "hclust", number.cex = .9,
         cl.lim= c(-0.5,1),is.corr = FALSE, addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.31, tl.offset = 0.4) # Text label color and rotation
