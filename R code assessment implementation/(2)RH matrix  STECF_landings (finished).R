library(lattice); library(rasterVis); library(latticeExtra);library(corrplot)
library(gridExtra); library(dplyr); library(tidyr);library("RColorBrewer")

rootPath <- "m:/robinhood/"

landings <- read.csv(file.path(rootPath,"data/stecf data/map__landings_by_rectangle_data.csv"), stringsAsFactors = F)

#select only relevant species from STECF landings data 
landings <- landings[landings$species %in% c("PLE","SOL", "WHG","COD","HAD","TUR","BLL", "LEM", "WIT", "DAB"),]

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

######################################################################
#MAKE cormatrix in species order for RH
##########################################################################
tempdat <- cortestyear_incl_rect_raw[,4:13]
tempdat <- tempdat[,c("TUR", "BLL","WIT","DAB","LEM","PLE","SOL","COD","WHG","HAD")]
corrmatrix <- cor(tempdat,use="complete.obs")
#rownames(corrmatrix)<-c("Turbot", "Brill", "Witch", "Dab", "Lemon sole", "Plaice", "Sole", "Cod", "Whiting", "Haddock")


#par(cex=1.5)
p2_j<-corrplot(corrmatrix,  col  =rev(adj_col2),
         diag=FALSE,
         method="color", type = "upper", order = "hclust",number.cex = 0.8,
         cl.lim= c(-0.5,1),is.corr = FALSE, addCoef.col = "black",  # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.1, tl.offset = 0.4) # Text label color and rotation


###convert data into dataframe for corrplot######here comes the second graph ##############-------

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

p2_j<-corrplot(corrmatrix,  col  =rev(adj_col2),
               diag=FALSE,
               method="color", type = "upper", order = "hclust",number.cex = 0.8,
               cl.lim= c(-0.5,1),is.corr = FALSE, addCoef.col = "black",  # Add coefficient of correlation
               tl.col = "black", tl.srt = 90, na.label = "x",tl.cex = 1.1, tl.offset = 0.4) # Text label color and rotation

mtext(text = "Species", side = 2, line = -4,  cex=1.3, at = 5.75 )
mtext(text = "Species", side = 1, cex=1.3,line = 3.5, at = 6 )

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





