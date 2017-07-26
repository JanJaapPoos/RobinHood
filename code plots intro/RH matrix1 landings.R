library(lattice); library(rasterVis); library(latticeExtra); library(gridExtra)

# this data is downloaded from STECF website

landings <- read.csv("w://imares/ijmuiden/afdeling/projecten/data poor mixed fisheries/robin hood/data/map__landings_by_rectangle_data.csv", stringsAsFactors = F)

#select only relevant species
landings <- landings[landings$species %in% c("PLE","SOL", "WHG","COD","HAD","TUR","BLL", "LEM"),]

#make all gears with small contributions into "OTHER"
landings[landings$regulated.gear %in% c("PEL_SEINE","PEL_TRAWL", "NONE","DREDGE","POTS", "TR3","OTTER","LL1","DEM_SEINE","BEAM"),]$regulated.gear <- "OTHER"

#make aggregations (summing)
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

my.theme <- BuRdTheme()

# Find the min and max values
my.min <- min(final$Landings.y)
my.max <- max(final$Landings.y)

# Customize the colorkey
my.at <- seq(my.min, my.max, length.out=length(my.theme$regions$col)-1)
my.ckey <- list(at=my.at, col=my.theme$regions$col)

levelplot(Landings.y ~  species + regulated.gear , data=final,at=my.at,par.settings=my.theme, colorkey=my.ckey)


# Find the min and max values
my.min <- min(final$frac)
my.max <- max(final$frac)

# Customize the colorkey
my.at <- seq(my.min, my.max, length.out=length(my.theme$regions$col)-1)
my.ckey <- list(at=my.at, col=my.theme$regions$col)

#customize panel
my.panel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,2), cex=0.8)
}

#levelplot(frac ~  species + regulated.gear , data=final,panel=my.panel,at=my.at,par.settings=my.theme, colorkey=my.ckey)
#levelplot(frac ~  species + regulated.gear , data=final,panel=my.panel,at=my.at,par.settings=my.theme, colorkey=my.ckey)



cortest                  <- tidyr::spread(final[,c("regulated.gear","species","frac")],key= species, value=frac)
cortestyear              <- tidyr::spread(finalyear[,c("regulated.gear","species","year","frac")],key= species, value=frac)
cortestyear_incl_veslen  <- tidyr::spread(finalyear_incl_veslen[,c("gear_vl","species","year","frac")],key= species, value=frac)



corall <- corall_incl_veslen <-  array(NA, dim = c(12,8,8))
ii <- 1
for (yrs in 2004:2015){
 corall[ii,,] <-  cor(cortestyear[cortestyear$year== yrs,3:10])
 corall_incl_veslen[ii,,] <-  cor(cortestyear_incl_veslen[cortestyear_incl_veslen$year== yrs,3:10], use="pairwise.complete.obs")
 ii <- ii+1
}

p1 <- levelplot(frac ~  as.factor(species) + as.factor(regulated.gear), data=final,panel=my.panel, aspect= "fill", col.regions = colorRampPalette(c('white','dark red'), space = "Lab"),            ylab= "Gear category", xlab="Species")

p2 <- levelplot(x=cor(cortest[,-1]), at=seq(-1,1,0.02),                  panel=my.panel, aspect= "fill", col.regions = colorRampPalette(c('dark blue','white','dark red'), space = "Lab"), ylab= "Species",      xlab="Species")

grid.arrange(p1,p2, ncol=2)


p3 <- levelplot(x=apply(corall,c(2,3),"mean"), at=seq(-1,1,0.02),         panel=my.panel, aspect= "fill", col.regions = colorRampPalette(c('dark blue','white','dark red'), space = "Lab"), ylab= "Species",      xlab="Species")
grid.arrange(p1,p3, ncol=2)


p4 <- levelplot(x=apply(corall_incl_veslen,c(2,3),"mean"), at=seq(-1,1,0.02),                          panel=my.panel, aspect= "fill", col.regions = colorRampPalette(c('dark blue','white','dark red'), space = "Lab"), ylab= "Species",      xlab="Species")

grid.arrange(p1,p4, ncol=2)
                
