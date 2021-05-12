setwd("M:\\RobinHood\\data\\landings and indices\\lemon sole index calcs")
lemIBTS <- read.csv("Lemon sole CPUE per length per area IBTS downloaded from DATRAS.csv")

lemIBTSagg <- aggregate(kgperhour~ Year + Quarter + Area, data=lemIBTS, FUN=sum) 

lemIBTSagg <- lemIBTSagg[!lemIBTSagg$Area==-9,] 
lemIBTSagg <- lemIBTSagg[!lemIBTSagg$Area==10,] 

lemIBTSagg <- lemIBTSagg[!lemIBTSagg$Quarter  %in% c(2,4),] 


xyplot(kgperhour~Year|Area + Quarter, data=lemIBTSagg)


lemIBTSindex <- aggregate( kgperhour~Year + Quarter, data= lemIBTSagg, FUN=mean)

xyplot(kgperhour~Year|Quarter, data=lemIBTSindex)



####
#per subarea

lemIBTS <- read.csv("Lemon sole CPUE per length per subarea IBTS downloaded from DATRAS.csv")

lemIBTSagg <- aggregate(kgperhour~ Year + Quarter + SubArea, data=lemIBTS, FUN=sum) 

########3
# select only rectangles that have been samples more than 10 times
######
rectcount <- aggregate(kgperhour~ Quarter + SubArea, data=lemIBTSagg, FUN="length")

lemIBTSagg <- merge(lemIBTSagg,rectcount, by=c("Quarter", "SubArea"))
lemIBTSagg <- lemIBTSagg[lemIBTSagg$kgperhour.y>25,]
lemIBTSagg

lemIBTSagg <- lemIBTSagg[!lemIBTSagg$Quarter  %in% c(2,4),] 


xyplot(kgperhour.x~Year|SubArea + Quarter, data=lemIBTSagg)


#check if there are rects without any lem sole in the entire time series: only one)
aggregate(kgperhour.x~ SubArea, data=lemIBTSagg, FUN="sum") 



lemIBTSindex <- aggregate( kgperhour.x~Year + Quarter, data= lemIBTSagg, FUN=mean)

xyplot(kgperhour.x~Year|Quarter, data=lemIBTSindex)

#select only data from 1980 onwards
lemIBTSindex <- lemIBTSindex[lemIBTSindex$Year >1979,]
write.csv(lemIBTSindex, file="Resulting lemon sole indices from IBTS.csv")

####################################################################
#BTS
###################################################################

lemBTS <- read.csv("Lemon sole CPUE per length per Hour and Swept Area BTS.csv")


lemBTS <- lemBTS[!lemBTS$Gear == "BT4AI",]


lemBTSagg <- aggregate(kgperhour~ Year + StatRec + Quarter + Gear, data=lemBTS, FUN=sum) 


lmindex <- lm(kgperhour~ as.factor(Year) + StatRec + Gear, data=lemBTSagg )


lemBTSindex <- aggregate( kgperhour~Year, data= lemIBTSagg, FUN=mean)
