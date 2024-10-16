##packages-------
library(tidyr); library(ggplot2); library(wesanderson); library(ggforce); 
library("cowplot")


####--------------------------------------------------#plotting the data#--------------------------------
##get the data from R_script(...insert_name...)(getting the data)----
##------preliminary plots(fast and ugly)-------

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  legend.title=element_text(size=14), 
  legend.text=element_text(size=13),
  plot.title = element_text(size=12))


#plot landings--------------
Lm_long<-gather(Lm,key=species, value="Landings", -years)

p <- ggplot(data = Lm_long, aes(x = years , y=as.numeric(Landings)/1000, colour=species))+
  ylab("Landings (1000 tonnes)") + xlab("Year") +
  geom_line(linetype="solid",size=0.8, alpha=0.4)+
  geom_point(size=3, aes(shape=species))+
  scale_shape_manual(values = c( 16, 15,17, 18,0))+
  scale_colour_manual(values = c("#3B9AB2","#9EBD91","#E8C520","#E29E00","#F21A00"))

p+ theme_light(base_size = 20)

#plot discards---------
Dm_long<-gather(Dm,key=species, value="Discards", -years)

Dm_long_cut<-Dm_long[Dm_long$species!="Dab",]

p <- ggplot(data = Dm_long_cut, aes(x = years ,y=as.numeric(Dm_long_cut$Discards)/1000,  colour=species))+
  ylab("Discards (1000 tonnes)") + xlab("Year") +
  geom_line(linetype="solid",size=0.8, alpha=0.4)+
  geom_point(size=3, aes(shape=species))+
  scale_shape_manual(values = c( 16, 15,17, 18,0))+
  scale_colour_manual(values = c("#EBCC2A", "#78B7C5", "#F21A00","#E1AF00","#3B9AB2"))#,"#3B9AB2" 
  #scale_y_continuous(trans = 'log10')
                       
p+ theme_light(base_size = 20)


#plot F's for actual visuals------
#actual lines of F
col<-wes_palette("Zissou1",10, type = "continuous")
col

Fm_long<-gather(Fm,key=species, value="val", -years)
colnames(Fm_long)[names(Fm_long) == "years"]<-"Year"
Fm_long[,1]<-as.numeric(Fm_long[,1])

#####################################################################
#here is the plot for the paper----------------
#####################################################################

ggplot(data = Fm_long,aes(x = Year, y = val ,colour= species, group=species))+
  ylab(expression(Fishing~mortality~(year^-1)))+
  geom_line(linetype="solid",size=1, alpha=0.4)+
  scale_colour_manual(values = c("#3B9AB2","#9EBD91","#E8C520","#E29E00","#F21A00"),labels=c("Cod (ages 2-4)","Haddock (ages 2-4)","Plaice (ages 2-6)","Sole (ages 2-6)","Whiting (ages 2-6)"))+
  #geom_point(size=2, aes(shape=species))+
  scale_shape_manual(values =rep(19,5),guide=F)+
  theme_bw()+
  #ggtitle("Fishing mortality etsimates from ICES")+
  geom_hline(aes(yintercept=0.152),colour="#E8C520",linetype="dashed", size=1, alpha=0.6 )+ #plaice
  geom_hline(aes(yintercept=0.157),colour="#E29E00",linetype="dashed", size=1, alpha=0.6 )+ #sole
  #geom_hline(aes(yintercept=0.31),colour="#3B9AB2",linetype="dashed" , size=1, alpha=0.6)+ #cod
  geom_hline(aes(yintercept=0.68),colour="#F21A00",linetype="dashed", size=1, alpha=0.6 )+ #whg
  geom_hline(aes(yintercept=0.174),colour= "#9EBD91",linetype="dashed", size=1, alpha=0.6 )+#haddock
  scale_x_continuous(limits = c(1975, 2025), expand = c(0, 0.0), breaks = c(seq(1980, 2030, by = 10)))+
  scale_y_continuous(limits = c(0.0, 1.3), expand = c(0, 0.0),breaks = c(0.1,seq(0.1, 1.3, by = 0.2)))+
  #My_Theme
  theme_light(base_size = 17)+
  theme(legend.position="bottom")+
  guides(col=guide_legend(nrow=2, byrow=TRUE))+
  theme(plot.margin=margin(t = 20, r = 20, b = 0, l = 5, unit = "pt"))


####################################################
#FIGURE WITH LANDINGS AND CATCHES FOR PAPER
####################################################

Cm_long<-gather(Cm,key=species, value="Catches", -years)


Lm_long$dsource <- "Landings"
names(Lm_long)[names(Lm_long)=="Landings"] <- "val"
Lm_long$csource <- "ICES"
Dm_long$dsource <- "Discards"
Cm_long$dsource <- "Catch"
Cm_long$Catches.1 <- NULL
names(Cm_long)[names(Cm_long)=="Catches"] <- "val"


#make csource of total catches "calculated" where we inferred discards 
Cm_long$csource <- "ICES"
Cm_long[Cm_long$species=="Brill" & Cm_long$years < 2014, ]$csource <- "inferred"
Cm_long[Cm_long$species=="Dab" & Cm_long$years < 2002, ]$csource <- "inferred"
Cm_long[Cm_long$species=="Dab" & Cm_long$years > 2016,  ]$csource <- "inferred.2"
Cm_long[Cm_long$species=="Lemon sole" & Cm_long$years < 2002,  ]$csource <- "inferred"
Cm_long[Cm_long$species=="Turbot" & Cm_long$years < 2013,  ]$csource <- "inferred"
Cm_long[Cm_long$species=="Witch" & Cm_long$years < 2002,  ]$csource <- "inferred"

Cm_comb <- rbind(Lm_long,Cm_long)
Cm_comb$species<-factor(Cm_comb$species, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))

ggplot(data = Cm_comb, aes(x = years ,y=as.numeric(Cm_comb$val)/1000, col=factor(dsource) ,shape=factor(csource) ) )+
  ylab("Catches (x 1000 tonnes)")+ xlab("Year")  + 
  geom_line(alpha=0.6) +
  geom_point(size=2.5)+
  scale_y_continuous(trans = "log10", breaks=c(1,5,10,50,100), limits=c(1,110))+
  scale_colour_manual(name="data type", values = c("#F8766D","#00BFC4"))+
  scale_shape_manual(name="data source",values=c(16,1,1), labels=c("ICES","inferred"), breaks=c("ICES","inferred"))+
  facet_grid(cols = vars(species)) +
  #My_Theme+
  #theme_bw()+
  theme_light(base_size = 17)+
  theme(axis.text.x = element_text(angle = 90))


#####################################################################################
#NOW INDICES: FIGURE FOR MANUSCRIPT 
####################################################################################

colnames(Im)[2:8]<-c("Turbot"  ,   "Brill"    ,  "Witch","Witch.1",      "Dab", "Lemon sole","Lemon sole.1")        

Im_long<-gather(Im,key=species, value="val", -years)
Im_long_simple<-Im_long

Im_long_simple$species_group<-Im_long_simple$species
Im_long_simple$species_group[Im_long_simple$species_group=="Witch.1"]<-"Witch"
Im_long_simple$species_group[Im_long_simple$species_group=="Lemon sole.1"]<-"Lemon sole"

Im_long_simple$quarter<-Im_long_simple$species
Im_long_simple$quarter<-1
Im_long_simple$quarter[Im_long_simple$species=="Dab"]<-2
Im_long_simple$quarter[Im_long_simple$species=="Witch"]<-2
Im_long_simple$quarter[Im_long_simple$species=="Lemon sole"]<-2

Im_long_simple$quarter[Im_long_simple$species=="Witch.1"]<-3
Im_long_simple$quarter[Im_long_simple$species=="Lemon sole.1"]<-3

Im_long_simple$species_group<-factor(Im_long_simple$species_group, levels=c("Dab","Brill","Lemon sole", "Turbot","Witch"))


ggplot(data = Im_long_simple, aes(x = years ,y=as.numeric(Im_long_simple$val),col=factor(quarter) ) )+
  ylab("Indices (kg/hour)")+ xlab("Year")  + 
  geom_line(size= 1,alpha=1) +
  scale_y_continuous(trans = "log10" ,breaks=c(0,0.05,0.1,0.5,1,5,10,20,40))+
  scale_colour_manual(name="index data origin",values = c("#FC717F","#6BB100","#00C0AF"),labels=c("Commercial","IBTS Q1","IBTS Q3"))+
  facet_grid(cols = vars(species_group))+
  My_Theme+
  #theme_bw(base_size = 16)+
  theme_light(base_size = 17)+
  theme(axis.text.x = element_text(angle = 90),plot.margin=unit(c(5.5,-2.5,5.5,5.5),"pt"))+
  annotation_logticks() 




####################################################
#JJP Data Rich F
####################################################
FMSY <- data.frame(species=unique((Fm_long$species)),Fmsy=c(0.21, 0.202,0.31,0.172,0.194))

ggplot(data = Fm_long,aes(x = Year, y = Fm_long$val, colour= species, shape=species))+
  ylab(expression(Fishing~mortality~(year^-1))) + xlab("Year") +
  geom_line(linetype="solid", alpha=0.4)+
  scale_y_continuous(limits=c(0,1.1),expand = c(0, 0.0))+
  geom_point(size=2.5)+
  geom_hline(data=FMSY, linetype = "dashed", mapping = aes(yintercept = Fmsy, color = species)) +
  theme_light(base_size = 15)




  
  
  
  
