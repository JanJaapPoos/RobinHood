#loop for making the final plot for differnt rh's
#getting results out of the objects-----
for(i in seq(1,1,0.05)){
  ### estimates for parameters in object------
  #let op, nog zelf de gewenste corrmult invoegen waarvan je de grafieken wilt kies uit:
  #(rep_rh0,rep_rh0.05,rep_rh0.1,rep_rh0.15,rep_rh0.2,rep_rh0.25,rep_rh0.3,rep_rh0.35,rep_rh0.4,rep_rh0.45,rep_rh0.5,rep_rh0.55,rep_rh0.6,rep_rh0.65,rep_rh0.7,rep_rh0.75,rep_rh0.8,rep_rh0.85,rep_rh0.9)
  rep.summ_rh <- summary(get(paste0("rep_rh",i)))
  
  #rep.summ_rh <- summary(rep_rh0.45)
  #print next line if you want to know the parameter estimates for different corrmults
  #rep.summ_rh[1:37,1:2]  
  
  
  ### estimates for calculated values--------------
  #B,BoBMSY,FoFMSY,Ihat,Sigma
  
  
  Bhat_rh   <- (matrix( rep.summ_rh[rownames(rep.summ_rh) == "B", "Estimate"], nc = mc)/1000)
  Blow_rh   <- Bhat_rh - 2 * (matrix(rep.summ_rh[rownames(rep.summ_rh) == "B", "Std. Error"],  nc = mc)/1000)
 # Blow_rh   <- apply(Blow_rh,c(1,2),function(x)if(x<0){0}else{x})
  Bhigh_rh  <- Bhat_rh + 2 * (matrix(rep.summ_rh[rownames(rep.summ_rh) == "B", "Std. Error"], nc = mc)/1000)
  
  BoBMSYhat_rh  <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "BoBMSY", "Estimate"], nc = mc)
  BoBMSYlow_rh  <- BoBMSYhat_rh - 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "BoBMSY", "Std. Error"],  nc = mc)
  #BoBMSYlow_rh <- apply(BoBMSYlow_rh,c(1,2),function(x)if(x<0){0}else{x})
  BoBMSYhigh_rh <- BoBMSYhat_rh + 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "BoBMSY", "Std. Error"], nc = mc)
  
  FoFMSYhat_rh  <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "FoFMSY", "Estimate"], nc = mc)
  FoFMSYlow_rh  <- FoFMSYhat_rh - 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "FoFMSY", "Std. Error"],  nc = mc)
  #FoFMSYlow_rh <- apply(FoFMSYlow_rh,c(1,2),function(x)if(x<0){0}else{x})
  FoFMSYhigh_rh <- FoFMSYhat_rh + 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "FoFMSY", "Std. Error"], nc = mc)
  
  logitUhat_rh  <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "logitU", "Estimate"], nc = mc+mcdr)
  Uhat_rh <- invcuslogit(logitUhat_rh)
  Ulow_rh  <- invcuslogit(logitUhat_rh - 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "logitU", "Std. Error"], nc = mc+mcdr))
  Ulow_rh <- apply(Ulow_rh,c(1,2),function(x)if(x<0){0}else{x})
  Uhigh_rh <- invcuslogit(logitUhat_rh + 2 * matrix(rep.summ_rh[rownames(rep.summ_rh) == "logitU", "Std. Error"], nc = mc+mcdr))
  
  BMSY_rh          <-  0.5* exp( matrix( rep.summ_rh[rownames(rep.summ_rh) == "logK", "Estimate"], nc=mc))
  
  #check fit to observations-----------
  Chat_rh <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "Chat", "Estimate"], nc = mc)
  Chat.se_rh <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "Chat", "Std. Error"], nc = mc)
  Clow_rh<-Chat_rh-2*Chat.se_rh
  Clow_rh <- apply(Clow_rh,c(1,2),function(x)if(x<0){0}else{x})
  Chigh_rh<-Chat_rh+2*Chat.se_rh
  
  Ihat_rh <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "Ihat", "Estimate"], nc = mi)
  Ihat.se_rh <- matrix( rep.summ_rh[rownames(rep.summ_rh) == "Ihat", "Std. Error"], nc = mi)
  Ilow_rh<-Ihat_rh-2*Ihat.se_rh
  Ilow_rh <- apply(Ilow_rh,c(1,2),function(x)if(x<0){0}else{x})
  Ihigh_rh<-Ihat_rh+2*Ihat.se_rh
  
  
  #plot estimates of Catches robin hood----
  #restructure the data
  colnames(Cm) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Cm_long<-gather(Cm,key=species, value="Catches", -years)
  mis_data_Chat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA)))
  colnames(mis_data_Chat) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  
  years_hat<-1978:2018
  
  
  Chat_incl_rh<-cbind(years_hat,Chat_rh)
  colnames(Chat_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Chat_incl_dataframe_rh<-as.data.frame(rbind(mis_data_Chat,Chat_incl_rh))
  C_hat_long_rh<-gather(Chat_incl_dataframe_rh, key=species,value="Est_Catches", -years)
  
  Chat.se_min_rh<- Chat_rh - (2*Chat.se_rh)
  Chat.se_min_incl_rh<-cbind(years_hat,Chat.se_min_rh)
  colnames(Chat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Chat.se_min_dataframe_rh<-as.data.frame(rbind(mis_data_Chat,Chat.se_min_incl_rh))
  Chat.se_min_long_rh<-gather(Chat.se_min_dataframe_rh, key=species,value="SE.min", -years)
  
  Chat.se_max_rh<- Chat_rh + (2*Chat.se_rh)
  Chat.se_max_incl_rh<-cbind(years_hat,Chat.se_max_rh)
  colnames(Chat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Chat.se_max_dataframe_rh<-as.data.frame(rbind(mis_data_Chat,Chat.se_max_incl_rh))
  Chat.se_max_long_rh<-gather(Chat.se_max_dataframe_rh, key=species,value="SE.max", -years)
  
  #------merge all indices into single dataframe for ggplot
  merge_catches_rh<-as.data.frame(cbind(Cm_long,C_hat_long_rh$Est_Catches,Chat.se_min_long_rh$SE.min, Chat.se_max_long_rh$SE.max ))
  merge_catches_rh[merge_catches_rh$years==2018,4]<- NA #c(4,5,6,9,10,11)
  merge_catches_rh[merge_catches_rh$years==2018,5]<- NA 
  merge_catches_rh[merge_catches_rh$years==2018,6]<- NA 
  colnames(merge_catches_rh) <- c("years", "species", "Catches", "est_Catches","Min","Max")
  #merge_catches_rh$species.1<-merge_catches_rh$species
  merge_catches_rh$species.1<- factor(merge_catches$species, levels=c("Dab","Brill", "Turbot","Lemon sole","Witch"))
  merge_catches_rh$Catches.1<-as.numeric(merge_catches_rh$Catches)/1000
  merge_catches_rh$est_catches.1<-as.numeric(merge_catches_rh$est_Catches)/1000
  merge_catches_rh$Min.1<-as.numeric(merge_catches_rh$Min)/1000
  merge_catches_rh$Max.1<-as.numeric(merge_catches_rh$Max)/1000
  
  C_plot_polite_rh<-
    ggplot(data = merge_catches_rh) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    #scale_y_continuous(trans = "log10" ,breaks=c(0,0.05,0.1,0.5,1,5,10,20,40), sec.axis = sec_axis(~., labels = NULL))+ #expand=c(0,0) ,limits=c(0,10))+
    geom_point(aes(x=years, y=as.numeric(merge_catches_rh$Catches.1), shape="16", colour=species.1),show.legend = F,size=1.2,alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow=1)+
    #scale_y_log10(limits = c(0.01,1e2))+
    #maybe change the scales manually
    theme_bw()+
    #5 COLOURS
    scale_colour_manual(values = col_values_5,guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(merge_catches_rh$est_catches.1),linetype="dashed",colour=species),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(merge_catches_rh$Min.1), ymax=as.numeric(merge_catches_rh$Max.1),fill=species.1), linetype="22",alpha=0.4) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    ggtitle("Estimates of Catches")
  
  
  #plot estimates of Indices robin hood----
  #restructure the data
  colnames(Im) <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
  Im_long<-gather(Im,key=species, value="Indices", -years)
  mis_data_Ihat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA,NA,NA)))
  colnames(mis_data_Ihat) <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
  
  years_hat<-1978:2018
  
  Ihat_incl_rh<-cbind(years_hat,Ihat_rh)
  colnames(Ihat_incl_rh) <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
  Ihat_incl_dataframe_rh<-as.data.frame(rbind(mis_data_Ihat,Ihat_incl_rh))
  I_hat_long_rh<-gather(Ihat_incl_dataframe_rh, key=species,value="Est_indices", -years)
  
  
  Ihat.se_min_rh<- Ihat_rh - (2*Ihat.se_rh)
  Ihat.se_min_incl_rh<-cbind(years_hat,Ihat.se_min_rh)
  colnames(Ihat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
  Ihat.se_min_dataframe_rh<-as.data.frame(rbind(mis_data_Ihat,Ihat.se_min_incl_rh))
  Ihat.se_min_long_rh<-gather(Ihat.se_min_dataframe_rh, key=species,value="SE.min", -years)
  
  Ihat.se_max_rh<- Ihat_rh + (2*Ihat.se_rh)
  Ihat.se_max_incl_rh<-cbind(years_hat,Ihat.se_max_rh)
  colnames(Ihat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Witch.1","Dab","Lemon sole","Lemon sole.1")
  Ihat.se_max_dataframe_rh<-as.data.frame(rbind(mis_data_Ihat,Ihat.se_max_incl_rh))
  Ihat.se_max_long_rh<-gather(Ihat.se_max_dataframe_rh, key=species,value="SE.max", -years)
  
  #------merge all indices into single dataframe for ggplot
  merge_indices_rh<-as.data.frame(cbind(Im_long,I_hat_long_rh$Est_indices,Ihat.se_min_long_rh$SE.min, Ihat.se_max_long_rh$SE.max,Ihat.se_max_long_rh$species ))
  colnames(merge_indices_rh) <- c("years", "species", "indices", "est_indices","Min","Max","species.1")
  merge_indices_rh$species.1<- factor(merge_indices_rh$species, levels=c("Dab","Brill","Turbot","Lemon sole","Lemon sole.1", "Witch","Witch.1"))
  
  merge_indices_rh$species_group<- merge_indices_rh$species
  merge_indices_rh$species_group[merge_indices_rh$species_group=="Witch.1"]<-"Witch"
  merge_indices_rh$species_group[merge_indices_rh$species_group=="Lemon sole.1"]<-"Lemon sole"
  
  #plot the data
  I_plot_polite_rh<-
    ggplot(data = merge_indices_rh) +
    ylab("Indices (kg/hour)")+
    xlab("Years")+
    geom_point(aes(x = years, y = as.numeric(merge_indices_rh$indices),shape="16",colour=species.1),show.legend = F,size=1.2, alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow=1)+
    #scale_y_log10(limits = c(0.01,1e2))+
    #maybe change the scales manually
    theme_bw()+
    scale_colour_manual(values = col_values_7,guide=F)+
    geom_line(aes(x=years,y=as.numeric(merge_indices_rh$est_indices),linetype="22", colour=species.1),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(merge_indices_rh$Min), ymax=as.numeric(merge_indices_rh$Max),fill=species.1),linetype="solid",alpha=0.4 ,show.legend = F)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_7,guide=F)+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    ggtitle("Estimates of Indices")
    #coord_cartesian(ylim = c(0, 7.5))
  
  
  #plot estimates of Uhat robin hood -----
  #restructure the data
  Uhat_rh<-cbind(years_hat,Uhat_rh)
  colnames(Uhat_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  Uhat_rh<-as.data.frame(Uhat_rh)
  Uhat_long_rh<-gather(Uhat_rh,key=species, value=U, -years)
  mis_data_Uhat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA)))
  mis_data_Uhat<-data.frame(cbind(mis_data_Uhat,Fm[1:3,2:6]))
  colnames(mis_data_Uhat) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  
  years_hat<-1978:2018
  
  
  Uhat_incl_dataframe_rh<-as.data.frame(rbind(mis_data_Uhat,Uhat_rh))
  Uhat_long_rh<-gather(Uhat_incl_dataframe_rh, key=species,value=U, -years)
  
  Uhat.se_min_rh<-Ulow_rh
  Uhat.se_min_incl_rh<-cbind(years_hat,Uhat.se_min_rh)
  colnames(Uhat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  mis_data_Uhat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)))
  colnames(mis_data_Uhat) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  Uhat.se_min_dataframe_rh<-as.data.frame(rbind(mis_data_Uhat,Uhat.se_min_incl_rh))
  Uhat.se_min_long_rh<-gather(Uhat.se_min_dataframe_rh, key=species,value="SE.min", -years)
  
  Uhat.se_max_rh<-Uhigh_rh
  Uhat.se_max_incl_rh<-cbind(years_hat,Uhat.se_max_rh)
  colnames(Uhat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  mis_data_Uhat<-data.frame(rbind(c(1975,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),c(1976,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),c(1977,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)))
  colnames(mis_data_Uhat) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole","Plaice","Sole","Cod","Whiting","Haddock")
  Uhat.se_max_dataframe_rh<-as.data.frame(rbind(mis_data_Uhat,Uhat.se_max_incl_rh))
  Uhat.se_max_long_rh<-gather(Uhat.se_max_dataframe_rh, key=species,value="SE.max", -years)
  
  #------merge all indices into single dataframe for ggplot
  merge_Uhat_rh<-as.data.frame(cbind(Uhat_long_rh,Uhat.se_min_long_rh$SE.min, Uhat.se_max_long_rh$SE.max,Uhat.se_max_long_rh$species ))
  colnames(merge_Uhat_rh) <- c("years", "species", "U","Min","Max","species.1")
  merge_Uhat_rh$species.1<- factor(merge_Uhat_rh$species, levels=c("Dab","Brill", "Turbot","Lemon sole","Witch","Cod","Plaice","Sole","Whiting"))
  
  merge_Uhat_cut_for_pres_rh<-merge_Uhat_rh[1:220,]
  
  #plot the data
  #U_plot_polite_rh<-
  ggplot(data=merge_Uhat_rh)+
    ylab("U (/year)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(merge_Uhat_rh$U),linetype="22",colour=species),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow =1)+
    theme_bw()+
    scale_colour_manual(values = wes_palette("Zissou1",9, type = "continuous"),guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(merge_Uhat_rh$Min), ymax=as.numeric(merge_Uhat_rh$Max),fill=species),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = wes_palette("Zissou1", 9, type = "continuous"),guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of Fishing mortality")
  
  
  U_plot_polite_rh<-
    ggplot(data=merge_Uhat_cut_for_pres_rh)+
    ylab("U (/year)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(merge_Uhat_cut_for_pres_rh$U),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow =1)+
    theme_bw()+
    #5 COLOURS
    scale_colour_manual(values = col_values_5, guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(merge_Uhat_cut_for_pres_rh$Min), ymax=as.numeric(merge_Uhat_cut_for_pres_rh$Max),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of Fishing mortality")
  
  
  #plot estimates of Bhat robin hood-----
  #restructure the data
  Bhat_rh<-cbind(years_hat,Bhat_rh)
  colnames(Bhat_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Bhat_rh<-as.data.frame(Bhat_rh)
  Bhat_long_rh<-gather(Bhat_rh,key=species, value=B, -years)
  
  Bhat.se_min_rh<-Blow_rh
  Bhat.se_min_incl_rh<-cbind(years_hat,Bhat.se_min_rh)
  colnames(Bhat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Bhat.se_min_incl_rh<-as.data.frame(Bhat.se_min_incl_rh)
  Bhat.se_min_long_rh<-gather(Bhat.se_min_incl_rh, key=species,value="SE.min", -years)
  
  Bhat.se_max_rh<-Bhigh_rh 
  Bhat.se_max_incl_rh<-cbind(years_hat,Bhat.se_max_rh)
  colnames(Bhat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  Bhat.se_max_incl_rh<-as.data.frame(Bhat.se_max_incl_rh)
  Bhat.se_max_long_rh<-gather(Bhat.se_max_incl_rh, key=species,value="SE.max", -years)
  
  
  #------merge all Biomass estimates into single dataframe for ggplot
  merge_Bhat_rh<-as.data.frame(cbind(Bhat_long_rh,Bhat.se_min_long_rh$SE.min, Bhat.se_max_long_rh$SE.max,Bhat.se_max_long_rh$species ))
  colnames(merge_Bhat_rh) <- c("years", "species", "B","Min","Max","species.1")
  merge_Bhat_rh$species.1<- factor(merge_Bhat_rh$species.1, levels=c("Dab","Brill","Turbot","Lemon sole","Witch"))
  
  #plot the data
  B_plot_polite_rh<-
    ggplot(data=merge_Bhat_rh)+
    ylab("Biomass (x1000 tonnes)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(merge_Bhat_rh$B),linetype="dashed",colour=species.1),linetype="22",size=0.75, alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow=1)+
    theme_bw()+
    #5 COLOURS
    scale_colour_manual(values = col_values_5, guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(merge_Bhat_rh$Min), ymax=as.numeric(merge_Bhat_rh$Max),fill=species.1),linetype="22",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of Biomasses")
  
  
  
  
  #plot estimates of BoBMSYhat robin hood-------
  #restructure the data
  
  BoBMSYhat_rh<-cbind(years_hat,BoBMSYhat_rh)
  colnames(BoBMSYhat_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  BoBMSYhat_rh<-as.data.frame(BoBMSYhat_rh)
  BoBMSYhat_long_rh<-gather(BoBMSYhat_rh,key=species, value=BoBMSY, -years)
  
  
  BoBMSYhat.se_min_rh<-BoBMSYlow_rh
  BoBMSYhat.se_min_incl_rh<-cbind(years_hat,BoBMSYhat.se_min_rh)
  colnames(BoBMSYhat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  BoBMSYhat.se_min_incl_rh<-as.data.frame(BoBMSYhat.se_min_incl_rh)
  BoBMSYhat.se_min_long_rh<-gather(BoBMSYhat.se_min_incl_rh, key=species,value="SE.min", -years)
  
  BoBMSYhat.se_max_rh<-BoBMSYhigh_rh 
  BoBMSYhat.se_max_incl_rh<-cbind(years_hat,BoBMSYhat.se_max_rh)
  colnames(BoBMSYhat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  BoBMSYhat.se_max_incl_rh<-as.data.frame(BoBMSYhat.se_max_incl_rh)
  BoBMSYhat.se_max_long_rh<-gather(BoBMSYhat.se_max_incl_rh, key=species,value="SE.max", -years)
  
  
  #------merge all BoBMSY estimates into single dataframe for ggplot
  merge_BoBMSYhat_rh<-as.data.frame(cbind(BoBMSYhat_long_rh,BoBMSYhat.se_min_long_rh$SE.min, BoBMSYhat.se_max_long_rh$SE.max,BoBMSYhat.se_max_long_rh$species ))
  colnames(merge_BoBMSYhat_rh) <- c("years", "species", "BoBMSY","Min","Max","species.1")
  merge_BoBMSYhat_rh$species.1<- factor(merge_BoBMSYhat_rh$species, levels=c("Dab","Brill", "Turbot","Lemon sole","Witch"))
  
  #plot the data
  
  BoBMSY_plot_polite_rh<-
    ggplot(data=merge_BoBMSYhat_rh)+
    ylab("B/BMSY")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(merge_BoBMSYhat_rh$B),linetype="dashed",colour=species),linetype="22",size=0.75, alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),nrow=1)+#,scales="free_y"
    theme_bw()+
    geom_hline(aes(yintercept=1,colour=species),linetype="twodash", size=0.8, alpha=0.4 )+
    scale_colour_manual(values = col_values_5,guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(merge_BoBMSYhat_rh$Min), ymax=as.numeric(merge_BoBMSYhat_rh$Max),fill=species),linetype="22",alpha=0.2)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of B/BMSY")+
    theme(plot.margin = margin(10,15 ,0, 10))+
    theme(panel.spacing = unit(1, "lines"))
  
  
  
  
  #plot estimates of FoFMSYhat robin hood-------
  #restructure the data
  
  FoFMSYhat_rh<-cbind(years_hat,FoFMSYhat_rh)
  colnames(FoFMSYhat_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  FoFMSYhat_rh<-as.data.frame(FoFMSYhat_rh)
  FoFMSYhat_long_rh<-gather(FoFMSYhat_rh,key=species, value=FoFMSY, -years)
  
  
  FoFMSYhat.se_min_rh<-FoFMSYlow_rh
  FoFMSYhat.se_min_incl_rh<-cbind(years_hat,FoFMSYhat.se_min_rh)
  colnames(FoFMSYhat.se_min_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  FoFMSYhat.se_min_incl_rh<-as.data.frame(FoFMSYhat.se_min_incl_rh)
  FoFMSYhat.se_min_long_rh<-gather(FoFMSYhat.se_min_incl_rh, key=species,value="SE.min", -years)
  
  FoFMSYhat.se_max_rh<-FoFMSYhigh_rh 
  FoFMSYhat.se_max_incl_rh<-cbind(years_hat,FoFMSYhat.se_max_rh)
  colnames(FoFMSYhat.se_max_incl_rh) <- c("years", "Turbot","Brill","Witch","Dab","Lemon sole")
  FoFMSYhat.se_max_incl_rh<-as.data.frame(FoFMSYhat.se_max_incl_rh)
  FoFMSYhat.se_max_long_rh<-gather(FoFMSYhat.se_max_incl_rh, key=species,value="SE.max", -years)
  
  #------merge all BoBMSY estimates into single dataframe for ggplot
  merge_FoFMSYhat_rh<-as.data.frame(cbind(FoFMSYhat_long_rh,FoFMSYhat.se_min_long_rh$SE.min, FoFMSYhat.se_max_long_rh$SE.max,FoFMSYhat.se_max_long_rh$species ))
  colnames(merge_FoFMSYhat_rh) <- c("years", "species", "FoFMSY","Min","Max","species.1")
  merge_FoFMSYhat_rh$species.1<- factor(merge_FoFMSYhat_rh$species, levels=c("Dab","Brill", "Turbot","Lemon sole","Witch"))
  
  #plot the data
  FoFMSY_plot_polite_rh<-
    ggplot(data=merge_FoFMSYhat_rh)+
    ylab("F/FMSY")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(merge_FoFMSYhat_rh$FoFMSY),linetype="dashed",colour=species),linetype="22",size=0.75, alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),nrow=1)+#,scales="free_y"
    theme_bw()+
    geom_hline(aes(yintercept=1,colour=species),linetype="twodash", size=0.8, alpha=0.4 )+
    scale_colour_manual(values = col_values_5,guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(merge_FoFMSYhat_rh$Min), ymax=as.numeric(merge_FoFMSYhat_rh$Max),fill=species),linetype="22",alpha=0.2)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of F/FMSY")+
    theme(plot.margin = margin(10,15 ,0, 10))+
    theme(panel.spacing = unit(1, "lines"))
  
  
  
  
  
  
  #combination plot of U, B, C, I robin hood--------------
  
  endplot<-
    
    ggdraw() +
    draw_plot(U_plot_polite_rh, y = 0.5, x = 0, width = 1, height = 0.5) +
    draw_plot(B_plot_polite_rh, y= 0, x = 0, width = 1, height = 0.5)
  
    ggdraw() +
    draw_plot(C_plot_polite_rh, y = 0.5, x = 0,width = 0.975, height = 0.5)+
    draw_plot(I_plot_polite_rh, y = 0, x = 0, width = 0.975, height = 0.5) 
    
  
  #give new plot a name and save as object in environment----------
  #might wanna change the name of "combiplot..."
  assign(paste0("combiplot_corrmult_", i),endplot)
  
  
  #save--------------
  #save.image("M:/robinhood/Code/R work/rh_UBIC_plotloop.RData")
  
  
  ############################################plot the combination plots of unrestricted and restricted data---------
  #extend the dataframes for ggplots
  #-----------------------------------------------------------catches dataframe:-----
  catch_combined<-merge_catches
  catch_combined$est_catches.1_rh<- merge_catches_rh$est_catches.1
  catch_combined$Min.1_rh<-merge_catches_rh$Min.1
  catch_combined$Max.1_rh<-merge_catches_rh$Max.1
  
  
  
  
  #plot combined data
ggplot(data = catch_combined) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Year")+
    geom_point(aes(x=years, y=as.numeric(catch_combined$Catches.1), shape="16", colour=species.1),show.legend = F,size=1.2,alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+
    #scale_y_log10(limits = c(0.01,1e2))+
    #maybe change the scales manually
    theme_bw()+
    scale_colour_manual(values = col_values_5,guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined$est_catches.1),linetype="dashed",colour=species.1),linetype="solid",alpha=0.8, size=0.75,show.legend = F)+
    geom_line(aes(x=years,y=as.numeric(catch_combined$est_catches.1_rh),linetype="dashed",colour=species.1),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined$Min.1), ymax=as.numeric(catch_combined$Max.1),fill=species.1), linetype="22",alpha=0.2) +
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined$Min.1_rh), ymax=as.numeric(catch_combined$Max.1_rh),fill=species.1), linetype="22",alpha=0.4) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Catches")
  C_plot_polite_combi
  
  C_plot_polite_combi<-ggplot(data = catch_combined) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Year")+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow=1)+
    geom_point(aes(x=years, y=as.numeric(catch_combined$Catches.1), shape="16"),show.legend = F,size=1.2,alpha=1)+
    theme_bw()+
    geom_line(aes(x=years,y=as.numeric(catch_combined$est_catches.1_rh)),colour="#00BFC4",alpha=0.8, size=0.75,show.legend = F)+
    geom_line(aes(x=years,y=as.numeric(catch_combined$est_catches.1)),linetype="22",colour="#F8766D",linetype="solid",alpha=0.6, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined$Min.1_rh), ymax=as.numeric(catch_combined$Max.1_rh),fill="00BFC4"),alpha=0.4)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined$Min.1), ymax=as.numeric(catch_combined$Max.1),fill= "F8766D"),alpha=0.2) +
    scale_fill_manual(values = c("#00BFC4","#F8766D"),name="data type",labels = c("Constrained", "Unconstrained"))+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    My_Theme+
    theme_light(base_size = 17)+
    theme(axis.text.x = element_text(angle = 90))+
    theme(plot.margin=unit(c(5.5,40.5,5.5,5.5),"pt"))
    #scale_x_continuous( limits=c(1975, 2020), expand = c(0, 0))+
    #guides(fill=guide_legend(nrow=2, byrow=TRUE))
  #scale_linetype_manual(values=c("solid","dashed"), name="data type")+
  
  
  #------------------------------------------------------------indices dataframe:----
  
  indices_combined<-merge_indices
  indices_combined$est_indices_rh<- merge_indices_rh$est_indices
  indices_combined$Min_rh<-merge_indices_rh$Min
  indices_combined$Max_rh<-merge_indices_rh$Max
  
  

  
    data_group<-c()
      for(i in as.vector(unique(indices_combined$species.1))){
        list<- c(rep(1,(length(indices_combined$species.1[indices_combined$species.1==i])/2)),rep(2,(length(indices_combined$species.1[indices_combined$species.1==i])/2)))
      data_group<-c(data_group,list)
        print(data_group)
      }
  length(data_group)
  

  indices_combined$data_group<-as.vector(data_group)
  indices_combined$data_group<-indices_combined$species
  indices_combined$data_group[indices_combined$data_group=="Witch.1"]<-"Witch"
  indices_combined$data_group[indices_combined$data_group=="Lemon sole.1"]<-"Lemon sole"
  
  
  
      
  
  #plot combined data
    ggplot(data = indices_combined, aes(colour=data_group)) +
    ylab("Indices (kg/hour)")+
    xlab("Years")+
    geom_point(aes(x = years, y = as.numeric(merge_indices_rh$indices),shape="16"),show.legend = F,size=1.2, alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",nrow=1)+
    #scale_y_log10(limits = c(0.01,1e2))+
    #maybe change the scales manually
    theme_bw()+
    #scale_colour_manual(values = col_values_7)+#,guide=F
    #scale_colour_manual(values = col_values_7)+#,guide=F
    geom_line(aes(x=years,y=as.numeric(indices_combined$est_indices_rh),linetype="22"),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
    geom_line(aes(x=years,y=as.numeric(indices_combined$est_indices),linetype="22"),linetype="solid",alpha=0.8, size=0.75,show.legend = F)+
    #geom_ribbon(aes(x=years,ymin=as.numeric(indices_combined$Min_rh), ymax=as.numeric(indices_combined$Max_rh)),alpha=0.4)+
    #geom_ribbon(aes(x=years,ymin=as.numeric(indices_combined$Min), ymax=as.numeric(indices_combined$Max)),alpha=0.2)+
    #scale_alpha_manual(guide=F)+
    #scale_fill_manual(values = rep("00BFC4",7),name="Species",guide=F)+ 
    #scale_fill_manual(values = col_values_7,name="Species",guide=F)+ 
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    #theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Indices")+
    theme(legend.direction  = "horizontal")
  
  
  
  
  I_plot_polite_combi<-ggplot(data=indices_combined)+
  ylab("Indices (kg/hour)")+
  xlab("Year")+
  facet_wrap(facets = vars(species.1),scales="free_y", nrow=1)+
  geom_point(aes(x = years, y = as.numeric(merge_indices_rh$indices),shape="16"),show.legend = F,size=1.2, alpha=1)+  
  theme_bw()+
  geom_line(aes(x=years,y=as.numeric(indices_combined$est_indices_rh)),col="#00BFC4",alpha=0.8, size=0.75, show.legend = F)+
  geom_line(aes(x=years,y=as.numeric(indices_combined$est_indices)),linetype="22", col= "#F8766D",alpha=0.6, size=0.75,show.legend=F)+
  geom_ribbon(aes(x=years,ymin=as.numeric(indices_combined$Min_rh), ymax=as.numeric(indices_combined$Max_rh), fill="00BFC4"),alpha=0.4,show.legend = F)+
  geom_ribbon(aes(x=years,ymin=as.numeric(indices_combined$Min), ymax=as.numeric(indices_combined$Max),fill= "F8766D"),alpha=0.2,show.legend = F)+
  scale_fill_manual(values = c("#00BFC4","#F8766D"),name="data type",labels = c("Constrained", "Unconstrained"))+
  My_Theme+
  theme_light(base_size = 17)+
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position="bottom",plot.margin=unit(c(5.5,12.5,5.5,5.5),"pt"),axis.title.x=element_blank())+
  scale_x_continuous( limits=c(1975, 2020), expand = c(0, 0))
  #scale_linetype_manual(values=c("solid","dashed"), name="data type")+
  
  #in case we want the guide, it is hidden in this graph here!!!
    
  
  C_I_plot<-ggdraw() +
    draw_plot(I_plot_polite_combi, x = 0, y = 0.5, width = 1, height = .5) +
    draw_plot(C_plot_polite_combi, x = 0, y = 0, width = 1, height = .5)
  C_I_plot
  
  
  assign(paste0("C_I_combined", i),C_I_plot)
  
  
  #------------------------------------------------------------Uhat dataframe:-----
  
  Uhat_combined<-merge_Uhat
  Uhat_combined$U_rh<- merge_Uhat_rh$U
  Uhat_combined$Min_rh<-merge_Uhat_rh$Min
  Uhat_combined$Max_rh<-merge_Uhat_rh$Max
  
  Uhat_combined_cut_for_pres<-Uhat_combined[1:220,]
  
  
  
  ggplot(data=Uhat_combined)+
    ylab("U")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined$U),linetype="22",colour=species),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined$U_rh),linetype="22",colour=species),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6)]),guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined$Min), ymax=as.numeric(Uhat_combined$Max),fill=species),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined$Min_rh), ymax=as.numeric(Uhat_combined$Max_rh),fill=species),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(1,3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6)]),guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of Fishing mortality")
  
  U_plot_polite_combi_cut<-ggplot(data=Uhat_combined_cut_for_pres)+
    ylab("U")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined_cut_for_pres$U),linetype="22",colour=species.1),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined_cut_for_pres$U_rh),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = col_values_5,guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined_cut_for_pres$Min), ymax=as.numeric(Uhat_combined_cut_for_pres$Max),fill=species.1),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined_cut_for_pres$Min_rh), ymax=as.numeric(Uhat_combined_cut_for_pres$Max_rh),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Fishing mortality")
  
  
  
  #------------------------------------------------------------Bhat dataframe:------
  Bhat_combined<-merge_Bhat
  Bhat_combined$B_rh<- merge_Bhat_rh$B
  Bhat_combined$Min_rh<-merge_Bhat_rh$Min
  Bhat_combined$Max_rh<-merge_Bhat_rh$Max
  
  
  
  B_plot_polite_combi<-ggplot(data=Bhat_combined)+
    ylab("Biomass (x 1000 tonnes)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(Bhat_combined$B.1),linetype="22",colour=species.1),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(Bhat_combined$B_rh),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = col_values_5,guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Bhat_combined$Min.1), ymax=as.numeric(Bhat_combined$Max.1),fill=species.1),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Bhat_combined$Min_rh), ymax=as.numeric(Bhat_combined$Max_rh),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = col_values_5,guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Biomasses")
  
  
  #------------------------------------------------------------BoBMSYhat dataframe:------
  BoBMSYhat_combined<-merge_BoBMSYhat
  BoBMSYhat_combined$BoBMSY_rh<- merge_BoBMSYhat_rh$BoBMSY
  BoBMSYhat_combined$Min_rh<-merge_BoBMSYhat_rh$Min
  BoBMSYhat_combined$Max_rh<-merge_BoBMSYhat_rh$Max
  
  
  BoBMSY_plot_polite_combi<-ggplot(data=BoBMSYhat_combined)+
    ylab("B/BMSY")+
    xlab("Year")+
    facet_wrap(facets = vars(species.1),nrow=1)+#,scales="free_y"
    theme_bw()+
    geom_hline(aes(yintercept=1),linetype="twodash", size=0.8, alpha=0.4 )+
    geom_line(aes(x = years, y = as.numeric(BoBMSYhat_combined$BoBMSY_rh),colour=species.1),col="#00BFC4",size=0.75, alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(BoBMSYhat_combined$BoBMSY),colour=species.1),linetype="22",col="#F8766D",size=0.75, alpha=0.6,show.legend=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(BoBMSYhat_combined$Min_rh), ymax=as.numeric(BoBMSYhat_combined$Max_rh),fill="00BFC4"),alpha=0.4,show.legend = F)+  
    geom_ribbon(aes(x=years, ymin=as.numeric(BoBMSYhat_combined$Min), ymax=as.numeric(BoBMSYhat_combined$Max),fill="F8766D"),alpha=0.2,show.legend = F)+
    scale_y_continuous( breaks=c(seq(-1.5,3.5,0.5)), limits = c(-1.25, 3.2), expand=c(0,0))+
    scale_fill_manual(values = c("#00BFC4","#F8766D"),name="data type",labels = c("Constrained", "Unconstrained"))+
    scale_alpha_manual(guide=F)+
    My_Theme+
    theme_light(base_size = 17)+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    theme(panel.spacing = unit(1, "lines"))+
    theme(legend.position="bottom",plot.margin=unit(c(-0.5,20.5,5.5,5.5),"pt"))
    #                                   margin(t = 0, r = 8, b = 0, l = 0)
    #scale_y_continuous(breaks = function(x) seq(from = 0, 
    #                                            to = 4, 
    #                                            length.out = 9))+
    #if (max(pretty(BoBMSYhat_combined$Max_rh))>3){coord_cartesian(ylim=c(0,3))}
  
  
  #------------------------------------------------------------FoFMSYhat dataframe:------
  FoFMSYhat_combined<-merge_FoFMSYhat
  FoFMSYhat_combined$FoFMSY_rh<- merge_FoFMSYhat_rh$FoFMSY
  FoFMSYhat_combined$Min_rh<-merge_FoFMSYhat_rh$Min
  FoFMSYhat_combined$Max_rh<-merge_FoFMSYhat_rh$Max
  
  
  FoFMSY_plot_polite_combi<-ggplot(data=FoFMSYhat_combined)+
    ylab("F/FMSY")+
    xlab("Year")+
    facet_wrap(facets = vars(species.1),nrow =1)+#,scales="free_y"
    theme_bw()+
    geom_hline(aes(yintercept=1),colour="grey",linetype="twodash", size=0.8, alpha=0.6 )+
    geom_line(aes(x = years, y = as.numeric(FoFMSYhat_combined$FoFMSY_rh)),col="#00BFC4",size=0.75, alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(FoFMSYhat_combined$FoFMSY)),col= "#F8766D",linetype="22",size=0.75, alpha=0.6,show.legend=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(FoFMSYhat_combined$Min_rh), ymax=as.numeric(FoFMSYhat_combined$Max_rh),fill="00BFC4"),alpha=0.4, show.legend = F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(FoFMSYhat_combined$Min), ymax=as.numeric(FoFMSYhat_combined$Max),fill="F8766D"),alpha=0.2, show.legend = F)+
    scale_fill_manual(values = c("#00BFC4","#F8766D"),name="data type",labels = c("Constrained", "Unconstrained"))+
    scale_alpha_manual(guide=F)+
    My_Theme+
    scale_y_continuous( breaks=c(seq(-1.5,3.5,0.5)), limits = c(-1.25, 3.0), expand=c(0,0))+
    theme_light(base_size = 17)+
    theme(legend.position="bottom",plot.margin=unit(c(5.5,20.5,0,5.5),"pt"), panel.spacing.x = unit(1,"lines"))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank())
    #scale_y_continuous(breaks = function(x) seq(from = 0, 
    #                                            to = 4, 
    #                                            length.out = 9) 
    #)+
    #if (max(pretty(FoFMSYhat_combined$Max_rh))>3.5){coord_cartesian(ylim=c(0,3.5))}
  
  
  
  #plot the combi plots-----------
  
  B_F_oMSY_plot<-ggdraw() +
    draw_plot(FoFMSY_plot_polite_combi, x = 0, y = 0.5, width = 1, height = .5)+
    draw_plot(BoBMSY_plot_polite_combi, x = 0, y = 0, width = 1, height = .5)
     
  B_F_oMSY_plot
  
  assign(paste0("overMSY_combined", i),B_F_oMSY_plot)
  
  #loop-part for the extended frames with fakewitch-----  
  #-----create split plot for biomass Dab---------
  cut_dab<-Bhat_combined[(Bhat_combined$species=="Dab"),]  
  #cut_dab$years<- c("'78",	"'79",	"'80",	"'81",	"'82",	"'83",	"'84",	"'85",	"'86",	"'87",	"'88",	"'89",	"'90",	"'91",	"'92",	"'93",	"'94",	"'95",	"'96",	"'97",	"'98",	"'99",	"'00",	"'01",	"'02",	"'03",	"'04",	"'05",	"'06",	"'07",	"'08",	"'09",	"'10",	"'11",	"'12",	"'13",	"'14",	"'15",	"'16",	"'17",	"'18")
  #cut_dab$years<-factor(cut_dab$years,levels=c("'78",	"'79",	"'80",	"'81",	"'82",	"'83",	"'84",	"'85",	"'86",	"'87",	"'88",	"'89",	"'90",	"'91",	"'92",	"'93",	"'94",	"'95",	"'96",	"'97",	"'98",	"'99",	"'00",	"'01",	"'02",	"'03",	"'04",	"'05",	"'06",	"'07",	"'08",	"'09",	"'10",	"'11",	"'12",	"'13",	"'14",	"'15",	"'16",	"'17",	"'18"))
  
  not_dab<-Bhat_combined[(Bhat_combined$species!="Dab"),]  
  
  B_plot_combi_cut_dab<-ggplot(data=cut_dab)+
    ylab("")+
    geom_line(aes(x = years, y = as.numeric(cut_dab$B.1),linetype="22",colour=species),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(cut_dab$B_rh),linetype="22",colour=species),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    #facet_wrap(facets = vars(species),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = wes_palette("Zissou1",6, type = "continuous")[c(1,2,3,4,6)],guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(cut_dab$Min.1), ymax=as.numeric(cut_dab$Max.1),fill=species),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(cut_dab$Min_rh), ymax=as.numeric(cut_dab$Max_rh),fill=species),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = wes_palette("Zissou1", 6, type = "continuous")[c(1)], guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    ggtitle("Estimates of Biomasses")+
    #facet_zoom(ylim=c(min(cut_dab$Min_rh), max(cut_dab$Max_rh)), zoom.size=1)+
    theme(axis.title.x=element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(plot.margin = margin(10,15 ,0, 10))+
    coord_cartesian(ylim=c(0,1000))
  
  #library(tidyr)
  
  
  B_plot_combi_not_dab<-ggplot(data=not_dab)+
    ylab("Biomass (x 1000 tonnes)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(not_dab$B.1),linetype="22",colour=species.1),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(not_dab$B_rh),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6)] ),guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(not_dab$Min.1), ymax=as.numeric(not_dab$Max.1),fill=species.1),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(not_dab$Min_rh), ymax=as.numeric(not_dab$Max_rh),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(wes_palette("Zissou1",6, type = "continuous")[c(3)],wes_palette("Zissou1",7, type = "continuous")[c(4,5,6)] ),guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    theme(plot.margin = margin(10,15 ,0, 10))
  #facet_zoom(ylim=c(min(cut_dab$Min_rh), max(cut_dab$Max_rh)), zoom.size=1)+
  #geom_point(data = not_dab %>% 
  #            group_by(species.1) %>% #group by facet variable
  #           summarise(y.min = if(pretty(Min.1)[1]<0){
  #                                                            0
  #                                                                } else{
  #                                                                        pretty(not_dab$Min.1)[1]
  #                                                            },
  #                     y.max = if (250<max(c(pretty(Max_rh)[length(pretty( Max_rh))],pretty(Max.1)[length(pretty(Max.1))]))){
  #                       (pretty(B_rh)[length(pretty(B_rh))]+100)
  #                                                                } else{
  #                       max(c(pretty(Max_rh)[length(pretty(Max_rh))],pretty(Max.1)[length(pretty(Max.1))]))
  #                          }
  #           )%>%
  #           gather(key, value, -species.1), 
  #         aes(x = 1, y = value),
  #         inherit.aes = FALSE, alpha = 0) +
  #scale_y_continuous(breaks = function(x) seq(from = x[1], 
  #                                           to = x[2], 
  #                                            length.out = 3), 
  #                   expand = c(0, 0))+
  #coord_cartesian(xlim=c(1975,2020))
  
  
  my_breaks <- function(x) { if (max(x) < 6000) seq(0, 5000, 1000) else seq(0, 15000, 5000) }
  scale_y_continuous(breaks = my_breaks)
  
  #if (test_expression) {
  #  statement1
  #} else {
  #  statement2
  #}
  
  #------------------extend plots of biomass, fishing mortality and catches---------
  #---for biomass
  fake_witch.1<-not_dab[1:41,]
  fake_lem.1<-not_dab[1:41,]
  fake_witch.1$species<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  fake_witch.1$species.1<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species.1<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  fake_witch.1[,3:5]<-NA
  fake_lem.1[,3:5]<-NA
  fake_witch.1[,7:9]<-NA
  fake_lem.1[,7:9]<-NA
  
  not_dab_extended<-as.data.frame(rbind(cut_dab,not_dab,fake_witch.1,fake_lem.1))#cut_dab can be removed if that helps to create better visual
  
  not_dab_extended$species.1<- factor(not_dab_extended$species.1, levels=c("Dab", "Brill","Turbot","Lemon sole", "Witch","Witch.1","Lemon sole.1"))
  
  #----for catches
  fake_witch.1<-catch_combined[1:41,]
  fake_lem.1<-catch_combined[1:41,]
  fake_witch.1$species<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  fake_witch.1$species.1<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species.1<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  head(fake_witch.1)
  fake_witch.1[,3:6]<-NA
  fake_lem.1[,3:6]<-NA
  fake_witch.1[,8:14]<-NA
  fake_lem.1[,8:13]<-NA
  
  catch_combined_extended<-as.data.frame(rbind(catch_combined,fake_witch.1,fake_lem.1))
  catch_combined_extended$species.1<- factor(catch_combined_extended$species.1, levels=c("Dab","Brill","Turbot","Lemon sole","Lemon sole.1", "Witch","Witch.1"))
  #----for fishing mortality
  fake_witch.1<-Uhat_combined_cut_for_pres[1:44,]
  fake_lem.1<-Uhat_combined_cut_for_pres[1:44,]
  fake_witch.1$species<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  fake_witch.1$species.1<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species.1<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  head(fake_witch.1)
  fake_witch.1[,3:5]<-NA
  fake_lem.1[,3:5]<-NA
  fake_witch.1[,7:9]<-NA
  fake_lem.1[,7:9]<-NA
  
  Uhat_combined_cut_for_pres_extended<-as.data.frame(rbind(Uhat_combined_cut_for_pres,fake_witch.1,fake_lem.1))
  Uhat_combined_cut_for_pres_extended$species.1<- factor(Uhat_combined_cut_for_pres_extended$species.1, levels=c("Dab","Brill","Turbot","Lemon sole", "Witch","Witch.1","Lemon sole.1"))
  
  #----for catch
  fake_witch.1<-catch_combined[1:44,]
  fake_lem.1<-catch_combined[1:44,]
  fake_witch.1$species<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  fake_witch.1$species.1<-rep("Witch.1",length(fake_witch.1[,1]))
  fake_lem.1$species.1<-rep("Lemon sole.1",length(fake_lem.1[,1]))
  head(fake_witch.1)
  fake_witch.1[,3:6]<-NA
  fake_lem.1[,3:6]<-NA
  fake_witch.1[,8:14]<-NA
  fake_lem.1[,8:14]<-NA
  
  catch_combined_extended<-as.data.frame(rbind(catch_combined,fake_witch.1,fake_lem.1))
  catch_combined_extended$species.1<- factor(catch_combined_extended$species.1, levels=c("Dab","Brill","Turbot","Lemon sole", "Witch","Witch.1","Lemon sole.1"))
  
  
  #make the plot objects---------
  #-------for biomass
  B_plot_combi_not_dab_ext<-ggplot(data=not_dab_extended)+
    ylab("Biomass (x 1000 tonnes)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(not_dab_extended$B.1),linetype="22",colour=species.1),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(not_dab_extended$B_rh),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(not_dab_extended$Min.1), ymax=as.numeric(not_dab_extended$Max.1),fill=species.1),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(not_dab_extended$Min_rh), ymax=as.numeric(not_dab_extended$Max_rh),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    theme(plot.margin = margin(10,15 ,0, 10))
  #facet_zoom(ylim=c(min(cut_dab$Min_rh), max(cut_dab$Max_rh)), zoom.size=1)
  
  
  B_plot_combi_not_dab_ext
  #library(gridExtra)
  #library(grid)
  something <- ggplotGrob(B_plot_combi_not_dab_ext)
  ## remove empty panels
  rm_grobs <- something$layout$name %in% c("panel-1-1","panel-1-6","panel-1-7","xlab-b","axis-l-1-1","axis-l-6-1","axis-l-7-1")
  # remove grobs
  something$grobs[rm_grobs] <- NULL
  something$layout <- something$layout[!rm_grobs, ]
  ## move axis-title closer to panel
  something$layout[something$layout$name == "axis-b-1-7",c("t", "b")]=c(30,30) #if not dab, then 24,24
  #something$layout[something$layout$name == "ylab-l", c("t","b")] = c(17,14)
  
  
  #install.packages("ggplotify")
  #library(ggplotify)
  #install.packages("ggpubr")
  #library(ggpubr)
  #B_plot_combi_not_dab_ext<-grid.draw(something)
  B_plot_combi_not_dab_ext<-as_ggplot(something)
  #B_plot_combi_not_dab_ext
  
  
  ###-----do the same for U and C as for B 
  
  U_plot_polite_combi_cut_ext<-ggplot(data=Uhat_combined_cut_for_pres_extended)+
    ylab("U (/year)")+
    xlab("Years")+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined_cut_for_pres_extended$U),linetype="22",colour=species.1),size=0.75 ,linetype="solid",alpha=0.8,show.legend=F)+
    geom_line(aes(x = years, y = as.numeric(Uhat_combined_cut_for_pres_extended$U_rh),linetype="22",colour=species.1),size=0.75 ,linetype="22",alpha=0.8,show.legend=F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol =1)+
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined_cut_for_pres_extended$Min), ymax=as.numeric(Uhat_combined_cut_for_pres_extended$Max),fill=species.1),linetype="solid",alpha=0.2)+
    geom_ribbon(aes(x=years, ymin=as.numeric(Uhat_combined_cut_for_pres_extended$Min_rh), ymax=as.numeric(Uhat_combined_cut_for_pres_extended$Max_rh),fill=species.1),linetype="solid",alpha=0.4)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #scale_linetype_manual(name="data type", values=c("22",""), labels=c("estimates",""))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Fishing mortality")
  
  
  #library(gridExtra)
  #library(grid)
  something <- ggplotGrob(U_plot_polite_combi_cut_ext)
  ## remove empty panels
  rm_grobs <- something$layout$name %in% c("panel-1-6","panel-1-7","xlab-b")
  # remove grobs
  something$grobs[rm_grobs] <- NULL
  something$layout <- something$layout[!rm_grobs, ]
  ## move axis-title closer to panel
  #something$layout[something$layout$name == , c("t","b")] = c(34,34)
  something$layout[something$layout$name == "axis-b-1-7",c("t", "b")]=c(30,30)
  
  #install.packages("ggplotify")
  #library(ggplotify)
  #B_plot_combi_not_dab_ext<-grid.draw(something)
  U_plot_polite_combi_cut_ext<-as_ggplot(something)
  #U_plot_polite_combi_cut_ext
  
  #---for catch
  C_plot_polite_combi_ext<-ggplot(data = catch_combined_extended) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    geom_point(aes(x=years, y=as.numeric(catch_combined_extended$Catches.1), shape="16", colour=species.1),show.legend = F,size=1.2,alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+
    #scale_y_log10(limits = c(0.01,1e2))+
    #maybe change the scales manually
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),linetype="dashed",colour=species.1),linetype="solid",alpha=0.8, size=0.75,show.legend = F)+
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1_rh),linetype="dashed",colour=species.1),linetype="22",alpha=0.8, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1), ymax=as.numeric(catch_combined_extended$Max.1),fill=species.1), linetype="22",alpha=0.2) +
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1_rh), ymax=as.numeric(catch_combined_extended$Max.1_rh),fill=species.1), linetype="22",alpha=0.4) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    ggtitle("Estimates of Catches")
  
  something <- ggplotGrob(C_plot_polite_combi_ext)
  ## remove empty panels
  rm_grobs <- something$layout$name %in% c("panel-1-6","panel-1-7","xlab-b")
  # remove grobs
  something$grobs[rm_grobs] <- NULL
  something$layout <- something$layout[!rm_grobs, ]
  ## move axis-title closer to panel
  something$layout[something$layout$name == "axis-b-1-7",c("t", "b")]=c(30,30)
  
  C_plot_polite_combi_ext<-as_ggplot(something)
  #C_plot_polite_combi_ext
  
  
  
  #------------draw the plots-----
  
  
  B_plot_polite_combi_alt<-ggdraw()+
    draw_plot(B_plot_combi_not_dab_ext, x = 0, y = 0.01, width = 1, height = 0.965) +
    draw_plot(B_plot_combi_cut_dab,     x = 0, y = 0.833, width = 1, height = 0.172)
  #B_plot_polite_combi_alt
  
  U_plot_polite_combi_alt<-ggdraw()+
    draw_plot(U_plot_polite_combi_cut_ext, x = 0, y = 0, width = 1, height = 1)
  
  #U_plot_polite_combi_alt
  
  C_plot_polite_combi_ext<-ggdraw()+
    draw_plot(C_plot_polite_combi_ext,x = 0, y = 0, width = 1, height = 1)
  #C_plot_polite_combi_ext
  
  
  I_plot_polite_combi
  
  #make some proper legends-----------------------------------------
  #--name_legend-------
  C_plot_polite_legend<-ggplot(data = catch_combined_extended) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    geom_point(aes(x=years, y=as.numeric(catch_combined_extended$Catches.1), shape="16", colour=species.1),size=1.2,alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),colour=species.1,linetype="dashed"),alpha=1, size=0.75, show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1), ymax=as.numeric(catch_combined_extended$Max.1),fill=species.1), linetype="22",alpha=0.4,show.legend = F) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    scale_shape_manual(name="data type", values=c(""), labels=c(""))+
    #scale_shape_manual(name="data type", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    ggtitle("Estimates of Catches")  
  C_plot_polite_legend
  
  
  
  legend_name<-get_legend(C_plot_polite_legend)
  as_ggplot(legend_name)
  #--observation_legend-------
  C_plot_polite_legend<-ggplot(data = catch_combined_extended) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    geom_point(aes(x=years, y=as.numeric(catch_combined_extended$Catches.1), shape="16", colour=species.1),size=1.2,alpha=1)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),colour=species.1,linetype="dashed"),alpha=1, size=0.75, show.legend = F)+
    #geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1), ymax=as.numeric(merge_catches_extended$Max.1),fill=species.1), linetype="22",alpha=0.4,show.legend = F) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #scale_shape_manual(name="data type", values=c(""), labels=c(""))+
    scale_shape_manual(name="", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    ggtitle("Esxtimates of Catches")  
  C_plot_polite_legend
  legend_obs<- get_legend(C_plot_polite_legend)
  as_ggplot(legend_obs)
  #--estimate_legend-------
  C_plot_polite_legend<-ggplot(data = catch_combined_extended) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    geom_point(aes(x=years, y=as.numeric(catch_combined_extended$Catches.1), shape="16", colour=species.1),size=1.2,alpha=1, show.legend = F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
    theme_bw()+
    scale_colour_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),colour=species.1,linetype="dashed"),alpha=1, size=0.75)+
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),colour=species.1,linetype="solid"),alpha=1, size=0.75)+
    #geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1), ymax=as.numeric(merge_catches_extended$Max.1),fill=species.1), linetype="22",alpha=0.4,show.legend = F) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values = c(col_values_5,col_values_5[1],col_values_5[1]),guide=F)+
    #scale_shape_manual(name="data type", values=c(""), labels=c(""))+
    #scale_shape_manual(name="", values=c(16), labels=c("observations"))+
    scale_linetype_manual(name="", values=c("22","solid"), labels=c("estimates(constrained)","estimates(unconstrained)"))+
    ggtitle("Estimates of Catches")  
  legend_est<-get_legend(C_plot_polite_legend)
  as_ggplot(legend_est)
  #species legend-----
  
  catch_combined_extended$species.1<- factor(catch_combined_extended$species.1, levels=c("Dab","Brill","Turbot","Lemon sole","Lemon sole.1", "Witch","Witch.1"))
  
  C_plot_polite_legend<-ggplot(data = catch_combined_extended) +
    ylab("Catch (x 1000 tonnes)")+
    xlab("Years")+
    geom_point(aes(x=years, y=as.numeric(catch_combined_extended$Catches.1), shape="16", colour=species.1),size=1.2,alpha=1, show.legend = F)+
    facet_wrap(facets = vars(species.1),scales="free_y",ncol=1)+#,ncol=1
    theme_bw()+
    scale_colour_manual(values = c(col_values_7))+
    #theme(legend.position = c(0, 1),legend.justification = c(0, 1))+
    #scale_shape_manual(values = c(rep(16,5)))
    #theme(legend.position = "none")
    geom_line(aes(x=years,y=as.numeric(catch_combined_extended$est_catches.1),colour=species.1,linetype="dashed"),alpha=1, size=0.75,show.legend = F)+
    geom_ribbon(aes(x=years,ymin=as.numeric(catch_combined_extended$Min.1), ymax=as.numeric(catch_combined_extended$Max.1),fill=species.1), linetype="22",alpha=1) +
    scale_alpha_manual(guide=F)+
    scale_fill_manual(name="species",values = c(col_values_7))+
    #scale_shape_manual(name="data type", values=c(""), labels=c(""))+
    #scale_shape_manual(name="", values=c(16), labels=c("observations"))+
    #scale_linetype_manual(name="", values=c("22"), labels=c("estimates"))+
    ggtitle("Estimates of Catches")+
    theme(legend.position = "bottom")
  C_plot_polite_legend
  
  catch_combined_extended$species.1<- factor(catch_combined_extended$species.1, levels=c("Dab","Brill","Turbot","Lemon sole", "Witch","Witch.1","Lemon sole.1"))  
  
  legend_species<-get_legend(C_plot_polite_legend)
  as_ggplot(legend_species)
  
  #--combine_legend-------
  
  
  legends<-ggdraw()+
    draw_plot(legend_est,                      y = 0.10, x = 0    , width = 1, height = 1)+
    draw_plot(legend_obs,                      y = 0.142, x = 0    , width = 1, height = 1)
  #draw_plot(legend_name,                      y = 0.18, x = 0    , width = 0.3, height = 1)+
  
  
  #legends
  
  #--make the plot of UBIC------------------------------
#  endplot_combined<-ggdraw() +
#    draw_plot(U_plot_polite_combi_alt, y = 0, x = 0.245, width = 0.245, height = 1) +
#    draw_plot(B_plot_polite_combi_alt, y= 0, x = 0, width = 0.245, height = 1)+
#    draw_label("Years",y = 0.28, x = 0.385, size=11 )+
#    draw_label("Years",y = 0.28, x = 0.14, size=11 )+
#    draw_plot(C_plot_polite_combi_ext, y = 0, x = 0.495,width = 0.245, height = 1)+
#    draw_plot(I_plot_polite_combi, y = 0, x = 0.745, width = 0.245, height = 1)+
#    draw_plot(legends,                      x = 0.08, y = -0.51    , width = 1, height = 1)+
#    #draw_label("Draft", colour = "#80404080", size = 120, angle = 45)
#  draw_plot(legend_name,                      x = -0.045, y = -.425    , width = 1, height = 1)+
#   draw_label("Years",y = 0.28, x = 0.64, size=11 )+
#    draw_plot(legend_species,                      x = -0.3, y = -.405    , width = 1, height = 1)

  
  #give new plot a name and save as object in environment----------
  assign(paste0("endplot_combi_corrmult_", i),endplot_combined)
}


#save.image("M:/copy of work thesis at 13.03.20/how to come to robin hood/saved trails/loop_results&plots_out_of_sd-reports_corrmult_sequence.RData")

#save.image("C:/Users/jappi/Bureaublad/robinhood/Code/Jasper final/rh_almost_good_to_go.RData")
#save.image("C:/Users/jappi/Bureaublad/robinhood/Code/Jasper final/all_data_for_figures.RData")

