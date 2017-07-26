## indices first is tur, second is brill (one more I obs than C obs per spec )
Im <-   read.csv("W://IMARES/IJmuiden/Afdeling/Projecten/data poor mixed fisheries/Robin Hood/data/Itur_bll.csv")
# Catches
Cm <- read.csv("W://IMARES/IJmuiden/Afdeling/Projecten/data poor mixed fisheries/Robin Hood/data/Cple_sol_tur_bll_lem.csv")
## F estimates from data rich assessments (ple, sol),same length as I )
Fm <- read.csv("W://IMARES/IJmuiden/Afdeling/Projecten/data poor mixed fisheries/Robin Hood/data/Fple_sol_cod.csv")

## print data for checks
merge(merge(Fm,Cm),Im, all=T)

#plot Fs
plot(x=Fm[,1], y=Fm[,2], ylim=c(0,1.3), xlab="Year", ylab="Fishing mortality (year-1)", type="b", col="black",bg="black", pch=21, yaxs="i", las=1, mgp=c(3,0.7,0), panel.first=grid() )
lines(x=Fm[,1], y=Fm[,3],  type="b", col="red", bg= "red", pch=22)
lines(x=Fm[,1], y=Fm[,4],  type="b", col="orange", bg= "orange", pch=0)
abline(h=0.21, col= "black", lty=2) #plaice Fmsy
abline(h=0.20, col= "red", lty=2) #sole Fmsy
abline(h=0.31, col= "orange", lty=2) #cod Fmsy

legend("topright",legend=c("European plaice", "Sole", "Cod"), col=c("black","red", "orange"), cex=0.8, pch=c(21,22,0), pt.bg=c("black","red", "white"), bty="n")
box()

#plot catches
options(scipen=8)
plot(x=Cm[,1],  y=(Cm[,2] + Cm[,3])/1000, pch=21, type="b", col="black",  bg="black", log="y", ylim=c(0.6,1800 ),  xlab="Year", ylab="Catch (1000 tonnes)", yaxs="i", las=1, mgp=c(3,0.7,0), panel.first=grid(equilogs=F) )

lines(x=Cm[,1], y=(Cm[,4] + Cm[,5])/1000, pch=22, type="b", col="red",    bg="red")
lines(x=Cm[,1], y=Cm[,6]/1000,            pch=23, type="b", col="blue",   bg="blue")
lines(x=Cm[,1], y=Cm[,7]/1000,            pch=24, type="b", col="magenta",bg="magenta")
lines(x=Cm[,1], y=(Cm[,8] + Cm[,9])/1000, pch=25, type="b", col="green",  bg="green")
legend("topleft",legend=c("European plaice", "Sole", "Turbot", "Brill", "Lemon sole"), col=c("black","red", "blue","magenta", "green"), pch=c(21,22,23,24, 25),pt.bg=c("black","red", "blue","magenta", "green"), bty="n", ncol=2)
box()

#plot indices
plot(x=Fm[,1],  y=Fm[,2], pch=21, type="b", col="black",  bg="black", log="y", ylim=c(0.1,1),  xlab="Year", ylab="I", yaxs="i", las=1, mgp=c(3,0.7,0), panel.first=grid(equilogs=F) )



