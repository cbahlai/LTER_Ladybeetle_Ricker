######################################################
# Code for figures generated in "Regime shifts in Harmonia" 
# by Bahlai, van der Werf and Landis
######################################################

#first, Harmonia population data
#pre-process data as it was used in analysis
#read in raw data
datahaxyweeklyunculled<-read.csv(file="C:/Rdata/harmoniaweekly.csv", header=TRUE, na.strings="")

#####################################
#     variables and descriptions
# Data has five columns
#
# Year	- year sample was taken
# Ordinal_date- day of year sample was taken
# Captures- Total number of Harmonia axyridis observed
# Traps	- Total number of traps reporting on that sampling date
# Per_trap- Average number of Harmonia per trap
######################################


library(reshape)
library(ggplot2)
library(gridExtra)
# for graphics, we will use the color palette GrandBudapest
# from Karthik Ram's wesanderson package
library(wesanderson)

#create dataset culled to a standard sampling date
cullpoint=241
datahaxyweekly<-datahaxyweeklyunculled[which(datahaxyweeklyunculled$Ordinal_date<cullpoint),]


#reshape the weekly observations to provide yearly total captures and traps
datahaxymelt<-melt(datahaxyweekly, id=1:2)
datahaxyraw<-cast(datahaxymelt, Year~variable, sum)

#compute Nt and Nt+1 in ladybeetles per trap based on yearly totals
datahaxyraw$Nt<-datahaxyraw$Captures/datahaxyraw$Traps
Nt1 = c()
for (i in 1:(length(datahaxyraw$Nt)-1)) {
  Nt1 = c(Nt1,datahaxyraw$Nt[i+1])
}

#assign sampling year a phase, based on output of model selection

phase = c()
for (i in 1:(length(datahaxyraw$Year))) {
	if(datahaxyraw$Year[i]<2001){
  		phase = c(phase, "A")
		}
	else if (datahaxyraw$Year[i]>2000& datahaxyraw$Year[i]<2006){
		phase = c(phase, "B")
		}
	else {
		phase = c(phase, "C")
		}
}
datahaxyraw$phase<-phase

#phases for the lines that are to join time series points
phasea = c()
for (i in 1:(length(datahaxyraw$Year))) {
  if(datahaxyraw$Year[i]<2001.1){
    phasea = c(phasea, "A")
  }
  else if (datahaxyraw$Year[i]>2000& datahaxyraw$Year[i]<2006.1){
    phasea = c(phasea, "B")
  }
  else if (datahaxyraw$Year[i]>2005.9){
    phasea = c(phasea, "C")
  }
}
datahaxyraw$phasea<-phasea

#cut out last sampling year, because there is no Nt+1 for that year
datahaxy<-datahaxyraw[which(datahaxyraw$Year<max(datahaxyraw$Year)),]

#append Nt+1 column to dataset
datahaxy$Nt1<-Nt1

######################################
#
# Generate time series figure
#
######################################
pal<-wes.palette(3, "GrandBudapest")

harmonia.timeseries<-ggplot(datahaxyraw, aes(Year, Nt, colour=phase, cex=1))+geom_point(size=4)+scale_color_manual(values = pal)+geom_line(data=datahaxyraw, aes(x=Year, y=Nt, group=phasea), size=1)+geom_line(size=1)+xlab("Year")+ylab("Average captures per trap")+theme_bw()+coord_equal(ratio=8)+geom_vline(xintercept=c(2000.9, 2005.9), colour="blue", linetype="longdash")+ theme(legend.key = element_blank())
harmonia.timeseries

######################################
#
# Generate Ricker model figure
#
######################################
phase.a<-function(x){x*exp(1.28*(1- x/0.33))}
phase.b<-function(x){x*exp(2.17*(1- x/0.47))}
phase.c<-function(x){x*exp(1.54*(1- x/0.30))}
phase.ac<-function(x){x*exp(1.47*(1- x/0.31))}



harmonia.ricker<-ggplot(datahaxy, aes(Nt, Nt1, colour=phase))+geom_point(size=4)+scale_color_manual(values = wes.palette(3, "GrandBudapest"))+xlab("N(t)")+ylab("N(t+1)")+theme_bw()+ theme(legend.key = element_blank())+stat_function(fun=phase.a, colour=pal[1], size=1)+stat_function(fun=phase.b, colour=pal[2], size=1)+stat_function(fun=phase.c, colour=pal[3], size=1)+stat_function(fun=phase.ac, colour="black", size=1, linetype="longdash")+coord_equal(ratio=1)
harmonia.ricker

######################################
#
# Pesticide by state data
#
######################################

#read in raw data
pesticide<-read.csv(file="C:/Rdata/midwest_soybean_pesticide_use.csv", header=TRUE, na.strings="")

#####################################
#     variables and descriptions
# Data has five columns
#
# Year	- year sample was taken
# State	- US state, where
#		IA= Iowa
#		IL= Illinois
#		MI= Michigan
#		WI= Wisconsin
# Compound- Insecticide active ingredient
# est_statewide_usage_kg- Statewide usage of this active 
#		ingredient in soybean, as estimated by
#		Stone et al of the United States Geological Survey
# soy_planted_ha- area of soybean planted in this state, as reported
#		by the United States Department of Agriculture
#		National Agricultural Statistics Service
#
######################################

pesticide$rate_ha<-pesticide$est_statewide_usage_kg/pesticide$soy_planted_ha

#remove states we were unable to get consistent scouting data for
pesticide<-pesticide[which(pesticide$State=="IA" | pesticide$State=="IL"| pesticide$State=="MI"| pesticide$State=="WI" ),]
#create data subsets for each active ingredient
cyhalo<-pesticide[which(pesticide$Compound=="CYHALOTHRINLAMBDA"),]
esfen<-pesticide[which(pesticide$Compound=="ESFENVALERATE"),]
imid<-pesticide[which(pesticide$Compound=="IMIDACLOPRID"),]
thiam<-pesticide[which(pesticide$Compound=="THIAMETHOXAM"),]

#choose colour and shape palettes
pal1<-c((wes.palette(5, "Zissou"))[c(1, 3, 5)],(wes.palette(4, "Rushmore"))[4])
shapepal<-c(15,17,19,8)

#create base figure for legend
pesticide.timeseries.leg<-ggplot(cyhalo, aes(Year, rate_ha, colour=State,shape=State))+geom_point(size=4)+scale_color_manual(values = pal1)+geom_line(size=1)+xlab("Year")+ylab("kg active ingredient per hectare")+theme_bw()+ theme(legend.key = element_blank())+scale_shape_manual(values = shapepal)

#plots for each individual active ingredient
cyhalo.timeseries<-ggplot(cyhalo, aes(Year, rate_ha, colour=State, shape=State))+geom_point(size=4)+scale_color_manual(values = pal1)+geom_line(size=1)+xlab(NULL)+ylab(NULL)+theme_bw()+ theme(legend.position="none")+scale_shape_manual(values = shapepal)+geom_vline(xintercept=2005.5,colour="blue", linetype="longdash")+ggtitle("Cyhalothrin-lambda")

esfen.timeseries<-ggplot(esfen, aes(Year, rate_ha, colour=State, shape=State))+geom_point(size=4) +scale_color_manual(values = pal1)+geom_line(size=1)+xlab(NULL)+ylab(NULL)+theme_bw()+ theme(legend.position="none")+scale_shape_manual(values = shapepal)+geom_vline(xintercept=2005.5,colour="blue", linetype="longdash")+ggtitle("Esfenvalerate")

imid.timeseries<-ggplot(imid, aes(Year, rate_ha, colour=State, shape=State))+geom_point(size=4)+scale_color_manual(values = pal1)+geom_line(size=1)+xlab(NULL)+ylab(NULL)+theme_bw()+ theme(legend.position="none")+scale_shape_manual(values = shapepal)+geom_vline(xintercept=2005.5,colour="blue", linetype="longdash")+ggtitle("Imidacloprid")

thiam.timeseries<-ggplot(thiam, aes(Year, rate_ha, colour=State, shape=State))+geom_point(size=4)+scale_color_manual(values = pal1)+geom_line(size=1)+xlab("Year")+ylab(NULL)+theme_bw()+ theme(legend.position="none")+scale_shape_manual(values = shapepal)+geom_vline(xintercept=2005.5,colour="blue", linetype="longdash")+ggtitle("Thiamethoxam")


######################################
#
# Stack pesticide plots together
#
######################################

#pull legend out of plot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

leg<-g_legend(pesticide.timeseries.leg)

#stick plots together in a vertical stack with a legend on the right

grid.arrange(arrangeGrob(arrangeGrob(cyhalo.timeseries, esfen.timeseries, imid.timeseries, thiam.timeseries, ncol=1),leg,ncol=2,widths=c(5/6,1/6), left="kg active ingredient per hectare") )

######################################
#
# boxplot of pesticide use by aphid infestation level
#
######################################
aphids<-read.csv(file="C:/Rdata/soybean_aphid_ratings_midwest.csv", header=TRUE, na.strings="")

aphidmelt<-melt(aphids, id=1)
names(aphidmelt)[names(aphidmelt)=="variable"]<-"Year"
names(aphidmelt)[names(aphidmelt)=="State."]<-"State"
aphidmelt$Year<-gsub("X","", aphidmelt$Year)
aphidmelt<-aphidmelt[which(aphidmelt$Year<2012),]

#merge infestation data with pesticide use data
aphid.pest<-merge(pesticide, aphidmelt, by=c("Year","State"))
names(aphid.pest)[names(aphid.pest)=="value"]<-"infestation"
#cull out data from before 2005, when pesticides weren't widely used
aphid.pest<-aphid.pest[which(aphid.pest$Year>2004),]


aphid.pest$infestation<-factor(aphid.pest$infestation, levels=c("Low","Spotty","Moderate","High"))


aphid.cyhalo<-aphid.pest[which(aphid.pest$Compound=="CYHALOTHRINLAMBDA"),]
aphid.esfen<-aphid.pest[which(aphid.pest$Compound=="ESFENVALERATE"),]
aphid.imid<-aphid.pest[which(aphid.pest$Compound=="IMIDACLOPRID"),]
aphid.thiam<-aphid.pest[which(aphid.pest$Compound=="THIAMETHOXAM"),]

#create new colour palatte
pal2<-c(wes.palette(4, "Royal2"))

#create plots of estimated pesticide use by aphid infestation
cyhalo.boxplot<-ggplot(aphid.cyhalo, aes(x=infestation, y=rate_ha, fill=Compound))+geom_boxplot()+xlab(NULL)+ylab(NULL)+ggtitle("Cyhalothrin-lambda")+theme_bw()+scale_fill_manual(values=pal2[1])+stat_summary(fun.y=median, geom="line", aes(group=1), cex=1, linetype="longdash")+ theme(legend.position="none")
esfen.boxplot<-ggplot(aphid.esfen, aes(x=infestation, y=rate_ha, fill=Compound))+geom_boxplot()+xlab(NULL)+ylab(NULL)+ggtitle("Esfenvalerate")+theme_bw()+scale_fill_manual(values=pal2[2])+stat_summary(fun.y=median, geom="line", aes(group=1), cex=1, linetype="longdash")+ theme(legend.position="none")
imid.boxplot<-ggplot(aphid.imid, aes(x=infestation, y=rate_ha, fill=Compound))+geom_boxplot()+xlab(NULL)+ylab(NULL)+ggtitle("Imidacloprid")+theme_bw()+scale_fill_manual(values=pal2[3])+stat_summary(fun.y=median, geom="line", aes(group=1), cex=1, linetype="longdash")+ theme(legend.position="none")
thiam.boxplot<-ggplot(aphid.thiam, aes(x=infestation, y=rate_ha, fill=Compound))+geom_boxplot()+xlab(NULL)+ylab(NULL)+ggtitle("Thiamethoxam")+theme_bw()+scale_fill_manual(values=pal2[4]) +stat_summary(fun.y=median, geom="line", aes(group=1), cex=1, linetype="longdash")+ theme(legend.position="none")
#stack plots together

grid.arrange(arrangeGrob(cyhalo.boxplot, esfen.boxplot, imid.boxplot, thiam.boxplot, ncol=2), left="kg active ingredient per hectare", sub="Aphid infestation")

