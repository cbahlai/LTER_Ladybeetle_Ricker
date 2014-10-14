######################################################
# Code for analysis performed in "Regime shifts in Harmonia" 
# by Bahlai, van der Werf and Landis
######################################################


#read in raw data
datahaxyweeklyunculled<-read.csv(file="https://raw.githubusercontent.com/cbahlai/LTER_Ladybeetle_Ricker/master/harmoniaweekly.csv", header=TRUE, na.strings="")

#####################################
# Data extracted from relational database based on 
# http://lter.kbs.msu.edu/datatables/67
# Only  observations of H. axyridis extracted
# Extracted data summed over sampling week, total number
# traps collected over that week from all treatments 
# in original data set
#
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

#cut out last sampling year, because there is no Nt+1 for that year
datahaxy<-datahaxyraw[which(datahaxyraw$Year<max(datahaxyraw$Year)),]

#append Nt+1 column to dataset
datahaxy$Nt1<-Nt1

#run model selection on average captures per trap by year 
library(minpack.lm)

#models with no breaks
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=datahaxy)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=datahaxy)
cat("AIC for Logistic model with no breaks is ", AIC(logistic.model.1),  "\n")
cat("AIC for Ricker model with no breaks is ", AIC(ricker.model.1),  "\n")
summary(logistic.model.1)
summary(ricker.model.1)


#models with one break
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move break point one unit at a time
while (Break1< (max(datahaxy$Year))) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	#if statement to discount data subsets with 3 or fewer observations
	if(nrow(part1)>3 & nrow(part2)>3){
		#run models on two subsets
		logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
		ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
		logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)
		ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)
		#output combined AIC for each model type, noting the maximum year in the data subset 
		combAIC.log<-AIC(logistic.model.1)+AIC(logistic.model.2)
		cat("AIC for Logistic model with 1 break at ", max(part1$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(part2),  "\n")
		combAIC.rick<-AIC(ricker.model.1)+AIC(ricker.model.2)
		cat("AIC for Ricker model with 1 break at ", max(part1$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(part2),  "\n")
		if (logisticaicbest>combAIC.log){
			logisticaicbest=combAIC.log
			Break1best.l=Break1
			}
		if (rickeraicbest>combAIC.rick){
			rickeraicbest=combAIC.rick
			Break1best.r=Break1
			}
		}
	Break1<-Break1+1
	}
#compute model parameters for best models with one break
#best logistic model with one break


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
#Breakpoint for logistic is at:
max(part1$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)	
summary(logistic.model.1)
summary(logistic.model.2)

#best ricker model with one break

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
#Breakpoint for ricker is at:
max(part1$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part2)
summary(ricker.model.1)
summary(ricker.model.2)


#models with two breaks
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set initial value of Break 2 to 3 year later than Break 1
Break2<-Break1+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move first break point one unit at a time, stopping four years before end of study
while (Break1< (max(datahaxy$Year)-4)) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	Break2<-Break1+3
	#while loop to move second break point one unit at a time
	while (Break2< (max(datahaxy$Year))) {
		partA<-part2[which(part2$Year<Break2),]
		partB<-part2[which(part2$Year>(Break2-1)),]
		#if statement to discount data subsets with 3 or fewer observations
		if(nrow(part1)>3 & nrow(partA)>3 &nrow(partB)>3){
			#run models on three subsets
			logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
			ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
			logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			logistic.model.3<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
			ricker.model.3<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
			#output combined AIC for each model type, noting the maximum year in the data subset 
			combAIC.log<-AIC(logistic.model.1)+AIC(logistic.model.2)+AIC(logistic.model.3)
			cat("AIC for Logistic model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),  "\n")
			combAIC.rick<-(AIC(ricker.model.1)+AIC(ricker.model.2)+AIC(ricker.model.3))
			cat("AIC for Ricker model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),  "\n")
			if (logisticaicbest>combAIC.log){
				logisticaicbest=combAIC.log
				Break1best.l=Break1
				Break2best.l=Break2
				}
			if (rickeraicbest>combAIC.rick){
				rickeraicbest=combAIC.rick
				Break1best.r=Break1
				Break2best.r=Break2
				}
			}
		Break2<-Break2+1
		}
	Break1<-Break1+1
	}


#compute model parameters for best models with two breaks
#best logistic model with two breaks


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
partA<-part2[which(part2$Year<Break2best.l),]
partB<-part2[which(part2$Year>(Break2best.l-1)),]
#Breakpoint for logistic is at:
max(part1$Year)
max(partA$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partA)
logistic.model.3<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partB)	
summary(logistic.model.1)
summary(logistic.model.2)
summary(logistic.model.3)

#best ricker model with two breaks

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
partA<-part2[which(part2$Year<Break2best.r),]
partB<-part2[which(part2$Year>(Break2best.r-1)),]
#Breakpoint for ricker is at:
max(part1$Year)
max(partA$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
ricker.model.3<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
summary(ricker.model.1)
summary(ricker.model.2)
summary(ricker.model.3)






#models with two breaks constraining parts 1 and B to same constants
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set initial value of Break 2 to 3 year later than Break 1
Break2<-Break1+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move first break point one unit at a time, stopping four years before end of study
while (Break1< (max(datahaxy$Year)-4)) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	Break2<-Break1+3
	#while loop to move second break point one unit at a time
	while (Break2< (max(datahaxy$Year))) {
		partA<-part2[which(part2$Year<Break2),]
		partB<-part2[which(part2$Year>(Break2-1)),]
		part1B<-rbind(part1, partB)
		#if statement to discount data subsets with 3 or fewer observations
		if(nrow(part1)>3 & nrow(partA)>3 &nrow(partB)>3){
			#run models on three subsets
			logistic.model.A<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			ricker.model.A<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			logistic.model.1B<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
			ricker.model.1B<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
			#output combined AIC for each model type, noting the maximum year in the data subset 
			combAIC.log<-AIC(logistic.model.A)+AIC(logistic.model.1B)
			cat("AIC for constrained Logistic model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),".\n")
			combAIC.rick<-AIC(ricker.model.A)+AIC(ricker.model.1B)
			cat("AIC for constrained Ricker model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),"\n")
			if (logisticaicbest>combAIC.log){
				logisticaicbest=combAIC.log
				Break1best.l=Break1
				Break2best.l=Break2
				}
			if (rickeraicbest>combAIC.rick){
				rickeraicbest=combAIC.rick
				Break1best.r=Break1
				Break2best.r=Break2
				}
			}
		Break2<-Break2+1
		}
	Break1<-Break1+1
	}



#compute model parameters for best constrained models with two breaks
#best constrained logistic model with two breaks


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
partA<-part2[which(part2$Year<Break2best.l),]
partB<-part2[which(part2$Year>(Break2best.l-1)),]
part1B<-rbind(part1, partB)
#Breakpoint for logistic is at:
max(part1$Year)
max(partA$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1B)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partA)

summary(logistic.model.1)
summary(logistic.model.2)


#best constrained ricker model with two breaks

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
partA<-part2[which(part2$Year<Break2best.r),]
partB<-part2[which(part2$Year>(Break2best.r-1)),]
part1B<-rbind(part1, partB)
#Breakpoint for ricker is at:
max(part1$Year)
max(partA$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
summary(ricker.model.1)
summary(ricker.model.2)





#calculate area under curve by year- use that instead of trap averages
datahaxyweekly <- datahaxyweekly[order(datahaxyweekly$Year, datahaxyweekly$Ordinal_date),]
avct = c()
avtrap = c()
diffsample=c()
year=c()
ord_date=c()
for (i in 1:length(datahaxyweekly$Year)){
	year=c(year, datahaxyweekly$Year[i])
	ord_date=c(ord_date, datahaxyweekly$Ordinal_date[i])
	avct=c(avct, (datahaxyweekly$Captures[i]+datahaxyweekly$Captures[i+1])/2)
	avtrap=c(avtrap, (datahaxyweekly$Traps[i]+datahaxyweekly$Traps[i+1])/2)
	diffsample=c(diffsample, (datahaxyweekly$Ordinal_date[i+1]-datahaxyweekly$Ordinal_date[i]))
}

datahaxyweekly$avct<-avct
datahaxyweekly$avtrap<-avtrap
datahaxyweekly$diffsample<-diffsample

#remove sample dates at the end of the season, which provide nonsensical averages 
#(flagged with negative number of days between samples)
datahaxyweekly<-datahaxyweekly[which(datahaxyweekly$diffsample>0),]
#calculate area under curve for each week for both captures and traps
datahaxyweekly$countarea<-datahaxyweekly$avct*datahaxyweekly$diffsample
datahaxyweekly$traparea<-datahaxyweekly$avtrap*datahaxyweekly$diffsample


#reshape the weekly observations to provide yearly total capture area and trap area


datahaxymelt<-melt(datahaxyweekly, id=1:2)
datahaxyraw<-cast(datahaxymelt, Year~variable, sum)
#compute Nt and Nt+1 in ladybeetles per trap based on yearly total area under curve
datahaxyraw$Nt<-datahaxyraw$countarea/datahaxyraw$traparea
Nt1 = c()
for (i in 1:(length(datahaxyraw$Nt)-1)) {
  Nt1 = c(Nt1,datahaxyraw$Nt[i+1])
}
#cut out last sampling year, because there is no Nt+1 for that year
datahaxy<-datahaxyraw[which(datahaxyraw$Year<max(datahaxyraw$Year)),]
#append Nt+1 column to dataset
datahaxy$Nt1<-Nt1
#repeat model selection with area under curve data




#models with no breaks
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=datahaxy)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=datahaxy)
cat("AIC for Logistic model with no breaks is ", AIC(logistic.model.1),  "\n")
cat("AIC for Ricker model with no breaks is ", AIC(ricker.model.1),  "\n")
summary(logistic.model.1)
summary(ricker.model.1)
#models with one break
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move break point one unit at a time
while (Break1< (max(datahaxy$Year))) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	#if statement to discount data subsets with 3 or fewer observations
	if(nrow(part1)>3 & nrow(part2)>3){
		#run models on two subsets
		logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
		ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
		logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)
		ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)
		#output combined AIC for each model type, noting the maximum year in the data subset 
		combAIC.log<-AIC(logistic.model.1)+AIC(logistic.model.2)
		cat("AIC for Logistic model with 1 break at ", max(part1$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(part2),  "\n")
		combAIC.rick<-AIC(ricker.model.1)+AIC(ricker.model.2)
		cat("AIC for Ricker model with 1 break at ", max(part1$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(part2),  "\n")
		if (logisticaicbest>combAIC.log){
			logisticaicbest=combAIC.log
			Break1best.l=Break1
			}
		if (rickeraicbest>combAIC.rick){
			rickeraicbest=combAIC.rick
			Break1best.r=Break1
			}
		}
	Break1<-Break1+1
	}
#compute model parameters for best models with one break
#best logistic model with one break


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
#Breakpoint for logistic is at:
max(part1$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part2)	
summary(logistic.model.1)
summary(logistic.model.2)

#best ricker model with one break

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
#Breakpoint for ricker is at:
max(part1$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part2)
summary(ricker.model.1)
summary(ricker.model.2)


#models with two breaks
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set initial value of Break 2 to 3 year later than Break 1
Break2<-Break1+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move first break point one unit at a time, stopping four years before end of study
while (Break1< (max(datahaxy$Year)-4)) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	Break2<-Break1+3
	#while loop to move second break point one unit at a time
	while (Break2< (max(datahaxy$Year))) {
		partA<-part2[which(part2$Year<Break2),]
		partB<-part2[which(part2$Year>(Break2-1)),]
		#if statement to discount data subsets with 3 or fewer observations
		if(nrow(part1)>3 & nrow(partA)>3 &nrow(partB)>3){
			#run models on three subsets
			logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
			ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
			logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			logistic.model.3<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
			ricker.model.3<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
			#output combined AIC for each model type, noting the maximum year in the data subset 
			combAIC.log<-AIC(logistic.model.1)+AIC(logistic.model.2)+AIC(logistic.model.3)
			cat("AIC for Logistic model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),  "\n")
			combAIC.rick<-(AIC(ricker.model.1)+AIC(ricker.model.2)+AIC(ricker.model.3))
			cat("AIC for Ricker model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),  "\n")
			if (logisticaicbest>combAIC.log){
				logisticaicbest=combAIC.log
				Break1best.l=Break1
				Break2best.l=Break2
				}
			if (rickeraicbest>combAIC.rick){
				rickeraicbest=combAIC.rick
				Break1best.r=Break1
				Break2best.r=Break2
				}
			}
		Break2<-Break2+1
		}
	Break1<-Break1+1
	}


#compute model parameters for best models with two breaks
#best logistic model with two breaks


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
partA<-part2[which(part2$Year<Break2best.l),]
partB<-part2[which(part2$Year>(Break2best.l-1)),]
#Breakpoint for logistic is at:
max(part1$Year)
max(partA$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partA)
logistic.model.3<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partB)	
summary(logistic.model.1)
summary(logistic.model.2)
summary(logistic.model.3)

#best ricker model with two breaks

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
partA<-part2[which(part2$Year<Break2best.r),]
partB<-part2[which(part2$Year>(Break2best.r-1)),]
#Breakpoint for ricker is at:
max(part1$Year)
max(partA$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
ricker.model.3<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partB)
summary(ricker.model.1)
summary(ricker.model.2)
summary(ricker.model.3)






#models with two breaks constraining parts 1 and B to same constants
#prevent data subsetting of <3 years from start of study
Break1=min(datahaxy$Year)+3
#set initial value of Break 2 to 3 year later than Break 1
Break2<-Break1+3
#set an arbitrarily high AIC to start comparing model AIcs to and best break point variables
logisticaicbest=2000
rickeraicbest=2000
Break1best.l=0
Break2best.l=0
Break1best.r=0
Break2best.r=0
#while loop to move first break point one unit at a time, stopping four years before end of study
while (Break1< (max(datahaxy$Year)-4)) {
	#subset data at Break1
	part1<-datahaxy[which(datahaxy$Year<Break1),]
	part2<-datahaxy[which(datahaxy$Year>(Break1-1)),]
	Break2<-Break1+3
	#while loop to move second break point one unit at a time
	while (Break2< (max(datahaxy$Year))) {
		partA<-part2[which(part2$Year<Break2),]
		partB<-part2[which(part2$Year>(Break2-1)),]
		part1B<-rbind(part1, partB)
		#if statement to discount data subsets with 3 or fewer observations
		if(nrow(part1)>3 & nrow(partA)>3 &nrow(partB)>3){
			#run models on three subsets
			logistic.model.A<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			ricker.model.A<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
			logistic.model.1B<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
			ricker.model.1B<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
			#output combined AIC for each model type, noting the maximum year in the data subset 
			combAIC.log<-AIC(logistic.model.A)+AIC(logistic.model.1B)
			cat("AIC for constrained Logistic model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.log,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),".\n")
			combAIC.rick<-AIC(ricker.model.A)+AIC(ricker.model.1B)
			cat("AIC for constrained Ricker model with 2 breaks at ", max(part1$Year), " and ", max(partA$Year), " is ", combAIC.rick,", # data points:", nrow(part1), " , ", nrow(partA), " , ", nrow(partB),"\n")
			if (logisticaicbest>combAIC.log){
				logisticaicbest=combAIC.log
				Break1best.l=Break1
				Break2best.l=Break2
				}
			if (rickeraicbest>combAIC.rick){
				rickeraicbest=combAIC.rick
				Break1best.r=Break1
				Break2best.r=Break2
				}
			}
		Break2<-Break2+1
		}
	Break1<-Break1+1
	}



#compute model parameters for best constrained models with two breaks
#best constrained logistic model with two breaks


part1<-datahaxy[which(datahaxy$Year<Break1best.l),]
part2<-datahaxy[which(datahaxy$Year>Break1best.l-1),]
partA<-part2[which(part2$Year<Break2best.l),]
partB<-part2[which(part2$Year>(Break2best.l-1)),]
part1B<-rbind(part1, partB)
#Breakpoint for logistic is at:
max(part1$Year)
max(partA$Year)
logistic.model.1<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=part1B)
logistic.model.2<-nlsLM(Nt1~ Nt*(1+r*(1- Nt/k)), start=list(r=1.5, k=0.5), data=partA)

summary(logistic.model.1)
summary(logistic.model.2)


#best constrained ricker model with two breaks

part1<-datahaxy[which(datahaxy$Year<Break1best.r),]
part2<-datahaxy[which(datahaxy$Year>Break1best.r-1),]	
partA<-part2[which(part2$Year<Break2best.r),]
partB<-part2[which(part2$Year>(Break2best.r-1)),]
part1B<-rbind(part1, partB)
#Breakpoint for ricker is at:
max(part1$Year)
max(partA$Year)
ricker.model.1<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=part1B)
ricker.model.2<-nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=3, k=0.5), data=partA)
summary(ricker.model.1)
summary(ricker.model.2)
