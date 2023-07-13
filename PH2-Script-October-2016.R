# Script for PH2 indicator calculation (Oct. 2016)

# The script is illustrated with the example of L4 Phyto Time-series (Chl a)

# The scv.file used is "PH2-PHYTO-L4-TOTAL-MARIE.csv", in general you should have your data exactly presented this way! (and the columns names should be EXACTLY similar to run the script without problems)
# Really important concerning your initial file: if you have several dates for a month, you should, in your exel file, allocate the same day to all the dates within a month. For ease, allocate the day "15" to all your dates.
# For instance, 04/05/2015 should be transformed into 15/05/2015
# Then save you file as a csv file wich contain two columns: first columns with dates (DATE) and second one with Chla 

#1 Set your directories where you have your file.
setwd("C:\\....)
setwd("C:\\Users\\Guy\\Desktop")


#Download packages (so install them first if they are not installed yet on your R)
library(boot)
library(pastecs) 

#2. Read your file
PH2L4=read.table("PH2-PHYTO-Exple.csv",sep=";", header=TRUE,na.string="-9")

# 3. calculate your monthly means 
date=PH2L4[,1] #rename the dates 
MONTHdata=by(PH2L4[,2],date,mean) # calculate your monthly means
monthdate=as.vector(names(MONTHdata))
monthCHLATOT=data.frame(matrix(unlist(MONTHdata),nrow=66,byrow=T),stringsAsFactors=FALSE) # For PHYTO, the number x put for "nrow=x" is the number of observations you have in your time-serie
data=as.matrix(monthCHLATOT)
date=as.matrix(monthdate)
datecor=as.Date(date,"%d/%m/%Y")# convert dates as text to R date format
DATA=cbind.data.frame(datecor,data) # put all your data into one data frame
DATA2=DATA[order(datecor),] # classify according to dates in a chronological order

# Now, you need to have a look at your times-serie. In order to continue the analysis, you will have to ensure
# that you do not have too many data are missing. Every year which have more than 4 months missing should be deleted.
#The time-series also needs to be consistent in time. The following analysis do not consider a time-serie with missing year. The fraction of the time-serie to use for the analysis 
#has thus to fullfill two conditions: no years missing within the time-series, no more than 4 months missing within a year. Also, the first year should start in January.
# If your first year with lot of data starts in february or march, it is a pitty to loose this year for the analysis. An idea is then to create "fake" data, based on the mean-month average for these months, based
#on the whole time-series


# If there are many years that need to be removed from the time-serie, you have to just take the part fo the time-serie of interest. For this you have to choose a part of yout time-series as folow:
DATA3<-DATA2[12:66,] # select the part of the time-serie which is relevant for the analysis, for instance here, you pick up data from 2009 instead of 2008

# Then after you have to use DATA3 instead of DATA2 in case you choose only a part of your time-serie
# Here, I already added a value at the start, and no year has less than 4 months missing. So we are using DATA2!!! (it is really important you understand well these steps!)

DATAOK=log10(DATA2[,2]+1) # normalisation of the data before start of analysis 

x=as.numeric(julian(DATA2[,1]))# convert dates as text to R date format
newdata=regul(x, DATAOK, units="daystoyears",frequency=12,
              methods="area",tol=3,n=100) # create a time series object from regul() output, you have to adapt the n to your data set which corresponds to the number of observations of your time-serie. You have to run it a first time with putting a higher number than n. 
# For instance, here n=66, so put n=100
newdata # you will see how many data are interpolated. Here, 6 data have been interpolated. So you run again the line, but with the exact number of lines you need, wether 66+6=72. It is necessary that there is no NA allocated at the end of he time-series
newdata=regul(x, DATAOK, units="daystoyears",frequency=12,
              methods="area",tol=3,n=72)
newdata  # to check that no NA have been padded at the end
rdata=tseries(newdata) # Now we have a regularized time-serie. You can just type "rdata" in R to vizualise your regularized time-serie and check everything looks fine.
# check that your time series is correct and that you don't have any hole in it!!! (otherwise you should not continue with the rest of the analysis, you have to fix the problem!)
rdata
      
#Start of the Decomposition
# decompose time series using seasonal difference method
# 1. compute mean annual cycle
ab=as.matrix(rdata)
mmean=rep(0,12) # monthly mean array
for (m in seq(1,12)){
aa=seq(m,nrow(ab),12)
mmean[m]=mean(ab[aa,])
}
      
# 2. remove seasonal cycle
nt=dim(ab)
nt=nt[1]
nyear=floor(nt/12)
repcycle=rep(mmean,nyear) # repeat mean cycle
repcycle=repcycle[1:nt]
ano=ab-repcycle # monthly anomalies
ano=as.ts(ano)
      
# save mensual anomalies
write.table(ano, file = "PH2PHYTO-anomonthlyL4.txt")
      
# 3. compute mean annual anomalies
anoyear=rep(0,nyear)
for (m in seq(1,nyear)){
aa=seq(12*(m-1)+1,12*m)
anoyear[m]=mean(ano[aa])
}
# convert series to time series1
repcycle=ts(repcycle,start=c(2008,1),frequency=12) # Here you have to adapt the start according to the first year of your time-serie of course
ano=ts(ano,start=c(2008,1),frequency=12) # same here, adapt to the starting year of your time-series
# remove singleton dimension of ano
ano=drop(ano)
anoyear=ts(anoyear,start=c(2008,1,1),frequency=1)
# save annual anomalies
write.table(anoyear, file = "PH2PHYTO-Example.txt") 

# Now, let's produce the graphs related to this analysis!
# Produce the grapg of the time-serie decomposition
x11()
par(mfrow=c(3,1))
plot(rdata,col=4)
title("Time series decomposition PH2-Phyto example") # for phyto
plot(repcycle,type="l",col=4)
plot(ano,col=4)
# You have to save the graph produced as a picture (jpeg or tiff) in the right file!
      
# Produce the grapg of mensual anomalies
windows()
year=2008:2014
plot(ano,type="h",xaxt="n",ylab="anomalies",col=4)
axis(1,year,las=1.9)
title("Phytoplankton PH2 Monthly anomalies Example") 
# You have to save the graph produced as a picture or pdf in the right file (in the following format: PH2-ANO-MENSUAL-PHYTO-L4)
      
# Produce the grapg of annual anomalies
windows()
year=2008:2013 
plot(anoyear,type="h",xaxt="n",ylab="anomalies",col=4)
axis(1,year,las=1.9)
title("PH2 Phyto Yearly anomalies example") 
# You have to save the graph produced as a picture in the right file (in the following format: PH2-ANO-ANUAL-PHYTO-L4)
      
Here we are, you have your graph of annual anomalies (and monthly ones as well)

# Create a data frame of your yearly anomalies
str(anoyear)
anoyearDF=as.data.frame(anoyear)
anoyearDF$Year= c(2008:2013) # you have to adapt this line with the years of your time-serie


# Create your final graph with the categorization
#We fix the limits based on quantiles
LCS = quantile(anoyearDF[,1],c(0.975)) # uppest superior limit
LCI = quantile(anoyearDF[,1],c(0.025)) # lowest inferior limit
LSS = quantile(anoyearDF[,1],c(0.75)) # 1st superior limit
LSI = quantile(anoyearDF[,1],c(0.25)) # 1st inferior limit
lim_sup = quantile(anoyearDF[,1],c(1))
lim_inf = quantile(anoyearDF[,1],c(0))
# Make the final result graph
x11()
plot(anoyearDF[,2],anoyearDF[,1] , ylim = c(lim_inf,lim_sup),xlim = c(2008, 2014),xaxt="n",las=2,xlab="",lab=c(8, 7, 2),pch=4,cex=1.2,lwd=3,ylab="yearly anomalies") # Of course, you have to change your xlim according to the start and the stop of your time-serie and adapt the numbers in the "lab" in order to fix correctly your x axis stickers
axis(1,xlim = c(2008, 2014),xlab="")
# To trace polygons zones
polygon(x=c(min(anoyearDF[,2])-20000,max(anoyearDF[,2])+20000,max(anoyearDF[,2])+20000,min(anoyearDF[,2])-20000),y=c(lim_inf,lim_inf,lim_sup,lim_sup),col="#6495ED",border=NA)
polygon(x=c(min(anoyearDF[,2])-20000,max(anoyearDF[,1])+20000,max(anoyearDF[,2])+20000,min(anoyearDF[,2])-20000),y=c(LCI,LCI,LCS,LCS),col="#87CEFA",border=NA)
polygon(x=c(min(anoyearDF[,2])-20000,max(anoyearDF[,2])+20000,max(anoyearDF[,2])+20000,min(anoyearDF[,2])-20000),y=c(LSI,LSI,LSS,LSS),col="#E0FFFF",border=NA)
points(anoyearDF[,2],anoyearDF[,1],pch=4,cex=1.2,lwd=2) # Replace the points above the polygons
# Add the lines which limit the different zones
abline(h=c(LSI,LSS), col="black", lty=3)
abline(h=c(LCI,LCS), col="black", lty=2)
col=c("#6495ED","#87CEFA","#E0FFFF")
title("PH2 plankton indice example") #. Of course adapt the name according to your data set name !!!!
legend("topright", xpd=TRUE,inset = c(-0.47,-0.16),cex=1,legend = c("extreme change (97.5%<x<2.5%)","important change (97.5%<x<75% and 25%<x<2.5%)","small change (25%<x<75%)"),fill=col,horiz=F,box.lty=0)# enlarge your graph in order to visualize the legend correctly (no interfering with the title)
# Keep the graph produced as a jpeg or TIFF picture with a name in this format: PH2ZOO-BOB-ALL.jpeg



