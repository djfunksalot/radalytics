# returns NA for logs without throwing a warning even if you give it a negative number 
 
"log10Nice"<-
function(myList)
{
    logList<-NULL
        for (val in myList) {
            if (val > 0) {
                logged<-log10(val)
            } else {
                logged<-NA
            }
            logList<-c(logList,logged)
        }
    return(logList)
}

"formatStandards"<-
function(data,standards) 
{
 i<-0;(while (i<len) {i<-i+1;if(data[i,]['type']=='standard') {print(data[i,]['name'])}})
}


"applyStandards"<-
function(data,standards) 
{
	i<-0
	while (i<nrow(standards)) {
		i<-i+1
		standard<-standards[i,]
		data[data$type == 'standard' & data$name == standard$name,]$conc<-standard$value
        }
	return(data)
}


##################################################START 4-PL FUNCTION#########################################


#  Defines the function with the following variables
#  run - this is the output datafile from the plate reader via SPS
#  ul - the upper limit of the assay as specified by the manufacturer
#  hlimL - horizonal limit low. This is the CV threshold for unknowns that have a concentration above the vertical threshold (vlim).
#  hlimH - horizontal limit high. This is the CV threshold for unknowns that have a concentration below the vertical threshold and above the limit of detection
#  vlim - vertical limit. This is selected based on a number of semi-arbitrary criteria. Examples; the second lowest calibrator, the upper limit divided by 10, or just by eyeballing the concentration vs. CV plot 

"fourpl"<-
function(data,ul,hlimL,hlimH,vlim,name='X') {
	# calls package 'dose response curve' which contains the drm() function used for non-linear regression
	require('drc')
	#check to see if we are running on a server
	is_server<-Sys.getenv("is_server", unset = FALSE)

	#split data array into samples and standards
	samples<-rbind(split(data,data$type)$sample,split(data,data$type)$control)
	standards<-split(data,data$type)$standard

	#create separate array for background standards
	#bgs<-standards[standards$conc==0,]

	#average background od
	#bg_mean<-mean(bgs$od)

	#remove background from rest of standards
	#standards<-standards[standards$conc!=0,]

	#Creates a scaling factor which will set the lowest calibator value to 1 in the case that it is less than one.
	#This is done to remove negative log vales which cannot be used as starting parameters in the drm() function later on.       
	#low<-ifelse(min(standards$conc)<1,low<-min(standards$conc),low<-1)

    
	#correct for background
	# odb = optical density background (corrected) not the rapper
	# samples$odb<-samples$od-bg_mean
	# standards$odb<-standards$od-bg_mean
	# bgs$odb<-bgs$od-bg_mean                                                                            

	#log transform of the calibrator dataset.
  	#standards$logodb<-log10Nice(standards$odb)
	#logconc<-log10(standards$conc/low)

	#log transform of backgrounds
	#bgs$logodb<-log10Nice(bgs$odb)
	#bgs$logconc<-log10Nice(bgs$conc)

	#log transform of unknown OD concs
	#samples$logodb<-log10Nice(samples$odb)

	#instantiate plot file
	if (is_server) {
		jpeg('rplot.jpg',width = 960, height = 960, units = "px", quality = 100)
	}

	#Dose response function which fits a 4pl model to the calibrator dataset
	fit<-drm(data.frame(standards$od,standards$conc),fct=LL.4())

#set IDs for the coefficients calculated by the 4-pl fit.
	b=hill_slope=round(fit$coef[1],6)
	cc=min_asymptote=round(fit$coef[2],6)
	d=max_asymptote=round(fit$coef[3],6)
	e=inflection_point=round(fit$coef[4],6)                            

	#Calculates a linear r^2 conc, rounded to six decimal places. NOTE this is only an estimated r^2 conc since a 4-pl curve fit doesn't normally contain an r^2 conc.
	rsquared<-round(cor(standards$od,standards$conc),6)
                      
	#Calculates the corresponding concentration from the 4pl equation. 
	samples$dilution[samples$type=='control']<-6
	samples$dilution[samples$type=='standard']<-6
	samples$calc_conc=e*(((-d+samples$od)/(cc-samples$od))^(1/b))
	samples$calc_conc_dil=samples$calc_conc*samples$dilution
	standards$calc_conc=e*(((-d+standards$od)/(cc-standards$od))^(1/b))
	standards$calc_conc_dil=standards$calc_conc
  
  
	#Removes NaN values from the logconc columns and replaces them with NA  
	samples$calc_conc[is.na(samples$calc_conc)]<-NA
	standards$calc_conc[is.na(standards$calc_conc)]<-NA
  
	#Creates a vector which contains the sample log concentration values which do not equal NaN  
	p2<-samples$calc_conc[!is.na(samples$calc_conc)]
  
  
#Set any remaining NaN concs to a concentration of 0
	#samples$calc_conc[is.na(samples$calc_conc)]<-0

 #creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.                      
	min<-min(p2)
	max<-max(p2)


#  x <-seq(0,max,length=10000)
#	y <-(cc+((d-cc)/(1+exp(b*log(x)-b*log(e)))))
#
#Creates a dataset of the log blank OD concs and the calculated log concentration
	unknowns<-data.frame(abs=c(samples$od),conc=c(samples$calc_conc))

#partitions the graph into 1 row and 2 columns
#  par(mfrow=c(2,1))

#set the title for the graph
  main=paste("Assay ",name,":  units")

#plots the log10 concentration vs. log10 OD concs for the calibrators. The x domain is set to include +/- 1 dilution for the calibrator to catch above and below range concs.
	plot(standards$calc_conc, main=main, standards$od, xlim=c(0,max), xlab="Concentration", ylab="Absorbance",pch=16,col="red",cex=1.5)
   
#plots the best-fit line                       
#	lines(x,y, lty="dotted", col="red")

#Superimposes the points for controls and unknowns on the best-fit line. Unknowns are blue and controls are green. The color has some transparency so that overlapping points can be visualized
#	lines(unknowns$calc_conc,unknowns$abs, type="points", pch=16,cex=1.5)

#Lays a background grid on the plot
	grid(nx=NULL,ny=NULL,col="gray")

#Creates two legends
	legend("topleft",title="4-PL Curve Fit",cex=0.75,pch=16,col=c("red","green","blue"),legend=c("Standard","Control","Unknown"),ncol=1)
	legend("bottomright",cex=0.75,pch=16,c(paste("r^2=",rsquared),paste("Min Asympt=",cc),paste("Max Asympt=",d),paste("Curvature=",b),paste("Inflection=",e)),ncol=1)
 
#Merges the standards, background and samples back into one dataset.           
  combined<-rbind(standards,samples)

#splits consecutive duplicate well values into two datasets, odd and even. This is performed so that paired well values may be aligned onto the same row.
  odd<-split(combined,combined$id%%2==0)$'FALSE'
  even<-split(combined,combined$id%%2==0)$'TRUE'

#creates a table containing sample type, name and the paired well values for each specimen
  pair<-data.frame(type=odd$type,name=odd$name,well_1=odd$calc_conc,well_2=even$calc_conc)

#calculates the mean paired well value and adds a column
  pair$conc_mean<-rowMeans(pair[,3:4],na.rm=FALSE,dim=1)
              
#calculates the paired well CV% and adds a column              
  pair$conc_cv<-(apply(pair[,3:4],1,sd)/pair$conc_mean)*100
  pair$conc_cv[is.na(pair$conc_cv)]<-NA

#generates a plot of the paired well concentration vs. paired well CV% with corresponding thresholds. Specimens that fail the criteria are colored "RED", those which pass are colored "BLUE".
  plot(pair$conc_mean,pair$conc_cv,pch=16,cex=0.7,main="Assay X: Concentration vs. CV Plot",xlab="units",ylab="CV%",col=ifelse((pair$conc_mean>vlim)&(pair$conc_cv>hlimL),"red",ifelse(pair$conc_cv>hlimH,"red","blue")),xlim=c(0,ul))

#thresholds
  abline(h=hlimL)
  abline(h=hlimH)
  abline(v=vlim)

#lays a grid in the background
  grid(nx=NULL,ny=NULL,col="gray")

#creates a column with the results of whether the specimen should be considered for retest based on the criteria. TRUE=review and potentially restest, FALSE=pass
  pair$retest<-ifelse((pair$conc_mean>vlim&pair$conc_cv>hlimL|pair$conc_cv>hlimH),"TRUE","FALSE")
  pair$retest[is.na(pair$retest)]<-FALSE

#assigns paired mean, paired CV and retest result columns to the odd and even data sets generated earlier
  odd$conc_mean<-pair$conc_mean
  odd$conc_cv<-pair$conc_cv
  odd$retest<-pair$retest
  even$conc_mean<-pair$conc_mean
  even$conc_cv<-pair$conc_cv
  even$retest<-pair$retest

#combines the odd and even datasets
  results<-rbind(odd,even)

#sorts results back into the original well position order
results<-results[with(results,order(results$id)),]

#adds TRUE or FALSE to whether the mean concentration is out of range and needs to be retested with a dilution
results$dilution<-ifelse(((results$conc_mean>ul)&(results$type=="sample")),"TRUE","FALSE")


   if (is_server) {
       dev.off()
   }                                               

#results<-results[with(results,order(results$name)),]
curve<-list(rsquared=rsquared,cc=cc,d=d,b=b,e=e,min=min,max=max)
return(list(curve=curve,results=results,standards=standards,pair=pair))
          
}

#########################################################END OF 4-PL FUNCTION######################################################################

##################################################START 4-PL FUNCTION#########################################


#  Defines the function with the following variables
#  run - this is the output datafile from the plate reader via SPS
#  ul - the upper limit of the assay as specified by the manufacturer
#  hlimL - horizonal limit low. This is the CV threshold for unknowns that have a concentration above the vertical threshold (vlim).
#  hlimH - horizontal limit high. This is the CV threshold for unknowns that have a concentration below the vertical threshold and above the limit of detection
#  vlim - vertical limit. This is selected based on a number of semi-arbitrary criteria. Examples; the second lowest calibrator, the upper limit divided by 10, or just by eyeballing the concentration vs. CV plot 

"fourplloglog"<-
function(data,ul,hlimL,hlimH,vlim,name='X')
{ 
# calls package 'dose response curve' which contains the drm() function used for non-linear regression
    #library('drc')
    require('drc')

#check to see if we are running on a server
    is_server<-Sys.getenv("is_server", unset = FALSE)

#split data array into samples and standards
    samples<-rbind(split(data,data$type)$sample,split(data,data$type)$control)
    standards<-split(data,data$type)$standard

#create separate array for background standards
    bgs<-standards[standards$conc==0,]

#average background od
    bg_mean<-mean(bgs$od)

#remove background from rest of standards
    standards<-standards[standards$conc!=0,]

#Creates a scaling factor which will set the lowest calibator value to 1 in the case that it is less than one.
#This is done to remove negative log vales which cannot be used as starting parameters in the drm() function later on.       
    low<-ifelse(min(standards$conc)<1,low<-min(standards$conc),low<-1)

    
#correct for background
# odb = optical density background (corrected) not the rapper
	 samples$odb<-samples$od-bg_mean
	 standards$odb<-standards$od-bg_mean
	 bgs$odb<-bgs$od-bg_mean                                                                            

#log transform of the calibrator dataset.
  standards$logodb<-log10Nice(standards$odb)
  logconc<-log10(standards$conc/low)

#log transform of backgrounds
	bgs$logodb<-log10Nice(bgs$odb)
	bgs$logconc<-log10Nice(bgs$conc)

#log transform of unknown OD concs
	samples$logodb<-log10Nice(samples$odb)

#instantiate plot file
  if (is_server) {
       jpeg('rplot.jpg',width = 960, height = 960, units = "px", quality = 100)
     }

#Dose response function which fits a 4pl model to the calibrator dataset
	fit<-drm(data.frame(standards$logodb,logconc),fct=LL.4())

#set IDs for the coefficients calculated by the 4-pl fit.
	b=hill_slope=round(fit$coef[1],6)
	cc=min_asymptote=round(fit$coef[2],6)
	d=max_asymptote=round(fit$coef[3],6)
	e=inflection_point=round(fit$coef[4],6)                            

#Calculates a linear r^2 conc, rounded to six decimal places. NOTE this is only an estimated r^2 conc since a 4-pl curve fit doesn't normally contain an r^2 conc.
	rsquared<-round(cor(standards$logodb,logconc),6)
                      
#creates a subset of the log OD concs which contain only those concs which do not contain NaN concs  
	p1<-samples$logodb[!is.na(samples$logodb)]
                                      
#sets all NaN concs equal to the minimum conc contained within both columns of OD concs which doesn't contain an NaN conc.
#NOTE: does this included standards and samples?
	samples$logodb[is.na(samples$logodb)]<-min(p1)

#Calculates the corresponding log concentration from the 4pl equation. 
	samples$logconc=e*(((-d+samples$logodb)/(cc-samples$logodb))^(1/b))
	standards$logconc=e*(((-d+standards$logodb)/(cc-standards$logodb))^(1/b))
  
  
#Removes NaN values from the logconc columns and replaces them with NA  
  samples$logconc[is.na(samples$logconc)]<-NA
  standards$logconc[is.na(standards$logconc)]<-NA
  
#Creates a vector which contains the sample log concentration values which do not equal NaN  
  p2<-samples$logconc[!is.na(samples$logconc)]
  
#Converts log concentration concs back into standard concentration
	samples$calc_conc<-(10^(samples$logconc))*low
  standards$calc_conc<-(10^(standards$logconc))*low
  bgs$calc_conc<-bgs$conc
  
#Set any remaining NaN concs to a concentration of 0
	samples$calc_conc[is.na(samples$calc_conc)]<-0

 #creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.                      
 min<-min(p2)
 max<-max(p2)


  x <-seq(min,max,length=10000)
	y <-(cc+((d-cc)/(1+exp(b*log(x)-b*log(e)))))

#Creates a dataset of the log blank OD concs and the calculated log concentration
	unknowns<-data.frame(abs=c(samples$od),conc=c(samples$calc_conc/samples$dilution))

#partitions the graph into 1 row and 2 columns
  par(mfrow=c(2,1))

#set the title for the graph
  main=paste("Assay ",name,": log-log, units")

#plots the log10 concentration vs. log10 OD concs for the calibrators. The x domain is set to include +/- 1 dilution for the calibrator to catch above and below range concs.
	plot(standards$calc_conc, main=main, standards$od, xlim=c(0,max), xlab="Concentration", ylab="Absorbance",pch=16,col="red",cex=1.5)
   
#plots the best-fit line                       
	lines(x,y, lty="dotted", col="red")

#Superimposes the points for controls and unknowns on the best-fit line. Unknowns are blue and controls are green. The color has some transparency so that overlapping points can be visualized
	lines(unknowns$conc,unknowns$abs, type="points", pch=16,cex=1.5)

#Lays a background grid on the plot
	grid(nx=NULL,ny=NULL,col="gray")

#Creates two legends
	legend("topleft",title="4-PL Curve Fit",cex=0.75,pch=16,col=c("red","green","blue"),legend=c("Standard","Control","Unknown"),ncol=1)
	legend("bottomright",cex=0.75,pch=16,c(paste("r^2=",rsquared),paste("Min Asympt=",cc),paste("Max Asympt=",d),paste("Curvature=",b),paste("Inflection=",e)),ncol=1)
 
#Merges the standards, background and samples back into one dataset.           
  combined<-rbind(standards,samples)

#splits consecutive duplicate well values into two datasets, odd and even. This is performed so that paired well values may be aligned onto the same row.
  odd<-split(combined,combined$id%%2==0)$'FALSE'
  even<-split(combined,combined$id%%2==0)$'TRUE'

#creates a table containing sample type, name and the paired well values for each specimen
  pair<-data.frame(type=odd$type,name=odd$name,well_1=odd$calc_conc,well_2=even$calc_conc)

#calculates the mean paired well value and adds a column
  pair$conc_mean<-rowMeans(pair[,3:4],na.rm=FALSE,dim=1)
              
#calculates the paired well CV% and adds a column              
  pair$conc_cv<-(apply(pair[,3:4],1,sd)/pair$conc_mean)*100
  pair$conc_cv[is.na(pair$conc_cv)]<-NA

#generates a plot of the paired well concentration vs. paired well CV% with corresponding thresholds. Specimens that fail the criteria are colored "RED", those which pass are colored "BLUE".
  plot(pair$conc_mean,pair$conc_cv,pch=16,cex=0.7,main="Assay X: Concentration vs. CV Plot",xlab="units",ylab="CV%",col=ifelse((pair$conc_mean>vlim)&(pair$conc_cv>hlimL),"red",ifelse(pair$conc_cv>hlimH,"red","blue")),xlim=c(0,ul))

#thresholds
  abline(h=hlimL)
  abline(h=hlimH)
  abline(v=vlim)

#lays a grid in the background
  grid(nx=NULL,ny=NULL,col="gray")

#creates a column with the results of whether the specimen should be considered for retest based on the criteria. TRUE=review and potentially restest, FALSE=pass
  pair$retest<-ifelse((pair$conc_mean>vlim&pair$conc_cv>hlimL|pair$conc_cv>hlimH),"TRUE","FALSE")
  pair$retest[is.na(pair$retest)]<-FALSE

#assigns paired mean, paired CV and retest result columns to the odd and even data sets generated earlier
  odd$conc_mean<-pair$conc_mean
  odd$conc_cv<-pair$conc_cv
  odd$retest<-pair$retest
  even$conc_mean<-pair$conc_mean
  even$conc_cv<-pair$conc_cv
  even$retest<-pair$retest

#combines the odd and even datasets
  results<-rbind(odd,even)

#sorts results back into the original well position order
results<-results[with(results,order(results$id)),]

#adds TRUE or FALSE to whether the mean concentration is out of range and needs to be retested with a dilution
results$dilution<-ifelse(((results$conc_mean>ul)&(results$type=="sample")),"TRUE","FALSE")


   if (is_server) {
       dev.off()
   }                                               

curve<-list(rsquared=rsquared,cc=cc,d=d,b=b,e=e,low=low,min=min,max=max)
return(list(curve=curve,results=results,standards=standards))
          
}

#########################################################END OF 4-PL FUNCTION######################################################################





########################################################START LINEAR FUNCTION########################################################################

#Defines the function with the following variables
# run - this is the output datafile from the plate reader via SPS
# ul - the upper limit of the assay as specified by the manufacturer
# hlimL - horizonal limit low. This is the CV threshold for unknowns that have a concentration above the vertical threshold (vlim).
# hlimH - horizontal limit high. This is the CV threshold for unknowns that have a concentration below the vertical threshold and above the limit of detection
# vlim - vertical limit. This is selected based on a number of semi-arbitrary criteria. Examples; the second lowest calibrator, the upper limit divided by 10, or just by eyeballing the concentration vs. CV plot 


"linearfit"<-
function(data,ul,hlimL,hlimH,vlim,name='X')
{

#check to see if we are running on a server
  is_server<-Sys.getenv("is_server", unset = FALSE)

#split data array into samples, standards and bg
        samples<-rbind(split(data,data$type)$sample,split(data,data$type)$control)
	standards<-split(data,data$type)$standard

#create seperate array for background standards
	bgs<-standards[standards$conc==0,]

#average background od
	bg_mean<-mean(bgs$od)

#remove background from rest of standards
	standards<-standards[standards$conc!=0,]

#correct for background
#odb = optical density bacground (corrected) not the rapper
	samples$odb<-samples$od-bg_mean
	standards$odb<-standards$od-bg_mean
	bgs$odb<-bgs$od-bg_mean

#log transform of the calibrator dataset.
	standards$logodb<-log10Nice(standards$odb)
	standards$logconc<-log10Nice(standards$conc)

#log transform of backgrounds
	bgs$logodb<-log10Nice(bgs$odb)
	bgs$logconc<-log10Nice(bgs$conc)

#log transform of unknown OD concs
	samples$logodb<-log10Nice(samples$odb)

	
#instantiate plot file
  if (is_server) {
    jpeg('rplot.jpg',width = 960, height = 960, units = "px", quality = 100)
    }

#Dose response function which fits a 4pl model to the calibrator dataset
	fit<-lm(data.frame(standards$logodb,standards$logconc))

#assigns letter IDs to the coefficients calculated by the linear fit.
  a=slope=round(fit$coef[2],6)
  b=yint=round(fit$coef[1],6)
                     
#creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.                      
  x <-seq(log2(ul/256),log2(ul*2), length=1000)
  y <-(fit$coef[2]*x)+fit$coef[1]

#partitions the graph into 1 row and 2 columns
  par(mfrow=c(2,1))

#set the title for the graph
  main=paste("Assay ",name,": log-log, units")

#plots the log10 concentration vs. log10 OD concs for the calibrators. The x domain is set to include +/- 1 dilution for the calibrator to catch above and below range concs.
	plot(standards$logconc, standards$logodb, main="Assay X: log-log, units", xlim=c(log10(ul/256),log10(ul*2)), xlab="x=log10(Concentration) units", ylab="y=log10(Absorbance)",pch=16,col="red",cex=1.5)
   
#plots the best-fit line                       
	lines(x,y, lty="dotted", col="red")

#Calculates a linear r^2 conc, rounded to six decimal places. NOTE this is only an estimated r^2 conc since a 4-pl curve fit doesn't normally contain an r^2 conc.
	rsquared<-round(cor(standards$logodb,standards$logconc),6)
                       
#creates a subset of the log OD concs which contain only those concs which do not contain NaN concs  
	p1<-samples$logodb[!is.na(samples$logodb)]

#sets all NaN concs equal to the minimum conc contained within both columns of OD concs which doesn't contain an NaN conc.
#NOTE: does this included standards and samples?

	samples$logodb[is.na(samples$logodb)]<-min(p1)

#Calculates the corresponding log concentration from the 4pl equation. 
	samples$logconc=(samples$logodb-b)/a 
  standards$logconc=(standards$logodb-b)/a
  
#Converts log concentration concs back into standard concentration
	samples$calc_conc<-10^(samples$logconc)
  standards$calc_conc<-10^(standards$logconc)
  bgs$calc_conc<-bgs$conc
  
#Set any remaining NaN concs to a concentration of 0
	samples$conc[is.na(samples$conc)]<-0
	samples$logconc[is.na(samples$logconc)]<-0

#Creates a dataset of the log blank OD concs and the calculated log concentration
	unknowns<-data.frame(abs=c(samples$logodb),conc=c(samples$logconc))

#Superimposes the points for controls and unknowns on the best-fit line. Unknowns are blue and controls are green. The color has some transparency so that overlapping points can be visualized

	lines(unknowns$conc,unknowns$abs, type="points", pch=16,cex=1.5)

#Lays a background grid on the plot
	grid(nx=NULL,ny=NULL,col="gray")

#Creates two legends
	legend("topleft",title="Linear Curve Fit",cex=0.75,pch=16,col=c("red","green","blue"),legend=c("Standard","Control","Unknown"),ncol=1)
	legend("bottomright",cex=0.75,pch=16,c(paste("r^2=",rsquared),paste("Y-int=",b),paste("Slope",a)),ncol=1)

#merges the standards, background and samples back into one dataset.           
  combined<-rbind(bgs,standards,samples)

#splits consecutive duplicate well values into two datasets, odd and even. This is performed so that paired well values may be aligned onto the same row.
  odd<-split(combined,combined$id%%2==0)$'FALSE'
  even<-split(combined,combined$id%%2==0)$'TRUE'

#creates a table containing sample type, name and the paired well values for each specimen
  pair<-data.frame(type=odd$type,name=odd$name,well_1=odd$calc_conc,well_2=even$calc_conc)

#calculates the mean paired well value and adds a column
  pair$conc_mean<-rowMeans(pair[,3:4],na.rm=FALSE,dim=1)
              
#calculates the paired well CV% and adds a column              
  pair$conc_cv<-(apply(pair[,3:4],1,sd)/pair$conc_mean)*100

#generates a plot of the paired well concentration vs. paired well CV% with corresponding thresholds. Specimens that fail the criteria are colored "RED", those which pass are colored "BLUE".
  plot(pair$conc_mean,pair$conc_cv,pch=16,cex=0.7,main="Assay X: Concentration vs. CV Plot",xlab="units",ylab="CV%",col=ifelse((pair$conc_mean>vlim)&(pair$conc_cv>hlimL),"red",ifelse(pair$conc_cv>hlimH,"red","blue")),xlim=c(0,ul))

#thresholds
  abline(h=hlimL)
  abline(h=hlimH)
  abline(v=vlim)

#lays a grid in the background
  grid(nx=NULL,ny=NULL,col="gray")

#creates a column with the results of whether the specimen should be considered for retest based on the criteria. TRUE=review and potentially restest, FALSE=pass
  pair$retest<-ifelse((pair$conc_mean>vlim&pair$conc_cv>hlimL|pair$conc_cv>hlimH),"TRUE","FALSE")

#assigns paired mean, paired CV and retest result columns to the odd and even data sets generated earlier
  odd$conc_mean<-pair$conc_mean
  odd$conc_cv<-pair$conc_cv
  odd$retest<-pair$retest
  even$conc_mean<-pair$conc_mean
  even$conc_cv<-pair$conc_cv
  even$retest<-pair$retest

#combines the odd and even datasets
  results<-rbind(odd,even)

#sorts results back into the original well position order
results<-results[with(results,order(results$id)),]

#adds TRUE or FALSE to whether the mean concentration is out of range and needs to be retested with a dilution
results$dilution<-results$conc_mean>250


  if (is_server) {
    dev.off()
  }                                               

#presents table with updated results
legend<-c('NA')
return(list(legend,results))
          
}


######################################################END OF LINEAR FUNCTION###################################################





######################################################START OF 5-PL FUNCTION############################################################


#  Defines the function with the following variables
#  run - this is the output datafile from the plate reader via SPS
#  ul - the upper limit of the assay as specified by the manufacturer
#  hlimL - horizonal limit low. This is the CV threshold for unknowns that have a concentration above the vertical threshold (vlim).
#  hlimH - horizontal limit high. This is the CV threshold for unknowns that have a concentration below the vertical threshold and above the limit of detection
#  vlim - vertical limit. This is selected based on a number of semi-arbitrary criteria. Examples; the second lowest calibrator, the upper limit divided by 10, or just by eyeballing the concentration vs. CV plot 

"fiveplfit"<-
function(data,ul,hlimL,hlimH,vlim,name='X')
{ 
# calls package 'dose response curve' which contains the drm() function used for non-linear regression
    #library('drc')
    require('drc')

#check to see if we are running on a server
    is_server<-Sys.getenv("is_server", unset = FALSE)


#split data array into samples, standards and bg
  samples<-rbind(split(data,data$type)$sample,split(data,data$type)$control)
  standards<-split(data,data$type)$standard


#create separate array for background standards
    bgs<-standards[standards$conc==0,]

#average background od
    bg_mean<-mean(bgs$od)

#remove background from rest of standards
    standards<-standards[standards$conc!=0,]

#Creates a scaling factor which will set the lowest calibator value to 1 in the case that it is less than one.
#This is done to remove negative log vales which cannot be used as starting parameters in the drm() function later on.       
    low<-ifelse(min(standards$conc)<1,low<-min(standards$conc),low<-1)

    
#correct for background
# odb = optical density background (corrected) not the rapper
	 samples$odb<-samples$od-bg_mean
	 standards$odb<-standards$od-bg_mean
	 bgs$odb<-bgs$od-bg_mean                                                                            

#log transform of the calibrator dataset.
  standards$logodb<-log10Nice(standards$odb)
  logconc<-log10(standards$conc/low)

#log transform of backgrounds
	bgs$logodb<-log10Nice(bgs$odb)
	bgs$logconc<-log10Nice(bgs$conc)

#log transform of unknown OD concs
	samples$logodb<-log10Nice(samples$odb)

#instantiate plot file
  if (is_server) {
     jpeg('rplot.jpg',width = 960, height = 960, units = "px", quality = 100)
     }

#Dose response function which fits a 4pl model to the calibrator dataset
	fit<-drm(data.frame(standards$logodb,logconc),fct=LL.5())

#set IDs for the coefficients calculated by the 4-pl fit.
	b=hill_slope=round(fit$coef[1],6)
	cc=min_asymptote=round(fit$coef[2],6)
	d=max_asymptote=round(fit$coef[3],6)
	e=inflection_point=round(fit$coef[4],6)       
  f=symmetry=round(fit$coef[5],6)                              

#Calculates a linear r^2 conc, rounded to six decimal places. NOTE this is only an estimated r^2 conc since a 4-pl curve fit doesn't normally contain an r^2 conc.
	rsquared<-round(cor(standards$logodb,logconc),6)
                      
#creates a subset of the log OD concs which contain only those concs which do not contain NaN concs  
	p1<-samples$logodb[!is.na(samples$logodb)]
                                      
#sets all NaN concs equal to the minimum conc contained within both columns of OD concs which doesn't contain an NaN conc.
#NOTE: does this included standards and samples?
	samples$logodb[is.na(samples$logodb)]<-min(p1)

#Calculates the corresponding log concentration from the 4pl equation. 
	samples$logconc=e*(((((d-cc)/(samples$logodb-cc))^(1/f))-1)^(1/b))
  standards$logconc=e*(((((d-cc)/(standards$logodb-cc))^(1/f))-1)^(1/b))

  
#Removes NaN values from the logconc columns and replaces them with NA  
  samples$logconc[is.na(samples$logconc)]<-NA
  standards$logconc[is.na(standards$logconc)]<-NA
  
#Creates a vector which contains the sample log concentration values which do not equal NaN  
  p2<-samples$logconc[!is.na(samples$logconc)]
  
#Converts log concentration concs back into standard concentration
	samples$calc_conc<-(10^(samples$logconc))*low
  standards$calc_conc<-(10^(standards$logconc))*low
  bgs$calc_conc<-bgs$conc
  
#Set any remaining NaN concs to a concentration of 0
	samples$calc_conc[is.na(samples$calc_conc)]<-0

 #creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.                      
  x <-seq(min(p2),max(p2),length=10000)
	y <-cc+((d-cc)/((1+(x/e)^b)^f))

#Creates a dataset of the log blank OD concs and the calculated log concentration
	unknowns<-data.frame(abs=c(samples$logodb),conc=c(samples$logconc))

#partitions the graph into 1 row and 2 columns
  par(mfrow=c(2,1))

#set the title for the graph
  main=paste("Assay ",name,": log-log, units")

#plots the log10 concentration vs. log10 OD concs for the calibrators. The x domain is set to include +/- 1 dilution for the calibrator to catch above and below range concs.
	plot(logconc, main=main, standards$logodb, xlim=c(min(p2),max(p2)), xlab="x=log10(Concentration) units", ylab="y=log10(Absorbance)",pch=16,col="red",cex=1.5)
   
#plots the best-fit line                       
	lines(x,y, lty="dotted", col="red")

#Superimposes the points for controls and unknowns on the best-fit line. Unknowns are blue and controls are green. The color has some transparency so that overlapping points can be visualized
	lines(unknowns$conc,unknowns$abs, type="points", pch=16,cex=1.5)

#Lays a background grid on the plot
	grid(nx=NULL,ny=NULL,col="gray")

#Creates two legends
	legend("topleft",title="5-PL Curve Fit",cex=0.75,pch=16,col=c("red","green","blue"),legend=c("Standard","Control","Unknown"),ncol=1)
	legend("bottomright",cex=0.75,pch=16,c(paste("r^2=",rsquared),paste("Min Asympt=",cc),paste("Max Asympt=",d),paste("Curvature=",b),paste("Inflection=",e),paste("Symmetry=",f),ncol=1))
 
#Merges the standards, background and samples back into one dataset.           
  combined<-rbind(bgs,standards,samples)

#splits consecutive duplicate well values into two datasets, odd and even. This is performed so that paired well values may be aligned onto the same row.
  odd<-split(combined,combined$id%%2==0)$'FALSE'
  even<-split(combined,combined$id%%2==0)$'TRUE'

#creates a table containing sample type, name and the paired well values for each specimen
  pair<-data.frame(type=odd$type,name=odd$name,well_1=odd$calc_conc,well_2=even$calc_conc)

#calculates the mean paired well value and adds a column
  pair$conc_mean<-rowMeans(pair[,3:4],na.rm=FALSE,dim=1)
              
#calculates the paired well CV% and adds a column              
  pair$conc_cv<-(apply(pair[,3:4],1,sd)/pair$conc_mean)*100
  pair$conc_cv[is.na(pair$conc_cv)]<-NA

#generates a plot of the paired well concentration vs. paired well CV% with corresponding thresholds. Specimens that fail the criteria are colored "RED", those which pass are colored "BLUE".
  plot(pair$conc_mean,pair$conc_cv,pch=16,cex=0.7,main="Assay X: Concentration vs. CV Plot",xlab="units",ylab="CV%",col=ifelse((pair$conc_mean>vlim)&(pair$conc_cv>hlimL),"red",ifelse(pair$conc_cv>hlimH,"red","blue")),xlim=c(0,ul))

#thresholds
  abline(h=hlimL)
  abline(h=hlimH)
  abline(v=vlim)

#lays a grid in the background
  grid(nx=NULL,ny=NULL,col="gray")

#creates a column with the results of whether the specimen should be considered for retest based on the criteria. TRUE=review and potentially restest, FALSE=pass
  pair$retest<-ifelse((pair$conc_mean>vlim&pair$conc_cv>hlimL|pair$conc_cv>hlimH),"TRUE","FALSE")
  pair$retest[is.na(pair$retest)]<-FALSE

#assigns paired mean, paired CV and retest result columns to the odd and even data sets generated earlier
  odd$conc_mean<-pair$conc_mean
  odd$conc_cv<-pair$conc_cv
  odd$retest<-pair$retest
  even$conc_mean<-pair$conc_mean
  even$conc_cv<-pair$conc_cv
  even$retest<-pair$retest

#combines the odd and even datasets
  results<-rbind(odd,even)

#sorts results back into the original well position order
results<-results[with(results,order(results$id)),]

#adds TRUE or FALSE to whether the mean concentration is out of range and needs to be retested with a dilution
results$dilution<-ifelse(((results$conc_mean>ul)&(results$type=="sample")),"TRUE","FALSE")


  if (is_server) {
    dev.off()
  }                                               
          
legend<-c('NA')
return(list(legend,results))

}

######################################################END OF 5-PL FUNCTION###################################################
westgard<-function(control,controls){
is_server<-Sys.getenv("is_server", unset = FALSE)
if (is_server) {
     jpeg('rplot.jpg',width = 960, height = 960, units = "px", quality = 100)
}
mean<-controls[controls$name==control[,'name'][1],2]
sd<-controls[controls$name==control[,'name'][1],3]


x<-control$conc_mean

#The following steps perform the run length encoding function (rle) on the data with different conditions on x.

res1<-rle(x>(3*sd+mean))
res1a<-rle(x>(3*sd+mean))
res1a$values=res1a$lengths>=1
res2<-rle((-3*sd+mean)>x)
res2a<-rle((-3*sd+mean)>x)
res2a$values=res2a$lengths>=1
res3<-rle(x>(2*sd+mean))
res3a<-rle(x>(2*sd+mean))
res3a$values=res3a$lengths>=2
res4<-rle((-2*sd+mean)>x)
res4a<-rle((-2*sd+mean)>x)
res4a$values=res4a$lengths>=2


#The following steps perform the inverse run length encoding function (inverse.rle) to give the indices where each of the runs are located in the vector. 

a=inverse.rle(res1)
a1=inverse.rle(res1a)
b=inverse.rle(res2)
b1=inverse.rle(res2a)
c=inverse.rle(res3)
c1=inverse.rle(res3a)
d=inverse.rle(res4)
d1=inverse.rle(res4a)


#The color value on the plot is set as an ifelse() function where, for each condition set on x (e.g. x>2 for length x>=2), if both are true, it will assign RED, if not, it will assign BLUE.

plot(x,pch=16,type="o",col=ifelse(((a==TRUE&a1==TRUE)|(b==TRUE&b1==TRUE)|(c==TRUE&c1==TRUE)|(d==TRUE&d1==TRUE)),"red","blue"),main="Assay Controls with Westgard Rules Applied, pass=blue  fail=red",xlab="",ylab="Concentration, units")


axis(1, 1:length(x),labels=control$plate,line=1,col="black",col.ticks="black",col.axis="black",cex.axis=0.5)
mtext("Plate",1,line=1,col="black",at=0.2,cex.axis=0.5)

axis(1, 1:length(x), labels=control$rungroup, line=3,cex.axis=0.5)
mtext("Run",1,line=3,at=0.2,cex.axis=0.5)

axis(1, 1:length(x), labels=control$rungroup,  line=5,cex.axis=0.5)
mtext("Plate",1,line=5,at=0.2,cex.axis=0.5)



#Lines for the mean and +/-1,2 and 3 SDs

abline(mean(x),0)
abline(mean(x)+sd(x),0)
abline(mean(x)-sd(x),0)
abline(mean(x)+2*sd(x),0)
abline(mean(x)-2*sd(x),0)
abline(mean(x)+3*sd(x),0)
abline(mean(x)-3*sd(x),0)

  if (is_server) {

    dev.off()
  }                                               
#return(which((a==TRUE&a1==TRUE)|(b==TRUE&b1==TRUE)|(c==TRUE&c1==TRUE)|(d==TRUE&d1==TRUE)))
return(control)

}



#title(main=paste("Assay X",
   # deparse(substitute(x)), "and", deparse(substitute(y)),
   # "standardized"), adj=".5")
