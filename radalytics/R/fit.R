# sets the ascribed value in a result panel
#
# rawresults - the raw results without ascribed values
# standards - standards by name with ascribed values
"applyStandards"<-
function()  {
data<-data.frame(id=rawresults$id,positionPlate=rawresults$positionPlate,objectType=rawresults$uqc,measuredValue=as.numeric(as.vector(rawresults$valueMeasured)),objectName=rawresults$layoutPlate,ascribedValue=0,wavelength=rawresults$wavelength,sequenceNum=rawresults$sequenceNum)
	i<-0
	while (i<nrow(standards)) {
		i<-i+1
		standard<-standards[i,]
                if(standard$ascribedValue) {
	 		data[data$objectType == 'standard' & data$objectName == standard$objectName,]$ascribedValue<-standard$ascribedValue
		}
        }
	return(data)
}
#########################################################START OF FIT FUNCTION######################################################################

#########################################################END OF FIT FUNCTION######################################################################
"fit"<-
  function() {
    if(experimenttype=='singlecurve' && fitmethod=='fourpl'){
      singlefourplfit()
    }
    if(experimenttype=='doublecurve' && fitmethod=='fivepl'){
      doublefiveplfit()
    }
    
  }


#########################################################START OF single 4-PL FUNCTION######################################################################


# data - this is the output data from plate reader reader via SPS
#
# ul - the upper limit of the assay as specified by the manufacturer
#
# vlim: (vertical limit) -  selected based on a number of semi-arbitrary criteria. 
#    Examples; the second lowest calibrator, the upper limit divided by 10
#    or just by eyeballing the concentration vs. CV plot 
#
# hlimL: (horizonal limit low) - cv threshold for unknowns that have a concentration above vlim
#  concentration above the vertical threshold (vlim).
#
# hlimH:  cv threshold for concentration below the vertical threshold and above the limit of detection
#
# name: this will be used only for r generated plot renderings (not normally used)
#

"singlefourplfit"<-
function() {

#required for drm() function used for non-linear regression
  require('drc')
  #check to see if we are running on a server
  hlimL<-as.numeric(hlimL)
  hlimH<-as.numeric(hlimH)
  vlim<-as.numeric(vlim)
  ul<-as.numeric(ul)
  
  #split data array into samples and standards, (controls are treated as samples)
  samples<-rbind(split(data,data$objectType)$sample,split(data,data$objectType)$control)
  standards<-split(data,data$objectType)$standard
  
  #creates a table containing sample type, name and the grouped well values for each specimen
  names<-unique(standards$objectName)
  groupedStandards<-data.frame()
  i<-0;while(i<=length(names)){
    if (toString(names[i])!='NA' && toString(names[i])!='') {
      objectName<-toString(names[i])
      groupedList<-subset(standards,standards$objectName==toString(names[i]) & standards$measuredValue >0)
      measuredMeanValue<-mean(groupedList$measuredValue)    
      ascribedValue<-groupedList[1,]$ascribedValue
      measuredSd<-sd(groupedList$measuredValue)        
      measuredSd[is.na(measuredSd)]<-NA
      measuredCv<-((measuredSd/measuredMeanValue)*100)           
      measuredCv[is.na(measuredCv)]<-NA
      #calculates the grouped well sd and adds a column 
      groupedStandards<-rbind(groupedStandards,data.frame(objectName,ascribedValue,measuredMeanValue,measuredSd,measuredCv))
    }
    i<-i+1
  }
  #print(groupedStandards)
  #Dose response function fits a 4pl model to the calibrator dataset
  fit<-drm(data.frame(groupedStandards$measuredMeanValue,groupedStandards$ascribedValue),fct=LL.4())
  
  #set IDs for the coefficients calculated by the 4-pl fit.
  b=hill_slope=round(fit$coef[1],6)
  cc=min_asymptote=round(fit$coef[2],6)
  d=max_asymptote=round(fit$coef[3],6)
  e=inflection_point=round(fit$coef[4],6)                            
  
  #Calculates a linear r^2 conc, rounded to six decimal places. 
  #NOTE this is only an estimated r^2 conc since a 4-pl curve fit
  #doesn't normally contain an r^2 conc.
  rsquared<-round(cor(groupedStandards$measuredMeanValue,groupedStandards$ascribedValue),6)
  #Calculates the corresponding concentration from the 4pl equation. 
  # set all to 1 for now
  samples$dilution[samples$objectType=='control']<-1
  samples$dilution[samples$objectType=='standard']<-1
  samples$dilution[samples$objectType=='sample']<-1
  #
  samples$calc_conc=e*(((-d+samples$measuredValue)/(cc-samples$measuredValue))^(1/b))
  samples$calc_conc_dil=samples$calc_conc*samples$dilution
  

  # standards are never diluted
  standards$dilution<-1
  # retrofit the standards
  standards$calc_conc=e*(((-d+standards$measuredValue)/(cc-standards$measuredValue))^(1/b))
  
  #Removes NaN values from the logconc columns and replaces them with NA  
  samples$calc_conc[is.na(samples$calc_conc)]<-NA
  standards$calc_conc[is.na(standards$calc_conc)]<-NA
  
  # they'll always be the same but we need this column later for merging with samples 
  standards$calc_conc_dil=standards$calc_conc
  
  #Creates a vector which contains the sample log concentration values which do not equal NaN  
  p2<-samples$calc_conc[!is.na(samples$calc_conc)]
  
  #creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.
  min<-min(p2)
  max<-max(p2)
  

  
  #Merges the standards, background and samples back into one dataset.           
  combined<-rbind(standards,samples)
  
  #creates a table containing sample type, name and the grouped well values for each specimen
  names<-unique(combined$objectName)
  group<-data.frame()
  i<-0;while(i<=length(names)){
    if (toString(names[i])!='NA' && toString(names[i])!='') {
      objectName<-toString(names[i])
      groupedList<-subset(combined,combined$objectName==toString(names[i]))
      interpretedMeanValue<-round(mean(groupedList$calc_conc),2)    
      interpretedSd<-round(sd(groupedList$calc_conc),2)       
      interpretedSd[is.na(interpretedSd)]<-NA
      interpretedCv<-round((interpretedSd/interpretedMeanValue)*100,2)           
      interpretedCv[is.na(interpretedCv)]<-NA
      #calculates the grouped well sd and adds a column 
      #print(objectName)
      group<-rbind(group,data.frame(objectName,interpretedMeanValue,interpretedSd,interpretedCv))
    }
    i<-i+1
  }
  
  results<-merge(combined,group,by.x=5,by.y=1)
  
  
  # generates a plot of the grouped well concentration vs. grouped well CV% 
  # with corresponding thresholds. Specimens that fail the criteria are colored "RED",
  # those which pass are colored "BLUE".
  
  # creates a column with the results of whether the specimen should be considered for retest
  # TRUE=review and potentially restest, FALSE=pass
  results$retest<-ifelse((results$interpretedMeanValue>vlim&results$interpretedCv>hlimL|results$interpretedCv>hlimH),"TRUE","FALSE")
  results$retest[is.na(results$retest)]<-TRUE

  #sorts results back into the original well position order
  results<-results[with(results,order(results$id)),]

  # if the mean concentration is out of range and needs to be retested with a dilution
   results$outOfRange<-ifelse(((results$interpretedMeanValue>ul)&(results$objectType=="sample")),"TRUE","FALSE")
   results$outOfRange[is.na(results$dilution)]<-FALSE
                                             

#results<-results[with(results,order(results$name)),]
curve<-list(rsquared=rsquared,cc=cc,d=d,b=b,e=e,min=min,max=max)
return(list(curve=curve,results=results,standards=standards,group=group))
          
}

#########################################################END OF single 4-PL FUNCTION######################################################################


#########################################################START OF double 5-PL FUNCTION######################################################################


# results - this is the output data from plate reader reader via SPS
#
# ul - the upper limit of the assay as specified by the manufacturer
#
# vlim: (vertical limit) -  selected based on a number of semi-arbitrary criteria. 
#    Examples; the second lowest calibrator, the upper limit divided by 10
#    or just by eyeballing the concentration vs. CV plot 
#
# hlimL: (horizonal limit low) - cv threshold for unknowns that have a concentration above vlim
#  concentration above the vertical threshold (vlim).
#
# hlimH:  cv threshold for concentration below the vertical threshold and above the limit of detection
#
# name: this will be used only for r generated plot renderings (not normally used)

"doublefiveplfit"<-function() {
  require('drc')
  data<-applyStandards()
  #split data array into samples and standards, (controls are treated as samples)
  samples<-rbind(split(data,data$objectType)$sample,split(data,data$objectType)$control)
  standards<-split(data,data$objectType)$standard
  #creates a table containing sample type, name and the grouped well values for each specimen
  names<-unique(standards$objectName)
  t<-1;while(t<=2){
    groupedStandards<-data.frame()
    w<-1;while(w<=2){
      primary_wavelength<-get(paste("read_",t,"_",1,sep=""))
      if(exists(paste("read_",t,"_",2,sep=""))){
        secondary_wavelength<-get(paste("read_",t,"_",2,sep=""))
      }
      assign('primary_wavelength', primary_wavelength, envir = .GlobalEnv)
      w<-w+1
    }
    varname<-paste("read_",t,"_",2,sep="")
    i<-0;while(i<=length(names)){
      if (toString(names[i])!='NA' && toString(names[i])!='') {
        objectName<-toString(names[i])
        groupedList<-subset(standards,standards$objectName==toString(names[i]) & standards$measuredValue >0)
        measuredMeanValue<-mean(groupedList$measuredValue)    
        ascribedValue<-groupedList[1,]$ascribedValue
        measuredSd<-sd(groupedList$measuredValue)        
        measuredSd[is.na(measuredSd)]<-NA
        measuredCv<-((measuredSd/measuredMeanValue)*100)           
        measuredCv[is.na(measuredCv)]<-NA
        #calculates the grouped well sd and adds a column 
        groupedStandards<-rbind(groupedStandards,data.frame(objectName,ascribedValue,measuredMeanValue,measuredSd,measuredCv))
      }
      i<-i+1
    }
    t<-t+1
  }
  
    #print(groupedStandards)
    #Dose response function fits a 4pl model to the calibrator dataset
    fit<-drm(data.frame(groupedStandards$measuredMeanValue,groupedStandards$ascribedValue),fct=LL.4())
    
    #set IDs for the coefficients calculated by the 4-pl fit.
    b=hill_slope=round(fit$coef[1],6)
    cc=min_asymptote=round(fit$coef[2],6)
    d=max_asymptote=round(fit$coef[3],6)
    e=inflection_point=round(fit$coef[4],6)                            
    
    #Calculates a linear r^2 conc, rounded to six decimal places. 
    #NOTE this is only an estimated r^2 conc since a 4-pl curve fit
    #doesn't normally contain an r^2 conc.
    rsquared<-round(cor(groupedStandards$measuredMeanValue,groupedStandards$ascribedValue),6)
    #Calculates the corresponding concentration from the 4pl equation. 
    # set all to 1 for now
    samples$dilution[samples$objectType=='control']<-1
    samples$dilution[samples$objectType=='standard']<-1
    samples$dilution[samples$objectType=='sample']<-1
    #
    samples$calc_conc=e*(((-d+samples$measuredValue)/(cc-samples$measuredValue))^(1/b))
    samples$calc_conc_dil=samples$calc_conc*samples$dilution
    
    
    # standards are never diluted
    standards$dilution<-1
    # retrofit the standards
    standards$calc_conc=e*(((-d+standards$measuredValue)/(cc-standards$measuredValue))^(1/b))
    
    #Removes NaN values from the logconc columns and replaces them with NA  
    samples$calc_conc[is.na(samples$calc_conc)]<-NA
    standards$calc_conc[is.na(standards$calc_conc)]<-NA
    
    # they'll always be the same but we need this column later for merging with samples 
    standards$calc_conc_dil=standards$calc_conc
    
    #Creates a vector which contains the sample log concentration values which do not equal NaN  
    p2<-samples$calc_conc[!is.na(samples$calc_conc)]
    
    #creates a sequence of concs which take on the concs of the best-fit line for a 4-pl curve.
    min<-min(p2)
    max<-max(p2)
    
    #Merges the standards, background and samples back into one dataset.           
    combined<-rbind(standards,samples)
    
    #creates a table containing sample type, name and the grouped well values for each specimen
    names<-unique(combined$objectName)
    group<-data.frame()
    i<-0;while(i<=length(names)){
      if (toString(names[i])!='NA' && toString(names[i])!='') {
        objectName<-toString(names[i])
        groupedList<-subset(combined,combined$objectName==toString(names[i]))
        interpretedMeanValue<-mean(groupedList$calc_conc)    
        interpretedSd<-sd(groupedList$calc_conc)        
        interpretedSd[is.na(interpretedSd)]<-NA
        interpretedCv<-((interpretedSd/interpretedMeanValue)*100)           
        interpretedCv[is.na(interpretedCv)]<-NA
        #calculates the grouped well sd and adds a column 
        #print(objectName)
        group<-rbind(group,data.frame(objectName,interpretedMeanValue,interpretedSd,interpretedCv))
      }
      i<-i+1
    }
    
    results<-merge(combined,group,by.x=5,by.y=1)
    
    
    # generates a plot of the grouped well concentration vs. grouped well CV% 
    # with corresponding thresholds. Specimens that fail the criteria are colored "RED",
    # those which pass are colored "BLUE".
    
    # creates a column with the results of whether the specimen should be considered for retest
    # TRUE=review and potentially restest, FALSE=pass
    results$retest<-ifelse((results$interpretedMeanValue>vlim&results$interpretedCv>hlimL|results$interpretedCv>hlimH),"TRUE","FALSE")
    results$retest[is.na(results$retest)]<-TRUE
    
    #sorts results back into the original well position order
    results<-results[with(results,order(results$id)),]
    
    # if the mean concentration is out of range and needs to be retested with a dilution
    results$outOfRange<-ifelse(((results$interpretedMeanValue>ul)&(results$objectType=="sample")),"TRUE","FALSE")
    results$outOfRange[is.na(results$dilution)]<-FALSE
    
    #results<-results[with(results,order(results$name)),]
    curve<-list(rsquared=rsquared,cc=cc,d=d,b=b,e=e,min=min,max=max)
    return(list(curve=curve,results=results,standards=standards,group=group))
    
  }

#########################################################END OF double 5-PL FUNCTION######################################################################
