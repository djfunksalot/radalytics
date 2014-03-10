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
