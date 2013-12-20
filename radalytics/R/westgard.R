westward<-function(control,mean,sd){

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

axis(1, 1:length(x), labels=control$run, line=3,cex.axis=0.5)
mtext("Run",1,line=3,at=0.2,cex.axis=0.5)

axis(1, 1:length(x), labels=control$day, adj="vertical", line=5,cex.axis=0.5)
mtext("Plate",1,line=5,at=0.2,cex.axis=0.5)



#Lines for the mean and +/-1,2 and 3 SDs

abline(mean(x),0)
abline(mean(x)+sd(x),0)
abline(mean(x)-sd(x),0)
abline(mean(x)+2*sd(x),0)
abline(mean(x)-2*sd(x),0)
abline(mean(x)+3*sd(x),0)
abline(mean(x)-3*sd(x),0)

print(which((a==TRUE&a1==TRUE)|(b==TRUE&b1==TRUE)|(c==TRUE&c1==TRUE)|(d==TRUE&d1==TRUE)))

}



#title(main=paste("Assay X",
   # deparse(substitute(x)), "and", deparse(substitute(y)),
   # "standardized"), adj=".5")
