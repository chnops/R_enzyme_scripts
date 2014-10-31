# These are functions for enzyme calculations based on the enzyme protocol last modified by SB on 07/2014
# Run all these first
MUB_blank_slope<-function(blank.plate, weighed.out.mub){
require(ggplot2)
	# mub.range is the range of values for your standard curve if there is any issue in the actual weighing out of MUB
	mub.range<-as.vector(c(weighed.out.mub*100/.881,weighed.out.mub*50/.881,weighed.out.mub*25/.881,weighed.out.mub*12.5/.881,weighed.out.mub*6.25/.881,weighed.out.mub*3.125/.881,weighed.out.mub*1.5625/.881,weighed.out.mub*.78125/.881))
	
	#making vectors of values for a regression for standard curve
	mub.flour<-rbind(c(blank.plate[,1],blank.plate[,2]))
	mub.conc<-rep(mub.range, 2)
	mub.st.curve.data<-data.frame(t(mub.flour),mub.conc)
	names(mub.st.curve.data)<-c("flour","conc")
	
	#calculating the slope
	mub.plate.slope<-coef(lm(flour~conc,data=mub.st.curve.data))[[2]]
p<-ggplot(mub.st.curve.data)+geom_point(aes(x=flour, y=conc))+geom_smooth(method="lm", aes(x=flour,y=conc))+theme_bw()+labs(x="Flourescence",y="Concentration")
	return(list(p,mub.plate.slope))
	
	
}
MUC_blank_slope<-function(blank.plate, weighed.out.muc){
require(ggplot2)
	# muc.range is the range of values for your standard curve if there is any issue in the actual weighing out of muc
	muc.range<-as.vector(c(weighed.out.muc*100/.876,weighed.out.muc*50/.876,weighed.out.muc*25/.876,weighed.out.muc*12.5/.876,weighed.out.muc*6.25/.876,weighed.out.muc*3.125/.876,weighed.out.muc*1.5625/.876,weighed.out.muc*.78125/.876))
	
	#making vectors of values for a regression for standard curve
	muc.flour<-rbind(c(blank.plate[,3],blank.plate[,4]))
	muc.conc<-rep(muc.range, 2)
	muc.st.curve.data<-data.frame(t(muc.flour),muc.conc)
	names(muc.st.curve.data)<-c("flour","conc")
	
	#calculating the slope
	muc.plate.slope<-coef(lm(flour~conc,data=muc.st.curve.data))[[2]]
p<-ggplot(muc.st.curve.data)+geom_point(aes(x=flour, y=conc))+geom_smooth(method="lm", aes(x=flour,y=conc))+theme_bw()+labs(x="Flourescence",y="Concentration")
	return(list(p,muc.plate.slope))
	
	
}
MUB_assay_slope<-function(assay.plate, weighed.out.mub){
	
	mub.range<-as.vector(c(weighed.out.mub*100/.881,weighed.out.mub*50/.881,weighed.out.mub*25/.881,weighed.out.mub*12.5/.881,weighed.out.mub*6.25/.881,weighed.out.mub*3.125/.881,weighed.out.mub*1.5625/.881,weighed.out.mub*.78125/.881))
	names(assay.plate)
	assayplate.mub.flour<-(c(assay.plate[,13],assay.plate[,14]))
	assayplate.mub.conc<-rep(mub.range,2)
	
	assayplate.mub.st.curve.data<-data.frame(assayplate.mub.flour,assayplate.mub.conc)
	names(assayplate.mub.st.curve.data)<-c("flour","conc")
	
	assayplate.mub.slope<-coef(lm(flour~conc,data=assayplate.mub.st.curve.data))[[2]]
	
	return(assayplate.mub.slope)
}
MUC_assay_slope<-function(assay.plate, weighed.out.muc){
	
	muc.range<-as.vector(c(weighed.out.muc*100/.876,weighed.out.muc*50/.876,weighed.out.muc*25/.876,weighed.out.muc*12.5/.876,weighed.out.muc*6.25/.876,weighed.out.muc*3.125/.876,weighed.out.muc*1.5625/.876,weighed.out.muc*.78125/.876))
	names(assay.plate)
	assayplate.muc.flour<-(c(assay.plate[,15],assay.plate[,16]))
	assayplate.muc.conc<-rep(muc.range,2)
	
	assayplate.muc.st.curve.data<-data.frame(assayplate.muc.flour,assayplate.muc.conc)
	names(assayplate.muc.st.curve.data)<-c("flour","conc")
	
	assayplate.muc.slope<-coef(lm(flour~conc,data=assayplate.muc.st.curve.data))[[2]]
	
	return(assayplate.muc.slope)
}

net_flourescence<-function(assay.plate,column.number,quench.coeff){
	flour<-(mean(assay.plate[,column.number+5],na.rm=TRUE)-mean(assay.plate[,17],na.rm=TRUE))/quench.coeff
	return(flour)
}
activity<-function(assay.plate,column.number,quench.coeff,emission.coefficient,Time,Mass){
	flour<-net_flourescence(assay.plate,column.number,quench.coeff)
	act<-(flour*.25)/(emission.coefficient*.25*Time*Mass)*100
	return(act)
}

enzyme_loop<-function(assay.plate){
	sample_ids<-as.vector(unique(assay.plate$sample.ID))

results<-data.frame()
for(i in 1:length(sample_ids)){
	#i<-1
	temp<-subset(assay.plate, sample.ID==sample_ids[i])
	
	mub_assay_slope<-MUB_assay_slope(temp,act_mub)
	muc_assay_slope<-MUC_assay_slope(temp,act_muc)
	
	mub_quench<-mub_assay_slope/mub_blank_slope
	muc_quench<-muc_assay_slope/muc_blank_slope
	
	new_row<-matrix(ncol=8,nrow=1)
	new_row[1,1]<-sample_ids[i]
	for(j in 1:7){
		if(j < 6){
			act<-activity(temp,j,mub_quench, mub_emission,mean(temp[,3]),mean(temp[,5]))
			new_row[1,1+j]<-act
			}
		if(j > 5){
			act<-activity(temp,j,muc_quench, muc_emission,mean(temp[,3]),mean(temp[,5]))
			new_row[1,1+j]<-act
			
		}
	}
	results<-rbind(results,new_row)
		
}
names(results)<-c("Sample_ID","BG","BX","CB","AP","NAG","Leu","Alal")

return(results)
}