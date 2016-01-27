#findpeaks.lowmem2.r
#by Nick Altemose
#20 Feb 2013
######

##step1: initialise inputs and outputs
library("parallel")
options(scipen=20)
btpath = "/data/smew2/altemose/software/bedtools2/bin/bedtools"
rep1=3
rep2=4
genomic=5

args=commandArgs(TRUE)
datapath = args[2]
sample = args[3]
wide=args[4]
slide=args[5]
rep1suffix = args[6]
rep2suffix = args[7]
genomicsuffix = args[8]
qual = as.integer(args[9])



#example hardwired input
#datapath="/data/smew2/altemose/MouseH3K4me3/PWDB6F1"
#sample = "Infertile"
#wide=1000
#slide=100
#rep1suffix = 243
#rep2suffix = 245
#genomicsuffix = 241
#qual=1





chrs=seq(1,19,1)
chrs=c(chrs,"X")

outfile1 = paste("FinalOutput/Constants.q",qual,".",sample,".",wide,"wide.",slide,"slide.bed",sep="")




getConstants=function(chr){

posfileA = paste(datapath,"/bychr/FragPos.q",qual,".",rep1suffix,".chr.chr",chr,".bed",sep="")
posfileB = paste(datapath,"/bychr/FragPos.q",qual,".",rep2suffix,".chr.chr",chr,".bed",sep="")
posfileG= paste(datapath,"/bychr/FragPos.q",qual,".",genomicsuffix,".chr.chr",chr,".bed",sep="")

windowfile = paste("bychr/chr",chr,".windows.",wide,"wide.",slide,"slide.bed",sep="")

infile1=paste("tmp/Temp0.FragDepth.",sample,".q",qual,".",rep1suffix,".chr",chr,".",wide,"wide.",slide,"slide.bed",sep="")
infile2=paste("tmp/Temp0.FragDepth.",sample,".q",qual,".",rep2suffix,".chr",chr,".",wide,"wide.",slide,"slide.bed",sep="")
infile3=paste("tmp/Temp0.FragDepth.",sample,".q",qual,".",genomicsuffix,".chr",chr,".",wide,"wide.",slide,"slide.bed",sep="")


##step2: load data, estimate constants alpha1, alpha2, and beta


system(paste(btpath," coverage -a ",posfileA," -b ",windowfile," -counts >",infile1,sep=""))
system(paste(btpath," coverage -a ",posfileB," -b ",windowfile," -counts >",infile2,sep=""))
system(paste(btpath," coverage -a ",posfileG," -b ",windowfile," -counts >",infile3,sep=""))


counts = read.table(infile1,header=FALSE,colClasses=c('NULL','integer','integer','integer'))

counts=counts[order(counts[,1]),]
countstemp = read.table(infile2,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
countstemp=countstemp[order(countstemp[,1]),]
counts[,4]=countstemp[,2]
countstemp = read.table(infile3,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
countstemp=countstemp[order(countstemp[,1]),]
counts[,5]=countstemp[,2]
rm(countstemp)

alpha1.est = sum(counts[(counts[,rep2]==0),rep1])/sum(counts[(counts[,rep2]==0),genomic])
alpha2.est = sum(counts[(counts[,rep1]==0),rep2])/sum(counts[(counts[,rep1]==0),genomic])
beta.est0 = (mean(counts[,rep2])-alpha2.est*mean(counts[,genomic]))/(mean(counts[,rep1])-alpha1.est*mean(counts[,genomic]))

zeroregions1=sum((counts[,rep2]==0) & (counts[,rep1] + counts[,genomic])>0)
zeroregions2=sum((counts[,rep1]==0) & (counts[,rep2] + counts[,genomic])>0)


#counts=counts[counts[,genomic]>0,] #remove regions with 0 genomic coverage
counts[(counts[,genomic]==0 & (counts[,rep1]+counts[,rep2])>0),genomic]=0.5 #pseudocount regions with 0 genomic coverage and >0 IP coverage to have genomic coverage of 0.5
counts=counts[counts[,genomic]>0,]


##step4: find initial set of p-values, find confident set of peaks, re-do estimate of beta and redo p-value calls


peaks=makepeaks(counts,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est0,r1=rep1,r2=rep2,g=genomic)

q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-10 & peaks[,"yhat_alt"]>5)
if(length(q)<100){
	q=which(!is.na(peaks[,"p-value"]) & peaks[,"p-value"]<1e-5 & peaks[,"yhat_alt"]>1)
}
rm(peaks)

beta.est = (mean(counts[q,rep2])-alpha2.est*mean(counts[q,genomic]))/(mean(counts[q,rep1])-alpha1.est*mean(counts[q,genomic]))


peaks=makepeaks(counts,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est,r1=rep1,r2=rep2,g=genomic)
gthresh = quantile(peaks[,"cov_g"],0.999)


#save peak positions and print constant values to terminal

r1comb=mean(peaks[,"bhat_alt"]*(peaks[,"yhat_alt"]+alpha1.est))
r1sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"])
r2comb=mean(peaks[,"bhat_alt"]*(beta.est*peaks[,"yhat_alt"]+alpha2.est))
r2sig=mean(peaks[,"bhat_alt"]*peaks[,"yhat_alt"]*beta.est)


return(c(chr,alpha1.est,alpha2.est,beta.est,mean(counts[,genomic]),mean(counts[,rep1]),mean(counts[,rep2]),as.integer(dim(counts)[1]),zeroregions1,zeroregions2,length(q),r1sig/r1comb,r2sig/r2comb,gthresh))

}


##step3: declare functions to find MLE values for each window

makepeaks=function(test=counts,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){

	sumcov = test[,r1]+test[,r2]+test[,g]

	term1=(sumcov)*(beta+1)
	term2=1+alpha1+alpha2
	term3=beta+1
	term4=beta*alpha1*test[,r2]+alpha2*test[,r1]
	term5=beta*(test[,r1]+test[,r2])
	term6=alpha1*alpha2
	term7=alpha1*beta+alpha2

	aterm=term1*beta-term3*term5
	bterm=term1*term7-term2*term5-term3*term4
	cterm=term1*term6-term2*term4
	
	rm(term1,term2,term3,term4,term5,term6,term7)

	yvals=(-bterm+sqrt(bterm^2-4*aterm*cterm))/2/aterm
	rm(aterm,bterm,cterm)
	
	yvals[yvals<0]=0
	yvals[yvals>1e9]=1e9

	bvals=(sumcov)/(1+alpha1+alpha2+(beta+1)*yvals)

	bvalsnull=(sumcov)/(1+alpha1+alpha2)
	yvalsnull=rep(0,length(bvalsnull))

	lhooddiff=2*(lhood(yvals,bvals,test,alpha1,alpha2,beta,r1,r2,g)-lhood(yvalsnull,bvalsnull,test,alpha1,alpha2,beta,r1,r2,g))
	
	rm(bvalsnull,yvalsnull)
	
	signif=pchisq(lhooddiff,df=1,lower.tail=F)

	results=cbind(test[,1],test[,2],yvals,signif,lhooddiff,test[,r1],test[,r2],test[,g],bvals)
	colnames(results)=c("start","stop","yhat_alt","p-value","Lhood_diff","cov_r1","cov_r2","cov_g","bhat_alt")
	return(results)
	
}
lhood=function(yhat,bhat,test,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){
	sumcov = test[,r1]+test[,r2]+test[,g]
	ourterm=sumcov*(log(bhat)-1)+test[,r1]*log(alpha1+yhat)+test[,r2]*log(alpha2+yhat*beta)
	return(ourterm)
}


coln= c("chr","alpha1","alpha2","beta","meancovgenomic","meancovrep1","meancovrep2","totalnonzerobins","alpha1trainingregions","alpha2trainingregions","betatrainingregions","rep1signal","rep2signal","genomiccov999thpctile")

print(date())
data=mclapply(chrs,getConstants,mc.preschedule=TRUE,mc.cores=20)
data2=t(simplify2array(data))
data2[,1]=unlist(lapply(data,function(x) as.character(x[[1]])))

data2=rbind(data2,c("autosomal",rep("NA",13)))
for(m in c(2,3,4,5,6,7,12,13,14)){
	data2[21,m]=weighted.mean(as.numeric(data2[1:19,m]),as.numeric(data2[1:19,8]))
}
for(m in c(8,9,10,11)){
	data2[21,m]=sum(as.numeric(data2[1:19,m]))
}

write.table(data2,file=outfile1,quote=FALSE,sep="\t",row.names=F,col.names=coln)

print(paste("printed results to",outfile1))
print(date())
quit(save="no",runLast=FALSE)

