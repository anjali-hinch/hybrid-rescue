#ForceCall.H3K4me3.Final.r
#by Nicolas Altemose
#10 November 2015
#implementation of an algorithm to perform LR testing of ChIP-seq data at a fixed set of sites
######

##Initialise inputs and outputs

library("parallel")
options(scipen=20)
btpath = "/data/smew2/altemose/software/bedtools2/bin/bedtools" #path to bedtools executable
rep1=3 #column in which fragment depth for IP replicate 1 is reported
rep2=4 #column in which fragment depth for IP replicate 2 is reported
genomic=5 #column in which fragment depth for genomic input is reported
chrs=c(seq(1,19,1),"X")
system("mkdir tmp",ignore.stdout = T, ignore.stderr = T)

args=commandArgs(TRUE)
constfile = args[2] #path to a file containing genome-wide estimates for constants alpha1/2 & beta (output of GetConstants.r)
bedfilebase = args[3] #path and base name for a 3-column bed file listing positions of windows in which to do force-calling, split by chromosome
posfileIP1base = args[4] #path and base name for files listing fragment positions for IP replicate 1 (ending in e.g. "chr1.bed")
posfileIP2base = args[5] #path and base name for files listing fragment positions for IP replicate 2 (ending in e.g. "chr1.bed")
posfileGbase = args[6] #path and base name for files listing fragment positions for genomic input (ending in e.g. "chr1.bed")
outfile = args[7] #name of outputfile

##example hardwired input
#constfile = "Constants.q1.Infertile.1000wide.100slide.bed"
#bedfilebase = "tmp/dmc1hotspots_PWDB6F1.PRDM9pb.txt.1kbslop"
#posfileIP1base = "/data/smew2/altemose/MouseH3K4me3/PWDB6F1/bychr/FragDepth.q0.243.chr"
#posfileIP2base = "/data/smew2/altemose/MouseH3K4me3/PWDB6F1/bychr/FragDepth.q0.245.chr"
#posfileGbase = "/data/smew2/altemose/MouseH3K4me3/PWDB6F1/bychr/FragDepth.q0.241.chr"
#outfile = "ForceCall.H3K4me3.1kb.PWDB6F1.PRDM9pb.Infertile.txt"

#read in constants
constdata=read.table(constfile,header=TRUE)

alpha1.est=constdata[dim(constdata)[1],2]
alpha2.est=constdata[dim(constdata)[1],3]
beta.est=constdata[dim(constdata)[1],4]


#declare function for each chromosome
getEnrichments=function(chr){

posfileIP1 = paste(posfileIP1base,".chr",chr,".bed",sep="")
posfileIP2 = paste(posfileIP2base,".chr",chr,".bed",sep="")
posfileG= paste(posfileGbase,".chr",chr,".bed",sep="")

bedfile= paste(bedfilebase,".chr",chr,".bed",sep="")

tempfile1=paste("/tmp/",outfile,".chr",chr,".temp1.bed",sep="")
tempfile2=paste("/tmp/",outfile,".chr",chr,".temp2.bed",sep="")
tempfileG=paste("/tmp/",outfile,".chr",chr,".tempG.bed",sep="")

outfiletemp = paste("tmp/",outfile,".chr",chr,".bed",sep="")


## now get coverage at windows centred on DSB midpoints

system(paste(btpath," coverage -a ",posfileIP1," -b ",bedfile," -counts >",tempfile1,sep=""))
system(paste(btpath," coverage -a ",posfileIP2," -b ",bedfile," -counts >",tempfile2,sep=""))
system(paste(btpath," coverage -a ",posfileG," -b ",bedfile," -counts >",tempfileG,sep=""))


#print(paste("reading in results and computing constants for DSB regions.",date()))
counts = read.table(tempfile1,header=FALSE,colClasses=c('character','integer','integer','integer'))

counts=counts[order(counts[,2]),]
countstemp = read.table(tempfile2,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
countstemp=countstemp[order(countstemp[,1]),]
counts[,5]=countstemp[,2]
countstemp = read.table(tempfileG,header=FALSE,colClasses=c('NULL','integer','NULL','integer'))
countstemp=countstemp[order(countstemp[,1]),]
counts[,6]=countstemp[,2]
rm(countstemp)

countsfilt=counts[,2:6]
countsfilt[(countsfilt[,genomic]==0 & (countsfilt[,rep1]+countsfilt[,rep2])>0),genomic]=0.5 #pseudocount regions with 0 genomic coverage and >0 IP coverage to have genomic coverage of 0.5
countsfilt=countsfilt[countsfilt[,genomic]>0,]

peaks=makepeaks(countsfilt,alpha1=alpha1.est,alpha2=alpha2.est,beta=beta.est,r1=rep1,r2=rep2,g=genomic)

validcount = 0
for (i in 1:dim(counts)[1]){
	startpos = counts[i,2]
	peakinfo = peaks[peaks[,1]==startpos,]
	if(length(peakinfo)>0){
		counts[i,7]=peakinfo[3]
		counts[i,8]=peakinfo[5]
		counts[i,9]=peakinfo[4]
		validcount=validcount+1
	}else{
		counts[i,7]="NA"
		counts[i,8]="NA"
		counts[i,9]="NA"
	}
}

colnames(counts)=c("chr","start","stop","cov_r1","cov_r2","cov_g","enrichment","Lhood_diff","p-value")


##write final output file

options(scipen=8)
write.table(counts,file=outfiletemp,quote=FALSE,sep="\t",row.names=F,col.names=T)

return(1)

}


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
	
	rm(term1,term4,term5,term6,term7)

	yvals=(-bterm+sqrt(bterm^2-4*aterm*cterm))/2/aterm
	rm(aterm,bterm,cterm)
	
	yvals[yvals<0]=0
	yvals[yvals>1e9]=1e9

	bvals=(sumcov)/(term2+term3*yvals)

	bvalsnull=(sumcov)/(term2)
	yvalsnull=rep(0,length(bvalsnull))

	lhooddiff=2*(lhood(yvals,bvals,test,alpha1,alpha2,beta,r1,r2,g)-lhood(yvalsnull,bvalsnull,test,alpha1,alpha2,beta,r1,r2,g))
	
	rm(bvalsnull,yvalsnull)
	
	signif=pchisq(lhooddiff,df=1,lower.tail=F)

	results=cbind(test[,1],test[,2],yvals,signif,lhooddiff,test[,r1],test[,r2],test[,g])
	colnames(results)=c("start","stop","yhat_alt","p-value","Lhood_diff","cov_r1","cov_r2","cov_g")
	return(results)
	
}
lhood=function(yhat,bhat,test,alpha1,alpha2,beta,r1=rep1,r2=rep2,g=genomic){
	sumcov = test[,r1]+test[,r2]+test[,g]
	ourterm=sumcov*(log(bhat)-1)+test[,r1]*log(alpha1+yhat)+test[,r2]*log(alpha2+yhat*beta)
	return(ourterm)
}


#####run all chromosomes in parallel

print(date())
funfunc = mclapply(chrs,getEnrichments,mc.preschedule=TRUE,mc.cores=20)
print(paste("done!:",date()))


finaltable=read.table(paste("tmp/",outfile,".chr1.bed",sep=""), header=T)
for(chr in chrs[2:length(chrs)]){
	outfiletemp = paste("tmp/",outfile,".chr",chr,".bed",sep="")
	tempdata1=read.table(outfiletemp, header=T)
	finaltable=rbind(finaltable,tempdata1)
}

write.table(finaltable,file=paste(outfile,".bed",sep=""),sep="\t",col.names=T,row.names=F,quote=F)

colnames(finaltable)=c("chr","start","stop","cov_r1","cov_r2","cov_g","enrichment","Lhood_diff","p-value")

quit(save="no",runLast=FALSE)

