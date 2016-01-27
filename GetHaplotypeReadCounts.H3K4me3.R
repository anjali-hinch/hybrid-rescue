###
### normal code here
###
library("Rsamtools")
library("inline")
library("Rcpp")
source("SNPExtractionFunctions.R")
library("parallel")
rm(reformatReads);reformatReads <- cxxfunction(signature(variablesIn="list"),cppReformatReads,plugin="RcppArmadillo")

options(scipen=20)
btpath = "/data/smew2/altemose/software/bedtools2/bin/bedtools"
bqFilter=20 # only use bases with bq greater than this number
iSizeUpperLimit=10000 # don't use reads with mate pairs this far apart
pthresh=0.01
mapqthresh=1

args=commandArgs(TRUE)
datapath = args[2]
sample = args[3]
wide = args[4]
slide = args[5]
rep1suffix = args[6]
rep2suffix = args[7]
genomicsuffix = args[8]
qual = args[9]


#datapath="/data/smew2/altemose/MouseH3K4me3/B6", datapath="/data/smew2/altemose/MouseH3K4me3/PWDB6F1"
#sample = "Reciprocal"
#wide=1000
#slide=100
#rep1suffix=285
#rep2suffix=285
#genomicsuffix=281
#qual=1





chrs=seq(1,19,1)
chrs=c(chrs,"X")
chrin=1


for(ch in 1:length(chrs)){

chrin=chrs[ch]

print(chrin)
print(date())

rep1file = paste(datapath,"/merged.rmdup.",rep1suffix,".bam",sep="")
rep2file = paste(datapath,"/merged.rmdup.",rep2suffix,".bam",sep="")
genomicfile = paste(datapath,"/merged.rmdup.",genomicsuffix,".bam",sep="")
vcffile = paste("bychr/PWDsnps.",sample,"DSBs.chr",chrin,".vcf",sep="")
regionfile = paste("FinalOutput/FinalDSBPeakRegions.q",qual,".",sample,".chr",chrin,".",wide,"wide.",slide,"slide.bed",sep="")
outfile = paste("FinalOutput/HaplotypesPlusEnrichmentDSBs.ndr.q",qual,".",sample,".chr",chrin,".",wide,"wide.",slide,"slide.bed",sep="")

chr=paste("chr",chrin,sep="") # chromosome of interest
regions = read.table(regionfile,header=TRUE)
vcf=read.table(vcffile,colClasses=c('factor','integer','NULL','factor','factor','NULL','NULL','NULL','NULL','character'),sep="\t")
vcf = subset(vcf,grepl("1/1",vcf[,5])==TRUE)

#chr=chr,regions=regions,vcf=vcf,genomicfile=genomicfile,rep1file=rep1file,rep2file=rep2file,mapqthresh=mapqthresh

combineCounts=function(m){
	result1=c(regions[m,],rep(NA,14))
	if(regions[m,9]<300){
		regionStart = regions[m,2]
		regionEnd = regions[m,3]
		
		#print(regionStart)
		L=as.integer(vcf[,2])
		which=L>regionStart & L<regionEnd & nchar(as.character(vcf[,3]))==1 # also only bi-allelic SNPs - c++ code doesn't work otherwise
		L=L[which]
		ref=as.character(vcf[which,3])
		alt=as.character(vcf[which,4])
		T=as.integer(length(L))
		
		midpoint=floor(regionStart+(regionEnd-regionStart)/2)
		ndrstart=midpoint-60
		ndrend=midpoint+60
		ndrL=as.integer(vcf[,2])
		ndrwhich=ndrL>ndrstart & ndrL<ndrend & nchar(as.character(vcf[,3]))==1 # also only bi-allelic SNPs - c++ code doesn't work otherwise
		ndrL=ndrL[ndrwhich]
		ndrnum=as.integer(length(ndrL))
		
		inputcounts=c(0,0,0,0)
		if(regions[m,9]>0){
			inputData=getReadInformation(bamName=genomicfile,ref=ref,alt=alt,L=L,T=T,mapqthresh=mapqthresh,chr=chr,regionStart=regionStart,regionEnd=regionEnd)
			inputcounts = getRegionCounts(inputData)
		}

		IP1counts=c(0,0,0,0)
		if(regions[m,7]>0){
			IP1Data=getReadInformation(bamName=rep1file,ref=ref,alt=alt,L=L,T=T,mapqthresh=mapqthresh,chr=chr,regionStart=regionStart,regionEnd=regionEnd)
			IP1counts = getRegionCounts(IP1Data)
		}
		
		IP2counts=c(0,0,0,0)
		if(regions[m,8]>0){
			IP2Data=getReadInformation(bamName=rep2file,ref=ref,alt=alt,L=L,T=T,mapqthresh=mapqthresh,chr=chr,regionStart=regionStart,regionEnd=regionEnd)
			IP2counts = getRegionCounts(IP2Data)
		}
		
		result1=c(regions[m,],IP1counts,IP2counts,inputcounts,T,ndrnum)
	}
	
	return(result1)
}

pos=as.list(seq(1,dim(regions)[1],1))

data=mclapply(pos,combineCounts,mc.preschedule=TRUE,mc.cores=16)
data2=t(simplify2array(data))
data2[,1]=unlist(lapply(data,function(x) as.character(x[[1]])))
print(date())

options(scipen=8)
write.table(data2,file=outfile,quote=FALSE,sep="\t",row.names=F,col.names=c(names(regions),"IP1.B6reads","IP1.PWDreads","IP1.noSNPreads","IP1.conflictreads","IP2.B6reads","IP2.PWDreads","IP2.noSNPreads","IP2.conflictreads","input.B6reads","input.PWDreads","input.noSNPreads","input.conflictreads","snps.region","snps.ndr"))
#write.table(data2,file=outfile,quote=FALSE,sep="\t",row.names=F,col.names=F)

}





quit(save="no",runLast=FALSE)
