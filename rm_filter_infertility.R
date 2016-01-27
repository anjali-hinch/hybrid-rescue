if(TRUE)
{
x=commandArgs(TRUE);
sample=x[1];
chr=x[2];
path=x[3];
} 

print(sample)
print(chr)

min_read_qual=20; 

#a=read.table(paste(path,mouse,"/chrs/ssDNA_type1.chr",chr,sep=""),as.is=TRUE);
a=read.table(paste(path,"ssDNA_",sample,"_type1.",chr,sep=""),as.is=TRUE);
quality=a[,4]; 
x=strsplit(quality,"_"); 
y=matrix(unlist(x),ncol=2,byrow=TRUE); 
leftqual=as.double(y[,1]); 
rightqual=as.double(y[,2]); 

good=which(leftqual >= min_read_qual | rightqual >= min_read_qual); 
a=a[good,]
leftqual=leftqual[good];
rightqual=rightqual[good]; 

type=a[,5]; 
x=strsplit(type,"_"); 
y=matrix(unlist(x),ncol=2,byrow=TRUE); 
itr=as.double(y[,1]); 
mho=as.double(y[,2]); 

readstart=a[,2]; 
readend=a[,3]; 
strand=a[,6]; 

good=1:nrow(a); 
res=cbind(readstart[good],readend[good],leftqual[good],rightqual[good],itr[good],mho[good],strand[good]); 
#write.table(res,file=paste(path,mouse,"/chrs/ssDNA_type1_filtered_only.chr",chr,sep=""),row.names=FALSE,col.names=FALSE); 
write.table(res,file=paste(path,"ssDNA_",sample,"_type1_filtered_only.",chr,sep=""),row.names=FALSE,col.names=FALSE); 
