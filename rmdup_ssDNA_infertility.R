x=commandArgs(TRUE);
sample=x[1];
chr=x[2];
path=x[3];


ssDNA=read.table(paste(path,"ssDNA_",sample,"_type1_filtered_only.",chr,sep=""))
watson=ssDNA[which(ssDNA[,7]=="+"),]
crick=ssDNA[which(ssDNA[,7]=="-"),]
keepwatson=rep(TRUE,length(watson[,1]))
keepcrick=rep(TRUE,length(crick[,1]))

watson=watson[order(watson[,1],watson[,2]),]
crick=crick[order(crick[,1],crick[,2]),]

watsonforloop=matrix(0,length(watson[,1]),2)
crickforloop=matrix(0,length(crick[,1]),2)

watsonforloop[,1]=watson[,1]
watsonforloop[,2]=watson[,2]
crickforloop[,1]=crick[,1]
crickforloop[,2]=crick[,2]

watsonlooptext<-'
Rcpp::NumericMatrix watson(watsonR);
Rcpp::LogicalVector keepwatson(keepwatsonR);
int nrows = watson.nrow();
for(int i = 0; i < nrows; i++){
    if(keepwatson[i]==TRUE){
        if(watson(i,0)==watson(i+1,0)){
            int j=i;
            while((j<=nrows-2) && (watson(j,0)==watson(j+1,0))){
                j=j+1;
            }
            int count=0;
            for(int l=i;l<=j;l++){
                if((watson(l,1)==watson(i,1)) && (keepwatson[l]==TRUE)){
                    count=count+1;
                }
            }
            int count2=0;
            if(count>1){
                for(int l=i;l<=j;l++){
                    if((watson(l,1)==watson(i,1)) && (keepwatson[l]==TRUE) && (count2<=count-2)){
                        keepwatson[l]=FALSE;
                        count2=count2+1;
                    }
                }
            }
        }
    }
}
return(keepwatson);
'
library("Rcpp")
library("inline")
watsonloop <- cxxfunction(signature(watsonR="matrix",keepwatsonR="vector"),watsonlooptext,plugin="Rcpp")


cricklooptext<-'
Rcpp::NumericMatrix crick(crickR);
Rcpp::LogicalVector keepcrick(keepcrickR);
int nrows = crick.nrow();
for(int i = 0; i < nrows; i++){
    if(keepcrick[i]==TRUE){
        if(crick(i,0)==crick(i+1,0)){
            int j=i;
            while((j<=nrows-2) && (crick(j,0)==crick(j+1,0))){
                j=j+1;
            }
            int count=0;
            for(int l=i;l<=j;l++){
                if((crick(l,1)==crick(i,1)) && (keepcrick[l]==TRUE)){
                    count=count+1;
                }
            }
            int count2=0;
            if(count>1){
                for(int l=i;l<=j;l++){
                    if((crick(l,1)==crick(i,1)) && (keepcrick[l]==TRUE) && (count2<=count-2)){
                        keepcrick[l]=FALSE;
                        count2=count2+1;
                    }
                }
            }
        }
    }
}
return(keepcrick);
'

crickloop <- cxxfunction(signature(crickR="matrix",keepcrickR="vector"),cricklooptext,plugin="Rcpp")


keepwatson=watsonloop(watsonR=watsonforloop,keepwatsonR=keepwatson)
keepcrick=crickloop(crickR=crickforloop,keepcrickR=keepcrick)

ssDNAunique=rbind(watson[keepwatson,],crick[keepcrick,])
ssDNAunique=ssDNAunique[order(ssDNAunique[,1],ssDNAunique[,2]),]

write.table(ssDNAunique,file=paste(path,"ssDNA_",sample,"_type1_filtered_only_rmdup.",chr,sep=""),quote=F,col.names=F,row.names=F)

