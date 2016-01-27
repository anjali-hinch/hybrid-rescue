x=commandArgs(TRUE);
sample=x[1];
chrname=x[2];
path=x[3]



pval_threshold=1e-4
hotspot_width=500
binding_width = 2*hotspot_width

getHotspots=function(path,chrname,sample,hotspot_width,binding_width){
    load(paste(path,chrname,"_rawhotspots_",sample,"_rmdup_",hotspot_width,"_",binding_width,sep=""))
    hotspot_list=matrix(0,length(res[,1]),6)
    include=FALSE
    hotspot_list_index=1
    start=0
    end=0
    maxpval=1
    heatmax=0
    centre=0
    background=0
    for(i in 1:dim(res)[1]){
        if(res[i,8]<pval_threshold){
            if(include==FALSE){
                start=res[i,1]-hotspot_width/2
                end=res[i,1]+hotspot_width/2
                maxpval=res[i,8]
                heatmax=res[i,3]
                background=res[i,4]
                centre=res[i,1]
                include=TRUE
            }
            else{
                end=res[i,1]+hotspot_width/2
                if(res[i,8]<maxpval){
                    maxpval=res[i,8]
                    centre=res[i,1]
                }
                if(res[i,3]>heatmax){
                    heatmax=res[i,3]
                    background=res[i,4]
                }
            }
        }
        else{
            if(include==TRUE){
                include=FALSE
                hotspot_list[hotspot_list_index,]=c(start,end,centre,maxpval,heatmax,background)
                hotspot_list_index=hotspot_list_index+1
            }
        }
    }
    if(include==TRUE){
        hotspot_list[hotspot_list_index,]=c(start,end,centre,maxpval,heatmax,background)
        hotspot_list_index=hotspot_list_index+1
    }
    hotspot_list=hotspot_list[1:(hotspot_list_index-1),]
    rm(res)
    return(hotspot_list)
}


hot=getHotspots(path,chrname,sample,hotspot_width,binding_width)

tmp=matrix(0,ncol=2,nrow=length(hot[,1]))
tmp[,1]=hot[,1]-400
tmp[,2]=hot[,2]+400

write.table(tmp,file=paste(path,"hotraw_toberefined_",sample,"_",chrname,sep=""),col.names=F,row.names=F,quote=F)

