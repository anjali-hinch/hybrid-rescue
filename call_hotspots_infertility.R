x=commandArgs(TRUE);
sample=x[1];
chrname=x[2];
path=x[3]

hotspot_width=500;

chrnamelist=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM", "chr1_GL456210_random", "chr1_GL456211_random", "chr1_GL456212_random", "chr1_GL456213_random", "chr1_GL456221_random", "chr4_GL456216_random", "chr4_GL456350_random", "chr4_JH584292_random", "chr4_JH584293_random", "chr4_JH584294_random", "chr4_JH584295_random", "chr5_GL456354_random", "chr5_JH584296_random", "chr5_JH584297_random", "chr5_JH584298_random", "chr5_JH584299_random", "chr7_GL456219_random", "chrX_GL456233_random", "chrY_JH584300_random", "chrY_JH584301_random", "chrY_JH584302_random", "chrY_JH584303_random", "chrUn_GL456239", "chrUn_GL456359", "chrUn_GL456360", "chrUn_GL456366", "chrUn_GL456367", "chrUn_GL456368", "chrUn_GL456370", "chrUn_GL456372", "chrUn_GL456378", "chrUn_GL456379", "chrUn_GL456381", "chrUn_GL456382", "chrUn_GL456383", "chrUn_GL456385", "chrUn_GL456387", "chrUn_GL456389", "chrUn_GL456390", "chrUn_GL456392", "chrUn_GL456393", "chrUn_GL456394", "chrUn_GL456396", "chrUn_JH584304")

chr=which(chrnamelist==chrname)


print(sample)
print(chr)


binding_width = 2*hotspot_width;
step=250;

#mm10!
chromLength = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698,16299,169725,241735,153618,39340,206961,66673,227966,14945,207968,191905,1976,195993,199368,205776,184189,953012,175968,336933,182347,259875,155838,158099,40056,22974,31704,47073,42057,20208,26764,28664,31602,72385,25871)

num_windows=ceiling(chromLength[chr]/step)

res=matrix(nrow=num_windows,ncol=8);

window_size=2*binding_width-hotspot_width



################ WITHOUT DUPLICATES

a=read.table(paste(path,"ssDNA_",sample,"_type1_filtered_only_rmdup.",chrname,sep=""),as.is=TRUE);
readstarts=as.double(a[,1]);
readends=as.double(a[,2]);
strand=a[,7];
rm(a)



bound_and_fraction=function(left_offset,size,a,b){
    index_left=ceiling((a-left_offset-size-1)/step)+1
    index_right=floor((b-left_offset-1)/step)+1
    if(index_left==index_right){
        return(c(1,index_left,index_right))
    }
    if(index_right==index_left+1){
        #fraction_left=(index_left*step+left_offset+1-a)/(b-a+1)
        #fraction_right=(b-(index_right-1)*step-left_offset)/(b-a+1)
        fraction_left=(min(b,(index_left-1)*step+left_offset+size)+1-a)/(b-a+1)
        fraction_right=(b-max(a-1,(index_right-1)*step+left_offset))/(b-a+1)
        if(fraction_left<0 || fraction_left>1 || fraction_right<0 || fraction_right>1)
        print(c(fraction_left,fraction_right))
        return(c(2,index_left,index_right,fraction_left,fraction_right))
    }
    if(index_right>index_left+1){
        fraction_left=(min(b,(index_left-1)*step+left_offset+size)+1-a)/(b-a+1)
        fraction_right=(b-max(a-1,(index_right-1)*step+left_offset))/(b-a+1)
        fraction_centre=min(size,b-a+1)/(b-a+1)
        if(fraction_left<0 || fraction_left>1 || fraction_centre<0 || fraction_centre>1 || fraction_right<0 || fraction_right>1)
        print(c(fraction_left,fraction_centre,fraction_right))
        return(c(3,index_left,index_right,fraction_left,fraction_right,fraction_centre))
    }
}

watson_left=vector("numeric",num_windows)
watson_centre=vector("numeric",num_windows)
watson_right=vector("numeric",num_windows)
crick_left=vector("numeric",num_windows)
crick_centre=vector("numeric",num_windows)
crick_right=vector("numeric",num_windows)

for(i in 1:length(readstarts)){
    if(i%%100000==0)
    print(i)
    if(strand[i]=="+"){
        #watson_left
        baf_res=bound_and_fraction(0,binding_width-hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            watson_left[baf_res[2]]=watson_left[baf_res[2]]+1
        }else if(baf_res[1]==2){
            watson_left[baf_res[2]]=watson_left[baf_res[2]]+baf_res[4]
            watson_left[baf_res[3]]=watson_left[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1)
            break
        }else if(baf_res[1]==3){
            watson_left[baf_res[2]]=watson_left[baf_res[2]]+baf_res[4]
            watson_left[(baf_res[2]+1):(baf_res[3]-1)]=watson_left[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            watson_left[baf_res[3]]=watson_left[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
        #watson_centre
        baf_res=bound_and_fraction(binding_width-hotspot_width,hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            watson_centre[baf_res[2]]=watson_centre[baf_res[2]]+1
        }else if(baf_res[1]==2){
            watson_centre[baf_res[2]]=watson_centre[baf_res[2]]+baf_res[4]
            watson_centre[baf_res[3]]=watson_centre[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1)
            break
        }else if(baf_res[1]==3){
            watson_centre[baf_res[2]]=watson_centre[baf_res[2]]+baf_res[4]
            watson_centre[(baf_res[2]+1):(baf_res[3]-1)]=watson_centre[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            watson_centre[baf_res[3]]=watson_centre[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
        #watson_right
        baf_res=bound_and_fraction(binding_width,binding_width-hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            watson_right[baf_res[2]]=watson_right[baf_res[2]]+1
        }else if(baf_res[1]==2){
            watson_right[baf_res[2]]=watson_right[baf_res[2]]+baf_res[4]
            watson_right[baf_res[3]]=watson_right[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 )
            break
        }else if(baf_res[1]==3){
            watson_right[baf_res[2]]=watson_right[baf_res[2]]+baf_res[4]
            watson_right[(baf_res[2]+1):(baf_res[3]-1)]=watson_right[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            watson_right[baf_res[3]]=watson_right[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
    } else if(strand[i]=="-"){
        #crick_left
        baf_res=bound_and_fraction(0,binding_width-hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            crick_left[baf_res[2]]=crick_left[baf_res[2]]+1
        }else if(baf_res[1]==2){
            crick_left[baf_res[2]]=crick_left[baf_res[2]]+baf_res[4]
            crick_left[baf_res[3]]=crick_left[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1)
            break
        }else if(baf_res[1]==3){
            crick_left[baf_res[2]]=crick_left[baf_res[2]]+baf_res[4]
            crick_left[(baf_res[2]+1):(baf_res[3]-1)]=crick_left[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            crick_left[baf_res[3]]=crick_left[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
        #crick_centre
        baf_res=bound_and_fraction(binding_width-hotspot_width,hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            crick_centre[baf_res[2]]=crick_centre[baf_res[2]]+1
        }else if(baf_res[1]==2){
            crick_centre[baf_res[2]]=crick_centre[baf_res[2]]+baf_res[4]
            crick_centre[baf_res[3]]=crick_centre[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 )
            break
        }else if(baf_res[1]==3){
            crick_centre[baf_res[2]]=crick_centre[baf_res[2]]+baf_res[4]
            crick_centre[(baf_res[2]+1):(baf_res[3]-1)]=crick_centre[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            crick_centre[baf_res[3]]=crick_centre[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
        #crick_right
        baf_res=bound_and_fraction(binding_width,binding_width-hotspot_width,readstarts[i],readends[i])
        if(baf_res[1]==1){
            crick_right[baf_res[2]]=crick_right[baf_res[2]]+1
        }else if(baf_res[1]==2){
            crick_right[baf_res[2]]=crick_right[baf_res[2]]+baf_res[4]
            crick_right[baf_res[3]]=crick_right[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 )
            break
        }else if(baf_res[1]==3){
            crick_right[baf_res[2]]=crick_right[baf_res[2]]+baf_res[4]
            crick_right[(baf_res[2]+1):(baf_res[3]-1)]=crick_right[(baf_res[2]+1):(baf_res[3]-1)]+baf_res[6]
            crick_right[baf_res[3]]=crick_right[baf_res[3]]+baf_res[5]
            if(baf_res[4]<0|| baf_res[4]>1 || baf_res[5]<0 || baf_res[5]>1 || baf_res[6]<0 || baf_res[6]>1)
            break
        }
    }
}



## NO CONTROL
res=matrix(nrow=num_windows,ncol=8);
for(i in 1:num_windows)
{
    if(i %% 10000 == 0)
    cat(paste("Tests completion:",as.integer(i/num_windows*100),"%\n"));
    
    #hotspot_centre = candidates[i];
    #left_left = hotspot_centre + hotspot_width/2 - binding_width;
    #left_centre = hotspot_centre - hotspot_width/2;
    #right_centre = hotspot_centre + hotspot_width/2;
    #right_right = hotspot_centre - hotspot_width/2 + binding_width;
    
    #watson_left = length(which(watson_starts <= left_centre & watson_ends >= left_left));
    #watson_centre = length(which(watson_starts <= right_centre & watson_ends >= left_centre));
    #watson_right = length(which(watson_starts <= right_right & watson_ends >= right_centre));
    #crick_left = length(which(crick_starts <= left_centre & crick_ends >= left_left));
    #crick_centre = length(which(crick_starts <= right_centre & crick_ends >= left_centre));
    #crick_right = length(which(crick_starts <= right_right & crick_ends >= right_centre));
    
    
    mle_background = max((watson_right[i] + crick_left[i]),1) / (2 * (binding_width - hotspot_width));
    mle_heat = (watson_left[i] + watson_centre[i] + crick_centre[i] + crick_right[i])/2 - (mle_background * binding_width);
    mle_background_null = max(watson_left[i] + watson_centre[i] + watson_right[i] + crick_left[i] + crick_centre[i] + crick_right[i],1)/ ( 2 * ( 2*binding_width - hotspot_width));
    if(mle_heat <= 0)
    {
        mle_heat = 0;
        mle_background=mle_background_null;
    }
    
    watson_left_model = crick_right_model = mle_heat/2 + mle_background * (binding_width - hotspot_width);
    watson_centre_model = crick_centre_model = mle_heat/2 + mle_background * hotspot_width;
    watson_right_model = crick_left_model = mle_background * (binding_width - hotspot_width);
    
    watson_left_null = crick_right_null = mle_background_null * (binding_width - hotspot_width);
    watson_centre_null = crick_centre_null = mle_background_null * hotspot_width;
    watson_right_null = crick_left_null = mle_background_null * (binding_width - hotspot_width);
    
    log_ll_model = watson_left[i] * log(watson_left_model) + watson_right[i] * log(watson_right_model) + crick_left[i] * log(crick_left_model) + crick_right[i] * log(crick_right_model) + watson_centre[i] * log(watson_centre_model) + crick_centre[i] * log(crick_centre_model)
    
    log_ll_null = watson_left[i] * log(watson_left_null) + watson_right[i] * log(watson_right_null) + crick_left[i] * log(crick_left_null) + crick_right[i] * log(crick_right_null) + watson_centre[i] * log(watson_centre_null) + crick_centre[i] * log(crick_centre_null);
    
    ll_ratio = 2 * (log_ll_model - log_ll_null);
    
    pvalue = pchisq(ll_ratio, df=1, lower.tail=FALSE);
    hotspot_centre=(i-1)*step+1+binding_width-hotspot_width/2
    res[i,]=c(hotspot_centre, mle_background, mle_heat, mle_background_null, log_ll_model, log_ll_null, ll_ratio, pvalue);
    
    #if(i %% 1e4 == 0)
    #{
    #print("saving work so far");
    #save(res,file=paste("~/chr",chr,"_algo_",mouse,"_",hotspot_width,"_",binding_width,sep=""));
    #}
}
save(res,file=paste(path,"/",chrname,"_rawhotspots_",sample,"_rmdup_",hotspot_width,"_",binding_width,sep=""));






