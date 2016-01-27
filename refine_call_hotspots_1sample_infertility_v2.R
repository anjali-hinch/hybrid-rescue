x=commandArgs(TRUE);
sample=x[1];
chrname=x[2];
path=x[3]

hotspot_width=300;
binding_width=700;

hotraw=read.table(paste(path,"hotraw_toberefined_",sample,"_",chrname,sep=""))

chrnamelist=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM", "chr1_GL456210_random", "chr1_GL456211_random", "chr1_GL456212_random", "chr1_GL456213_random", "chr1_GL456221_random", "chr4_GL456216_random", "chr4_GL456350_random", "chr4_JH584292_random", "chr4_JH584293_random", "chr4_JH584294_random", "chr4_JH584295_random", "chr5_GL456354_random", "chr5_JH584296_random", "chr5_JH584297_random", "chr5_JH584298_random", "chr5_JH584299_random", "chr7_GL456219_random", "chrX_GL456233_random", "chrY_JH584300_random", "chrY_JH584301_random", "chrY_JH584302_random", "chrY_JH584303_random", "chrUn_GL456239", "chrUn_GL456359", "chrUn_GL456360", "chrUn_GL456366", "chrUn_GL456367", "chrUn_GL456368", "chrUn_GL456370", "chrUn_GL456372", "chrUn_GL456378", "chrUn_GL456379", "chrUn_GL456381", "chrUn_GL456382", "chrUn_GL456383", "chrUn_GL456385", "chrUn_GL456387", "chrUn_GL456389", "chrUn_GL456390", "chrUn_GL456392", "chrUn_GL456393", "chrUn_GL456394", "chrUn_GL456396", "chrUn_JH584304")

chr=which(chrnamelist==chrname)

print(sample)


a=read.table(paste(path,"ssDNA_",sample,"_type1_filtered_only_rmdup.",chrname,sep=""),as.is=TRUE);
readstarts=as.double(a[,1]);
readends=as.double(a[,2]);
strand=a[,7];
rm(a)

#hotraw # matrix nx2 of start/end positions # mm10

chromLength = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698,16299,169725,241735,153618,39340,206961,66673,227966,14945,207968,191905,1976,195993,199368,205776,184189,953012,175968,336933,182347,259875,155838,158099,40056,22974,31704,47073,42057,20208,26764,28664,31602,72385,25871)

count_reads_in_bin=function(starts,ends,bin_left,bin_right){
    overlap_left=intersect(which(starts<bin_left),which(ends>=bin_left))
    fully_in=intersect(which(starts>=bin_left),which(ends<=bin_right))
    overlap_right=intersect(which(starts<=bin_right),which(ends>bin_right))
    tmpendsleft=ends[overlap_left]
    tmpendsleft[which(tmpendsleft>bin_right)]=bin_right
    tmpstartsright=starts[overlap_right]
    tmpstartsright[which(tmpstartsright<bin_left)]=bin_left
    count=length(fully_in)+sum((tmpendsleft-bin_left+1)/(ends[overlap_left]-starts[overlap_left]+1))+sum((bin_right-tmpstartsright+1)/(ends[overlap_right]-starts[overlap_right]+1))
    return(count)
}

get_likelihood_profile=function(hotindex,offset,hotspot_width=hotspot_width,binding_width=binding_width){
    readindexes=intersect(
        which(readstarts<=hotraw[hotindex,2]+offset),
        which(readends>=hotraw[hotindex,1]-offset))
    tmpstarts=readstarts[readindexes]
    tmpends=readends[readindexes]
    
    tmpindexwatson=which(strand[readindexes]=="+")
    tmpindexcrick=which(strand[readindexes]=="-")
    
    positions=seq(from=hotraw[hotindex,1]-offset,to=hotraw[hotindex,2]+offset,by=1)
    l=vector("numeric",length(positions))
    heat=vector("numeric",length(positions))
    
    for(i in 1:length(positions)){
        leftbigbin=positions[i]-(binding_width-hotspot_width)-ceiling(hotspot_width/2)

        watson_left=count_reads_in_bin(tmpstarts[tmpindexwatson],tmpends[tmpindexwatson],leftbigbin,leftbigbin+(binding_width-hotspot_width)-1)
        watson_centre=count_reads_in_bin(tmpstarts[tmpindexwatson],tmpends[tmpindexwatson],leftbigbin+(binding_width-hotspot_width),leftbigbin+binding_width-1)
        watson_right=count_reads_in_bin(tmpstarts[tmpindexwatson],tmpends[tmpindexwatson],leftbigbin+binding_width,leftbigbin+binding_width+(binding_width-hotspot_width)-1)
        crick_left=count_reads_in_bin(tmpstarts[tmpindexcrick],tmpends[tmpindexcrick],leftbigbin,leftbigbin+(binding_width-hotspot_width)-1)
        crick_centre=count_reads_in_bin(tmpstarts[tmpindexcrick],tmpends[tmpindexcrick],leftbigbin+(binding_width-hotspot_width),leftbigbin+binding_width-1)
        crick_right=count_reads_in_bin(tmpstarts[tmpindexcrick],tmpends[tmpindexcrick],leftbigbin+binding_width,leftbigbin+binding_width+(binding_width-hotspot_width)-1)
        
        
        mle_background = max((watson_right + crick_left),1) / (2 * (binding_width - hotspot_width));
        mle_heat = (watson_left + watson_centre + crick_centre + crick_right)/2 - (mle_background * binding_width);
        mle_background_null = max(watson_left + watson_centre + watson_right + crick_left + crick_centre + crick_right,1)/ ( 2 * ( 2*binding_width - hotspot_width));
        if(mle_heat <= 0)
        {
            mle_heat = 0;
            mle_background=mle_background_null;
        }
        
        watson_left_model = crick_right_model = mle_heat/binding_width*(binding_width - hotspot_width) + mle_background * (binding_width - hotspot_width);
        watson_centre_model = crick_centre_model = mle_heat/binding_width*(hotspot_width) + mle_background * hotspot_width;
        watson_right_model = crick_left_model = mle_background * (binding_width - hotspot_width);
        
        watson_left_null = crick_right_null = mle_background_null * (binding_width - hotspot_width);
        watson_centre_null = crick_centre_null = mle_background_null * hotspot_width;
        watson_right_null = crick_left_null = mle_background_null * (binding_width - hotspot_width);
        
        log_ll_model = watson_left * log(watson_left_model) + watson_right * log(watson_right_model) + crick_left * log(crick_left_model) + crick_right * log(crick_right_model) + watson_centre * log(watson_centre_model) + crick_centre * log(crick_centre_model)
        
        log_ll_null = watson_left * log(watson_left_null) + watson_right * log(watson_right_null) + crick_left * log(crick_left_null) + crick_right * log(crick_right_null) + watson_centre * log(watson_centre_null) + crick_centre * log(crick_centre_null);
        
        l[i] = 2 * (log_ll_model - log_ll_null);
        heat[i] = mle_heat
    }
    return(list(positions,l,heat))
}


    
centrelist=matrix(ncol=4)
print(c(hotspot_width,binding_width))
for(i in 1:length(hotraw[,1])){
    tmp=get_likelihood_profile(i,400,hotspot_width,binding_width)
    if(i==1){
        centrelist=matrix(c(tmp[[1]][which.max(tmp[[2]])],tmp[[3]][which.max(tmp[[2]])],tmp[[2]][which.max(tmp[[2]])],pchisq(tmp[[2]][which.max(tmp[[2]])], df=1, lower.tail=FALSE)),ncol=4)
    } else{
        centrelist=rbind(centrelist,matrix(c(tmp[[1]][which.max(tmp[[2]])],tmp[[3]][which.max(tmp[[2]])],tmp[[2]][which.max(tmp[[2]])],pchisq(tmp[[2]][which.max(tmp[[2]])], df=1, lower.tail=FALSE)),ncol=4))
    }
    print(i)
}
    
write.table(centrelist[which(centrelist[,4]<=1e-4),],file=paste(path,"refinedcentrev2_mm10_",sample,"_",chrname,sep=""),quote=F,col.names=F,row.names=F)

print(paste("****chr ",chr,"done"))
