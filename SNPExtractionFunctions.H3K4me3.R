###
### c++ function to be called within R
###

library("Rsamtools")
library("inline")
library("Rcpp")

###
### c++ function to be called within R
###
cppReformatReads <- '
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>
//
// only a single list with everything comes in - declare off this list
//
  Rcpp::List cvariablesIn(variablesIn);
//
// now subdivide into components in c++
//
  arma::ivec firstReadSpot = as<arma::ivec>(cvariablesIn["firstReadSpot"]);
  arma::ivec secondReadSpot = as<arma::ivec>(cvariablesIn["secondReadSpot"]);
  arma::ivec posRead = as<arma::ivec>(cvariablesIn["posRead"]);
  const int numberOfReads = as<int>(cvariablesIn["numberOfReads"]);
  const int T = as<int>(cvariablesIn["T"]);
  arma::ivec L = as<arma::ivec>(cvariablesIn["L"]);
  Rcpp::CharacterVector seqRead = as<Rcpp::CharacterVector>(cvariablesIn["seqRead"]);
  Rcpp::CharacterVector qualRead = as<Rcpp::CharacterVector>(cvariablesIn["qualRead"]);
  Rcpp::CharacterVector ref = as<Rcpp::CharacterVector>(cvariablesIn["ref"]);
  Rcpp::CharacterVector alt = as<Rcpp::CharacterVector>(cvariablesIn["alt"]);
  Rcpp::List splitCigarRead = as<Rcpp::List>(cvariablesIn["splitCigarRead"]);
  arma::ivec lengthOfSplitCigarRead = as<arma::ivec>(cvariablesIn["lengthOfSplitCigarRead"]);
  arma::ivec iSizeTooBigRead = as<arma::ivec>(cvariablesIn["iSizeTooBigRead"]);
  arma::ivec mapq = as<arma::ivec>(cvariablesIn["mapq"]);
  const int bqFilter = as<int>(cvariablesIn["bqFilter"]);
//
// new variables
//
  int x1, x2, y, i, iPair, iM, iRead, t;
  int nSNPInRead = 0;
  int nReadSpanningSNPs = 0;
  char s;
  double d, phiTemp, phiRef, phiAlt, phi;
  Rcpp::List sampleReads;
  arma::ivec seqLocal(200);
  arma::ivec qualLocal(200);
  arma::ivec posLocal(200); // there shouldnt be this many SNPs
  int refPosition, refOffset, strandOffset;
  int iNumOfMs, readPosOverall, loopEnd, cigarLength;
  Rcpp::CharacterVector cigarType(1);
  Rcpp::CharacterVector cigar2(1);
  int iU, nU, iAll, curPos, whileVar, localbq;
  int tMin=0;
  int tMax=0; // left and right boundaries of what SNPs to look at
//
// loop over each pair of reads
//
  for(iRead=0; iRead<=numberOfReads-1; iRead++)
//  for(iRead=16909; iRead<=16909; iRead++)
  {
//std::cout << "iRead " << iRead <<  "\\n";
    nSNPInRead=-1; // reset to 0
    // only use if the insert size is within acceptable margins
    if(iSizeTooBigRead(firstReadSpot(iRead))==0)
    {
      // determine whether there are 2 or 1 read
      loopEnd=1;
      if(secondReadSpot(iRead)==-1) // no second read
        loopEnd=0; // ergo only loop over one read
      // loop twice over the two parts of the read
      for(iPair=0; iPair<=loopEnd; iPair++)
      {
        // set position in reads of current read
        if(iPair==0)
          readPosOverall=firstReadSpot[iRead];
        if(iPair==1)
          readPosOverall=secondReadSpot[iRead];
        //
        // also as the GATK bounds BQ my MQ, and only uses BQ>17, only use mapq>17
        //
        if(mapq(readPosOverall)>=bqFilter)
        {
          //
          // for this read, find eligible SNPs (tStart, tEnd)
          //
          double curPos;
          if(loopEnd==0) // only one read
            curPos=posRead(firstReadSpot(iRead));
          if(loopEnd==1)
            curPos=(posRead(firstReadSpot(iRead))+posRead(secondReadSpot(iRead)))/2;
          // now move tMin forward until A)first within 10000 and B) not in front of SNP
          whileVar=0;
          while(whileVar==0)
          {
            // dont continue if too far
            if(tMin<(T-1))
            {
              // continue while no more than 10000 bp before
              if((2000+ L(tMin))<curPos)
              {
                tMin++;
              } else {
                whileVar=1; // break loop - done!
              }
            } else {
              whileVar=1; // break loop
            }
          }
          // now - go right until >1000 bp away
          whileVar=0;
          while(whileVar==0)
          {
            // dont continue if too far
            if(tMax<(T-1))
            {
              // continue while no more than 1000 bp after
              if(( L(tMax)-2000)<curPos)
              {
                tMax++;
              } else {
                whileVar=1; // break loop - done!
              }
            } else {
              whileVar=1; // break loop
            }
          }
          //  if(iRead % 10 == 0 && iRead<200)
          //
          //
          // now, for this read, calculate whether there are SNPs
          //
          refPosition=posRead(readPosOverall);
          Rcpp::List cigarReadInfo = as<Rcpp::List>(splitCigarRead(readPosOverall));
          iNumOfMs = lengthOfSplitCigarRead(readPosOverall);
          // set some things
          refOffset=0; // offset against the reference sequence
          strandOffset=0; // offset in the strand
          // get cigar info from the read
          arma::ivec cigarLengthVec = as<arma::ivec>(cigarReadInfo(0));
          Rcpp::CharacterVector cigarTypeVec = as<Rcpp::CharacterVector>(cigarReadInfo(1));
          //
          // now, loop over each part of the read (M, D=del, I=ins)
          //
          for(iM=0;iM<=iNumOfMs;iM++)
          {
            cigarLength=cigarLengthVec(iM);
            cigarType(0)=cigarTypeVec(iM);
            // if its an M - scan
            if(cigarType(0)=="M")
            {
              x1 = refPosition + refOffset; // left part of M
              x2 = refPosition + refOffset + cigarLength-1; // right part of M
              for(t=tMin; t<=tMax; t++) // determine whether that snps is spanned by the read
              {
                y = L[t];
                if(x1 <= y && y <= x2) // if this is true - have a SNP!
                {
                  s = seqRead[readPosOverall][y-refPosition-refOffset+strandOffset];
                  // check if ref or ALT - only keep if true
                  // also only use if BQ at least bqFilter (17) (as in 17 or greater)
                  localbq=int(qualRead[readPosOverall][y-refPosition-refOffset+strandOffset])-33;
                  // also bound BQ above by MQ
                  if(localbq>mapq(readPosOverall)) // if greater, than reduce
                    localbq=mapq(readPosOverall);
                  if((s==ref[t][0] || s==alt[t][0]) && (localbq>=bqFilter))
                  {
                    // is this the reference or alternate?
                    nSNPInRead = nSNPInRead+1;
                    if(s==ref[t][0])
                      seqLocal[nSNPInRead] = 0;
                    if(s==alt[t][0])
                      seqLocal[nSNPInRead] = 1;
                    qualLocal[nSNPInRead] = localbq;
                    posLocal[nSNPInRead] = t;
                  } // end of check if ref or alt
                } // end of whether this SNP intersects read
              } // end of loop on SNP
              // now, bump up ref and pos offset by x1
              refOffset=refOffset + cigarLength;
              strandOffset=strandOffset + cigarLength;
            } // end of if statement on whether cigar type is M
            // if it is an insertion - bump the strand offset
            if(cigarType(0)=="I")
              strandOffset=strandOffset+cigarLength;
            // if it is a deletion - bump the reference position
            if(cigarType(0)=="D")
              refOffset=refOffset+cigarLength;
          } // close for loop on each M type within read
        } // end of check on mapping quality of this read
      } // end of loop on 1 or two reads
      if(nSNPInRead > -1) // save result!
      {
        arma::ivec pR = posLocal.subvec(0,nSNPInRead); // position (R means Read)
        arma::ivec sR = seqLocal.subvec(0,nSNPInRead); // sequence
        arma::ivec qR = qualLocal.subvec(0,nSNPInRead); // quality
        // get average physical location
        //
        // turn this into one unique value per SNP
        //
        arma::ivec pRU = arma::unique(pR); // the U means unique
        nU = pRU.n_elem-1; // 0-based length of pRU - ie (n)umber of (U)nique SNPs in read
        arma::vec phiU(nU+1); // keep phis here
        for(iU=0; iU<=nU; iU++) // for each unique entry
        {
          // reset phis
          phiAlt=1;
          phiRef=1;
          // go through all elements looking for that SNP
          for(iAll=0; iAll<=nSNPInRead; iAll++)
          {
            if(pR(iAll)==pRU(iU)) // there is a match - consider
            {
              // turn into phi
              //   calculate probability from phred scale
              phiTemp= 1 - pow(10,-(double(qR(iAll))/10));
              if(qR(iAll)<=0)
                phiTemp=0.5; // BQ = 0 so no probability so make it 0.5
              // scale to appropriate base
              phi = (1-phiTemp) * ( 1-sR(iAll)) + phiTemp * sR(iAll);
              phiAlt=phiAlt * phi;
              phiRef=phiRef * (1-phi);
            } // end if statement on whether there is a match
          } // end of for loop going through all SNPs in the read
          // now calculate probability - calculate numerator and denominator
          // where for a= product_{i=0} P(alt,i)
          // where for b= product_{i=0} (1-P(alt,i))
          // phi = P(alt) = a/(a+b)
          phiU(iU)=phiAlt/(phiAlt+phiRef); // done!
        } // end of for loop for each unique element
        //
        // get physical position for the read
        // hmmm - right now in effect weighted by occurence
        //
        d=0;
        for(i=0; i<=nSNPInRead; i++)
          d = d + L(pR(i));
        d = d/(1+nSNPInRead);
        //
        // save results but dont label list elements to save space
        //
        // save a smaller version unless need to debug
        sampleReads.push_back(Rcpp::List::create(
          Rcpp::Named("nSNPs")= nU,
          Rcpp::Named("avPos")= d,
          Rcpp::Named("snpInL")= pRU,
          Rcpp::Named("snpProb")= phiU,
          Rcpp::Named("iRead")= iRead,
          Rcpp::Named("frs")= firstReadSpot[iRead],
          Rcpp::Named("srs")= secondReadSpot[iRead]));
        //sampleReads.push_back(Rcpp::List::create(iRead,nU,d,phiU,pRU,pR,sR,qR));
        //sampleReads.push_back(Rcpp::List::create(iRead,nSNPInRead,d,phiU,pRU,pR,sR,qR));
        // save number of SNPs as well
        nReadSpanningSNPs = nReadSpanningSNPs +1;
      }  // end of save result
    } // close if statement on whether or not the insert size is too big
  } // close for loop on read pair
  //
  // done
  //
  return(wrap(sampleReads));
'



###
### for a bam, for a set of SNPs, get information
###
getReadInformation=function(bamName,ref,alt,L,T,mapqthresh,chr=chr,regionStart=regionStart,regionEnd=regionEnd)
{
  ###
  ### set some flags and load in a single result
  ###
  flag=     scanBamFlag(isPaired = TRUE, isProperPair = NA, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, isMinusStrand = NA, isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA, isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
  ### okay actually start loading here
  what=c("qname","strand","pos","seq","qual","cigar","isize","mapq")
  ### set which region of the genome to interrogate
  eval(parse(text=  (paste("which = RangesList(\"",chr,"\"=IRanges(",regionStart,",",regionEnd,"))",sep=""))))
  #idx=paste(substr(bamName,1,nchar(bamName)-3),"bam.bai",sep="")
  idx1=paste(substr(bamName,1,nchar(bamName)-3),"bai",sep="")
  idx2=paste(bamName,".bai",sep="")
  if(file.exists(idx1)) idx=idx1
  if(file.exists(idx2)) idx=idx2
  param=ScanBamParam(flag=flag,which=which,what=what) # define parameters
  sampleData=scanBam(file=bamName,index=idx,param=param) # load the data
  ###
  ### reformat some things
  ###
  # super basic - read or pair
  qname=sampleData[[1]]$qname
  qnameUnique=unique(qname)
  if(length(qnameUnique)==0) {
    return(NA)
  }
  else{
  qnameInteger=match(qname,qnameUnique)
  # get positions - NOTE - THEY ARE 0 BASED
  firstReadSpot=as.integer(match(1:max(qnameInteger),qnameInteger)-1) # first instance
  y=qnameInteger
  y[firstReadSpot+1]=NA
  secondReadSpot=as.integer(match(1:max(qnameInteger),y)-1) # second read - MAY BE NA
  secondReadSpot[is.na(secondReadSpot)  ]=as.integer(-1) # Switch to -1 - skip over if -1 in c++ code
  numberOfReads=as.integer(length(firstReadSpot))
  # get more info read as well
  mapq=as.integer(sampleData[[1]]$mapq)
  posRead=as.integer(sampleData[[1]]$pos)
  cigarRead=sampleData[[1]]$cigar
  strandRead=sampleData[[1]]$strand
  seqRead=as.character(sampleData[[1]]$seq)
  qualRead=as.character(sampleData[[1]]$qual) # hmm
  qualRead=as.character(sampleData[[1]]$qual) # hmm
  iSizeTooBigRead=as.integer(abs(sampleData[[1]]$isize)>iSizeUpperLimit | as.integer(sampleData[[1]]$mapq)<mapqthresh)
  iSizeTooBigRead[is.na(iSizeTooBigRead)==TRUE]=0
  
  ###
  ### also, need to reconfigure cigar properly
  ###
  splitCigarRead=lapply(1:length(cigarRead),  function(x)  list(as.integer(100),"M"))
  # also lost of 51M - skip these
  which=  cigarRead=="51M"
  splitCigarRead[which]=lapply(1:sum(which),  function(x)  list(as.integer(51),"M"))
  which=cigarRead!="51M" & cigarRead!="100M"
  splitCigarRead[which]=   lapply(cigarRead[which],function(c) {
  y=unlist(strsplit(c,""))
  t1=is.na(as.numeric(y))
  t2=(1:length(t1))[t1]
  t3=c(1,t2[-length(t2)]+1)
  n=length(t2)
  t5=array(0,n)
  t6=array("",n)
  for(i in 1:n)
  {
    t5[i]=substr(c,t3[i],t2[i]-1)
    t6[i]=substr(c,t2[i],t2[i])
  }
       return(list(as.integer(t5),t6))
  })
  lengthOfSplitCigarRead=as.integer(unlist(lapply(splitCigarRead,function(x) length(x[[1]])-1))) # 0 BASED
  ###
  ### push through c++ function
  ###
  sampleReads=reformatReads(list(
    firstReadSpot=firstReadSpot,
    secondReadSpot=secondReadSpot,
    posRead=posRead,
    numberOfReads=numberOfReads,
    T=T,
    L=L,
    seqRead=seqRead,
    qualRead=qualRead,
    ref=ref,
    alt=alt,
    splitCigarRead=splitCigarRead,
    lengthOfSplitCigarRead=lengthOfSplitCigarRead,
    iSizeTooBigRead=iSizeTooBigRead,
    mapq=mapq,
    bqFilter=bqFilter))
  ###
  ### return something clean
  ###
  # quick check
  a=unlist(lapply(sampleReads,function(x) qname[x$frs+1]))
  b=unlist(lapply(sampleReads,function(x) {
    if(x$srs==-1) return(NA) else return(qname[x$srs+1]) }))
  if(sum(a!=b,na.rm=TRUE)>0) {
    print("ERROR - problem associating reads");     return(NA)
  }
  sampleReads=lapply(sampleReads,function(x) return(x[names(x)!="frs" & names(x)!="srs" & names(x)!="iRead"]))
  names(sampleReads)=a
  return(sampleReads)
  }
}





getRegionCounts=function(nData)
{
	b6count=0
	pwdcount=0
	nosnpcount=0
	contracount=0
	if(length(nData)>0){
		counts=sapply(nData, getSNPcounts,simplify="matrix")
		nosnpcount=sum(counts[1,]==0 & counts[2,]==0)
		b6count= sum(counts[1,]>0 & counts[2,]==0)
		pwdcount= sum(counts[1,]==0 & counts[2,]>0)
		contracount = sum(counts[1,]>0 & counts[2,]>0)
	}
	return(c(b6count,pwdcount,nosnpcount,contracount))
}

getSNPcounts = function(nlist){
	bcount=as.integer(sum(nlist$snpProb<=pthresh))
	pcount=as.integer(sum(nlist$snpProb>=(1-pthresh)))
	return (c(bcount,pcount))
}






rm(reformatReads);reformatReads <- cxxfunction(signature(variablesIn="list"),cppReformatReads,plugin="RcppArmadillo")
