#Extracts the sample size for nbdaData in order for aicc to be calculated. As Burnham and Anderson (2002) point out, sample size is not always a straithforward issue. Here we take it to be the number of acquisition events
sampSizeExtract<-function(nbdadata){

	if(class(nbdadata)=="nbdaData"|class(nbdadata)=="dTADAData"){

		return(sum(nbdadata@status))

	}

  if(is.list(nbdadata)){

    totalSS<-0;

    for(i in 1:length(nbdadata)){

      totalSS<-totalSS+sampSizeExtract(nbdadata[[i]]);

    }

    return(totalSS);
  }

}

