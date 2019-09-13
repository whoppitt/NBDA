plotConnections<-function(nbdadata,network=1,rescale=F,lty=1, symbol=NULL,xlab="Acquisition event", ylab=NULL,title=NULL,plotID=NULL,offset=c(0.1,0),xlim=NULL, ylim=NULL, titlePos=c(0,0),
                          averageLinelty=2,averageLinecol=4){

  connections<-nbdadata@stMetric[,network]
  if(rescale){

    if(is.null(ylab)){ylab="Relative connection to informed individuals"}

    connections[tapply(connections,nbdadata@time2,sum)[nbdadata@time2]>0]<-(connections/tapply(connections,nbdadata@time2,sum)[nbdadata@time2])[tapply(connections,nbdadata@time2,sum)[nbdadata@time2]>0]
    connections<-connections-tapply(connections,nbdadata@time2,mean)[nbdadata@time2]

  }else{
    if(is.null(ylab)){ylab="Total connection to informed individuals"}
  }

    if(is.null(symbol)){
      plot(nbdadata@time2,connections,col=nbdadata@status+1,xlab=xlab, ylab=ylab,main="",xlim=xlim, ylim=ylim);
      points(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=2);
      lines(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=2, lty=lty);

    }else{

      plot(nbdadata@time2,connections,col=nbdadata@status+1,pch=as.numeric(as.factor(symbol)),xlab=xlab, ylab=ylab, main="",xlim=xlim, ylim=ylim)
      points(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=2,pch=as.numeric(as.factor(symbol))[nbdadata@status==1]);
      lines(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=2, lty=lty);
    }
    if(!is.null(plotID)){
      text(nbdadata@time2[nbdadata@status==1]+offset[1],connections[nbdadata@status==1]+offset[2],col=2,labels=plotID);
    }
    text(x=titlePos[1],y=titlePos[2],labels=title)
    lines(unique(nbdadata@time2),tapply(connections,nbdadata@time2,mean),lty=averageLinelty,col=averageLinecol)
}
