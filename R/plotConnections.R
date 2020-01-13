#'Plot the connections to informed individuals against event order
#'
#'For each acquisition event, the total connection of each naive individual to informed individuals is plotted for a specified
#'network.
#'
#'The individual to learn each event is plotted in red and these are connected by a red line (by default). A dotted
#'line showing the average connection to informed individuals each event is plotted in blue (by default). If the order of
#'acquisition does not follow the network in question, we would expect the red points to be distributed evenly above and
#'below the blue line. The more closely the order of acquisition follows the network, the more points will be above the blue
#'line and the higher they will be above it.
#'
#'@seealso \code{\link{nbdaData}}
#'
#'@param nbdadata an object of class \code{\link{nbdaData}}.
#'@param network numeric giving the number of the network for which the connections are to be plotted, i.e. corresponding to
#'the third dimension of assMatrix in the \code{\link{nbdaData}} object.
#'@param rescale logical indicating whether the connections should be rescaled to be plotted relative to the average
#'connection to informed individuals for each event.
#'@param lty numeric or string determining the connecting line type, see (\code{\link[graphics]{par}}).
#'@param lwd numeric determining the connecting line width, see (\code{\link[graphics]{par}}).
#'@param pchVector an optional vector of plotting symbols to be used for each type of individual as defined by \code{symbol}
#'below. See \code{pch} in (\code{\link[graphics]{par}}).
#'@param xlab passed to the \code{\link[graphics]{plot}} function.
#'@param ylab passed to the \code{\link[graphics]{plot}} function.
#'@param title an optional string giving a title for the plot.
#'@param titlePos  an optional numerical vector of length two, specifying where to plot the title.
#'@param symbol an optional vector allowing different symbols to be plotted for each point. e.g. to plot different symbols for
#'learners an non learners, specify \code{symbol=<name of nbdadata>@@status}. To choose which symbol is used for each category,
#'use \code{symbol}, e.g. \code{pchVector = c(3,8)}. Alternatively, the user could plot points based on a specific ILV, e.g.
#'\code{symbol=<name of nbdadata>@@asocILVdata[,1]} to the first asoc_ilv.
#'@param plotID an optional vector of names or indentifiers to be plotted next to the symbol for each learner. If \code{idname} was
#'provided to the \code{\link{nbdaData}} function, then this can be done by specifiying
#'\code{plotID = <name of nbdadata>@@idname[<name of nbdadata>@@orderAcq]}
#'@param offset an optional numerical vector of length two, specifying where to plot id relative to its corresponding point.
#'@param xlim passed to the \code{\link[graphics]{plot}} function.
#'@param ylim passed to the \code{\link[graphics]{plot}} function.
#'@param averageLinelty  numeric or string determining the average line type, see (\code{\link[graphics]{par}}).
#'@param averageLinecol  numeric or string determining the average line colour, see (\code{\link[graphics]{par}}).
#'@param averageLinelwd  numeric or string determining the average line colour, see (\code{\link[graphics]{par}}).
#'@param plotConnectLine logical specifying whether a connecting line is to be plotted.
#'@param plotAverageLine logical specifying whether an average line is to be plotted.
#'@return NULL


plotConnections<-function(nbdadata,network=1,rescale=F,lty=1,lwd=1,col=2,pchVector=NULL,symbol=NULL,xlab="Acquisition event", ylab=NULL,title=NULL,plotID=NULL,offset=c(0.1,0),xlim=NULL, ylim=NULL, titlePos=c(0,0),
                          averageLinelty=2,averageLinecol=4,averageLinelwd=1,plotConnectLine=T,plotAverageLine=T){

  connections<-nbdadata@stMetric[,network]
  if(rescale){

    if(is.null(ylab)){ylab="Relative connection to informed individuals"}

    connections[tapply(connections,nbdadata@time2,sum)[nbdadata@time2]>0]<-(connections/tapply(connections,nbdadata@time2,sum)[nbdadata@time2])[tapply(connections,nbdadata@time2,sum)[nbdadata@time2]>0]
    connections<-connections-tapply(connections,nbdadata@time2,mean)[nbdadata@time2]

  }else{
    if(is.null(ylab)){ylab="Total connection to informed individuals"}
  }

  if(is.null(symbol)){symbol<-rep(1,length(connections))}

  if(is.null(pchVector)){
    pch<-as.numeric(as.factor(symbol))
  }else{
    pch<-pchVector[as.numeric(as.factor(symbol))]
  }

  plot(nbdadata@time2,connections,col=nbdadata@status+1,pch=pch,xlab=xlab, ylab=ylab, main="",xlim=xlim, ylim=ylim)
  points(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=col,pch=pch[nbdadata@status==1]);
  if(plotConnectLine) lines(nbdadata@time2[nbdadata@status==1],connections[nbdadata@status==1],col=col, lty=lty,lwd=lwd);

  if(!is.null(plotID)){
    text(nbdadata@time2[nbdadata@status==1]+offset[1],connections[nbdadata@status==1]+offset[2],col=col,labels=plotID);
  }
  text(x=titlePos[1],y=titlePos[2],labels=title)
  if(plotAverageLine) lines(unique(nbdadata@time2),tapply(connections,nbdadata@time2,mean),lty=averageLinelty,col=averageLinecol,lwd=averageLinelwd)
}
