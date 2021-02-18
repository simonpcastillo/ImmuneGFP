 # This function computes the Morisita-Horn index for two vector of abundances over discrete sampling units.
 # The length of both vectors MUST be equal.
 # The function returns the Morisita-Horn index
 
  .mh<-function(x1, x2){
    if(length(x1)!=length(x2)) stop('length of x1 and x2 differs')
    mhdf<-data.frame(x1=as.numeric(x1), x2=as.numeric(x2))
    mhdf$x1x2 = mhdf$x1*mhdf$x2
    mhdf$sqx1<-mhdf$x1*mhdf$x1 
    mhdf$sqx2<-mhdf$x2*mhdf$x2
    #sq$fx1<-sq$x1/sum(sq$x1)
    #sq$fx2<-sq$x2/sum(x2)
    
    mh<-2*sum(mhdf$x1x2)/((sum(mhdf$sqx1)/(sum(mhdf$x1))^2 + (sum(mhdf$sqx2)/(sum(mhdf$x2))^2))*sum(mhdf$x1)*sum(mhdf$x2))
    return(mh)
  }