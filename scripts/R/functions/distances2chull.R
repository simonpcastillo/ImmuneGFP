#For two vectors of coordinates x and y, this function computes the distance
# from each  cell to any of the edges of the convex hull enclosing the tumour.

#This function returns a list with several df. The first element, 'distances', is a matrix m x n; with m being the total number of cells and n the total number of convex hull edges
  
  
  
distances2chull<-function(x,y){
    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(grDevices)

    dfx<-data.frame(x=x, y=y)
    chu<-chull(dfx$x, dfx$y)
    cuh<-c(chu, chu[1])
    ver<-dfx[cuh,]
    
    for(n in 1:nrow(ver)){
      ver$x2[n]<- ver$x[n+1]
      ver$y2[n]<- ver$y[n+1]
    }
    ver<- ver[-nrow(ver),]
    for(m in 1:nrow(ver)){
      
      ver$m[m] <-(ver$y2[m] - ver$y[m])/(ver$x2[m] - ver$x[m])
      ver$b[m] <- -ver$x[m]*ver$m[m] + ver$y[m]
    }
    
    
    distances<-matrix(nrow = nrow(dfx), ncol=nrow(ver))
    x_int<-matrix(nrow = nrow(dfx), ncol=nrow(ver))
    y_int<-matrix(nrow = nrow(dfx), ncol=nrow(ver))
    slopes_m2 <- matrix(nrow = nrow(dfx), ncol=nrow(ver))
    intercept_b2<-matrix(nrow = nrow(dfx), ncol=nrow(ver))
    
    for(i in 1:nrow(dfx)){
      for(j in 1:nrow(ver)){
        
        if(ver$m[j] == 0){slopes_m2[i,j] <- 1} #slope eq point
        if(ver$m[j] != 0){slopes_m2[i,j] <- -1/ver$m[j]} #slope eq point
        
        intercept_b2[i,j]<-dfx$y[i] - dfx$x[i]*slopes_m2[i,j]
        
        
        xint <- x_int[i,j] <- (intercept_b2[i,j] - ver$b[j])/(ver$m[j] - slopes_m2[i,j])
        yint <- y_int[i,j]<- slopes_m2[i,j]*x_int[i,j] + intercept_b2[i,j]
        
        
        if(xint >= max(ver[j,c("x", "x2")]) & yint <= min(ver[j,c("y", "y2")])){
          x_int[i,j] <-max(ver[j,c("x", "x2")]); y_int[i,j] <-min(ver[j,c("y", "y2")])}
        
        if(xint >= max(ver[j,c("x", "x2")]) & yint >= max(ver[j,c("y", "y2")])){
          x_int[i,j] <-max(ver[j,c("x", "x2")]); y_int[i,j] <-max(ver[j,c("y", "y2")])}
        
        if(xint <= min(ver[j,c("x", "x2")]) & yint <= min(ver[j,c("y", "y2")])){
          x_int[i,j] <-min(ver[j,c("x", "x2")]); y_int[i,j] <-min(ver[j,c("y", "y2")])}
        
        if(xint <= min(ver[j,c("x", "x2")]) & yint >= min(ver[j,c("y", "y2")])){
          x_int[i,j] <-min(ver[j,c("x", "x2")]); y_int[i,j] <-max(ver[j,c("y", "y2")])}


        
        distances[i,j] = sqrt((x_int[i,j]-dfx$x[i])^2+(y_int[i,j]-dfx$y[i])^2 )
        
      }
      
    }
    listRet<- list(distances=distances,vertexchull=ver, x_int=x_int, y_int=y_int,
                  slopes_m2=slopes_m2,intercept_b2=intercept_b2)
    return(listRet)
  }#ElFin
