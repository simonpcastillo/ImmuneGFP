#This short script computes statistical descriptors of the Treg landscape at regional resolution (core/periphery).
# Note that the object all_2 is a dataframe with the following columns:
#- sample: sample id
#- class: cell's class Tumor or immune
#- region: cell's locaiton core or periphery
#- x, y: two columns for cell's location

#This script only uses GFP locations to draw a convexhull defining the tumour border


if (!require("pacman")) install.packages("pacman")
pacman::p_load(moments, spatstat, raster, sf)

kurtosis_Treg<-c()
skewness_Treg<-c()
median_Treg<-c()
mean_Treg<-c()

samples<-c('S174F3', 'S197M1', 'S199M1', 'S203M1')

for ( sam in samples){
  rast<- raster(kde2d(all_2[all_2$sample==sam & all_2$class=='immune','x'],
                             all_2[all_2$sample==sam & all_2$class=='immune','y'],
                             n=length(all_2[all_2$sample== sam & all_2$class=='immune','x'])))
  nameraster<- paste0('rimm_',sam)
  assign(nameraster, rast)
}


for ( sam in samples){

                      rimm<- raster(kde2d(all_2[all_2$sample==sam & all_2$class=='immune','x'],
                             all_2[all_2$sample==sam& all_2$class=='immune','y'],
                             n=length(all_2[all_2$sample==sam & all_2$class=='immune','x'])))
                      
                      
                      df<- all_2[all_2$sample==sam& all_2$class=='Tumor' & all_2$region=='core',c('x','y')]
                      chu<-chull(df$x, df$y)
                      cuh<-c(chu, chu[1])
                      ver<-df[cuh,]
                      
                      for(n in 1:nrow(ver)){
                        ver$x2[n]<- ver$x[n+1]
                        ver$y2[n]<- ver$y[n+1]
                      }
                      poly <- SpatialPolygons(list(Polygons(list(Polygon(ver[,c('x','y')])), ID='a')))
                      
                      df<- all_2[all_2$sample==sam & all_2$class=='Tumor' & all_2$region=='periphery',c('x','y')]
                      chu<-chull(df$x, df$y)
                      cuh<-c(chu, chu[1])
                      ver2<-df[cuh,]
                      
                      a<-st_polygon(list(as.matrix(ver[,c('x','y')])))
                      b<-st_polygon(list(as.matrix(ver2[,c('x','y')])))
                      resi<-st_difference(b, a)
                      plot(resi)
                      m<-as(resi, 'Spatial')
                      
                      envTreg_in<- raster::extract(rimm, poly) 
                      envTreg_out<- raster::extract(rimm, m)
                      
                      namout<- paste0(sam, '_envTreg_per')
                      namin<- paste0(sam, '_envTreg_core')
                      
                      assign(namout, envTreg_out)
                      assign(namin, envTreg_in)
                      
                      kurtosis_Treg<-append(kurtosis_Treg, c(kurtosis(envTreg_in[[1]]), kurtosis(envTreg_out[[1]])))
                      skewness_Treg<-append(skewness_Treg, c(skewness(envTreg_in[[1]]), skewness(envTreg_out[[1]])))
                      median_Treg<-append(median_Treg, c(median(envTreg_in[[1]]), median(envTreg_out[[1]])))
                      mean_Treg<-append(mean_Treg, c(mean(envTreg_in[[1]]), mean(envTreg_out[[1]])))
                      
    
}


distro_Treg<-data.frame(sample= rep(samples, each=2),
                       position= as.factor(rep(c('core', 'periphery'), 4)),
                       mean=mean_Treg,
                       median=median_Treg,
                       kurtosis=kurtosis_Treg,
                       skewness=skewness_Treg)




