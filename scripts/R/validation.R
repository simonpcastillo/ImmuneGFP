
if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster, caret)

#setup path for annotations
labels<- Sys.glob('.../annotations/*_C2-*_MIP.tsv')

#setup path for masks
masks<- Sys.glob('.../masks/C2*.tiff')

#IMPORTANT
#Before running check the match between labels and masks
#######


TP=c()
FP=c()
TN=c()
FN=c()
PRF=list()
CMat =list()

for (m in 1:length(masks)){
print(masks[m])
a<- raster(masks[m])
ann<- as(a, 'matrix')

points_bkg<- read.table(labels[m]);colnames(points_bkg)<-points_bkg[1,]; points_bkg <-points_bkg[-1,]
points_bkg$x<-round(as.numeric(points_bkg$x));points_bkg$y<-round(as.numeric(points_bkg$y))

detVect<- c( )
for(i in 1:nrow(points_bkg)){
  #print(i/nrow(points_bkg)*100)
  det <- ann[points_bkg$y[i],points_bkg$x[i]]
  detVect[i]<-det
}

points_bkg$found<- detVect
FP[m] = nrow(points_bkg[points_bkg$found==255,])
TN[m] = nrow(points_bkg[points_bkg$found==0,])

points_bkg$observed<-'bkg'
for(p in 1:nrow(points_bkg)){
  if(points_bkg$found[p]==255){points_bkg$predicted[p]<-'pos'}
  if(points_bkg$found[p]==0){points_bkg$predicted[p]<-'bkg'}
}

points_pos<- read.table(labels[(m+length(masks))]);colnames(points_pos)<-points_pos[1,]; points_pos <-points_pos[-1,]
points_pos$x<-round(as.numeric(points_pos$x));points_pos$y<-round(as.numeric(points_pos$y))

detVect<- c( )
for(i in 1:nrow(points_pos)){
  #print(i/nrow(points_pos)*100)
  det <- ann[points_pos$y[i],points_pos$x[i]]
  detVect[i]<-det
}

points_pos$found<- detVect
TP[m] = nrow(points_pos[points_pos$found==255,])
FN[m] = nrow(points_pos[points_pos$found==0,])
points_pos$observed<-'pos'
for(p in 1:nrow(points_pos)){
  if(points_pos$found[p]==255){points_pos$predicted[p]<-'pos'}
  if(points_pos$found[p]==0){points_pos$predicted[p]<-'bkg'}
}

df0 <- rbind(points_bkg, points_pos)
tab<-table(df0$observed, df0$predicted)
CMat[[m]]<- confusionMatrix(tab, positive = 'pos')

n = sum(tab) 
nc = nrow(tab) 
diag = diag(tab)  
rowsums = apply(tab, 1, sum)
colsums = apply(tab, 2, sum) 

accuracy = sum(diag) / n 
precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
PRF[[m]]=data.frame(precision, recall, f1) 

}

