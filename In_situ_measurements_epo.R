#In_situ measurements and Horizon classification#
library(tripack)
library(SDMTools)
library(rgl) 
library(spectroscopy)
library(pbapply)
library(pca3d)

load('huntervalley_2cm_ground.RDATA')
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('pr_varExp.RData')
load('check_plots.RData')
load('EPO_transformation_matrix.RDATA')


prev_dir <-getwd()

#####filter the new spectra dataset and EPO tronsform it####

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_hv_DATA<-pblapply(raw_hv_spectra,correct_step)    
trim_hv_DATA<- pblapply(cont_hv_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
abs_hv_DATA<-pblapply(trim_hv_DATA,function(x) absorbance<-log(1/x))  
filt_hv_DATA<- pblapply(abs_hv_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_hv_DATA<- pblapply(filt_hv_DATA,function(x) strip_spectra(x,c(500:2450),which=10)) 
snv_hv_DATA <- pblapply(strip_hv_DATA,snvBLC)
epo_hv_DATA <- pblapply(snv_hv_DATA,function(x) {
  epo_tmp<-as.matrix(x) %*% P
  epo_tmp<-data.frame(epo_tmp)
  colnames(epo_tmp)<-seq(500,2450,by =10)
  epo_tmp}
)

####then bring in and remove outliers from global dataset####

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(500:2450),which=10))
snv_DATA <- pblapply(strip_DATA,snvBLC)
epo_DATA <- pblapply(snv_DATA,function(x) {
   epo_tmp<-as.matrix(x) %*% P
   epo_tmp<-data.frame(epo_tmp)
   colnames(epo_tmp)<-seq(500,2450,by =10)
   epo_tmp}
)

pr_spectra<-prcomp(do.call(rbind,epo_DATA), center=T,scale=T) 
screeplot(pr_spectra)#visualize the PC
pr_varExp(do.call(rbind,epo_DATA))#check the acummulative variation on the data
pr_scores <- pr_spectra$x 

###do some convex.hull analysis###
rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2])
rand.ch <- convex.hull(rand_tr,plot.it=F)
pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))

plot(pr_scores[,1],pr_scores[,2])
lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
###outlier detection and removal procedure###
mean_pcaA<- (colMeans(pr_scores[,1:5])) # mean of each of the components
mean_pcaA
cov_pcaA<- cov(pr_scores[,1:5]) # covariance matrix of the componets
cov_pcaA
chiMat<- matrix(NA,ncol=3,nrow=nrow(pr_scores))
chiMat[,1]<-mahalanobis(pr_scores[,1:5], mean_pcaA, cov_pcaA) # calculate mahalanobis distance

# Fit chi squared distribution and determine which rows can be removed based on the pcrit value (generally a pcrit of 0.975 is acceptable)
chiMat[,2]<- pchisq(c(chiMat[,1]), df =  5) # chi squared distribution with 5 degrees of freedom (df = number of components)
pcrit<- 1-((0.24-0.003*5)/sqrt(nrow(pr_scores))) # critical probability for outlier removal
pcrit
for (i in 1: nrow(chiMat)){
  if (chiMat[i,2] >= pcrit)
  {chiMat[i,3]<-0}else{chiMat[i,3]<-1}} # which spectra are outliers
par(mfrow=c(1,2))
plot(chiMat[,1],chiMat[,2], xlab= "distance", ylab="cumlative prob")
points(chiMat[which(chiMat[,3]==0),1:2],pch='X',col='red')
abline(v=qchisq(pcrit,5),col="green")
plot(pr_scores[,1], pr_scores[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])),ylim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])))
points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')

###Remove outliers from original dataset
original_data <-do.call(rbind,epo_DATA)
original_details <- data.frame(do.call(rbind,sampleg_details))

new_data <- original_data[chiMat[,3] == 1,]
new_details<- original_details[chiMat[,3] == 1,]
new_pr_scores<- pr_scores[chiMat[,3] == 1,]

rand.tr<-tri.mesh(new_pr_scores[,1],new_pr_scores[,2])
rand.ch<-convex.hull(rand.tr, plot.it=F) #convex hull
pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y))

newg_spectra <- split.data.frame(new_data,new_details$Sample)
newg_sampleg_details <- split.data.frame(new_details,new_details$Sample)

#####create the horizon confidence spheres with global database####

pr_spectra<-prcomp(do.call(rbind,newg_spectra), center=T,scale=T) 
pr_scores <- pr_spectra$x 

sample3d_details<-data.frame(do.call(rbind,newg_sampleg_details))

#####visualize the different horizons on pca space####
open3d(cex=1)
pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.01,
      group=sample3d_details$b.master_hor,
      show.plane=F,
      show.group.labels=T,
)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  c1<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.4),add=T,alpha=.12)
  #segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col=palette()[i])
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
  #wire3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.5,alpha=.7))
}


#####prepare the in_situ samples #####
spec_hv_cont <- as.data.frame(do.call(rbind,epo_hv_DATA))

hv_projected<-predict(pr_spectra,spec_hv_cont) #projected scores
sample3d_hv_details<-data.frame(do.call(rbind,sample_hv_details))

specMatch <- pnt.in.poly(hv_projected[,1:2],pr_poly)
specMatch
sum(specMatch$pip)/nrow(specMatch)*100 #

number_of_classes_hv<-length(levels(sample3d_hv_details$b.master_hor))
palette(terrain.colors(number_of_classes_hv))

plot3d(hv_projected[,1:3],
       col=replace(sample3d_hv_details$b.master_hor,palette(),c("A1","B2","C" )),
       type='s',
       radius=.5,
       group=sample3d_hv_details$b.master_hor,
       add=T)

number_of_classes<-length(levels(sample3d_hv_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  
  c1<- mean(hv_projected[sample3d_hv_details$b.master_hor==levels(sample3d_hv_details$b.master_hor)[i],1])
  c2<- mean(hv_projected[sample3d_hv_details$b.master_hor==levels(sample3d_hv_details$b.master_hor)[i],2])
  c3<- mean(hv_projected[sample3d_hv_details$b.master_hor==levels(sample3d_hv_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(hv_projected[sample3d_hv_details$b.master_hor==levels(sample3d_hv_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.60),add=T,alpha=.1)
}


###plot the effect of EPO ###

open3d(cex=1)
pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.1,
      group=sample3d_details$b.master_hor,
      show.plane=F,
      show.group.labels=F,
)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
plot3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
       type='s',
       radius=.4,
       group=sample3d_details$b.master_hor,
       show.plane=F,
       show.group.labels=T,
       lit=F,
       add=T,
       
)

plot3d(hv_projected[,1:3],
       col='black',
       type='s',
       radius=.4,
       group=sample3d_hv_details$b.master_hor,
       lit=F,
       add=T)
