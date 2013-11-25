###Canonical Variates###

library(tripack)
library(SDMTools)
library(rgl) 
library(spectroscopy)
library(pbapply)
library(pca3d)
library(ggplot2)
library(ellipse)
library(plyr)
library(aqp)
library(lattice)


load('huntervalley_2cm_ground.RDATA')
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('pr_varExp.RData')
load('check_plots.RData')
load('EPO_transformation_matrix.RDATA')


prev_dir <-getwd()

####Load in and remove outliers from global dataset####

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

check_plots(rawg_spectra,'Raw Spectra','Reflectance')
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

check_plots(epo_DATA,'Processed Spectra','Epo_Units')

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

#condition
count <- table(original_details[chiMat[,3] == 0,]$Sample)
exclude <- names(count)[count > 5]
no_out_details <- new_details[!(new_details$Sample %in% exclude),]
no_out_details$Sample <- droplevels(no_out_details$Sample)

#new list
no_out_ground_DATA <- split.data.frame(cbind(no_out_details,new_data[!(new_details$Sample %in% exclude),]),no_out_details$Sample,drop=T)

#####fuzzy classification####
load('fanny_data_by_sample_pc_euc_no_out.RData')
num_clusters <-3
fuzzy_data <-list()

for (i in 1:length(no_out_ground_DATA)){
  x<-i
  
  fuzzy_data[[i]]<- fanny_data_by_sample_pc_euc_no_out[[x]][[num_clusters]]
}


fanny_by_sample_cont <- unlist(sapply(fuzzy_data,function(x) x$clustering))

#####Discriminant analysis in jmp####

tmp <- data.frame(cluster=fanny_by_sample_cont,do.call(rbind,no_out_ground_DATA))


# cluster1 <- replace(tmp[,1],tmp[,1]!=1,0)
# cluster2 <- replace(tmp[,1],tmp[,1]!=2,0)
# cluster3 <- replace(tmp[,1],tmp[,1]!=3,0)

write.csv(tmp,file='cda_test.csv')


responseY <- as.matrix(tmp[,1])
predictorX <- as.matrix(tmp[,-c(1:7)])


pca <- princomp(predictorX, cor=T) # principal components analysis using correlation matrix
pc.comp <- pca$scores
pc.comp1 <- -1*pc.comp[,1] # principal component 1 scores (negated for convenience)
pc.comp2 <- -1*pc.comp[,2] # principal component 2 scores (negated for convenience)

library(MASS)
model.lda <- MASS::(responseY ~ pc.comp1+pc.comp2)
predict(model.lda)$class
plot(model.lda)

#####Cluster means####

tmp <- data.frame(cluster=fanny_by_sample_cont,do.call(rbind,no_out_ground_DATA))
tmp <- split.data.frame(tmp,tmp$Sample)
tmp1 <- lapply(tmp,function (x) {
  tmp2 <- split(x,x$cluster)
  tmp3  <- lapply(tmp2,function (y){
    tmp3_det  <- y[,1:7] 
    tmp3_data <- y[,8:length(y)]
    tmp_final <- data.frame(Sample=tmp3_det[1,3],cluster=tmp3_det[1,1],t(colMeans(tmp3_data)))
  })
  tmp4<-do.call(rbind,tmp3)
})
clus_mean<-do.call(rbind,tmp1)  

View(clus_mean)

write.csv(clus_mean,file='cda_mean_clus_test.csv')

