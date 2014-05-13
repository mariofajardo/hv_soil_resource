###Soil description algorithm###
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


prev_dir <-getwd()
setwd('RData/')
load('huntervalley_2cm_ground.RDATA')
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('pr_varExp.RData')
load('check_plots.RData')
load('EPO_transformation_matrix.RDATA')

setwd(prev_dir)

####Load in Dataset A (laboratory dataset) and remove outliers####

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

# check_plots(rawg_spectra,'Raw Spectra','Reflectance')
# check_plots(trim_DATA,'Trimmed Spectra','Reflectance')
# check_plots(abs_DATA,'Abs Spectra','Reflectance')
# check_plots(filt_DATA,'Filt Spectra','Reflectance')
# check_plots(epo_DATA,'Processed Spectra','Epo_Units')

pr_spectra<-prcomp(do.call(rbind,epo_DATA), center=T,scale=T) 
screeplot(pr_spectra)#visualize the PC
pr_varExp(do.call(rbind,epo_DATA))#check the acummulative variation on the data
pr_scores <- pr_spectra$x 

#####convex.hull analysis#####
rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2])
rand.ch <- convex.hull(rand_tr,plot.it=F)
pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))

plot(pr_scores[,1],pr_scores[,2])
lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
#####outlier detection #####
mean_pcaA<- (colMeans(pr_scores[,1:5])) # mean of each of the components
mean_pcaA
cov_pcaA<- cov(pr_scores[,1:5])         # covariance matrix of the componets
cov_pcaA
chiMat<- matrix(NA,ncol=3,nrow=nrow(pr_scores)) 
chiMat[,1]<-mahalanobis(pr_scores[,1:5], mean_pcaA, cov_pcaA) # calculate mahalanobis distances

# Fit chi squared distribution and determine which rows can be removed based on the pcrit value (generally a pcrit of 0.975 is acceptable)
chiMat[,2]<- pchisq(c(chiMat[,1]), df =  5)       # chi squared distribution with 5 degrees of freedom (df = number of components)
pcrit<- 1-((0.24-0.003*5)/sqrt(nrow(pr_scores)))  # critical probability for outlier removal
pcrit
for (i in 1: nrow(chiMat)){
  if (chiMat[i,2] >= pcrit)
  {chiMat[i,3]<-0}else{chiMat[i,3]<-1}}           # which spectra are outliers
par(mfrow=c(1,2))
plot(chiMat[,1],chiMat[,2], xlab= "distance", ylab="cumlative prob")
points(chiMat[which(chiMat[,3]==0),1:2],pch='X',col='red')
abline(v=qchisq(pcrit,5),col="green")
plot(pr_scores[,1], pr_scores[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])),ylim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])))
points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')

####Remove outliers from original dataset####
original_data <-do.call(rbind,epo_DATA)
original_details <- data.frame(do.call(rbind,sampleg_details))

new_data <- original_data[chiMat[,3] == 1,]
new_details<- original_details[chiMat[,3] == 1,]
new_pr_scores<- pr_scores[chiMat[,3] == 1,]


####exclude the samples with too more than 10 outliers (20% of the the observations on a soil profile of 1m)####
count <- table(original_details[chiMat[,3] == 0,]$Sample)
exclude <- names(count)[count > 10]
no_out_details <- new_details[!(new_details$Sample %in% exclude),]
no_out_details$Sample <- droplevels(no_out_details$Sample)

######Dataset A filtered and transformed without outliers####
no_out_ground_DATA <- split.data.frame(cbind(no_out_details,new_data[!(new_details$Sample %in% exclude),]),no_out_details$Sample,drop=T)

length(epo_DATA)              # original samples
length(no_out_ground_DATA)    # no outliers samples

#Data preparation
no_out_data<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[7:ncol(x)]
  colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
  spectra})

no_out_details<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

no_out_pr_spectra<-prcomp(do.call(rbind,no_out_data), center=T,scale=T) 
no_out_pr_Exp <- pr_varExp(do.call(rbind,no_out_data))
no_out_pr_scores <- no_out_pr_spectra$x 

####convex hull ploygon####
rand.tr<-tri.mesh(no_out_pr_scores[,1],no_out_pr_scores[,2])
rand.ch<-convex.hull(rand.tr, plot.it=F) 
pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y))


####Visualize Dataset A without the outliers#### 


sample_details_pc<-data.frame(do.call(rbind,no_out_details))$b.master_hor
#####visualize the different horizons on pca space####

data <- data.frame(Soil_Types=sample_details_pc, no_out_pr_scores)
ggplot(data, aes_string(x='PC1', y='PC2', col='Soil_Types')) + geom_point(alpha=.9,size=5) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),title=paste0('Cumulative explanation :',round(no_out_pr_Exp[2,2],2),'%'))

####create confidence ellipses####
conf_ellipse <- data.frame()
for(i in levels(data$Soil_Types)){
  conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(data[data$Soil_Types==i,], 
                                                               ellipse(cor(PC1,PC2), 
                                                                       scale=c(sd(PC1),sd(PC2)), 
                                                                       centre=c(mean(PC1),mean(PC2)),
                                                                       level=.95))),Soil_Types=i))
}
#####plot confidence ellipses####
ggplot(data, aes_string(x='PC1', y='PC2', colour='Soil_Types')) + geom_point(alpha=.6,size=6) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
       y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
       title=paste0('Cumulative explanation :',
                    round(no_out_pr_Exp[2,2],2),'%')) +
  geom_path(data=conf_ellipse, aes_string(x='x', y='y',colour='Soil_Types'), size=1, linetype=1)


#####ploting convex hulls####

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(data,'Soil_Types', find_hull)

ggplot(data, aes_string(x='PC1', y='PC2', colour='Soil_Types')) + geom_point(alpha=.6,size=6) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
       y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
       title=paste0('Cumulative explanation :',
                    round(no_out_pr_Exp[2,2],2),'%')) +
  geom_polygon(data = hulls, alpha = 0.2)



####and in 3d pca space###
sample3d_details <- do.call(rbind,no_out_details)
open3d(cex=1)
pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.2,
      group=sample3d_details$b.master_hor,
      show.plane=F,
      show.group.labels=T,
)

#confidence spheres 40% conf#
number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  c1<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.4),add=T,alpha=.12)
  #segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col=palette()[i])
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
  #wire3d(ellipse3d(cov(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.5,alpha=.7))
}


#####Bring in the fuzzy cluster of Dataset A and do some processing#####
source('fuzzy_single_sample.R')

num_clusters <-3  #extract the number of clusters that we want to use from 1 to 6 possibles
fuzzy_data <-list()

for (i in 1:length(no_out_ground_DATA)){
  x<-i
  
  fuzzy_data[[i]]<- fanny_data_by_sample_pc_euc_no_out[[x]][[num_clusters]]
}

length(fuzzy_data)   #one fuzzy object by sample
# str(fuzzy_data[[1]]) #structure of the fuzzy object

fanny_by_sample_cont <- unlist(sapply(fuzzy_data,function(x) x$clustering))

#####Create the Dataset A dataframes with assigned clusters and details####
spec_cluster_cont <- as.data.frame(do.call(rbind,no_out_data))
sample3d_cluster_details<-data.frame(do.call(rbind,no_out_details),cluster=fanny_by_sample_cont)

# #####Save dataset without outlier samples and epo transformation for use later####
# no_out_epo_clus_DATA <- split.data.frame(cbind(sample3d_cluster_details,spec_cluster_cont),sample3d_cluster_details$Sample,drop=T)
# # save(no_out_epo_clus_DATA,file='no_out_epo_clus_DATA.RData')
# 
# #####Visualize Dataset A clusters on a  3PCA space#### 
palette(terrain.colors(num_clusters))
open3d(cex=1)
pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_cluster_details$cluster),
      radius=.1,
      group=sample3d_cluster_details$cluster,
      show.plane=F,
      show.group.labels=T)

for (i in 1:num_clusters){
  
  c1<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,1])
  c2<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,2])
  c3<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_cluster_details$cluster==i,1:3]),col=palette()[i],centre=centremean,level=.4),alpha=.12,add=T)
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
}

#####Visualize the same clusters with the traditional horizon dessignations#### 
# open3d(cex=1)
number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))

pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_cluster_details$b.master_hor),
      radius=.1,
      group=sample3d_cluster_details$b.master_hor,
      show.plane=F,
      show.group.labels=T)

palette(terrain.colors(num_clusters))

for (i in 1:num_clusters){
  
  c1<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,1])
  c2<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,2])
  c3<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_cluster_details$cluster==i,1:3]),col=palette()[i],centre=centremean,level=.7),alpha=.12,add=T)
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
}

#####evaluation in percentage of cluster and master horizons#####
#or how much of the traditional horizons is on each cluster

test_clas_details <-sample3d_cluster_details

par(mfcol=c(1,3))
clus_1<-test_clas_details[test_clas_details$cluster==1,] 
x<-clus_1$b.master_hor
clus_1_result<-nrow(clus_1[x=='A1' | x=='A2e' | x=='A2',])*100/nrow(clus_1) #percentage of 'A' horizons present in cluster 1

clus_2<-test_clas_details[test_clas_details$cluster==2,] 
y<-clus_2$b.master_hor
clus_2_result<-nrow(clus_2[y=='B2' | y=='B2w'| y=='B1',])*100/nrow(clus_2) #percentage of 'B' horizons present in cluster 1

clus_3<-test_clas_details[test_clas_details$cluster==3,] 
z<-clus_3$b.master_hor
clus_3_result<-nrow(clus_3[z=='B3' | z=='C',])*100/nrow(clus_3) #percentage of 'C' horizons present in cluster 1

cluster_perc_of_dataset<- c(nrow(clus_1)*100/nrow(test_clas_details),
  nrow(clus_2)*100/nrow(test_clas_details),
  nrow(clus_3)*100/nrow(test_clas_details))
cluster_perc_of_dataset
sum(cluster_perc_of_dataset)


abc_horizon_perc_on_cluster<-c(clus_1_result,clus_2_result,clus_3_result) #vector of percentages of each traditional horizon on each cluster
abc_horizon_perc_on_cluster

#####Visualize percentages of each traditional horizon on each cluster####
histogram(~test_clas_details$b.master_hor|test_clas_details$cluster)



#####filter Dataset B (in_situ spectra) and EPO tronsform it####
raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

# check_plots(raw_hv_spectra,'Raw Spectra_hv','Reflectance')
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
# check_plots(epo_hv_DATA,'EPO Spectra_hv','EPO units')


#####project samples in pca space  #####
spec_hv_cont <- as.data.frame(do.call(rbind,epo_hv_DATA))

hv_projected<-predict(no_out_pr_spectra,spec_hv_cont)      #projected scores
sample3d_hv_details<-data.frame(do.call(rbind,sample_hv_details))


#####check what percentage of Datset B samples are on Dataset A convex hull polygon####
specMatch <- pnt.in.poly(hv_projected[,1:2],pr_poly)
specMatch
sum(specMatch$pip)/nrow(specMatch)*100 ###percentage of points that are inside the convex hull


#####Visualize Dataset B horizons on Dataset A on a  3PCA space#### 
####With traditional horizons####
number_of_classes_hv<-length(levels(sample3d_hv_details$b.master_hor))
palette(terrain.colors(number_of_classes_hv))

 open3d(cex=1)
pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.01,
      group=sample3d_details$b.master_hor,
      show.plane=F,
      show.group.labels=T,
)

plot3d(hv_projected[,1:3],
       col=mapvalues(sample3d_hv_details$b.master_hor,from=levels(sample3d_hv_details$b.master_hor),to=palette()),
       type='s',
       radius=.5,
       group=sample3d_hv_details$b.master_hor,
       lit=F,
       add=T)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  c1<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.40),add=T,alpha=.12)
  #segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col=palette()[i])
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
  #wire3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.5,alpha=.7))
}

####Visualize Samples of Dataset B by horizon on clusters of Dataset A####
number_of_classes_hv<-length(levels(sample3d_cluster_details$cluster))
palette(terrain.colors(number_of_classes_hv))

#open3d(cex=1)
pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_cluster_details$cluster),
      radius=.01,
      group=sample3d_cluster_details$cluster,
      show.plane=F,
      show.group.labels=T,
)

number_of_classes_hv<-length(levels(sample3d_hv_details$b.master_hor))
palette(terrain.colors(number_of_classes_hv))

plot3d(hv_projected[,1:3],
       col=mapvalues(sample3d_hv_details$b.master_hor,from=levels(sample3d_hv_details$b.master_hor),to=palette()),
       type='s',
       radius=.5,
       group=sample3d_hv_details$b.master_hor,
       lit=F,
       add=T)

palette(terrain.colors(num_clusters))

for (i in 1:num_clusters){
  
  c1<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,1])
  c2<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,2])
  c3<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_cluster_details$cluster==i,1:3]),col=palette()[i],centre=centremean,level=.7),alpha=.12,add=T)
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
}

######Bring in the fuzzy cluster of Dataset B and do some processing####
source('fuzzy_single_sample_in_situ.R')

for (j in 1:length(epo_hv_DATA)){  
  
  num_clusters <-3
  
  fuzzy_hv_data <-list()
  
  for (i in 1:length(epo_hv_DATA)){
    fuzzy_hv_data[[i]]<- fanny_pits_pc_euc_EPO[[i]][[num_clusters]]
  }

}

fanny_by_sample_cont_hv <- unlist(sapply(fuzzy_hv_data,function(x) x$clustering))

sample3d_cluster_details_hv<-data.frame(sample3d_hv_details,cluster=fanny_by_sample_cont_hv)

#####Visualize Samples of Dataset B by cluster on clusters of Dataset A#####
palette(terrain.colors(num_clusters))
open3d(cex=1)
pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_cluster_details$cluster),
      radius=.2,
      group=sample3d_cluster_details$cluster,
      show.plane=F,
      show.group.labels=T)

for (i in 1:num_clusters){
  
  c1<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,1])
  c2<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,2])
  c3<- mean(no_out_pr_scores[sample3d_cluster_details$cluster==i,3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_cluster_details$cluster==i,1:3]),col=palette()[i],centre=centremean,level=.95),alpha=.12,add=T)
  segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
}

plot3d(hv_projected[,1:3],
       col=replace(sample3d_cluster_details_hv$cluster,palette(),levels(as.factor(sample3d_cluster_details_hv$cluster))),
       type='s',
       radius=.3,
       group=sample3d_cluster_details_hv$b.master_hor,
       add=T,
       lit=F)

####the same but with a convex hull####
####cluster1###
par(mfrow=c(1,3))
cluster<- 1
rand_tr_hv <-tri.mesh(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
                      no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2])
rand.ch_hv <- convex.hull(rand_tr_hv,plot.it=F)
pr_poly_hv <-cbind(x=c(rand.ch_hv$y),y=c(rand.ch_hv$y))

specMatch <- pnt.in.poly(hv_projected[sample3d_cluster_details_hv$cluster==cluster,1:2],pr_poly)
sum_points_in<-round(sum(specMatch$pip)/nrow(specMatch)*100,2) ##percentage of samples on the convex hull
sum_points_in

plot(no_out_pr_scores[,1],no_out_pr_scores[,2],col='white',
     main=paste0('% of points in convex hull','=',sum_points_in,'%'),
     xlab=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'),
     ylab=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'))

points(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
       no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2],col='gray') 
points(hv_projected[sample3d_cluster_details_hv$cluster==cluster,],col='red')           #Samples of Dataset B cluster 1
lines(c(rand.ch_hv$x,rand.ch_hv$x[1]),c(rand.ch_hv$y,rand.ch_hv$y[1]),col='blue') #convex hull of Dataset A cluster 1


####cluster 2###
cluster<- 2
rand_tr_hv <-tri.mesh(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
                      no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2])
rand.ch_hv <- convex.hull(rand_tr_hv,plot.it=F)
pr_poly_hv <-cbind(x=c(rand.ch_hv$y),y=c(rand.ch_hv$y))

specMatch <- pnt.in.poly(hv_projected[sample3d_cluster_details_hv$cluster==cluster,1:2],pr_poly)
sum_points_in<-round(sum(specMatch$pip)/nrow(specMatch)*100,2) ##percentage of samples on the convex hull
sum_points_in

plot(no_out_pr_scores[,1],no_out_pr_scores[,2],col='white',
     main=paste0('% of points in convex hull','=',sum_points_in,'%'),
     xlab=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'),
     ylab=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'))

points(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
       no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2],col='gray') 
points(hv_projected[sample3d_cluster_details_hv$cluster==cluster,],col='red')           #Samples of Dataset B cluster 1
lines(c(rand.ch_hv$x,rand.ch_hv$x[1]),c(rand.ch_hv$y,rand.ch_hv$y[1]),col='blue') #convex hull of Dataset A cluster 1


####cluster 3###
cluster<- 3
rand_tr_hv <-tri.mesh(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
                      no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2])
rand.ch_hv <- convex.hull(rand_tr_hv,plot.it=F)
pr_poly_hv <-cbind(x=c(rand.ch_hv$y),y=c(rand.ch_hv$y))

specMatch <- pnt.in.poly(hv_projected[sample3d_cluster_details_hv$cluster==cluster,1:2],pr_poly)
sum_points_in<-round(sum(specMatch$pip)/nrow(specMatch)*100,2) ##percentage of samples on the convex hull
sum_points_in

plot(no_out_pr_scores[,1],no_out_pr_scores[,2],col='white',
     main=paste0('% of points in convex hull','=',sum_points_in,'%'),
     xlab=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'),
     ylab=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'))

points(no_out_pr_scores[sample3d_cluster_details$cluster==cluster,1],
       no_out_pr_scores[sample3d_cluster_details$cluster==cluster,2],col='gray') 
points(hv_projected[sample3d_cluster_details_hv$cluster==cluster,],col='red')           #Samples of Dataset B cluster 1
lines(c(rand.ch_hv$x,rand.ch_hv$x[1]),c(rand.ch_hv$y,rand.ch_hv$y[1]),col='blue') #convex hull of Dataset A cluster 1


######evaluation in percentage of cluster and master horizons######
#or how much of the traditional horizons is on each cluster


test_clas_details_hv <-sample3d_cluster_details_hv

par(mfcol=c(1,3))
clus_1_hv<-test_clas_details_hv[test_clas_details_hv$cluster==1,] 
x<-clus_1_hv$b.master_hor
clus_1_result_hv<-nrow(clus_1_hv[x=='A1' | x=='A2e' | x=='A2',])*100/nrow(clus_1_hv) #percentage of 'A' horizons present in cluster 1

clus_2_hv<-test_clas_details_hv[test_clas_details_hv$cluster==2,] 
y<-clus_2_hv$b.master_hor
clus_2_result_hv<-nrow(clus_2_hv[y=='B2' | y=='B2w'| y=='B1',])*100/nrow(clus_2_hv) #percentage of 'B' horizons present in cluster 1

clus_3_hv<-test_clas_details_hv[test_clas_details_hv$cluster==3,] 
z<-clus_3_hv$b.master_hor
clus_3_result_hv<-nrow(clus_3_hv[z=='B3' | z=='C',])*100/nrow(clus_3_hv) #percentage of 'C' horizons present in cluster 1

cluster_perc_of_dataset_hv<- c(nrow(clus_1_hv)*100/nrow(test_clas_details_hv),
                            nrow(clus_2_hv)*100/nrow(test_clas_details_hv),
                            nrow(clus_3_hv)*100/nrow(test_clas_details_hv))
cluster_perc_of_dataset_hv
sum(cluster_perc_of_dataset_hv)


abc_horizon_perc_on_cluster_hv<-c(clus_1_result_hv,clus_2_result_hv,clus_3_result_hv) #vector of percentages of each traditional horizon on each cluster
abc_horizon_perc_on_cluster_hv

#####Visualize percentages of each traditional horizon on each cluster####
histogram(~test_clas_details_hv$b.master_hor|test_clas_details_hv$cluster)

####Horizon predictions####

####Dataset A###
setwd('plots/')
pdf('observed_vs_predicted_horizons_3clust.pdf',width=10,height=6)
input_data <-no_out_ground_DATA
par(mfrow=c(1,2))
for (j in 1:length(input_data)){  
num_clusters <-3 ##select the decided number of clusters (up to 6)
fuzzy_data <-list()
for (i in 1:length(no_out_ground_DATA)){
  fuzzy_data[[i]]<- fanny_data_by_sample_pc_euc_no_out[[i]][[num_clusters]]
}
z<-cbind(input_data[[j]][,1:6],cluster=fuzzy_data[[j]]$clustering)
z1 <- z2 <- z

z1_colors <- factor(z$b.master_hor)
levels(z1_colors) <- terrain.colors(length(unique(z1$b.master_hor)))
z1$soil_color <- as.character(z1_colors)
z1$Sample <- paste0(z1$Sample,'_Horizon')

z2_colors <- factor(z$cluster)
levels(z2_colors) <- terrain.colors(length(unique(z$cluster)))
z2$soil_color <- as.character(z2_colors)
z2$Sample <- paste0(z1$Sample,'_Cluster')

z_ <- rbind(z1,z2)
depths(z_)<- Sample ~ top + bottom
plot(z_)
}
dev.off()
shell.exec('observed_vs_predicted_horizons_3clust.pdf')
 
setwd(prev_dir)

##note : close the previous pdf to create the next one##

####Dataset B###
setwd('..//plots')
par(mfrow=c(1,1))
pdf('observed_vs_predicted_horizons_in_situ_3clust.pdf',width=10,height=6)
input_data <-hv_pits_DATA

 for (j in 1:length(input_data)){  
  num_clusters <-3 ##select the decided number of clusters (up to 6)
  fuzzy_hv_data <-list()
  for (i in 1:length(input_data)){
    fuzzy_hv_data[[i]]<- fanny_pits_pc_euc_EPO[[i]][[num_clusters]]
  }
  z<-cbind(input_data[[j]][,1:6],cluster=fuzzy_hv_data[[j]]$clustering)
  z1 <- z2 <- z

  z1_colors <- factor(z$b.master_hor)
  levels(z1_colors) <- terrain.colors(length(unique(z1$b.master_hor)))
  z1$soil_color <- as.character(z1_colors)
  z1$Sample <- paste0(z1$Sample,'_Horizon')
  
  z2_colors <- factor(z$cluster)
  levels(z2_colors) <- terrain.colors(length(unique(z$cluster)))
  z2$soil_color <- as.character(z2_colors)
  z2$Sample <- paste0(z1$Sample,'_Cluster')

  z_ <- rbind(z1,z2)
  depths(z_)<- Sample ~ top + bottom
  
  plot(z_)
  
 }
dev.off()
shell.exec('observed_vs_predicted_horizons_in_situ_3clust.pdf')

setwd(prev_dir)

####plotmemberships####
par(mfrow=c(2,1))

a<-fanny_data_by_sample_pc_euc_no_out[[1]][[3]] #sample 3 with 4 clusters
str(a)

plot(a$membership[,1],type='l')
lines(a$membership[,2],type='l',col='blue')
lines(a$membership[,3],type='l',col='red')

a<-fanny_data_by_sample_pc_euc_no_out[[1]][[2]] #sample 3 with 4 clusters

str(a)

plot(a$membership[,1],type='l')
lines(a$membership[,2],type='l',col='blue')
# lines(a$membership[,3],type='l',col='red')
# lines(a$membership[,4],type='l',col='green')


####plotmemberships pits####
par(mfrow=c(2,1))

a<-fanny_pits_pc_euc_EPO[[8]][[2]] #sample 3 with 4 clusters
str(a)

plot(a$membership[,1],type='l')
lines(a$membership[,2],type='l',col='blue')
lines(a$membership[,3],type='l',col='red')

a<-fanny_pits_pc_euc_EPO[[8]][[4]] #sample 3 with 4 clusters

str(a)

plot(a$membership[,1],type='l')
lines(a$membership[,2],type='l',col='blue')
lines(a$membership[,3],type='l',col='red')
lines(a$membership[,4],type='l',col='green')
par(mfrow=c(1,1))
####PREDICTIONS####
require(Cubist)
require(naturalsort)

#####Bring the Carbon data####
setwd('RData/')
load('lhc_samples.RDATA')
load('selected_rnd_samples_HVHR.RDATA')
setwd(prev_dir)
setwd('../../CNS')
lhc_calib_tmp <-read.csv('lhc_last_edition.csv')
lhc_calib_tmp <-data.frame(instrument_id=lhc_calib_tmp$Name,TOC=lhc_calib_tmp$Percent2)
lhc_calib_tmp <- lhc_calib_tmp[naturalorder(as.character(lhc_calib_tmp$instrument_id)),][1:100,]

rnd_valid_tmp <-read.csv('rnd_calibration_points.csv')
rnd_valid_tmp <-data.frame(instrument_id=rnd_valid_tmp$Name,TOC=rnd_valid_tmp$Percent2)
rnd_valid_tmp <- rnd_valid_tmp[naturalorder(as.character(rnd_valid_tmp$instrument_id)),][1:100,]
setwd(prev_dir)

data_a_details <- sample3d_cluster_details
data_a_spec<-spec_cluster_cont

#Calibration set
spec_c <-data_a_spec[data_a_details$File.Name%in%lhc_selected_details$File.Name,]            # spectra used for calibration
colnames(spec_c)<-seq(500,2450,by=10)
details<-data_a_details[data_a_details$File.Name%in%lhc_selected_details$File.Name,]
details$cluster<-as.factor(details$cluster)
soil_c <-lhc_calib_tmp[lhc_selected_details$File.Name%in%details$File.Name,]$TOC         # data used for calibration

#Validation set
setwd('../../../Hunter/pits_2013/')
pits_analysis <-read.csv(file='analysis.csv',dec='.')
pits_analysis$sample<-factor(pits_analysis$sample,levels=unique(pits_analysis$sample))
####fitting splines to horizon observations####
# Start
source("C:/Users/mfaj1435/Documents/University of Sydney/PhD/Toolbox_scripts_etc/ea_spline_f.r")

var <- 'TOC'
lam <- 0.5  # Lambda value
d <- t(seq(0,100,by=5)) # GlobalSoilMap.Net specifications or user defined depths
mxd <- 100 # max depth for spline
s <- 0.05*mean(pits_analysis$TOC)  # 5% of the standard deviation of the target attribute 

splines<- lapply(levels(pits_analysis$sample),function(x){
  ans <- ea_spline(mxd,pits_analysis[pits_analysis$sample==x,],var,lam,d,s)
})
names(splines) <-levels(pits_analysis$sample)

soil_v<-do.call(rbind,lapply(levels(pits_analysis$sample),function(x){
  rows<-nrow(pits_analysis[pits_analysis$sample==x,])
  tmp<-data.frame(sample=x,
                  top=seq(from=0,by=5,length.out=20),
                  bottom=seq(from=5,by=5,length.out=20),
                  TOC=t(splines[[x]][[2]]))
  }))

soil_v<-na.exclude(soil_v)

soil_v$cluster <- do.call(c,lapply(levels(pits_analysis$sample),function(x) sample3d_cluster_details_hv$cluster[sample3d_cluster_details_hv$Sample==x]))
soil_v$pred_TOC <-NA
                                                                        


setwd(prev_dir)


soil_c.cubist_models <-list()
soil_c.cubist_models<-lapply(levels(details$cluster),function(x) cubist(spec_c[details$cluster==x,],
                                                                        soil_c[details$cluster==x],
                                                                        control=cubistControl(rules=2),
                                                                        committees=100))     # fit cubist model
names(soil_c.cubist_models)<-paste0('cluster',levels(details$cluster))



for (n in levels(details$cluster)){
soil_v$pred_TOC[soil_v$cluster==n] <- predict(soil_c.cubist_models[[paste0('cluster',n)]],spec_hv_cont[sample3d_cluster_details_hv$cluster==n,])
}
test<-goof(soil_v$TOC,soil_v$pred_TOC,col=soil_v$cluster)

global_model<-cubist(spec_c,soil_c,rules=2,committees=100)
soil_v$pred_global_TOC <-predict(global_model,spec_hv_cont)

test2<-goof(soil_v$TOC,soil_v$pred_global_TOC,col=soil_v$cluster)

test
test2
par(mfrow=c(2,4))
data_validation <- list ()
for (n in levels(soil_v$sample)){
  
plot(soil_v$TOC[soil_v$sample==n],
     seq(from=5,by=5,length.out=nrow(soil_v[soil_v$sample==n,])),
     type='l',
     ylim=c(nrow(soil_v[soil_v$sample==n,])*5,0),
     xlim=c(0,10),main=n,
     ylab='Depth',
     xlab='TOC (%)')

lines(soil_v$pred_TOC[soil_v$sample==n],
      seq(from=5,by=5,length.out=nrow(soil_v[soil_v$sample==n,])),
      type='l',
      col='blue',
      ylim=c(nrow(soil_v[soil_v$sample==n,])*5,0),
      xlim=c(0,10))

points(pits_analysis$TOC[pits_analysis$sample==n],
       rowMeans(cbind(pits_analysis$top[pits_analysis$sample==n],pits_analysis$bottom[pits_analysis$sample==n])),
       col='red',
       ylim=c(nrow(soil_v[soil_v$sample==n,])*5,0),
       xlim=c(0,10))

lines(soil_v$pred_global_TOC[soil_v$sample==n],
      seq(from=5,by=5,length.out=nrow(soil_v[soil_v$sample==n,])),
      type='l',
      col='green',
      ylim=c(nrow(soil_v[soil_v$sample==n,])*5,0),
      xlim=c(0,10))

points(x=mapply(function(x,y) mean(soil_v$pred_TOC[soil_v$sample==n][which(soil_v$top[soil_v$sample==n] %in% x):
                                                   which(soil_v$bottom[soil_v$sample==n] %in% y)]),
                pits_analysis$top[pits_analysis$sample==n],
                pits_analysis$bottom[pits_analysis$sample==n],
                SIMPLIFY=T),
       y=rowMeans(cbind(pits_analysis$top[pits_analysis$sample==n],
                        pits_analysis$bottom[pits_analysis$sample==n]))
       )

data_validation[[n]] <- data.frame(sample=n,
                                   obs=pits_analysis$TOC[pits_analysis$sample==n],
                                   pred=mapply(function(x,y) mean(soil_v$pred_TOC[soil_v$sample==n][which(soil_v$top[soil_v$sample==n] %in% x):
                                                                                                 which(soil_v$bottom[soil_v$sample==n] %in% y)]),
                                            pits_analysis$top[pits_analysis$sample==n],
                                            pits_analysis$bottom[pits_analysis$sample==n],
                                            SIMPLIFY=T),
                                   cluster= sapply(pits_analysis$top[pits_analysis$sample==n],function(x){
                                               soil_v$cluster[soil_v$sample==n][which(soil_v$top[soil_v$sample==n] %in% x)]})
                              )

}

data_validation<-do.call(rbind,data_validation)

pit<-7 

goof(data_validation$obs[data_validation$sample==levels(data_validation$sample)[pit]],data_validation$pred[data_validation$sample==levels(data_validation$sample)[pit]],col=data_validation$cluster)
###end###

