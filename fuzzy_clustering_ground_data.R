#####Fuzzy Clustering_FANNY_Principal_components_whole_spectra####
library(spectroscopy)
library(pbapply)
library(tripack)
library(SDMTools)
library(cluster)
library(multicore)

prev_dir <-getwd()
setwd('RData/')
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('pr_varExp.RData')
load('EPO_transformation_matrix.RDATA')
#bring in and remove outliers from global dataset#
setwd(prev_dir)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(500:2300),which=2))
snv_DATA <- pblapply(strip_DATA,snvBLC)
# epo_DATA <- pblapply(snv_DATA,function(x) {
#   epo_tmp<-as.matrix(x) %*% P
#   epo_tmp<-data.frame(epo_tmp)
#   colnames(epo_tmp)<-seq(500,2300,by =2)
#   epo_tmp}
# )

epo_DATA<-snv_DATA


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

pca_data <-do.call(rbind,strip_DATA)

test_data <- round(as.data.frame(new_pr_scores[,1:26]),12)

num_clusters <- 1:8

fanny_data_pca_no_out <- mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean',maxit=200)
})

save(file='fanny_data_pca_no_out.RData',fanny_data_pca_no_out)

#####Fuzzy cluster by sample ####
library(spectroscopy)
library(pbapply)
library(tripack)
library(SDMTools)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('pr_varExp.RData')
#bring in and remove outliers from global dataset#

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


pr_spectra<-prcomp(do.call(rbind,snv_DATA), center=T,scale=T) 
screeplot(pr_spectra)#visualize the PC
pr_varExp(do.call(rbind,snv_DATA))#check the acummulative variation on the data
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
original_data <-do.call(rbind,snv_DATA)
original_details <- data.frame(do.call(rbind,sampleg_details))

new_data <- original_data[chiMat[,3] == 1,]
new_details<- original_details[chiMat[,3] == 1,]

#condition
count <- table(original_details[chiMat[,3] == 0,]$Sample)
exclude <- names(count)[count > 5]
no_out_details <- new_details[!(new_details$Sample %in% exclude),]
no_out_details$Sample <- droplevels(no_out_details$Sample)

#new lists
no_out_DATA <- split.data.frame(new_data[!(new_details$Sample %in% exclude),],no_out_details$Sample,drop=T)
details_no_out_DATA <- split.data.frame(no_out_details,no_out_details$Sample,drop=T)


pr_spectra_DATA <-lapply(no_out_DATA,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$Sample,z$x[,1:5]),details_no_out_DATA,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=F) # principal component scores
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- 1:6
#   i="gku17"
fanny_data_by_sample_pc_euc_no_out_snv <- lapply(levels(pr_scores_DATAFRAME$group),function(i){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
  mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) fanny(x=test_data,k=y,metric='euclidean'))
})





save(file='fanny_data_by_sample_pc_euc_no_out_snv.RData',fanny_data_by_sample_pc_euc_no_out_snv)




#####Fuzzy data by samples on pca sapce without outliers####
#####Fuzzy cluster by sample ####
library(spectroscopy)
library(pbapply)
library(tripack)
library(SDMTools)
library(cluster)
library(multicore)

prev_dir <-getwd()
setwd('RData/')
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('pr_varExp.RData')
load('EPO_transformation_matrix.RDATA')

setwd<-prev_dir
#####First remove outliers from global dataset####
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
# epo_DATA <- pblapply(snv_DATA,function(x) {
#   epo_tmp<-as.matrix(x) %*% P
#   epo_tmp<-data.frame(epo_tmp)
#   colnames(epo_tmp)<-seq(500,2450,by =10)
#   epo_tmp}
# )

epo_DATA<-snv_DATA
pr_spectra<-prcomp(do.call(rbind,epo_DATA), center=T,scale=T) 
screeplot(pr_spectra)               #visualize the PC
pr_varExp(do.call(rbind,epo_DATA))  #check the acummulative variation on the data
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

#####Remove samples that have more than 5 outliers####
count <- table(original_details[chiMat[,3] == 0,]$Sample)
exclude <- names(count)[count > 5]
no_out_details <- new_details[!(new_details$Sample %in% exclude),]
no_out_details$Sample <- droplevels(no_out_details$Sample)

#####new dataset without outliers####
no_out_DATA <- split.data.frame(new_data[!(new_details$Sample %in% exclude),],no_out_details$Sample,drop=T)
details_no_out_DATA <- split.data.frame(no_out_details,no_out_details$Sample,drop=T)


pr_spectra_DATA <-lapply(no_out_DATA,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$Sample,z$x[,1:5]),details_no_out_DATA,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=F) # principal component scores
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- 1:6
#   i="gku17"
fanny_data_by_sample_pc_euc_no_out_SNV <- lapply(levels(pr_scores_DATAFRAME$group),function(i){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
  mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) fanny(x=test_data,k=y,metric='euclidean'))
})


save(file='fanny_data_by_sample_pc_euc_no_out_SNV.RData',fanny_data_by_sample_pc_euc_no_out_SNV)