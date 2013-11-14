#library(plyr)
#Check_equal_variables <- function (x){
#  var <- 1:ncol(x)
#  pboptions(type='win',style=3,title='Checking equal variables')
#  pblapply(var,function(y){
#     browser()
#    var2 <-(y+1):ncol(x)
#    print(y)
#    if(length(var2)==1 | (y==ncol(x))){
#      NA
#    } else {
#      any(sapply(var2,function(z)identical(x[,y],x[,z])))
#    }
#  })
#}
#addNoise <- function(mtx) {
#  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
#  random.stuff <- matrix(runif(prod(dim(mtx)), min = 0.0001, max = 0.0009), nrow = dim(mtx)[1])
#  random.stuff + mtx
#}
#equal_values <- unlist(Check_equal_variables(do.call(rbind,filt_DATA)))
#data_fuzzy <- abs_DATA[,equal_values[1:(length(equal_values)-2)]]
###########Fuzzy_clustering_Extragrades#####################################
# library(spectroscopy)
# library(pbapply)
# load('mahaldist.RData')
# load('fuzzy_padarian.RData')
# load('run_multiple_fuzzy.RData')
# load('huntervalley_2cm_dry_cores.RDATA')
# load('correct_steps.RData')
# raw_spectra<-lapply(DATA,function(x) {
#   spectra<-x[7:2157]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_details<-lapply(DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# pboptions(type='txt',style=3)
# cont_DATA<-pblapply(raw_spectra,correct_step)    
# trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2300))) 
# abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
# strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(500:2300),which=10)) 
# filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# snVC_DATA<- pblapply(filt_DATA,snvBLC) 
# strip_DATA<- pblapply(snVC_DATA,function(x) strip_spectra(x,c(500:2300),which=2))
# 
# 
# 
# test_data <- round(do.call(rbind,strip_DATA),12)

#snVC_DATA_10_to_20_cluster<-run_multiple_fuzzy(test_data,phi=1.5,10:20,3,exp_eg=0.05)


# ###########Performance measurements#####################################
# ####FPI####
# 
# extract_FPI<-function(fuzzy_cluster_group) 
#  sapply(fuzzy_cluster_group@clusters,function(x){
#      memberships<-x@U
#      part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#      FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#    })
# 
# FPI<-c(extract_FPI(snVC_DATA_2_to_5_cluster),(extract_FPI(snVC_DATA_5_to_10_cluster)),(extract_FPI(snVC_DATA_10_to_20_cluster)))
# 
# plot(2:22,FPI,type='l')
# 
# abline(v=12)
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_cluster_group) 
#   sapply(fuzzy_cluster_group@clusters,function(x){
#     memberships<-x@U
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(snVC_DATA_2_to_5_cluster)),(extract_NCE(snVC_DATA_5_to_10_cluster)),(extract_NCE(snVC_DATA_10_to_20_cluster)))
# 
# plot(2:22,NCE,type='l',xlab='Number of fuzzy clusters')



#fuzzy_test<-fuzzy_extragrades(0.8,test_data,phi=1.5,3L,3,exp_eg=0.05,disttype=1)



#####Processing####
# setClass('FuzzyCluster',representation(data='matrix',U='matrix',W='matrix',centroids='matrix',phi='numeric',classes='integer',distance='character',
#                                        alpha='numeric',`Ue_mean - Ue_req`='numeric',iterations='integer'))
# setClass('FuzzyClusterGroup',representation(clusters='list'))
# 
# setMethod('length','FuzzyClusterGroup',function(x) length(x@clusters)) 
# 
# load('snVC_DATA_10_to_20_cluster.RData')
# 
# a<-snVC_DATA_10_to_20_cluster@clusters$'14'
# b<-attributes(a@U)$hard_clust

#####Fuzzy Clustering_FANNY####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_dry_cores.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(raw_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 


spec_cont <- as.data.frame(do.call(rbind,strip_DATA))
test_data <- spec_cont


num_clusters <- 1:36

fanny_data_whole_spectra <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

#save(file='fanny_data_whole_spectra.RData',fanny_data_whole_spectra)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#    sapply(fuzzy_data,function(x){
#         memberships<-x$membership
#         part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#         FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#       })
#    
# FPI<-c(extract_FPI(fanny_data_whole_spectra))
# plot(1:length(FPI),FPI,type='l')
# abline(v=8)
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#    sapply(fuzzy_data,function(x){
#      memberships<-x$membership
#      entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#      NCE<-entropy/log10(ncol(memberships)-1)
#    })
#  
# NCE<-c((extract_NCE(fanny_data_whole_spectra)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')
# abline(v=8)

#####Fuzzy Clustering_FANNY_Visible spectra####
library(cluster)
library(multicore)
library(pbapply)
library(spectroscopy)
load('huntervalley_2cm_dry_cores.RDATA')
load('correct_steps.RData')

raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(raw_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 


test_data <- round(do.call(rbind,strip_DATA),12)
num_clusters <- 1:36

fanny_data_visibe <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

#save(file='fanny_data_visibe.RData',fanny_data_visibe)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#     FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#   })
# 
# FPI<-c(extract_FPI(fanny_data_visibe))
# plot(1:length(FPI),FPI,type='l')
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(fanny_data_visibe)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')


#####Fuzzy Clustering_FANNY_Visible spectra and depth####
library(cluster)
library(multicore)
library(pbapply)
library(spectroscopy)
load('huntervalley_2cm_dry_cores.RDATA')
load('correct_steps.RData')

raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(raw_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 

spec_strip_depth <- as.data.frame(do.call(rbind,strip_DATA),depth=do.call(rbind,sample_details)$bottom)

test_data <- round(spec_strip_depth,12)
num_clusters <- 1:36

fanny_data_v_d <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

#save(file='fanny_data_v_d.RData',fanny_data_v_d)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#     FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#   })
# 
# FPI<-c(extract_FPI(fanny_data_v_d_2_1m))
# plot(1:length(FPI),FPI,type='l')
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(fanny_data_v_d_2_1m)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')


#####GROUND DATA####

library(spectroscopy)
library(pbapply)
load('mahaldist.RData')
load('fuzzy_padarian.RData')
load('run_multiple_fuzzy.RData')
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

pboptions(type='txt',style=3)
cont_DATA<-pblapply(raw_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(500:2300),which=10)) 
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
snVC_DATA<- pblapply(filt_DATA,snvBLC) 
strip_DATA<- pblapply(snVC_DATA,function(x) strip_spectra(x,c(500:2300),which=2))



test_data <- round(do.call(rbind,strip_DATA),12)

#snVC_DATA_10_to_20_cluster<-run_multiple_fuzzy(test_data,phi=1.5,10:20,3,exp_eg=0.05)


# ###########Performance measurements#####################################
# ####FPI####
# 
# extract_FPI<-function(fuzzy_cluster_group) 
#  sapply(fuzzy_cluster_group@clusters,function(x){
#      memberships<-x@U
#      part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#      FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#    })
# 
# FPI<-c(extract_FPI(snVC_DATA_2_to_5_cluster),(extract_FPI(snVC_DATA_5_to_10_cluster)),(extract_FPI(snVC_DATA_10_to_20_cluster)))
# 
# plot(2:22,FPI,type='l')
# 
# abline(v=12)
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_cluster_group) 
#   sapply(fuzzy_cluster_group@clusters,function(x){
#     memberships<-x@U
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(snVC_DATA_2_to_5_cluster)),(extract_NCE(snVC_DATA_5_to_10_cluster)),(extract_NCE(snVC_DATA_10_to_20_cluster)))
# 
# plot(2:22,NCE,type='l',xlab='Number of fuzzy clusters')



#fuzzy_test<-fuzzy_extragrades(0.8,test_data,phi=1.5,3L,3,exp_eg=0.05,disttype=1)



#####Processing####
# setClass('FuzzyCluster',representation(data='matrix',U='matrix',W='matrix',centroids='matrix',phi='numeric',classes='integer',distance='character',
#                                        alpha='numeric',`Ue_mean - Ue_req`='numeric',iterations='integer'))
# setClass('FuzzyClusterGroup',representation(clusters='list'))
# 
# setMethod('length','FuzzyClusterGroup',function(x) length(x@clusters)) 
# 
# load('snVC_DATA_10_to_20_cluster.RData')
# 
# a<-snVC_DATA_10_to_20_cluster@clusters$'14'
# b<-attributes(a@U)$hard_clust

#####Fuzzy Clustering_FANNY####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 


spec_cont <- as.data.frame(do.call(rbind,strip_DATA))
test_data <- spec_cont


num_clusters <- 1:36

fanny_data_g_whole_spectra <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

save(file='fanny_data_g_whole_spectra.RData',fanny_data_g_whole_spectra)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#    sapply(fuzzy_data,function(x){
#         memberships<-x$membership
#         part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#         FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#       })
#    
# FPI<-c(extract_FPI(fanny_data_whole_spectra))
# plot(1:length(FPI),FPI,type='l')
# abline(v=8)
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#    sapply(fuzzy_data,function(x){
#      memberships<-x$membership
#      entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#      NCE<-entropy/log10(ncol(memberships)-1)
#    })
#  
# NCE<-c((extract_NCE(fanny_data_whole_spectra)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')
# abline(v=8)

#####Fuzzy Clustering_FANNY_Visible spectra####
library(cluster)
library(multicore)
library(pbapply)
library(spectroscopy)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 


test_data <- round(do.call(rbind,strip_DATA),12)
num_clusters <- 1:36

fanny_data_g_visibe <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

save(file='fanny_data_g_visibe.RData',fanny_data_g_visibe)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#     FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#   })
# 
# FPI<-c(extract_FPI(fanny_data_visibe))
# plot(1:length(FPI),FPI,type='l')
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(fanny_data_visibe)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')


#####Fuzzy Clustering_FANNY_Visible spectra and depth####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 

spec_strip_depth <- as.data.frame(do.call(rbind,strip_DATA),depth=do.call(rbind,sampleg_details)$bottom)

test_data <- round(spec_strip_depth,12)
num_clusters <- 1:36

fanny_data_g_v_d <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

save(file='fanny_data_g_v_d.RData',fanny_data_g_v_d)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#     FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#   })
# 
# FPI<-c(extract_FPI(fanny_data_v_d_2_1m))
# plot(1:length(FPI),FPI,type='l')
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(fanny_data_v_d_2_1m)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')



#####Fuzzy Clustering_FANNY_Principal_components_visible spectra####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 
pca_data <-do.call(rbind,strip_DATA)

pr_spectra<-prcomp(pca_data, center=T,scale=T) 
pr_Exp<- pr_varExp(pca_data) # amount of vartiation explained by first 10 PCs
pr_scores <- pr_spectra$x # principal component scores

test_data <- pr_scores

num_clusters <- 1:12

fanny_data_visible_pca <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

save(file='fanny_data_pca_visible.RData',fanny_data_visible_pca)
###########Performance measurements#####################################
####FPI####
# extract_FPI<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
#     FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
#   })
# 
# FPI<-c(extract_FPI(fanny_data_v_d_2_1m))
# plot(1:length(FPI),FPI,type='l')
# 
# 
# ####NCE####
# extract_NCE<-function(fuzzy_data) 
#   sapply(fuzzy_data,function(x){
#     memberships<-x$membership
#     entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
#     NCE<-entropy/log10(ncol(memberships)-1)
#   })
# 
# NCE<-c((extract_NCE(fanny_data_v_d_2_1m)))
# 
# plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')



#####Fuzzy Clustering_FANNY_Principal_components_whole_spectra####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:2300)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:2300),which=2)) 
pca_data <-do.call(rbind,strip_DATA)

pr_spectra<-prcomp(pca_data, center=T,scale=T) 
pr_Exp<- pr_varExp(pca_data) # amount of vartiation explained by first 10 PCs
pr_scores <- pr_spectra$x # principal component scores

test_data <- round(as.data.frame(pr_scores[,1:6]),12)

num_clusters <- 1:8

fanny_data_pca <- mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean',maxit=200)
})

#save(file='fanny_data_pca.RData',fanny_data_pca)
###########Performance measurements#####################################
###FPI####
extract_FPI<-function(fuzzy_data) 
  sapply(fuzzy_data,function(x){
    memberships<-x$membership
    part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
    FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
  })

FPI<-c(extract_FPI(fanny_data_pca))
plot(1:length(FPI),FPI,type='l')


####NCE####
extract_NCE<-function(fuzzy_data) 
  sapply(fuzzy_data,function(x){
    memberships<-x$membership
    entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
    NCE<-entropy/log10(ncol(memberships)-1)
  })

NCE<-c((extract_NCE(fanny_data_pca)))

plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')

###########FuzMe data_spectra_whole #####################################

library(spectroscopy)
library(pbapply)
library(cluster)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:2300)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
snVC_DATA <- lapply(filt_DATA,function(x) snvBLC(x))
strip_DATA<- pblapply(snVC_DATA,function(x) strip_spectra(x,c(390:2300),which=20)) 

#  pca_data <-do.call(rbind,strip_DATA)
# # 
#  pr_spectra<-prcomp(pca_data, center=T,scale=T) 
#  pr_Exp<- pr_varExp(pca_data) # amount of vartiation explained by first 10 PCs
#  pr_scores <- pr_spectra$x # principal component scores

test_data <- as.data.frame(do.call(rbind,strip_DATA))
rownames(test_data) <- 1:3190
prev_dir <- getwd()
setwd('C:/fuzme/spectra_whole/')
write.table(test_data,file='spectra_whole.txt',sep='\t',quote=F,dec='.')
setwd('..')
shell('FuzME3_5.exe - Control.txt',intern=T,wait=T)

setwd <- prev_dir
setwd('fuzme_classes/spectra_whole/')

fuz_pca_summary<- read.table(file='summary.txt',header=T,sep='')


###########mclust data #####################################

library(spectroscopy)
library(pbapply)
library(cluster)
library(mclust)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:2300)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:2300),which=2)) 

test_data <- round(as.data.frame(do.call(rbind,strip_DATA),12))

test<-Mclust(data=test_data)

clas <- test$classification

#####Fuzzy Clustering_FANNY_PC_Euc_dist####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_dry_cores.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(raw_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 

spec_cont <- as.data.frame(do.call(rbind,strip_DATA))

pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
pr_scores <- pr_spectra$x # principal component scores

test_data <- as.data.frame(pr_scores[,1:26])


num_clusters <- 1:24

fanny_data_whole_spectra_pc_euc <- mclapply(mc.preschedule=T,mc.cores=8,X=num_clusters, FUN=function (y) {
  temp <- fanny(x=test_data,k=y,metric='euclidean')
})

save(file='fanny_data_whole_spectra_pc_euc.RData',fanny_data_whole_spectra_pc_euc)

##########Performance measurements#####################################
###FPI####

fanny_data <-fanny_data_whole_spectra_pc_euc

extract_FPI<-function(fuzzy_data) 
  sapply(fuzzy_data,function(x){
    memberships<-x$membership
    part_coef<-sum(memberships[,-ncol(memberships)]^2)/nrow(memberships)
    FPI<-1-((ncol(memberships)-1)*part_coef-1)/((ncol(memberships)-1)-1)
  })


FPI<-c(extract_FPI(fanny_data))
plot(1:length(FPI),FPI,type='l')


####NCE####
extract_NCE<-function(fuzzy_data) 
  sapply(fuzzy_data,function(x){
    memberships<-x$membership
    entropy<--1*(sum(memberships[,-ncol(memberships)]*log10(memberships[,-ncol(memberships)]))/nrow(memberships))
    NCE<-entropy/log10(ncol(memberships)-1)
  })

NCE<-c((extract_NCE(fanny_data)))

plot(1:length(NCE),NCE,type='l',xlab='Number of fuzzy clusters')

#####Fuzzy cluster by sample ####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('huntervalley_2cm_ground.RDATA')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

rawg_spectra<-lapply(ground_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sampleg_details<-lapply(ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(rawg_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 



pr_spectra_DATA <-lapply(strip_DATA,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$Sample,z$x[,1:5]),sampleg_details,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=T) # principal component scores
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- 1:6
i<-'gch15'
fanny_data_whole_spectra_pc_euc <- lapply(levels(pr_scores_DATAFRAME$group),function(i){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
  mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) fanny(x=test_data,k=y,metric='euclidean'))
})



# save(file='fanny_data_whole_spectra_pc_euc.RData',fanny_data_whole_spectra_pc_euc)



#####Fuzzy cluster by sample with EPO_trasnformation####
library(spectroscopy)
library(pbapply)
library(cluster)
library(multicore)
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('specFu.RData')
load('EPO_transformation_matrix.RDATA')
rm(trimSpec)
rm(filter_sg)

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

cont_DATA<-pblapply(raw_hv_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
strip_DATA<- pblapply(trim_DATA,function(x) strip_spectra(x,c(500:2450),which=10)) 
abs_DATA<-pblapply(strip_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
snVC_DATA<- pblapply(filt_DATA,snvBLC) 
epo_DATA <- pblapply(snVC_DATA,function(x) {
  epo_tmp<-as.matrix(x) %*% P
  epo_tmp<-data.frame(epo_tmp)
  colnames(epo_tmp)<-seq(500,2450,by =10)
  epo_tmp}
)



pr_spectra_DATA <-lapply(epo_DATA,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$a.Pit,z$x[,1:5]),sample_hv_details,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=T) # principal component scores
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- 1:6
fanny_data_whole_spectra_pc_euc_EPO <- lapply(levels(pr_scores_DATAFRAME$group),function(i){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
  mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) fanny(x=test_data,k=y,metric='euclidean'))
})



save(file='fanny_data_whole_spectra_pc_euc_EPO.RData',fanny_data_whole_spectra_pc_euc_EPO)