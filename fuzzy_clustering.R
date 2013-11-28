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
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$File.Name,z$x[,1:5]),sample_hv_details,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=T) # principal component scores
a<-pr_scores_DATA[[1]]

num_clusters <- 1:6
fanny_pits_pc_euc_EPO <- lapply(pr_scores_DATA,function(i){
  test_data<-i[,2:ncol(i)]
  mclapply(mc.preschedule=F,mc.cores=8,X=num_clusters, FUN=function (y) fanny(x=test_data,k=y,metric='euclidean'))
})



save(file='fanny_pits_pc_euc_EPO.RData',fanny_pits_pc_euc_EPO)