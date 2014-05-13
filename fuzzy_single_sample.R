#####Fuzzy cluster by sample EPO ####
set.seed(412)
require(doSNOW)
require(foreach)
library(cluster)

no_out_data<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[7:ncol(x)]
  colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
  spectra})

no_out_details<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

pr_spectra_DATA <-lapply(no_out_data,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$Sample,z$x[,1:5]),no_out_details,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=F) # principal component scores
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- 1:6

#   i="gku17"

cl <-makeCluster(8)
registerDoSNOW(cl)

fanny_data_by_sample_pc_euc_no_out <- lapply(levels(pr_scores_DATAFRAME$group),function(i){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
foreach(i=num_clusters,.packages='cluster') %dopar% fanny(x=test_data,k=i,metric='euclidean')
})
stopCluster(cl)
