#####Fuzzy cluster by in situ sample with EPO ####
set.seed(412)
require(doSNOW)
require(foreach)
library(cluster)

pr_spectra_hv_DATA <-lapply(epo_hv_DATA,function(x) prcomp(x, center=T,scale=T)) 
pr_scores_hv_DATA <- mapply(function(y,z) data.frame(group=y$File.Name,z$x[,1:5]),sample_hv_details,pr_spectra_hv_DATA,SIMPLIFY=F,USE.NAMES=T) # principal component scores
a<-pr_scores_hv_DATA[[1]]

num_clusters <- 1:6

cl <-makeCluster(8)
registerDoSNOW(cl)

fanny_pits_pc_euc_EPO <- lapply(pr_scores_hv_DATA,function(i){
  test_data<-test_data<-i[,2:ncol(i)]
  foreach(i=num_clusters,.packages='cluster') %dopar% fanny(x=test_data,k=i,metric='euclidean')
})
stopCluster(cl)
