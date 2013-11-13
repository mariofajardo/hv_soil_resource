require(rgl) 
require(spectroscopy)
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('check_plots.RData')
library(pbapply)
prev_dir <-getwd()

#####colors by horizon####
raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[6:2156]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:5]
  spectra})

pboptions(type='txt',style=3)
# check_plots(raw_hv_spectra,'Absorb','Absorb')
cont_DATA<-pblapply(raw_hv_spectra,correct_step)   
# check_plots(cont_DATA,'Absorbance','Absorbance')
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
# check_plots(abs_DATA,'Absorbance','Absorbance')
strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(500:2300),which=10)) 
#check_plots(strip_DATA,'Stripped spectra','Absorbance')
filt_DATA<- pblapply(strip_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# check_plots(filt_DATA,'S-G Filter','Absorbance')
snVC_DATA<- pblapply(filt_DATA,snvBLC) 

#zlim <- range(data3d)
#zlen <- zlim[2] - zlim[1] 
#colorlut <- terrain.colors(zlen,alpha=0) # height color lookup table
#colorlut <- colorRampPalette(c("blue", "orange"))( zlen )
#col <- colorlut[data3d-zlim[1]+1 ]
#new<-data3d[rep( 1:nrow(data3d),each=1),]

input_data<-filt_DATA

for (i in 1:length(input_data)){
  x<-i
  y<-input_data
  
  sample3d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  sample3d_details<-sample_hv_details[[x]][1:nrow(sample_hv_details[[x]]),]
  levels(sample3d_details$b.hor)<-terrain.colors(n=length(levels(sample_hv_details[[x]]$b.hor)),alpha=.75)
  col <- matrix(sample3d_details$b.hor,nrow=nrow(sample3d),ncol=ncol(sample3d))
  
  
  persp3d(x=1:nrow(sample3d)*2,y=as.numeric(colnames(sample3d)),z=sample3d,
          aspect=c(2, 1, .5),
          col=col,
          xlab='Depth',
          ylab='wavelength',
          main=paste0('sample: ', regmatches(as.character(sample3d_details[1,1]),(regexpr(pattern="(?!=[_])[[:alnum:]]+",as.character(sample3d_details[1,1]),perl=T))))
  )
  
  play3d(spin3d(axis=c(0,0,1), rpm=.001), duration=15)
}

#####colors by spectra####
input_data<-filt_DATA

for (i in 1:length(input_data)){
  x<-i
  y<-input_data
  
  sample3d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  sample3d_details<-sample_hv_details[[x]][1:nrow(sample_hv_details[[x]]),]
  colspectra <- spectra2colour(cont_DATA[[i]],colnames(cont_DATA[[i]]))
  col <- matrix(colspectra[,4],nrow=nrow(cont_DATA[[i]]),ncol=ncol(cont_DATA[[i]]))
  
  persp3d(x=1:nrow(sample3d),y=as.numeric(colnames(sample3d)),z=sample3d,
          aspect=c(2, 1, .5),
          col=col,
          xlab='Depth',
          ylab='wavelength',
          main=paste0('sample: ', regmatches(as.character(sample3d_details[1,1]),(regexpr(pattern="(?!=[_])[[:alnum:]]+",as.character(sample3d_details[1,1]),perl=T))))
  )
  
  play3d(spin3d(axis=c(0,0,1), rpm=0.001), duration=15)
}


#####colors by fanny clusters####

# load('fanny_data_visibe.RData')
# 
# raw_spectra<-lapply(DATA,function(x) {
#   spectra<-x[7:2157]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_details<-lapply(DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# 
# fuzzy_data_selected<- fuzzy_data[[4]]
# 
# input_data<-snVC_DATA
# sample_details <-data.frame(do.call(rbind,sample_details),cluster=fuzzy_data_selected$clustering)
# sample_details <-split(sample_details,f=sample_details$Sample)
# 
# 
# for (i in 1:length(input_data)){
#   x<-i
#   y<-input_data
#   
#   sample3d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
#   sample3d_details<-sample_details[[x]][1:nrow(sample_details[[x]]),]
#   levels(sample3d_details$cluster)<-terrain.colors(n=15,alpha=.75)
#   col <- matrix(sample3d_details$cluster,nrow=nrow(sample3d),ncol=ncol(sample3d))
#   
#   persp3d(x=1:nrow(sample3d)*2,y=as.numeric(colnames(sample3d)),z=sample3d,
#           aspect=c(2, 1, .5),
#           col=col,
#           xlab='Depth',
#           ylab='wavelength',
#           main=paste0('sample: ', regmatches(as.character(sample3d_details[1,1]),(regexpr(pattern="(?!=[_])[[:alnum:]]+",as.character(sample3d_details[1,1]),perl=T))))
#   )
#   
#   play3d(spin3d(axis=c(0,0,1),rpm=.001), duration=15)
# }

#####3d pca colors by horizon####
library(rgl)
library(pca3d)
library(spectroscopy)
library(pbapply)
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[6:2156]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:5]
  spectra})


cont_DATA<-pblapply(raw_hv_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 


spec_cont <- as.data.frame(do.call(rbind,strip_DATA))

pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
pr_scores <- pr_spectra$x # principal component scores

sample3d_details<-data.frame(do.call(rbind,sample_hv_details))

open3d()
pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.3,
      group=sample3d_details$b.master_hor,
      show.plane=F,
      show.group.labels=T)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  
  c1<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.5),add=T,alpha=.5)
}


#####3d pca visible spectra clusters ####
# library(rgl)
# library(pca3d)
# library(spectroscopy)
# library(pbapply)
# load('huntervalley_2cm_dry_cores.RDATA')
# load('correct_steps.RData')
# load('specFu.RData')
# load('fanny_data_visibe.RData')
# rm(trimSpec)
# rm(filter_sg)
# 
# raw_spectra<-lapply(DATA,function(x) {
#   spectra<-x[7:2157]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_details<-lapply(DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# 
# cont_DATA<-pblapply(raw_spectra,correct_step)    
# trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
# abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
# strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(400:2300),which=2)) 
# filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# spec_cont <- as.data.frame(do.call(rbind,filt_DATA))
# 
# pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
# pr_scores <- pr_spectra$x # principal component scores
# 
# number_of_classes<-3
# palette(terrain.colors(number_of_classes))
# fuzzy_data_selected<- fanny_data_visibe[[number_of_classes-2]]
# sample3d_details<-data.frame(do.call(rbind,sample_details),fuzzy_data_selected$clustering)
# 
# open3d()
# pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$fuzzy_data_selected.clustering),
#       radius=.3,
#       group=sample3d_details$fuzzy_data_selected.clustering,
#       show.plane=F,
#       show.group.labels=T)
# 
# for (i in 1:number_of_classes){
#   
#   c1<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1])
#   c2<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,2])
#   c3<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,3])
#   centremean <- c(c1,c2,c3)
#   plot3d(ellipse3d(cov(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1:3]),col=palette()[i],centre=centremean,level=.6),add=T,alpha=.5)
# }
#####3d pca whole spectra clusters####
# library(rgl)
# library(pca3d)
# library(spectroscopy)
# library(pbapply)
# load('hv_soil_resource_5cm_resolution_pits_2013.RData')
# load('correct_steps.RData')
# load('specFu.RData')
# load('fanny_data_whole_spectra.RData')
# rm(trimSpec)
# rm(filter_sg)
# 
# raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
#   spectra<-x[7:2157]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_hv_details<-lapply(hv_pits_DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# 
# cont_DATA<-pblapply(raw_spectra,correct_step)    
# trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(400:2300))) 
# abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
# filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(400:2300),which=10)) 
# 
# 
# spec_cont <- as.data.frame(do.call(rbind,strip_DATA))
# 
# pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
# pr_scores <- pr_spectra$x # principal component scores
# 
# number_of_classes<-4 ###add number
# palette(terrain.colors(number_of_classes))
# fuzzy_data_selected<- fanny_data_whole_spectra[[number_of_classes]]
# sample3d_details<- data.frame(do.call(rbind,sample_details),fuzzy_data_selected$clustering)
# 
# 
# open3d()
# pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$fuzzy_data_selected.clustering),
#       radius=.6,
#       group=sample3d_details$fuzzy_data_selected.clustering,
#       show.plane=F,
#       show.group.labels=T)
# 
# for (i in 1:number_of_classes){
#   
#   c1<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1])
#   c2<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,2])
#   c3<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,3])
#   centremean <- c(c1,c2,c3)
#   plot3d(ellipse3d(cov(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1:3]),col=palette()[i],centre=centremean,level=.6),alpha=.5,add=T)
# }
# 


#####3d pca colors by horizon on visible spectra pca space####
library(rgl)
library(pca3d)
library(spectroscopy)
library(pbapply)
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[6:2156]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:5]
  spectra})


cont_DATA<-pblapply(raw_hv_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(390:700),which=2)) 
filt_DATA<- pblapply(strip_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
snVC_DATA<- pblapply(filt_DATA,snvBLC) 
spec_cont <- as.data.frame(do.call(rbind,snVC_DATA))

pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
pr_scores <- pr_spectra$x # principal component scores

sample3d_details<-data.frame(do.call(rbind,sample_hv_details))

open3d()
pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.5,
      group=sample3d_details$b.master_hor,
      show.group.labels=T,
      show.plane=F)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  
  c1<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.4),alpha=.5,add=T)
}

#####3d pca colors by cluster on visible spectra pca space####
# library(rgl)
# library(pca3d)
# library(spectroscopy)
# library(pbapply)
# load('huntervalley_2cm_dry_cores.RDATA')
# load('correct_steps.RData')
# load('specFu.RData')
# load('fanny_data_visibe.RData')
# rm(trimSpec)
# rm(filter_sg)
# 
# raw_spectra<-lapply(DATA,function(x) {
#   spectra<-x[7:2157]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_details<-lapply(DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# 
# cont_DATA<-pblapply(raw_spectra,correct_step)    
# trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700))) 
# abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
# strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(390:700),which=2)) 
# filt_DATA<- pblapply(strip_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
# snVC_DATA<- pblapply(filt_DATA,snvBLC) 
# spec_cont <- as.data.frame(do.call(rbind,snVC_DATA))
# 
# pr_spectra<-prcomp(spec_cont, center=T,scale=T) 
# pr_scores <- pr_spectra$x # principal component scores
# 
# 
# number_of_classes<-5
# palette(terrain.colors(number_of_classes))
# 
# fuzzy_data_selected<- fanny_data_visibe[[number_of_classes-2]]
# 
# sample3d_details<-data.frame(do.call(rbind,sample_details),fuzzy_data_selected$clustering)
# 
# open3d()
# pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$fuzzy_data_selected.clustering),
#       radius=.6,
#       group=sample3d_details$fuzzy_data_selected.clustering,
#       show.plane=F,
#       show.group.labels=T)
# 
# for (i in 1:number_of_classes){
#   
#   cluster <- i
#   c1<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==cluster,1])
#   c2<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==cluster,2])
#   c3<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==cluster,3])
#   centremean <- c(c1,c2,c3)
#   plot3d(ellipse3d(cov(pr_scores[sample3d_details$fuzzy_data_selected.clustering==cluster,1:3]),col=palette()[i],centre=centremean,level=.6,lwd=5),alpha=.5,add=T)
# }


#####3d pca colors by horizon on visible spectra pca space plus depth####
library(rgl)
library(pca3d)
library(spectroscopy)
library(pbapply)
load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('specFu.RData')
rm(trimSpec)
rm(filter_sg)

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[6:2156]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:5]
  spectra})


cont_DATA<-pblapply(raw_hv_spectra,correct_step)  
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 


spec_strip_depth <- data.frame(as.matrix.data.frame(do.call(rbind,strip_DATA)),depth=do.call(rbind,sample_hv_details)$bottom)

pr_spectra<-prcomp(spec_strip_depth, center=T,scale=T) 
pr_scores <- pr_spectra$x # principal component scores


sample3d_details<-data.frame(do.call(rbind,sample_hv_details))

open3d()
pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
      radius=.3,
      group=sample3d_details$b.master_hor,
      show.group.labels=T,
      show.plane=F)

number_of_classes<-length(levels(sample3d_details$b.master_hor))
palette(terrain.colors(number_of_classes))
for (i in 1:number_of_classes){
  
  c1<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
  c2<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
  c3<- mean(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
  centremean <- c(c1,c2,c3)
  plot3d(ellipse3d(cov(pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.90),add=T,alpha=.4)
  #lines3d(x=c(0,c1),y=(0,c2),z=(0,c3))
}

#####3d pca colors by cluster on visible spectra pca space plus depth####
# library(rgl)
# library(pca3d)
# library(spectroscopy)
# library(pbapply)
# load('hv_soil_resource_5cm_resolution_pits_2013.RData')
# load('correct_steps.RData')
# load('specFu.RData')
# load('fanny_data_v_d_2_1m.RData')
# rm(trimSpec)
# rm(filter_sg)
# 
# raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
#   spectra<-x[6:2156]
#   names(spectra)<-as.numeric(350:2500)
#   spectra})
# 
# sample_hv_details<-lapply(hv_pits_DATA,function(x) {
#   spectra<-x[1:5]
#   spectra})
# 
# cont_DATA<-pblapply(raw_hv_spectra,correct_step)  
# trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(390:700)))
# abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))
# filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0))
# strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(390:700),which=5)) 
# 
# spec_strip_depth <- data.frame(do.call(rbind,strip_DATA),depth=do.call(rbind,sample_hv_details)$bottom)
# 
# pr_spectra<-prcomp(spec_strip_depth, center=T,scale=T) 
# pr_scores <- pr_spectra$x # principal component scores#check_plots(strip_DATA,'Stripped spectra','Absorbance')
# 
# 
# sample3d_details<-data.frame(do.call(rbind,sample_details))
# 
# number_of_classes<-3
# palette(terrain.colors(number_of_classes))
# 
# fuzzy_data_selected<- fanny_data_v_d_2_1m[[number_of_classes-2]]
# 
# sample3d_details<-data.frame(do.call(rbind,sample_hv_details),fuzzy_data_selected$clustering)
# 
# open3d()
# pca3d(pr_scores[,1:3],col=as.numeric(sample3d_details$fuzzy_data_selected.clustering),
#       radius=.4,
#       group=sample3d_details$fuzzy_data_selected.clustering,
#       show.plane=F,
#       show.group.labels=T)
# 
# for (i in 1:number_of_classes){
#   
#   
#   c1<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1])
#   c2<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,2])
#   c3<- mean(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,3])
#   centremean <- c(c1,c2,c3)
#   plot3d(ellipse3d(cov(pr_scores[sample3d_details$fuzzy_data_selected.clustering==i,1:3]),col=palette()[i],centre=centremean,level=.9,lwd=5),alpha=.4,add=T)
# }
# 

#####GROUND DATA####



