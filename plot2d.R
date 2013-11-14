#####DATA####
library(spectroscopy)
library(plyr)
library(pbapply)
library(ggplot2)
library(reshape2)

load('hv_soil_resource_5cm_resolution_pits_2013.RData')
load('correct_steps.RData')
load('fanny_data_whole_spectra_pc_euc_EPO.RData')


#####colors by horizon####
raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[7:length(x)]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

trim_hv_DATA<- pblapply(raw_hv_spectra,function(x) trimSpec(x, wavlimits=range(500:2450))) 

pboptions(type='txt',style=3)
abs_hv_DATA<-pblapply(trim_hv_DATA,function(x) absorbance<-log(1/x))  
strip_hv_DATA<- pblapply(abs_hv_DATA,function(x) strip_spectra(x,c(500:2451),which=10)) 
filt_hv_DATA<- pblapply(strip_hv_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# mscg_DATA <- pblapply(filtg_DATA,mscBLC)
snVC_hv_DATA<- pblapply(filt_hv_DATA,snvBLC) 
# chcg_DATA <- pblapply(filtg_DATA,chBLC)
inputg_data<-snVC_hv_DATA

#####accumulated absorbance variation ground spectra
spectrag <- list()
for (i in 1:length(inputg_data)){
  x<-i
  y<-inputg_data
  sample2d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  spectrum <- matrix(0,nrow=nrow(sample2d)+1,ncol=ncol(sample2d),dimnames=list(rownames(sample2d)+1,colnames(sample2d)))
  for (j in 2:nrow(sample2d)){
    spectrum[j,] <- spectrum[j-1,]+c(sample2d[j,]-sample2d[j-1,])
  }
  spectrum <- spectrum[-c(1,nrow(spectrum)),]
  names_var<-colnames(spectrum)
  spectrum<- matrix(rep(spectrum,each=2),nrow=nrow(spectrum)*2)
  colnames(spectrum)<-names_var
  spectrag[[i]] <- melt(spectrum,measure.vars=c(names(spectrum)),value.name='Acc_absorvance_variation')
}

#absorbance
# spectrag1 <- list()
# for (i in 1:length(inputg_data)){
#   x<-i
#   y<-inputg_data
#   sample2d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
#   sample2d <- melt(sample2d,measure.vars=c(names(sample2d)),value.name='Acc_absorvance_variation')
#   spectrag1[[i]]<-sample2d
# }

# ggplot(spectrag[[5]], aes(X2, X1, fill = value)) + geom_tile() +
#     labs(x='Wavelenght',y='Depth',title = "Acumulated residuals")

#####Air_dried_DATA####
library(spectroscopy)
library(plyr)
library(pbapply)
library(ggplot2)
library(reshape2)
load('huntervalley_2cm_dry_cores.RDATA')
load('correct_steps.RData')



#####colors by horizon####
raw_spectra<-lapply(DATA,function(x) {
  spectra<-x[7:length(x)]
  names(spectra)<-as.numeric(500:2451)
  spectra})

sample_details<-lapply(DATA,function(x) {
  spectra<-x[1:6]
  spectra})

pboptions(type='txt',style=3)
trim_DATA<- pblapply(raw_spectra,function(x) trimSpec(x, wavlimits=range(500:2451))) 

abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
strip_DATA<- pblapply(abs_DATA,function(x) strip_spectra(x,c(500:2451),which=10)) 
filt_DATA<- pblapply(strip_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
# msc_DATA <- pblapply(filt_DATA,mscBLC)
snVC_DATA<- pblapply(filt_DATA,snvBLC) 
# chc_DATA <- pblapply(filt_DATA,chBLC)
input_data<-snVC_DATA

#####accumulated absorbance variation dried cores
spectra <- list()
for (i in 1:length(input_data)){
  x<-i
  y<-input_data
  sample2d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  spectrum <- matrix(0,nrow=nrow(sample2d)+1,ncol=ncol(sample2d),dimnames=list(rownames(sample2d)+1,colnames(sample2d)))
  for (j in 2:nrow(sample2d)){
    spectrum[j,] <- spectrum[j-1,]+c(sample2d[j,]-sample2d[j-1,])
  }
  spectrum <- spectrum[-c(1,nrow(spectrum)),] 
  names_var<-colnames(spectrum)
  spectrum<- matrix(rep(spectrum,each=2),nrow=nrow(spectrum)*2)
  colnames(spectrum)<-names_var
  spectra[[i]] <- melt(spectrum,measure.vars=c(names(spectrum)),value.name='Acc_absorvance_variation')
}

#####absorbance  dried cores
spectra1 <- list()
for (i in 1:length(input_data)){
  x<-i
  y<-input_data
  sample2d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  sample2d <- sample2d[-1,] 
  names_var<-colnames(sample2d)
  sample2d<- matrix(rep(sample2d,each=2),nrow=nrow(sample2d)*2)
  colnames(sample2d)<-names_var
  spectra1[[i]] <- melt(sample2d,measure.vars=c(names(sample2d)),value.name='Acc_absorvance_variation')
}

#####Plot the ground effect####
# pdf('ground_effect.pdf',width=10,height=6)
# for (i in 1:length(input_data)){
# 
# x<-1
# a <-data.frame(type='air_dried',spectra[[x]])
# b <-data.frame(type='ground',spectrag[[x]])
# sample_spectra <- rbind(a,b)
# 
# print(ggplot(sample_spectra, aes(X2, X1, fill = value)) + 
#         geom_tile() +
#         labs(x='Wavelength',y='Depth',title = "Accumulated residuals air dried core") +
#         facet_wrap(~ type)+
#         scale_y_reverse())
# 
# }
# dev.off()

#####Filter the air dried sample row-wise####
myfun <- function(y) filter_sg(y,n = 11, p = 2, m = 0)
spectraf <- list()
for (i in 1:length(input_data)){
  x<-i
  y<-input_data
  sample2d<-as.matrix(y[[x]][1:nrow(y[[x]]),])
  names_var<-colnames(sample2d)
  filt_sample <- t(myfun(t(sample2d)))  
  spectrum <- matrix(0,nrow=nrow(filt_sample)+1,ncol=ncol(filt_sample),dimnames=list(rownames(filt_sample)+1,colnames(filt_sample)))
  for (j in 2:nrow(filt_sample)){
    spectrum[j,] <- spectrum[j-1,]+c(filt_sample[j,]-filt_sample[j-1,])
  }
  spectrum <- spectrum[-c(1,nrow(spectrum)),] 
  spectrum<- matrix(rep(spectrum,each=2),nrow=nrow(spectrum)*2)
  colnames(spectrum)<-names_var
  spectraf[[i]] <- melt(spectrum,measure.vars=c(names(spectrum)),value.name='Acc_absorvance_variation')
}

#####Plot the ground effect plus a S_V Filter####


pdf('ground_effect_filtered.pdf',width=10,height=6)
for (i in 1:length(input_data)){  
  
  x<-i
  
  a <-data.frame(type='air_dried',spectra[[x]])
  b <-data.frame(type='ground',spectrag[[x]])
  c <-data.frame(type='air_dried_filtered',spectraf[[x]])
  d <-data.frame(type='resid. ground-air',cbind(a[2:3],Acc_absorvance_variation=b$Acc_absorvance_variation-a$Acc_absorvance_variation))
  e <-data.frame(type='resid. ground-air_filt',cbind(a[2:3],Acc_absorvance_variation=b$Acc_absorvance_variation-c$Acc_absorvance_variation))
  
  sample_spectra <- rbind(a,b,c,d,e)
  
  print(ggplot(sample_spectra, aes(Var2, Var1, fill = Acc_absorvance_variation)) + 
          geom_tile() +
          labs(x='Wavelength',
               y='Depth',
               title = "Effect of grinding samples and a possible shortcut") +
          scale_fill_gradient2( low = "brown", high = "green",mid='white') +
          scale_x_continuous(breaks = c(700,1400,2100))+
          facet_wrap(~ type,ncol=5) +
          scale_y_reverse())
  
}
dev.off()
shell.exec('ground_effect_filtered.pdf')