##############HUNTER_VALLEY_2cm_RESOLUTION_PITS_2013####
#####################pit_samples_5_rep_5cm_resolution###
library(pbapply)
library(reshape)
library(spectroscopy)
setwd("..//Functions")#####    SETWD
load('correct_steps.RData')
load('check_plots.RData')
prev_dir <- getwd()
setwd("../..//Spectra//PITS_HV/")#####    SETWD


###########read the spectra and assign horizon names###########
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=F)
tmp_DATA <- pblapply(files,function(x){       ##tmp_data <-DATA                     
  tmp <- read.csv(x,header=T)
  data <- tmp[-1]
  rep <- as.factor(rep(1:nrow(data), each=5, length.out=nrow(data)))
  tmp2 <- as.data.frame(sapply(data, tapply, rep, mean))
  tmp2$File.Name <- paste(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),paste0('Spectrum',sprintf('%05i',c(1:nrow(tmp2)))),sep='_')
  Pit <- regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T))
  top <- seq(0,by=5,length.out=nrow(tmp2))
  bottom <- seq(5,by=5,length.out=nrow(tmp2))
  tmp3 <- cbind(File.Name=tmp2$File.Name,Pit=Pit,top=top,bottom=bottom,tmp2[-length(tmp2)])
  tmp3[order(tmp3$File.Name,decreasing=F),] 
  
}) 

hor <- sort_df(read.csv('NIR_hor.csv',header=T))     
list_df <- as.list.data.frame(split(hor,hor$id))
hor_names <- lapply(list_df,function(x){
 tmp <- data.frame(hor = rep(x$hor,times=x$bottom-x$top+1),master_hor = rep(x$master_hor,times=x$bottom-x$top+1))
})

insert_columns<- function(a,b){   
 data.frame(File.Name=a$File.Name,a$Pit,top=a$top,bottom=a$bottom,b$hor,b$master_hor,a[-(1:4)])
}

hv_pits_DATA <- mapply(insert_columns,tmp_DATA,hor_names,SIMPLIFY=F)

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
 spectra<-x[6:2156]
 names(spectra)<-as.numeric(350:2500)
 spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
 spectra<-x[1:5]
 spectra})


setwd('..//PITS_HV')
 
save(hv_pits_DATA,file='hv_soil_resource_5cm_resolution_pits_2013.RData')

hv_pits_DATA[[1]][1:2,1:10]
hv_pits_DATA[[2]][1:18,1:6]

