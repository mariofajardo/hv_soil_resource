#removing moisture effect from the spectra
library(signal)
library(Cubist)
library(pls)
library(spectroscopy)
library(pbapply)
load('RData/hv_soil_resource_5cm_resolution_pits_2013.RData')
load('RData/check_plots.RData')
load('RData/correct_steps.RData')
prev_dir<- getwd()
setwd('../../../Lab_work/Spectra/EPO')
# Data from Minasny et al. (2010), 100 soil samples under 3 different moisture conditions
spectra0<-read.csv('moisture_dry.txt',header=FALSE)
spectra1<-read.csv('moisture_wet.txt',header=FALSE)
spectra2<-read.csv('moisture_wet2.txt',header=FALSE)

wavelength<-seq(350,2500,by =1) 

# SAVITSKY-GOLAY SMOOTHING FILTER 
abs_filtered0<-filter_sg(spectra0, n = 11, p = 2, m = 0)
abs_filtered1<-filter_sg(spectra1, n = 11, p = 2, m = 0)
abs_filtered2<-filter_sg(spectra2, n = 11, p = 2, m = 0)

#CORRECT STEP OF NIR INSTRUMENT
cont_spectra0<-as.matrix(correct_step(nir_spectra=as.data.frame(abs_filtered0)))
cont_spectra1<-as.matrix(correct_step(nir_spectra=as.data.frame(abs_filtered1)))
cont_spectra2<-as.matrix(correct_step(nir_spectra=as.data.frame(abs_filtered2)))

spectra_trim0<-strip_spectra(cont_spectra0, wavelength,wavlimits = range(500:2450), which = 10)
spectra_trim1<-strip_spectra(cont_spectra1, wavelength,wavlimits = range(500:2450), which = 10)
spectra_trim2<-strip_spectra(cont_spectra2, wavelength,wavlimits = range(500:2450), which = 10)

wavelength10<-seq(500,2450,by =10)

# plot the absorbance spectra
plot(wavelength10,spectra_trim0[1,],type="l",ylim=c(0,1.5))
lines(wavelength10,spectra_trim1[1,],col="blue")
lines(wavelength10,spectra_trim2[1,],col="green")

# perform snv
spec_snvC0<- snvBLC(spectra_trim0)
spec_snvC1<- snvBLC(spectra_trim1)
spec_snvC2<- snvBLC(spectra_trim2)

# plot the absorbance-snv spectra
plot(wavelength10,spec_snvC0[1,],type="l")
lines(wavelength10,spec_snvC1[1,],col="blue")
lines(wavelength10,spec_snvC2[1,],col="green")


#####Calculate the EPO####
#D is the difference matrix (between dry and wet spectra)

D=as.matrix(spec_snvC0-spec_snvC2)
npc<-4  # define no. EPO factors
P<- epo(D,npc)


# save(P,file='EPO_transformation_matrix.RDATA')
setwd(prev_dir)
save(P,file='RData/EPO_transformation_matrix.RDATA')

# #Project the spectra
# Z0 <- as.matrix(spec_snvC0) %*% P    # EPO projected spectra of spec0
# Z1 <- as.matrix(spec_snvC1) %*% P    # EPO projected spectra of spec1
# Z2 <- as.matrix(spec_snvC2) %*% P    # EPO projected spectra of spec2
# 
# # plot  the epo transformed spectra of sample 1
# plot(wavelength10,Z0[1,],"l")
# lines(wavelength10,Z1[1,],"l",col="blue")
# lines(wavelength10,Z2[1,],"l",col="green")

setwd(prev_dir)
##Project the library spectra into transformed spectra

raw_hv_spectra<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[7:2157]
  names(spectra)<-as.numeric(350:2500)
  spectra})

sample_hv_details<-lapply(hv_pits_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

pboptions(type='txt',style=3)
# check_plots(raw_hv_spectra,'Absorb','Absorb')
cont_DATA<-pblapply(raw_hv_spectra,correct_step)   
abs_DATA<-pblapply(cont_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
trim_DATA<- pblapply(filt_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
strip_DATA<- pblapply(trim_DATA,function(x) strip_spectra(x,c(500:2450),which=10)) 
snVC_DATA<- pblapply(strip_DATA,snvBLC) 

epo_DATA <- pblapply(snVC_DATA,function(x) {
  epo_tmp<-as.matrix(x) %*% P
  epo_tmp<-data.frame(epo_tmp)
  colnames(epo_tmp)<-seq(500,2450,by =10)
  epo_tmp}
)

check_plots(epo_DATA,'EPO_transformed_data','Absorbance')


# # Make a Cubist model from EPO transformed dry spectra
# epo.cubist_model<-cubist(x= specZ, y=soilv)
# gf.epo_cubist<- goof(soilv,predict(epo.cubist_model,specZ))
# gf.epo_cubist
# 
# # Predict the values from spectra at different moisture content
# epo.cubist_predict.dry<-predict(epo.cubist_model, as.data.frame(Z0))
# epo.cubist_predict.wet<-predict(epo.cubist_model, as.data.frame(Z1))
# epo.cubist_predict.wet2<-predict(epo.cubist_model, as.data.frame(Z2))
# 
# gfdry.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.dry)
# gfwet.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.wet)
# gfwet2.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.wet2)
# 
# gfdry.epo_cubist
# gfwet.epo_cubist
# gfwet2.epo_cubist
# 
# plot(soilC$TotalC,epo.cubist_predict.dry,xlim=c(0,10), ylim=c(0,10),xlab="Observed",ylab="Predicted")
# points(soilC$TotalC, epo.cubist_predict.wet, col="red")
# points(soilC$TotalC, epo.cubist_predict.wet2,col="green")
# abline(a = 0, b = 1, col = "brown4")
# 
# # Plot dry vs wet predicted
# plot(epo.cubist_predict.dry, epo.cubist_predict.wet)
# abline(a = 0, b = 1, col = "brown4")
# 
# 
