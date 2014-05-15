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


save(P,file='EPO_transformation_matrix.RDATA')
setwd(prev_dir)
