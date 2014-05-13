#removing moisture effect from the spectra
library(signal)
library(Cubist)
library(pls)
library(spectroscopy)
library(pbapply)
load('RData/check_plots.RData')
load('RData/correct_steps.RData')
prev_dir<- getwd()
setwd('../../../Lab_work/Spectra/EPO')
# Data from Minasny et al. (2010), 100 soil samples under 3 different moisture conditions
spectra0<-read.csv('TS_dry_sub.txt',header=TRUE)[,-1]
spectra1<-read.csv('TS_wet_sub.txt',header=TRUE)[,-1]

wavelength<-seq(350,2500,by =1) 

# SAVITSKY-GOLAY SMOOTHING FILTER 
abs_filtered0<-filter_sg(spectra0, n = 11, p = 2, m = 0)
abs_filtered1<-filter_sg(spectra1, n = 11, p = 2, m = 0)

#CORRECT STEP OF NIR INSTRUMENT
cont_spectra0<-as.matrix(correct_step(nir_spectra=as.data.frame(abs_filtered0)))
cont_spectra1<-as.matrix(correct_step(nir_spectra=as.data.frame(abs_filtered1)))

spectra_trim0<-strip_spectra(cont_spectra0, wavelength,wavlimits = range(500:2450), which = 10)
spectra_trim1<-strip_spectra(cont_spectra1, wavelength,wavlimits = range(500:2450), which = 10)

wavelength10<-seq(500,2450,by =10)

# plot the absorbance spectra
plot(wavelength10,spectra_trim0[1,],type="l",ylim=c(0,1.5))
lines(wavelength10,spectra_trim1[1,],col="blue")

# perform snv
spec_snvC0<- snvBLC(spectra_trim0)
spec_snvC1<- snvBLC(spectra_trim1)

# plot the absorbance-snv spectra
plot(wavelength10,spec_snvC0[1,],type="l")
lines(wavelength10,spec_snvC1[1,],col="blue")


#####Calculate the EPO####
#D is the difference matrix (between dry and wet spectra)

D=as.matrix(spec_snvC0-spec_snvC1)
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
