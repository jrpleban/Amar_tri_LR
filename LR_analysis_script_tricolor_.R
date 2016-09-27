## Bayesian Estimate of PS paramters from  LR curves ###
# R 3.2.1
# 09_26_2016. Jonathan R Plean.
# Reviewed by: No one

### SCript depending on Light Responce Model  
# LR_model.R  quadratic function often using including in Text: Plant Physiological Ecology (Lamber, 2008) 
# model estimates Amax (Assimilation at Saturating Irradiance)
# Rd ( respiration rate Dark)
# phiCO2  (quantum yield based on gas exchange)


# Species Amaranthus tricolor  
# Tri_Tri are Tricolor plant exhibiting tricolor behavior (green,yellow, red)
# Tri are Tricolor plant not displaying tricolor beahvri (green and yellow only)
# HYP are Amaranthus genotype which exhbitis only green leaves

# Sources: collaborator Chris Mure & James Berry
#   UB Biology 

setwd("~/Desktop/C_Mure_data/")   ### Change 4 server

### package for Bayesian Sampleing
library("rjags", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library") ### Change 4 server
###  Bayesian LR model
source("Mod_scripts/LR_model.R")


#### Set up for rjags #########   ## some of the code is modified from the Puppy Book
###    Doing Bayesian data analysis: A tutorial introduction with R  (Kruschke, 2010)
parameters = c("Amax", "Rd",
               "phiCO2","tau") ### pars to be monitored
adaptSteps = 10000             # Number of steps to "tune" the samplers.
burnInSteps = 10000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 10000        # Total number of steps in chains to save.
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.datalist1<-list(N=N, An=A[,1], Ca=CP[,1], g=g[,1], Jf=Jf[,1], O=O[,1])



## data from all measuremnts days (7_28_16 and 8_04_16) compiled into CM_tricolor_data_clean_09_298_16.txt
LR_1 <- read.delim("CMure_tricolor_clean_compiled_data_092816.txt",sep="\t", header=TRUE)

### data 7_28-16 Tri Base REd curve has min value of -13 and two poor data points
## consider redoing that curve


names(LR_1)

### set up Identifter for Posterior outputs
IGT<-unique(LR_1[,2:6][c("Date","Plant", "Leaf_Section", "Color", "ID")]) 
Plant<-IGT$Plant; Leaf_Section<-IGT$Leaf_Section;Color<-IGT$Color; IDs<-IGT$ID; N<-length(IDs)
Date<-IGT$Date


datalist<-vector("list", N)
for(j in 1:N){
  datalist[[j]]<-list(N=length(LR_1[LR_1$ID==IDs[j],]$PARi),
                         An=c(LR_1[LR_1$ID==IDs[j],]$Photo), Q=c(LR_1[LR_1$ID==IDs[j],]$PARi))
}


### establishing all models & doing burn in
models<-vector("list", N)
for(j in 1:N){
models[[j]] <- jags.model(textConnection(LR_model), 
                     data = datalist[[j]], n.chains=nChains , n.adapt=adaptSteps)
      update(models[[j]], burnInSteps)
}


### sampling all models and calc sigma and converting to martrix form 
mcmcsamples<-vector("list", nChains); 
for(j in 1:N){
  mcmcsamples[[j]] <- coda.samples(models[[j]],
                                   variable.names=parameters,
                                   n.iter=nPerChain , thin=thinSteps )
}

mcmcChain<-vector("list", N); sigma<-vector("list", N)
  for(i in 1:N){
    mcmcChain[[i]] <- as.matrix( mcmcsamples[[i]])
    sigma[[i]] =1  / sqrt( mcmcChain[[i]][, "tau" ] )
    mcmcChain[[i]] = as.data.frame(cbind( mcmcChain[[i]], sigma[[i]] ))
  }

##### pull out medain values estimates and write to file
meds<-matrix(, nrow = N, ncol = 5)
for(j in 1:N){
meds[j,]<- apply(mcmcChain[[j]], 2, median)
}
Medians<-cbind(IGT[,1:5], meds)
colnames(Medians)<-c( "Date", "Plant", "Leaf_Section", "Color" ,"ID","Amax", "Rd", "phiCO2", "tau", "sigma")
### write table for Median values

write.table(Medians, "medians_LR_trait_estimates_tricolor_alldata_09_26_16.txt", sep="\t", col.name=TRUE,row.names=FALSE)



