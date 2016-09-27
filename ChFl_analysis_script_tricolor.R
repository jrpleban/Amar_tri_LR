## Bayesian Estimate of PS paramters from  LR curves ###

###   Chris Mure Amaranthus tricolor #####

setwd("~/Desktop/C_Mure_data/")   ### Change 4 server

### package for Bayesian Sampleing
library("rjags", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library") ### Change 4 server
###  Bayesian LR model
source("Mod_scripts/RdPhill_Ind_model.R")

######  traditional Light responce function requiring light (Q) and parameter (Pars)

#### Set up for rjags #########
parameters = c("s1", "Rd25","mphi",
               "phi2ll","tau_Rd", "tau_phi") ### pars to be monitored
adaptSteps = 10000             # Number of steps to "tune" the samplers.
burnInSteps = 10000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 10000        # Total number of steps in chains to save.
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.datalist1<-list(N=N, An=A[,1], Ca=CP[,1], g=g[,1], Jf=Jf[,1], O=O[,1])


######  load data
library(data.table)
## data from two measuremnts days (7_28_16 and 8_04_16) compiled into CM_tricolor_data_clean_08_12_16.txt
LR_1 <- read.delim("CMure_tricolor_clean_compiled_data_092816.txt",sep="\t", header=TRUE)

### data 7_28-16 Tri Base REd curve has min value of -13 and two poor data points
## consider redoing that curve


names(LR_1)

### set up Identifter for Posterior outputs
IGT<-unique(LR_1[,2:6][c("Date","Plant", "Leaf_Section", "Color", "ID")]) 
Plant<-IGT$Plant; Leaf_Section<-IGT$Leaf_Section;Color<-IGT$Color; IDs<-IGT$ID; N<-length(IDs)
Date<-IGT$Date

###  FOR Plotting Raw data 
### PARi rounded to produce mean values for Plant types
LR_1$PARi_grp<-(round(LR_1$PARi, -1))
###  Plant group, combination of Species, Sectino and color
LR_1$Plant_grp<-paste(LR_1$Plant, LR_1$Leaf_Section,LR_1$Color, sep='_')
#### mean values based on plant/leaf type and PARi_group
LR_means<-aggregate(Photo ~ PARi_grp * Plant_grp, data=LR_1, FUN=mean)
LR_means$SD<-aggregate(Photo ~ PARi_grp * Plant_grp, data=LR_1, FUN=sd)[,3]



datalist<-vector("list", N)
for(j in 1:N){
  datalist[[j]]<-list(N_ll=length(LR_1[LR_1$ID==IDs[j],]$PARi),
                         A_ll=c(LR_1[LR_1$ID==IDs[j],]$Photo), Inc_ll=c(LR_1[LR_1$ID==IDs[j],]$PARi),
                         phi2_ll=c(LR_1[LR_1$ID==IDs[j],]$PhiPS2))
}


### establishing all models & doing burn in
models<-vector("list", N)
for(j in 1:N){
models[[j]] <- jags.model(textConnection(RdPhi2ll_mod), 
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



mcmcChain<-vector("list", N); sigma_Rd<-vector("list", N);  sigma_phi2<-vector("list", N)
  for(i in 1:N){
    mcmcChain[[i]] <- as.matrix( mcmcsamples[[i]])
    sigma_Rd[[i]] =1  / sqrt( mcmcChain[[i]][, "tau_Rd" ] )
    sigma_phi2[[i]] =1  / sqrt( mcmcChain[[i]][, "tau_phi" ] )
    mcmcChain[[i]] = as.data.frame(cbind( mcmcChain[[i]], sigma_Rd[[i]], sigma_phi2[[i]] ))
  }

##### pull out medain values estimates and write to file
meds<-matrix(, nrow = N, ncol = 8)
for(j in 1:N){
meds[j,]<- apply(mcmcChain[[j]], 2, median)
}
Medians<-cbind(IGT[,1:3], meds)
colnames(Medians)<-c( "Date", "Plant", "Leaf_Section", "Color" ,"ID", "Rd", "mphi","phi2ll", "s1",  "tau_Rd", "tau_phi",
                     "sigma_Rd",  "sigma_phi")
### write table for Median values

write.table(Medians, "medians_LR_Chl_flo_trait_estimates_tricolor_09_26_16.txt", sep="\t", col.name=TRUE,row.names=FALSE)


