## Bayesian Estimate of chlorophyll fluorescence paramters from  LR curves ###
# R 3.2.1
# 05_31_2017. Jonathan R Plean.
# fixed to estimate with only low light portion of the curve (> 250 PARi)
# Reviewed by: No one

### Script depends on fluorescence Light Responce Model  see Yin(2009) PCE for model details
# RdPhill_Ind_model.R   quadratic function often using including in Text: Plant Physiological Ecology (Lamber, 2008) 
# model estimates s1 (lumped parameter described in Yin(2009) PCE))
# mphi (lowlight slope parameter described in Yin(2009) PCE))
# Rd (respiration rate Dark, y intercept)
# phi2ll  (quantum yield based on fluorescence, initial slope)


# Species Amaranthus tricolor  
# Tri_Tri are Tricolor plant exhibiting tricolor behavior (green,yellow, red)
# Tri are Tricolor plant not displaying tricolor beahvri (green and yellow only)
# HYP are Amaranthus genotype which exhbitis only green leaves

# Sources: collaborator Chris Mure & James Berry
#   UB Biology 


setwd("/Volumes/My Passport for Mac/C_Mure_data") ### Change 4 server

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

## data from all measuremnts days (7_28_16 and 9_29_16) compiled into CMure_tricolor_clean_compiled_data_093016_clean
LR <- read.delim("CMure_tricolor_clean_compiled_data_093016_clean.txt",sep="\t", header=TRUE)
#   LR_1 removes curve for ID 13 after CM and JRP evaluated data on 11-16-19
LR_1<-LR[LR$ID!="13",]

### clean version has ID 3 and 15 removed
### ID 3   Tri Base Red curve has min value of -13 
###  ID 15 has high variance at light level > 500


names(LR_1)

### set up Identifter for Posterior outputs
IGT<-unique(LR_1[,2:6][c("Date","Plant", "Leaf_Section", "Color", "ID")]) 
Plant<-IGT$Plant; Leaf_Section<-IGT$Leaf_Section;Color<-IGT$Color; IDs<-IGT$ID; N<-length(IDs)
Date<-IGT$Date

LR_low<-LR_1[LR_1$PARi<=400,]
hist(LR_low$PARi)
datalist<-vector("list", 7)
for(j in 1:21){
  datalist[[j]]<-list(N=length(LR_low[LR_low$ID==IDs[j],]$PARi), N_ll=length(LR_low[LR_low$ID==IDs[j],]$PARi),
                      A_ll=c(LR_low[LR_low$ID==IDs[j],]$Photo) , Inc_ll=c(LR_low[LR_low$ID==IDs[j],]$PARi) ,phi2_ll=c(LR_low[LR_low$ID==IDs[j],]$PhiPS2))
}
datalist[[3]]

### establishing all models & doing burn in
models<-vector("list", 21)
for(j in 1:21){
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

plot(mcmcsamples[[16]])
gelman.diag(mcmcsamples[[6]])

names(mcmcChain) <- paste(IGT$Plant, IGT$Leaf_Section,IGT$Color,IGT$ID, sep='_')
IGT$Rec<- paste(IGT$Plant, IGT$Leaf_Section,IGT$Color,IGT$ID, sep='_')

##### pull out medain values estimates and write to file
meds<-matrix(, nrow = N, ncol = 8)
for(j in 1:N){
meds[j,]<- apply(mcmcChain[[j]], 2, median)
}
Medians_CF<-cbind(IGT[,1:5], meds)
colnames(Medians_CF)<-c( "Date", "Plant", "Leaf_Section", "Color" ,"ID", "Rd_CF", "mphi","phi2ll", "s1",  "tau_Rd", "tau_phi",
                     "sigma_Rd",  "sigma_phi")
### write table for Median values

write.table(Medians_CF, "medians_LR_Chl_flo_trait_estimates_tricolor_05_31_17.txt", sep="\t", col.name=TRUE,row.names=FALSE)



### save mcmcChains for later analysis as individual files
setwd("/Volumes/My Passport for Mac/C_Mure_data/mcmc_data_CF")

for(k in 1:N){
  write.table( mcmcChain[[k]], c(IGT$Rec[k]), sep="\t", col.name=TRUE,row.names=FALSE )
}



######  95% CI's based on pooled MCMC samples

source("/Volumes/My Passport for Mac/ProgramsDoingBayesianDataAnalysis-2/HDIofMCMC.R")


### leaf section-level posterior distributions
Ahy_b<-rbind(mcmcChain[[1]],mcmcChain[[2]],mcmcChain[[3]])
Ahy_t<-rbind(mcmcChain[[4]],mcmcChain[[5]],mcmcChain[[6]])
At_b<-rbind(mcmcChain[[7]],mcmcChain[[8]],mcmcChain[[9]])
At_t<-rbind(mcmcChain[[10]],mcmcChain[[11]],mcmcChain[[12]])
Att_br<-rbind(mcmcChain[[13]],mcmcChain[[14]],mcmcChain[[15]])
Att_my<-rbind(mcmcChain[[16]],mcmcChain[[17]],mcmcChain[[18]])
Att_tg<-rbind(mcmcChain[[19]],mcmcChain[[20]],mcmcChain[[21]])
###   60% HDR for each geno
l1<-apply(Ahy_b,2 ,HDIofMCMC)
l2<-apply(Ahy_t,2 ,HDIofMCMC) 
l3<-apply(At_b,2 ,HDIofMCMC) 
l4<-apply(At_t,2 ,HDIofMCMC) 
l5<-apply(Att_br,2 ,HDIofMCMC) 
l6<-apply(Att_my,2 ,HDIofMCMC) 
l7<-apply(Att_tg,2 ,HDIofMCMC) 


CIs<-cbind(rep(c("Ahy_b", "Ahy_t", "At_b", "At_t", "Att_br", "Att_my", "Att_tg"), each=2),rbind(l1,l2,l3,l4,l5,l6,l7))

write.table(CIs, "CIs_LR_Chl_flo_trait_estimates_tricolor_05_31_17.txt", sep="\t", col.name=TRUE,row.names=FALSE)




