## Bayesian Estimate of PS paramters from  LR curves ###
# R 3.2.1
# Updated: 05_31_2076. Jonathan R Plean.
# Reviewed by: No one

### Script depends on Light Responce Model  
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

setwd("/Volumes/My Passport for Mac/C_Mure_data")   ### Change 4 server

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



## data from all measuremnts days (7_28_16 and 9_29_16) compiled into CMure_tricolor_clean_compiled_data_093016_clean
LR <- read.delim("CMure_tricolor_clean_compiled_data_093016_clean.txt",sep="\t", header=TRUE)
#   LR_1 removes curve for ID 13 after CM and JRP evaluated data on 11-16-19
LR_1<-LR[LR$ID!="13",]
### clean version has ID 3 and 15 removed
### ID 3   Tri Base Red curve has min value of -13 
###  ID 15 has high variance at light level > 500
## consider redoing that curve


names(LR_1)

### set up Identifter for Posterior outputs
IGT<-unique(LR_1[,2:6][c("Date","Plant", "Leaf_Section", "Color", "ID")]) 
Plant<-IGT$Plant; Leaf_Section<-IGT$Leaf_Section;Color<-IGT$Color; IDs<-IGT$ID; N<-length(IDs)
Date<-IGT$Date

length(sort(unique(IGT$ID)))
IGT[order(Plant),] 
datalist<-vector("list", 21)
for(j in 1:21){
  datalist[[j]]<-list(N=length(LR_1[LR_1$ID==IDs[j],]$PARi),
                         An=c(LR_1[LR_1$ID==IDs[j],]$Photo), Q=c(LR_1[LR_1$ID==IDs[j],]$PARi))
}


### establishing all models & doing burn in
models<-vector("list", 21)
for(j in 1:21){
models[[j]] <- jags.model(textConnection(LR_model), 
                     data = datalist[[j]], n.chains=nChains , n.adapt=adaptSteps)
      update(models[[j]], burnInSteps)
}


### sampling all models and calc sigma and converting to martrix form 
mcmcsamples<-vector("list", nChains); 
for(j in 1:21){
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

names(mcmcChain) <- paste(IGT$Plant, IGT$Leaf_Section,IGT$Color,IGT$ID, sep='_')
IGT$Rec<- paste(IGT$Plant, IGT$Leaf_Section,IGT$Color,IGT$ID, sep='_')






##### pull out medain values estimates and write to file
meds<-matrix(, nrow = N, ncol = 5)
for(j in 1:N){
meds[j,]<- apply(mcmcChain[[j]], 2, median)
}
MediansGE<-cbind(IGT[,1:5], meds)
colnames(MediansGE)<-c( "Date", "Plant", "Leaf_Section", "Color" ,"ID","Amax", "Rd_GE", "phiCO2", "tau", "sigma")
### write table for Median values

write.table(MediansGE, "medians_LR_trait_estimates_tricolor_alldata_05_31_17.txt", sep="\t", col.name=TRUE,row.names=FALSE)



### save mcmcChains for later analysis as individual files
setwd("/Volumes/My Passport for Mac/C_Mure_data/mcmc_data_GE")

for(k in 1:N){
write.table( mcmcChain[[k]], c(IGT$Rec[k]), sep="\t", col.name=TRUE,row.names=FALSE )
}


###   95% credible intervals using function from John K. Kruschke Doing Bayesian Data Analysis in R (2014)
source("/Volumes/My Passport for Mac/ProgramsDoingBayesianDataAnalysis-2/HDIofMCMC.R")


### leaf section-level posterior distributions
Ahy_b<-rbind(mcmcChain[[1]],mcmcChain[[2]],mcmcChain[[3]])
Ahy_t<-rbind(mcmcChain[[4]],mcmcChain[[5]],mcmcChain[[6]])
At_b<-rbind(mcmcChain[[7]],mcmcChain[[8]],mcmcChain[[9]])
At_t<-rbind(mcmcChain[[10]],mcmcChain[[11]],mcmcChain[[12]])
Att_br<-rbind(mcmcChain[[13]],mcmcChain[[14]],mcmcChain[[15]])
Att_my<-rbind(mcmcChain[[16]],mcmcChain[[17]],mcmcChain[[18]])
Att_tg<-rbind(mcmcChain[[19]],mcmcChain[[20]],mcmcChain[[21]])
###   95% HDR for each geno
l1<-apply(Ahy_b,2 ,HDIofMCMC)
l2<-apply(Ahy_t,2 ,HDIofMCMC) 
l3<-apply(At_b,2 ,HDIofMCMC) 
l4<-apply(At_t,2 ,HDIofMCMC) 
l5<-apply(Att_br,2 ,HDIofMCMC) 
l6<-apply(Att_my,2 ,HDIofMCMC) 
l7<-apply(Att_tg,2 ,HDIofMCMC) 


CIs<-cbind(rep(c("Ahy_b", "Ahy_t", "At_b", "At_t", "Att_br", "Att_my", "Att_tg"), each=2),rbind(l1,l2,l3,l4,l5,l6,l7))

write.table(CIs, "CIs_LR_GS_trait_estimates_tricolor_05_31_17.txt", sep="\t", col.name=TRUE,row.names=FALSE)








