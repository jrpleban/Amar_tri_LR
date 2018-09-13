### Photosynthesis Research using MultisynQ
##  modeling fitting exponential decay Phi2 vs light intensity
# R 3.2.1
# Updated: 05_12_2017. Jonathan R Pleban UB Geography  
# Reviewed by: No one

### Script depends on  model called Phi2_decay_model_hier
## currently a 1 level hierarchy can be set by species, treatment, leaf type esc.
################################################
# Developed using data from PhotosynQ project:
#   A. tricolor leaf photosynthesis by leaf type and section
# Hierarch based on leaf type and section 
# Collaborators: J. Berry & C. Mure UB Biology
#                 D.S Mackay UB Geography
###  packages used 
require("stringi")
require(plyr)
require("coda")
require("rjags")
require(ggmcmc)

## set as source directory
setwd("~/Desktop/Amaranth_MS/PhotosynQ_modeling_R")
#Dat <- as.data.frame(read.delim("data/data_05_09_17.txt"))
#names(Dat)
#length(Dat[,1])
#Dat$PARi<-Dat[,10]
#names(Dat)
### sort by decreasing light intesity
#D<-arrange(Dat,desc(PARi))
#plot(D$Light.Intensity..PAR.)
#str(D)
#colN<-names(D)
#D$Leaf_Sec<-paste(D[,20],D[,18], sep="_")
#unique(D$Leaf_Sec)
## simpify name coding for leaf type based on earlier versin of photsynQ questioning
#D$L_sec1<-mapvalues(D$Leaf_Sec, from=c("Tri Tri_Green", "Tri Tri_Red", "Tri Uniform_Green", "Tri Uniform_Red",  
 #                               "Tri Tri_Yellow"), to=c("TTG","TTR", "TUG","TUR","TTY" ))
#D$L_sec<-stri_sub(D$L_sec,1,3)
#write.csv(D,"UBdata_05_09_17_cleaned.csv")


## LOADING Combined WY and UB datasets from Base MultispeQ protocol . ## Leaf Photosynthesis MultispeQ V1.0
D1<-as.data.frame(read.csv("data/WY_MulitspeQ_base_protocol_cleaned_09022018.csv"))
D1$PARi<-D1$Light.Intensity..PAR.
D2<-as.data.frame(read.csv("data/UBdata_05_09_17_cleaned.csv"))
names(D)
plot(D1$PARi,D1$Phi2, xlim=c(0,2000))
points(D2$PARi,D2$Phi2, pch=3)


#####   Expondential decay function used to make median predictions
###  Phi2= alpha* (1 - exp(-beta * PARi)   
###alpha (y-intercept) is P[1], beta (decay exp) P[2]
ExD_Fun<-function(P, Dat){
  Phi2P<-(P[1]-P[3]) * (exp(P[2] * Dat))+P[3]
  return(Phi2P)
}



######## Begin Bayesian evaluation using exp decay model
######## 
#### Set up for rjags #########  this sets up how many interation to run jags sampler
parameters = c("alpha", "beta","kappa",
               "mu.alpha", "mu.beta","mu.kappa",
               "tau.alpha", "tau.beta","tau.kappa","tau")### pars to be monitored (could add mu.alpha and mu.beta esc if you want to monitor that part of the hierarchy)
adaptSteps = 50000             # Number of steps to "tune" the samplers.
burnInSteps = 50000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 5000       # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.



###############################
###  Hierarchical Bayesian Model    ####
###       Expondenial decap model phi2 vs light        ####
###############################
source("Modeling_scripts/Phi2_exp_decay_model2.R")   
### feel free to play around with model 
###   such as changing prior distributins to see differences in posteriors

### list is used to pass data to model
datalist1<-list(N=length(c(D1[,1], D2[,1])), N_GT=5, Phi2=c(D1$Phi2,D2$Phi2) ,GT=c(D1$L_secF,D2$L_secF),
                PARi=c(D1$PARi,D2$PARi))
datalist1  ### N= sample size, N_L_T_S = types of leaf sections
           ###  Phi2 is observatin of phi2, PARi light intesity (umol m-2 s-1)
           #####   L_T_S is coding for each leaf type
model <- jags.model(textConnection(Phi2_decay_model_hier),
                    data = datalist1, n.chains=nChains , n.adapt=adaptSteps)
#### some burn in period is necessary to modle begins to converge at posterior distribution before actual sampling
update(model, burnInSteps)
### generating sampling distribtion this will take a while maybe 5-10 miniutes depending on system and sample size
mcmc_samples<- coda.samples(model,
                            variable.names=parameters,
                            n.iter=nPerChain , thin=thinSteps )
### convert rjags mcmc object to matrix
mcmcChain = as.matrix( mcmc_samples)
### convert tau to sigma and add to matrix easier interpretatin of error.
sigma =1  / sqrt(mcmcChain[, "tau" ] )
mcmcChain = as.data.frame(cbind( mcmcChain, sigma ))
### quick diagnistics
## error plt
hist(mcmcChain$sigma)
## gelmen and rubin convergence diagnostics (should be very close to 1 for all paramters (1.1 max)
## if too high look at trace plots below consider running longer chains or changing priors
gelman.diag(mcmc_samples)
### quick look at median estimates
apply(mcmcChain,2,median)
hist(mcmcChain$"tau.alpha")
#### visualize samples 
### with this data set and this set up notice poor converge on L_T_S 2   TTR
### you can see it has poor Gelman diagnostics also  
#ggs_traceplot(ggs(mcmc_samples))
## NOTE:density plots will not have uniform x axis
#ggs_density(ggs(mcmc_samples))

names(mcmcChain)
### if satified save median parameter estimates for validatin plots 
Pars<-cbind(apply(mcmcChain,2,median)[1:5],betas<-apply(mcmcChain,2,median)[6:10],kappas<-apply(mcmcChain,2,median)[11:15])
### Hyp PAR data for validation plotting
PAR<-seq(10,2000,by=10)
### medians predictions for PAR between 10-2000 using Exp decay func
PREDS <- d <- as.data.frame(matrix(nrow=length(PAR), ncol=6))
for (i in 1:5) {PREDS[,i] <- ExD_Fun(c(Pars[i,]), PAR)
                PREDS[,6] <-  PAR
}

setwd("~/Desktop/Amaranth_MS/PhotosynQ_modeling_R")
D<-as.data.frame(read.csv("data/Combined_data_09_11_18_cleaned.csv"))
colnames(D)<-c("Phi2","PARi","L_sec")
str(D)
### setting up for validaiton plots
Leafs<-c("TUG","TUR", "TTG","TTR","TTY")
colnames(PREDS)<-c(Leafs, "PAR")



Acols<-c("green1", "red1", "seagreen2", "deeppink4" ,"yellow")

pdf("figs/Fig_MultispeQ_Combined_random_data_phi2_vs_PARi_Amaranth_Tricolor.pdf", width=6, height=6)



plot(D[D$L_sec=="TUG",]$PARi, D[D$L_sec=="TUG",]$Phi2, col=Acols[1],ylim=c(0,0.83),
     xlim=c(0,2005),xlab="" , ylab="",main="Combined Data"  ,pch=1, cex=.3)
for (i in 1:4){
  points(D[D$L_sec==Leafs[i+1],]$PARi, D[D$L_sec==Leafs[i+1],]$Phi2, col=Acols[i+1],pch=i+1, cex=.3)
}
for (i in 1:5){
  lines(PREDS[,6], PREDS[,i], col=Acols[i], lwd=3, lty=3)
}
mtext(expression(~Phi~""[2]~"(mol e"^{-1}~"mol photon"^{-1}~")"), side=2.4, cex=1, line=2)
mtext(expression("PAR"[i]~"("~mu~"mol m"^{-2}~"s"^{-1}~")"), side=1, line=2.5,cex=1)
#text(1200, .4, c("data fit to exponential function"))
#text(1200, .36, c("Phi2= (alpha-kappa) * (exp(beta) * PARi))+kappa"))
legend("topright", Leafs, col=Acols, pch=15)

dev.off()


plot(D[D$L_sec=="TTY",]$PARi, D[D$L_sec=="TTY",]$Phi2, col=Acols[5],ylim=c(0,0.83),
     xlim=c(0,2005),xlab="" , ylab="",main="Combined Data"  )
lines(PREDS[,6], PREDS[,5], col=Acols[5])


pdf("figs/Fig_MultispeQ_combined_data_density_plot_phi2_intercept_Amaranth_Tricolor.pdf", width=6, height=6)



hist(mcmcChain[,1], prob=TRUE, breaks=10,xlim=c(0.05,.75),ylim=c(0,55),
     col=Acols[1] , xlab="",main="")
for (i in 1:4){
hist(mcmcChain[,i+1], prob=TRUE, breaks=10,col=Acols[i+1], add=TRUE)
}
legend("topleft", Leafs, col=Acols, pch=15)
mtext(expression(~Phi~""[PSII]~" intercept (mol e- mol photon"^{-1}~")"), side=1, cex=1, line=2.4)

dev.off()

pdf("figs/Fig_MultispeQ_combined_dat_density_plot_beta_decay_par_Amaranth_Tricolor.pdf", width=6, height=6)



hist(mcmcChain[,6], prob=TRUE, breaks=20,ylim=c(0,2600),xlim=c(-.01, 0),
     col=Acols[1] , xlab="", main="PhotosynQ data")
for (i in 1:4){
  hist(mcmcChain[,i+6], prob=TRUE, breaks=20,col=Acols[i+1], add=TRUE)
}
legend("topleft", Leafs, col=Acols, pch=15)
mtext(expression(~beta~"decay parameter"), side=1, cex=1, line=2.4)

dev.off()

##########   
##########   
##########   


source("Modeling_scripts/HDIofMCMC.R") 


HDI<-apply(mcmcChain,2,HDIofMCMC)
meds<-apply(mcmcChain,2,median)

plot(meds[1:5], ylim=c(0.1,0.7), cex=.1)

for (i in 1:5){
arrows(i,HDI[1,i],
       i,HDI[2,i],  lwd=1, length=0.05, angle=90, code=3, col=Acols[i])
}

for (i in 1:5){
  points(c(i) ,meds[i], col=Acols[i], pch=15)
}



plot(meds[6:10], cex=.1, ylim=c(-0.0001, -0.008))

for (i in 1:5){
  arrows(i,HDI[1,i+5],
         i,HDI[2,i+5],  lwd=1, length=0.05, angle=90, code=3, col=Acols[i])
}

for (i in 1:5){
  points(c(i) ,meds[i], col=Acols[i], pch=15)
}
##########  (STOP) some kinks still need to be worked out here (STOP)  ############

##########   
##########   
###########  
##########   

##########   HDI posterior predictive check, much more powerful than median estimates as shows uncertainty
source("Modeling_scripts/HDIofMCMC.R")  ### this function is from John Kruschke's book
                    ####### Doing Bayesian Data Analysis, Second Edition: 
                          ########    A Tutorial with R, JAGS, and Stan.   


chainLength=nrow(mcmcChain)
###"Kc"     "Ko"     "Rd"     "Vcmax"  "gammaS" "gm"     "tau"    "sigma"
xPostPred = as.vector(D[D$L_sec=="TTR",]$PARi)
# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.
yPostPred = matrix( 0 , nrow=length(xPostPred), ncol=chainLength )
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
# Generate posterior predicted y values.
for ( chainIdx in 1:chainLength ) {
  yPostPred[,chainIdx] = rnorm( length(D[D$L_sec=="TTR",]$PARi) ,
                                mean = ExD_Fun(as.numeric(matrix(cbind(mcmcChain[chainIdx,1],mcmcChain[chainIdx,5] ))),
                                                as.data.frame(D[D$L_sec=="TTR",]$PARi) )[,1],
                                sd = rep( mcmcChain$sigma[chainIdx] , length(xPostPred) ) )
}

for ( xIdx in 1:length(D[D$L_sec=="TTR",]) ) {
  yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}





chainLength=nrow(mcmcChain)
str(mcmcChain)
for(i in 1:5){
  
  chainLength=nrow(mcmcChain)
  ###"Kc"     "Ko"     "Rd"     "Vcmax"  "gammaS" "gm"     "tau"    "sigma"
  xPostPred = as.vector(D[D$L_sec==Leafs[i],]$PARi)
  # Define matrix for recording posterior predicted y values at each x value.
  # One row per x value, with each row holding random predicted y values.
  yPostPred = matrix( 0 , nrow=length(xPostPred), ncol=chainLength )
  # Define matrix for recording HDI limits of posterior predicted y values:
  yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
  # Generate posterior predicted y values.
  for ( chainIdx in 1:chainLength ) {
    yPostPred[,chainIdx] = rnorm( length(D[D$L_sec==Leafs[i],]$PARi) ,
                                  mean = ExD_Fun(as.numeric(matrix(cbind(mcmcChain[chainIdx,1],mcmcChain[chainIdx,i+5] ))),
                                                 as.data.frame(D[D$L_sec=="i",]$PARi) )[,1],
                                  sd = rep( mcmcChain$sigma[chainIdx] , length(xPostPred) ) )
  }
  
  for ( xIdx in 1:length(D[D$L_sec==Leafs[i],]) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
  }
}
  
  # Plot data values:
  plot( ACi[ACi$geno==geno[i],]$PARi , ACi[ACi$geno==geno[i],]$ETR ,lwd=1.5, cex=.8, col="goldenrod", ylim=c(-10, 350) , xlim=c(0, 2100) ,
        xlab="" , ylab="" , xaxt='n',  yaxt='n' )
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=c(genocol[i]))
  
  axis(2,at=seq(0,500, by=50), tck=-0.02, cex.axis=0.7,  las=2, lwd=1.2)
  axis(2,at=seq(0,500, by=50), tck=0.02, labels=F, lwd=1.2)
  axis(1,at=seq(0,2500, by=250), tck=0.02, cex.axis=0.7, lwd=1.2)
  axis(1,at=seq(0,2500, by=250), tck=-0.02, labels=F, lwd = 1.2)
  
  mtext(expression("ETR"~mu~"mol m"^{-2}~"s"^{-1}~")"), side=2, line=1.8,cex=.8)
  mtext(expression("PAR"[i]~ "("~mu~"mol m"^{-2}~"s"^{-1}~")"), side=1, line=2.4,cex=.8)
  lo<-loess(yHDIlim[,1] ~ xPostPred)
  lines(sort(xPostPred), sort(predict(lo)), lwd=2, col="goldenrod" )
  lo2<-loess(yHDIlim[,2] ~ xPostPred)
  lines( sort(xPostPred),sort(predict(lo2)), lwd=2, col="goldenrod"  )
  #########  Add confidence intercal for Jm
  
  # change baed on # of parameters varying by genotpye
  nvar<-8;nstatic<-10; ntv= c(nvar+nstatic)
  chainLength=nrow(mcmcChain_Jm)
  ###"Kc"     "Ko"     "Rd"     "Vcmax"  "gammaS" "gm"     "tau"    "sigma"
  xPostPred = as.vector(ACi[ACi$geno==geno[i],]$PARi)
  # Define matrix for recording posterior predicted y values at each x value.
  # One row per x value, with each row holding random predicted y values.
  yPostPred2 = matrix( 0 , nrow=length(xPostPred), ncol=chainLength )
  # Define matrix for recording HDI limits of posterior predicted y values:
  yHDIlim2 = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
  # Generate posterior predicted y values.
  for ( chainIdx in 1:chainLength ) {
    yPostPred2[,chainIdx] = rnorm( length(dat_Jm[[i]][,1]) ,
                                   mean = farQ_CiCc_Jm_temp(as.data.frame(dat_Jm[[i]]),
                                                            as.numeric(as.data.frame(gList_Jm[[i]])[chainIdx,1:c(ntv-1)]))[,4],
                                   sd = rep( as.data.frame(gList_Jm[[i]])[chainIdx,ntv], length(xPostPred) ) )
  }
  
  for ( xIdx in 1:length(dat_Jf[[i]][,1]) ) {
    yHDIlim2[xIdx,] = HDIofMCMC( yPostPred2[xIdx,] )
  }
  
  # ADD Aj simulation points with CIs:
  
  lo_2<-loess(yHDIlim2[,1] ~ xPostPred)
  lines(sort(xPostPred), sort(predict(lo_2)), lwd=2, lty=2, col="darkseagreen" )
  lo2_2<-loess(yHDIlim2[,2] ~ xPostPred)
  lines( sort(xPostPred),sort(predict(lo2_2)), lwd=2, lty=2, col="darkseagreen"  )
  
  text(200,325, geno[i])
  
  
  
