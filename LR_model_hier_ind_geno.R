###############################
###      Bayesian Model    ####
### For light response curve ####
#### with hierarchical structure
######  & no temp dependency
#    works for two level hierarch ind and genotpye ###
###############################
LR_model <- "model {
for (i in 1:N){
An[i] ~ dnorm( mu.A[i] , tau )
#######  light responce model based on Ogren PLant
mu.A[i] <-  (Q[i]*phiCO2[I[i]] + Amax[I[i]] - sqrt((phiCO2[I[i]]*Q[i] +Amax[I[i]])^2-4*Q[i]*phiCO2[I[i]]*thetaJ*Amax[I[i]]))/(2*thetaJ)-Rd[I[i]]
}

#mu.A[i] <-  ((Q[i]*phiCO2 + Amax[I[i]] - sqrt((phiCO2*Q[i] +Amax[I[i]])^2-4*Q[i]*phiCO2*thetaJ*Amax[I[i]]))/(2*thetaJ))-Rd
#}
for (k in 1:N_I){
    Amax[k] ~ dnorm(AmaxG[G[k]], tau.Amax)
    Rd[k] ~ dnorm(RdG[G[k]], tau.Rd)
    phiCO2[k] ~ dnorm(phiCO2G[G[k]], tau.phiCO2)
}

for (g in 1:N_G){
AmaxG[g] ~ dnorm(mu.Amax, tau.Amax)
RdG[g] ~ dnorm(mu.Rd, tau.Rd)
phiCO2G[g] ~ dnorm(mu.phiCO2, tau.phiCO2)
#thetaJ[g] ~ dnorm(mu.thetaJ, tau.thetaJ)
}

mu.Amax ~ dnorm(18, .57)T(0,1000)
tau.Amax ~ dgamma(.001 , .001)#dnorm(5, .1)T(0,1000)
mu.Rd ~ dnorm(2, 10)
tau.Rd ~ dgamma(.1 , 1)
mu.phiCO2 ~ dnorm(0.044, 6.714897e+04)T(0,1)
tau.phiCO2 ~ dnorm(1.489226e-05, 1000000)
#mu.thetaJ ~ dunif(0,1)
#tau.thetaJ ~ dnorm(1.489226e-05, 1000000)
#phiCO2 ~ dgamma(.5 , 10)
#thetaJ ~ dnorm(.62, 1.416861e+02)T(0,1)
#Rd ~ dnorm(2, 10)
#Amax ~ dnorm(28, .57)T(0,1000)
#Rd ~ dnorm(3.1, 2.9)T(0,100)

tau ~ dgamma(.001 , .001)
}
" 