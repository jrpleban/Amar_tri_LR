Phi2_decay_model_hier <-
  "model {
for (i in 1:N){
Phi2[i] ~ dnorm( mu.Phi2[i] , tau )
mu.Phi2[i] <- (alphaP[Plant[i]]-kappaP[Plant[i]]) * (exp(betaP[Plant[i]] * (PARi[i])))+kappaP[Plant[i]]
}
for (g in 1:N_GT){
    alpha[g] ~ dnorm(mu.alpha,tau.alpha)
    beta[g] ~ dnorm(mu.beta,tau.beta) 
    kappa[g] ~ dnorm(mu.kappa,tau.kappa)T(0,200)
}
## plant level paramters
for (p in 1:Nplant){
alphaP[p]~ dnorm(alpha[GT[p]], tau.alphaP)
betaP[p] ~ dnorm(beta[GT[p]], tau.betaP) 
kappaP[p]~ dnorm(kappa[GT[p]], tau.kappaP)T(0,200)
}
mu.alpha ~ dnorm(0.75, 1/(0.05^2))
mu.beta~ dnorm(-0.002,1/(0.001^2))#Uninformative
mu.kappa ~ dnorm(0.07, 1/(.01^2))T(0,200)
tau.alpha <- pow(sigmaA,-2)
sigmaA ~ dnorm(0.05,1/(0.01^2))T(0,10)
tau.beta <-  pow(sigmaB,-2)
sigmaB ~ dnorm(0.001,1/(0.0005^2))T(0,10)
tau.kappa <-  pow(sigmaK,-2)
sigmaK ~ dnorm(0.1,1/(0.1^2))T(0,10)
tau.alphaP <- pow(sigmaAP,-2)
sigmaAP ~ dnorm(0.05,1/(0.01^2))T(0,10)
tau.betaP <-  pow(sigmaBP,-2)
sigmaBP ~ dnorm(0.001,1/(0.0005^2))T(0,10)
tau.kappaP <- pow(sigmaKP,-2)
sigmaKP ~ dnorm(0.1,1/(0.1^2))T(0,10)
tau ~ dgamma(0.001, 0.001)
}
"

