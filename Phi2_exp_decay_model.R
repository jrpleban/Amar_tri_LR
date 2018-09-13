Phi2_decay_model_hier <-
  "model {
for (i in 1:N){

Phi2[i] ~ dnorm( mu.Phi2[i] , tau )
mu.Phi2[i] <- alpha[L_T_S[i]] * (exp(beta[L_T_S[i]] * (PARi[i])))

}
for (g in 1:N_L_T_S){
beta[g]~ dnorm(mu.beta, tau.beta)
alpha[g] ~ dnorm(mu.alpha, tau.alpha)T(0,1)
}

mu.alpha ~ dnorm(0.4, 0.1)T(0,1)
mu.beta~ dnorm(-0.003,20)#T(0,-0.01) #Uninformative
tau.alpha ~ dnorm(.2, 1)
tau.beta~ dnorm(1.0,1) #Uninformative
tau ~ dgamma(0.001, 0.001)
}
"

