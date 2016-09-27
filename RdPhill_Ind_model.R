RdPhi2ll_mod <- "model
{
    for (j in 1:N_ll){
    A_ll[j] ~ dnorm( mu.A_ll[j] , tau_Rd )
    mu.A_ll[j] <- s1 * (Inc_ll[j]*phi2_ll[j])/4  - Rd25

}
  
for (k in 1:N_ll){
 phi2_ll[k] ~ dnorm(mu.phi2_ll[k] , tau_phi )
 mu.phi2_ll[k] <- mphi * Inc_ll[k]  + phi2ll
}

 
 
 Rd25 ~ dunif(1,5)
 s1 ~ dunif(.2,.7)
 phi2ll ~dunif(.1,.9)
 mphi ~ dunif(-1,1)

 ### precision for Rd from low light
 tau_Rd ~ dgamma(.001 , .001)
 tau_phi~ dgamma(.001 , .001)

} "