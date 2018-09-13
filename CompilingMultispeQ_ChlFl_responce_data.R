### Photosynthesis Research using MultispeQ
##  compiling data set
# R 3.2.1
# Updated: 09_11_2018. Jonathan R Pleban UB Geography  
# Reviewed by: No one yet

### This script compiles data from MultispeQ rapid light responce protocols
################################################
# 
# Collaborators: J Berry and C Mure UB Biology
#                 L. Guadagno & B. Ewers U Wyoming Botany
#                 D.S Mackay UB Geography

###  packages used 
require("stringr")
require(plyr)

## set as source directory
setwd("~/Desktop/Amaranth_MS/PhotosynQ_modeling_R")

####  loading data from saved files 
# load low light portion of curve
P1 <- read.csv("data/Rapid_WY_0_125_raw.csv")
names(P1)
str(P1)
### duplicate each obs 5x and add PAR column
P1p<-P1[rep(1:nrow(P1),each=5),] 
P1p$PARi<-rep(c(0,15,30,60,125),435)
# load high light portion of curve
P2 <- read.csv("data/Rapid_WY_250_2000_raw.csv")
names(P2)
str(P2)
### duplicate each obs 5x and add PAR column
P2p<-P2[rep(1:nrow(P2),each=5),] 
P2p$PARi<-rep(c(250,500,1000, 1500,2000),435)

## merge two sections of LR
fastmerge <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)
  
  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)
  
  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)
  
  # add blank columns to d2
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }
  
  # add blank columns to d1
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }
  
  return(rbind(d1, d2))
}

## merge all three protocals into one data frame
P<-fastmerge(P1p, P2p)

plot(P$PARi)
names(P)


## process variables of interest against PAR
P$PHI2 <- ifelse(P$PARi == 0, P$Phi2_0,
                 ifelse(P$PARi == 15, P$Phi2_15,
                        ifelse(P$PARi == 30, P$Phi2_30,
                               ifelse(P$PARi == 60, P$Phi2_60,
                                      ifelse(P$PARi == 125, P$Phi2_125,
                                             ifelse(P$PARi == 250, P$Phi2_250,
                                                    ifelse(P$PARi == 500, P$Phi2_500,
                                                           ifelse(P$PARi == 1000, P$Phi2_1000,
                                                                  ifelse(P$PARi == 1500, P$Phi2_1500,
                                                                         ifelse(P$PARi == 2000, P$Phi2_2000, NA))))))))))
P$PHINPQt <- ifelse(P$PARi == 0, P$PhiNPQ_0,
                    ifelse(P$PARi == 15, P$PhiNPQ_15,
                           ifelse(P$PARi == 30, P$PhiNPQ_30,
                                  ifelse(P$PARi == 60, P$PhiNPQ_60,
                                         ifelse(P$PARi == 125, P$PhiNPQ_125,
                                                ifelse(P$PARi == 250, P$PhiNPQ_250,
                                                       ifelse(P$PARi == 500, P$PhiNPQ_500,
                                                              ifelse(P$PARi == 1000, P$PhiNPQ_1000,
                                                                     ifelse(P$PARi == 1500, P$PhiNPQ_1500,
                                                                            ifelse(P$PARi == 2000, P$PhiNPQ_2000, NA))))))))))


plot(P$PARi, P$PHINPQt)


P$PHINOt <- ifelse(P$PARi == 0, P$PhiNO_0,
                   ifelse(P$PARi == 15, P$PhiNO_15,
                          ifelse(P$PARi == 30, P$PhiNO_30,
                                 ifelse(P$PARi == 60, P$PhiNO_60,
                                        ifelse(P$PARi == 125, P$PhiNO_125,
                                               ifelse(P$PARi == 250, P$PhiNO_250,
                                                      ifelse(P$PARi == 500, P$PhiNO_500,
                                                             ifelse(P$PARi == 1000, P$PhiNO_1000,
                                                                    ifelse(P$PARi == 1500, P$PhiNO_1500,
                                                                           ifelse(P$PARi == 2000, P$PhiNO_2000, NA))))))))))


plot(P$PARi, P$PHINOt)


P$NPQt <- ifelse(P$PARi == 0, P$NPQt_0,
                 ifelse(P$PARi == 15, P$NPQt_15,
                        ifelse(P$PARi == 30, P$NPQt_30,
                               ifelse(P$PARi == 60, P$NPQt_60,
                                      ifelse(P$PARi == 125, P$NPQt_125,
                                             ifelse(P$PARi == 250, P$NPQt_250,
                                                    ifelse(P$PARi == 500, P$NPQt_500,
                                                           ifelse(P$PARi == 1000, P$NPQt_1000,
                                                                  ifelse(P$PARi == 1500, P$NPQt_1500,
                                                                         ifelse(P$PARi == 2000, P$NPQt_2000, NA))))))))))


plot(P$PARi, P$NPQt)


P$qL <- ifelse(P$PARi == 0, P$qL_0,
               ifelse(P$PARi == 15, P$qL_15,
                      ifelse(P$PARi == 30, P$qL_30,
                             ifelse(P$PARi == 60, P$qL_60,
                                    ifelse(P$PARi == 125, P$qL_125,
                                           ifelse(P$PARi == 250, P$qL_250,
                                                  ifelse(P$PARi == 500, P$qL_500,
                                                         ifelse(P$PARi == 1000, P$qL_1000,
                                                                ifelse(P$PARi == 1500, P$qL_1500,
                                                                       ifelse(P$PARi == 2000, P$qL_2000, NA))))))))))


plot(P$PARi, P$qL)

P$qP <- ifelse(P$PARi == 0, P$qP_0,
               ifelse(P$PARi == 15, P$qP_15,
                      ifelse(P$PARi == 30, P$qP_30,
                             ifelse(P$PARi == 60, P$qP_60,
                                    ifelse(P$PARi == 125, P$qP_125,
                                           ifelse(P$PARi == 250, P$qP_250,
                                                  ifelse(P$PARi == 500, P$qP_500,
                                                         ifelse(P$PARi == 1000, P$qP_1000,
                                                                ifelse(P$PARi == 1500, P$qP_1500,
                                                                       ifelse(P$PARi == 2000, P$qP_2000, NA))))))))))


plot(P$PARi, P$qP)

### might be able to simplify this ???
P$FvP.Fm <- ifelse(P$PARi == 0, P$FvP.Fm_0_P,
               ifelse(P$PARi == 15, P$FvP.Fm_15_P,
                      ifelse(P$PARi == 30, P$FvP.Fm_30_P,
                             ifelse(P$PARi == 60, P$FvP.Fm_60_P,
                                    ifelse(P$PARi == 125, P$FvP.Fm_125_P,
                                           ifelse(P$PARi == 250, P$FvP.Fm_250_P,
                                                  ifelse(P$PARi == 500, P$FvP.Fm_500_P,
                                                         ifelse(P$PARi == 1000, P$FvP.Fm_1000_P,
                                                                ifelse(P$PARi == 1500, P$FvP.Fm_1500_P,
                                                                       ifelse(P$PARi == 2000, P$FvP.Fm_2000_P, NA))))))))))
plot(P$PARi, P$FvP.Fm)
names(P)



P$LEF<-P$PARi*P$PHI2
plot(P$PARi, P$NPQt)

plot(P$PARi, P$LEF)



write.csv(P,"data/Full_set_MultispeQ_Phi2_PARLresponce_Guadagno_Amaranth_Rapid.csv")




Leafs<-c("TUG","TUR", "TTG","TTR","TTY")
colnames(PREDS)<-c(Leafs, "PAR")

names(P)

Acols<-c("green1", "red1", "seagreen2", "deeppink4" ,"yellow")

plot(P[P$Leaf_Section_Type==Leafs[1],]$PARi, P[P$Leaf_Section_Type==Leafs[1],]$LEF, col=Acols[1], ylim=c(0,450))
plot(P[P$Leaf_Section_Type==Leafs[2],]$PARi, P[P$Leaf_Section_Type==Leafs[2],]$LEF, col=Acols[2], ylim=c(0,450))
plot(P[P$Leaf_Section_Type==Leafs[3],]$PARi, P[P$Leaf_Section_Type==Leafs[3],]$LEF, col=Acols[3], ylim=c(0,450))
plot(P[P$Leaf_Section_Type==Leafs[4],]$PARi, P[P$Leaf_Section_Type==Leafs[4],]$LEF, col=Acols[4], ylim=c(0,450))
plot(P[P$Leaf_Section_Type==Leafs[5],]$PARi, P[P$Leaf_Section_Type==Leafs[5],]$LEF, col=Acols[5], ylim=c(0,450))

plot(P[P$Leaf_Section_Type==Leafs[1],]$PARi, P[P$Leaf_Section_Type==Leafs[1],]$LEF, col=Acols[1])
for (i in 1:4){
  points(P[P$Leaf_Section_Type==Leafs[i+1],]$PARi, P[P$Leaf_Section_Type==Leafs[i+1],]$LEF, col=Acols[i+1])
}
