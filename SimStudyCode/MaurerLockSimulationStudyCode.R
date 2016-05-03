### Full Simulation ###
# Violate independence in Bivariate MANCOVA
# assess effect on type I error rate

## Preliminaries
setwd("C:\\Users\\Karsten\\Dropbox\\Dissertation\\CurriculumStudy\\NewSimulationStudy")
library(ggplot2)
library(reshape2)
library(MASS)
library(Matrix)
library(plyr)
library(gridExtra)

dat <- read.csv("DataForStudyAnalysis.csv",header=T)
#change factor level order to put sweeney as default control group
dat$room <- factor(dat$room.x, levels=c("sweeney", "kildee"))
dat$lab5perc <- dat$lab5perc * 100
dat$hw2perc <- dat$hw2perc * 100
dat$treatment <- "Simulation-Based"
dat$treatment[which(dat$room=="sweeney")] <- "Traditional"
n= nrow(dat)



### 
# Cohorts in real data
table(dat$Section, dat$treatment)



# fit model to data using all effects except treatment to get 
# estimates for "true" effects of midterm and lab5 to put into "truth" for simulation
truemod <-  lm(cbind(ConfMC,HypMC)~ midterm + lab5perc, data=dat)

### Create function to simulate M times, with 

# Covariance within students
# cov(y_ijkl, y_ijkl') = sigma_ll' = sigma_ll' * (1 + z + g + d)
# sigma_ll' is estimated from independence model

# sigma_ll' * z is cov for all students on the l and l' pair of responses 
# sigma_ll' * g is added if students are classmates
# sigma_ll' * d is added if students are labmates
# sigma_ll' is added if the covariance is within a student

SimViolationNoTrtEffect <- function(M, z, g, d, truemod){
  # Use fitted coefficients, b, as "true" betas with beta for treatment forced to 0
  # then for any given set of possible errors Y = X'b + errors
  truebetas <-  rbind(truemod$coef, c(0,0))
  xb <- truemod$fitted.values
  
  # variance covariance estimates from model assuming indepence to get plug-in
  #   estimators for within student var/cov elements sigma_ll' to use in sim
  trueresids <- manova(truemod)$residuals
  trueSig <- t(trueresids)%*%trueresids / (101-3-1)
  
  # Simulation Procedure:
  # randomly assign student to cohort (classroom and labroom)
  # Use room structure to set up independence violated variance structure
  # take draw error (errm) from violated variance structure
  # take Ym = xb + errm
  # fit model assuming independence to Ym's
  # record betam's, sigmam's and SE[betam's] into keepers dataframe
  # Repeat Simulation M times
  
  keepers <- data.frame( betamidterm1s = rep(NA, M), betamidterm2s = rep(NA, M),
                         SEbetamidterm1s =  rep(NA, M), SEbetamidterm2s =  rep(NA, M),       
                         pvalmidterm1s =  rep(NA, M), pvalmidterm2s =  rep(NA, M),
                         betalab1s = rep(NA, M), betalab2s = rep(NA, M),
                         SEbetalab1s =  rep(NA, M), SEbetalab2s =  rep(NA, M),                       
                         pvallab1s =  rep(NA, M), pvallab2s =  rep(NA, M),
                         beta1s = rep(NA, M), beta2s = rep(NA, M),
                         SEbeta1s =  rep(NA, M), SEbeta2s =  rep(NA, M),
                         pval1s =  rep(NA, M), pval2s =  rep(NA, M),
                         pillais = rep(NA,M), pillaipvals = rep(NA,M),
                         sig11 = rep(NA,M) , sig12  = rep(NA,M), sig22= rep(NA,M))
  smalldat = dat[,c("midterm","lab5perc")]
  for (m in 1:M){
    # Simulate student assignment to section / treatment cohorts
    #   Real data has corhort counts that we will mimic 
    smalldat$cohortm <- sample(c(rep("A",22),rep("B",28),rep("C",24),rep("D",27)), 101, replace=FALSE)
    smalldat$room <- "Simulation-Based" ; smalldat$room[smalldat$cohortm %in% c("C","D")] <- "Traditional" 
    
    # Set var/cov matrix using nested variability
    #   build row by row based on if row and column share classroom, cohort and student
    roomsm <- rep(smalldat$room,2) 
    cohortm <- rep(smalldat$cohortm,2) 
    student <- rep(1:n,2)
    n=nrow(dat)
    #matrix of additive multiplier for sigma_ll'
    Varm <- matrix( rep(1, 4*n^2) , nrow=2*n, byrow=TRUE )
    for (i in 1:(2*n)) {
      Varm[i,which(roomsm==roomsm[i])] <- Varm[i,which(roomsm==roomsm[i])] + g
      Varm[i,which(cohortm==cohortm[i])] <- Varm[i,which(cohortm==cohortm[i])] + d
      Varm[i,which(student==student[i])] <- Varm[i,which(student==student[i])] + 1
    }
    # scale each cov by appropriate sigma_ll'
    for (i in 1:n) {
      Varm[i,1:n] <- Varm[i,1:n] * trueSig[1,1]
      Varm[i,(n+1):(2*n)] <- Varm[i,(n+1):(2*n)] * trueSig[1,2]
      Varm[(i+n),1:n] <- Varm[(i+n),1:n] * trueSig[2,1]
      Varm[(i+n),(n+1):(2*n)] <- Varm[(i+n),(n+1):(2*n)] * trueSig[2,2]        
    }
    errm <-  mvrnorm(1, rep(0, 2*n), Varm)
    ym <- xb + matrix(errm,ncol=2,byrow=F)
    smalldat$ym1 <- ym[,1]
    smalldat$ym2 <- ym[,2]
    modm <- lm(cbind(ym1,ym2)~ midterm + lab5perc + room, data=smalldat)
    modmresids <- manova(modm)$residuals
    modmSig <- t(modmresids)%*%modmresids / (101-3-1)
    keepers[m,c("betamidterm1s","SEbetamidterm1s","pvalmidterm1s")] <- summary(modm)[[1]][[4]][2,c(1,2,4)]
    keepers[m,c("betamidterm2s","SEbetamidterm2s","pvalmidterm2s")] <- summary(modm)[[2]][[4]][2,c(1,2,4)]
    keepers[m,c("betalab1s","SEbetalab1s","pvallab1s")] <- summary(modm)[[1]][[4]][3,c(1,2,4)]
    keepers[m,c("betalab2s","SEbetalab2s","pvallab2s")] <- summary(modm)[[2]][[4]][3,c(1,2,4)]
    keepers[m,c("beta1s","SEbeta1s","pval1s")] <- summary(modm)[[1]][[4]][4,c(1,2,4)]
    keepers[m,c("beta2s","SEbeta2s","pval2s")] <- summary(modm)[[2]][[4]][4,c(1,2,4)]
    keepers[m,c("pillais","pillaipvals")] <- summary(manova(modm))[[4]][3,c(2,6)]
    keepers[m,c("sig11","sig12","sig22")] <- c(modmSig[1,1],modmSig[1,2],modmSig[2,2])
  } 
  return(keepers)
}

### Try Simulations for  tau=.5, gamma=.05, delta=.05, eta=.4
timer <- Sys.time()
set.seed(62668)
M= 100
truemod <-  lm(cbind(ConfMC,HypMC)~ midterm + lab5perc, data=dat)
simvalues <- SimViolationNoTrtEffect(M, g=0, d=.1, z=0, truemod)
Sys.time()-timer

length(which(simvalues$pillaipvals < 0.05))/M



########################################################################################

#Create a set of variance parameter values to simulate over
simseq <- data.frame(g=c(0,rep(seq(0,.1,by=.02),6)), d=c(0,rep(seq(0,.1,by=.02),each=6)), z=0) 
head(simseq)

M= 20000
L = nrow(simseq)
shortsim <- data.frame(g=rep(NA,L),d=rep(NA,L),z=rep(NA,L),
                       T1CI=rep(NA,L),T1HT=rep(NA,L),T1pillai=rep(NA,L)) 
resultlist <- vector("list",L)
counter <- 1
timer <- Sys.time()
for (i in 1:nrow(simseq)){
  set.seed(62668)
  truemod <-  lm(cbind(ConfMC,HypMC)~ midterm + lab5perc, data=dat)
  simvalues <- SimViolationNoTrtEffect(M, g=simseq$g[i], d=simseq$d[i], z=simseq$z[i], truemod)
  print(Sys.time()-timer)  # 
  print(i)
  shortsim[counter, ] <- c(simseq$g[i],  simseq$d[i], simseq$z[i],
                           length(which(simvalues$pval1s <= .05))/M,
                           length(which(simvalues$pval2s <= .05))/M,
                           length(which(simvalues$pillaipvals <= .05))/M)
  resultlist[[counter]] <- simvalues
  counter <- counter+1
}


# setwd(C:\\Users\\Karsten\\Dropbox\\Dissertation\\CurriculumStudy\\NewSimulationStudy)
# save(resultlist, file = "NewSimulationResultsReparam.RData")
# write.csv(shortsim,"NewSimShortResultsReparam.csv",row.names=F)

qplot(eta, T1pillai, data=shortsim)





#### Begin Initial Data Visualization Code ####
# setwd("C:\\Users\\Karsten\\Dropbox\\Dissertation\\CurriculumStudy\\NewSimulationStudy")
load("NewSimulationResultsReparam.Rdata")
#str(resultlist[[1]])
shortsim <- read.csv("Sim5000NewResults.csv")
#head(shortsim)

L <- length(resultlist)
M = nrow(resultlist[[1]]) 
errorsummary2 <- data.frame(T2ErrorRateBetaMidterm1 = rep(NA, L),
                            T2ErrorRateBetaMidterm2 = rep(NA, L),
                            T2ErrorRateBetaLab1 = rep(NA, L),
                            T2ErrorRateBetaLab2 = rep(NA, L),
                            T1ErrorRateBeta1 = rep(NA, L),
                            T1ErrorRateBeta2 = rep(NA, L),
                            T1ErrorRatePillai = rep(NA, L))

alpha <- 0.05
for (l in 1:L){
  errorsummary2[l,1] <- length(which(resultlist[[l]]$pvalmidterm1s > alpha)) / M
  errorsummary2[l,2] <- length(which(resultlist[[l]]$pvalmidterm2s > alpha)) / M
  errorsummary2[l,3] <- length(which(resultlist[[l]]$pvallab1s > alpha)) / M
  errorsummary2[l,4] <- length(which(resultlist[[l]]$pvallab2s > alpha)) / M
  errorsummary2[l,5] <- length(which(resultlist[[l]]$pval1s <= alpha)) / M
  errorsummary2[l,6] <- length(which(resultlist[[l]]$pval2s <= alpha)) / M
  errorsummary2[l,7] <- length(which(resultlist[[l]]$pillaipvals <= alpha)) / M
}
simseq <- data.frame(g=c(0,rep(seq(0,.1,by=.02),6)), d=c(0,rep(seq(0,.1,by=.02),each=6)), z=0) 
errorsummary2 <- cbind(errorsummary2,simseq)
#head(errorsummary2,10)
#write.csv(errorsummary2, "SimulationStudyErrorSummary.csv", row.names=FALSE)
errorsummary2 <- read.csv("SimulationStudyErrorSummary.csv")

### Type 1 error plots: 
# error increases as delta/gamma (class/labmate effects) increase in
# relation to the within student eta.  Tau held to zero for all
p1 <- qplot(g,T1ErrorRatePillai, data=errorsummary2, geom="line", group=d, 
            color=d, size=I(.7)) + geom_hline(yintercept=alpha)+
  xlab("g") + theme_bw() + ylim(c(0,1)) 

p2 <- qplot(g,T1ErrorRateBeta1, data=errorsummary2, geom="line", group=d, 
            color=d, size=I(.7)) + geom_hline(yintercept=alpha)+
  xlab("g") + theme_bw() + ylim(c(0,1)) + 
  theme(legend.position="none")

p3 <- qplot(g,T1ErrorRateBeta2, data=errorsummary2, geom="line", group=d, 
            color=d, size=I(.7)) + geom_hline(yintercept=alpha)+
  xlab("g") + theme_bw() + ylim(c(0,1)) + 
  theme(legend.position="none")

grid.arrange(p2,p3,p1, nrow=1, widths=c(1,1,1.5))

### beta estimates expectation (check for bias)
# true betas in generative model
truebetas <- rbind(truemod$coef, c(0,0))
trueresids <- manova(truemod)$residuals
trueSig <- t(trueresids)%*%trueresids / (101-3-1)

# average fitted value in simulations under each variance parameterization
head(resultlist[[1]])
betaMeans <- ldply(resultlist, summarize,
                   betamidterm1 = mean(betamidterm1s),
                   betamidterm2 = mean(betamidterm2s),
                   betalab1 = mean(betalab1s),
                   betalab2 = mean(betalab2s),
                   beta1 = mean(beta1s),
                   beta2 = mean(beta2s))
# Compare true to average fitted
summary(betaMeans)
truebetas
head(betaMeans,10)

### fitted coefficient variabilities decrease as the eta decreases
betaSds <- ldply(resultlist, summarize,
                 SDbetamidterm1 = sd(betamidterm1s),
                 SDbetamidterm2 = sd(betamidterm2s),
                 SDbetalab1 = sd(betalab1s),
                 SDbetalab2 = sd(betalab2s),
                 SDbeta1 = sd(beta1s),
                 SDbeta2 = sd(beta2s), 
                 sig1 = mean(sig11))
betaSds <- cbind(betaSds,simseq)
qplot(1-eta, betaSds$SDbetamidterm1, data=simseq, color=LabPropOfBetweenVar ) + 
  geom_smooth() + theme_bw()

qplot(d+g , sig1 , data=betaSds, color=d, size=g)

# how often was my simulation getting an effect size as large as what I saw?
errorsummary3 <- data.frame(EffectAsBig1 = rep(NA, L),
                            EffectAsBig2 = rep(NA, L))

for (l in 1:L){
  errorsummary3[l,1] <- length(which(resultlist[[l]]$pval1s < 0.03714 )) / M
  errorsummary3[l,2] <- length(which(resultlist[[l]]$pvalmidterm2s < 0.38998)) / M
}

head(cbind(errorsummary3,simseq),10)
