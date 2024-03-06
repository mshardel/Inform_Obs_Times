###########################################################
#Simulation study for 
#"Time-Varying Latent Eect Models for Repeated Measurements to Address
#Informative Observation Times in the U.S. Medicare Minimum Data Set"   
#By Shardell M, Chen C, and Falvey J
#University of Maryland School of Medicine

#Three outcome model covariates of interest: X1 and X2 are time-invariant, X3 is time-varying
#Four observation model covariates: Z1=X1, Z2=X2, Z3=X3, and Z4 is time-varying
#Comparing proposed inverse conditional intensity rate ratio (ICIRR)
#method to Buzkova and Lumley (2009) (BL)

###########################################################

##directory of simulation results
res.dir <- "./wt_sims_results/"

## Load packages
library(plyr)
library(nleqslv)
library(MASS)


###number of simulations:
n.sim <- 500

######################vector of possible sample sizes
nn.vec <- c(250, 500)

#####################specify vector of possible times
tau <-2 #max follow-up time
time_discrete <- seq(0.02,tau,0.02)
n.times <- length(time_discrete)

###################specify alpha.10(t)
#alpha.10(t) = 1 + sin(t)
alpha.10.t <- 1 + sin(time_discrete)

##################possible functions for alpha.20(t)
alpha.20.t.mat <- cbind((1/(1+time_discrete)), sin(time_discrete), rep(0, length(time_discrete)) )
alpha.20.t.text.vec <- c("1.over.(1+t)", "sin(t)","0")

##################possible values for delta4 for Z4, (auxiliary observation-model covariate)
delta4.vec <- c(-0.5, 0, 0.5)

##################iterate over different sample size
for(ii.nn in 1:length(nn.vec)){

#sample size per simulation:
  nn <- nn.vec[ii.nn]

################iterate over alpha.20(t) and gamma3)

  for(ii.alpha2 in 1:ncol(alpha.20.t.mat)){

    alpha.20.t <- alpha.20.t.mat[,ii.alpha2]
    alpha.20.t.text <- alpha.20.t.text.vec[ii.alpha2]

    for(ii.delta4 in 1:length(delta4.vec)){

      delta4 <- delta4.vec[ii.delta4]


############storage matrix for beta (3 covariates) (ICIRR)
      beta.ICIRR.mat <- matrix(c(0), nrow=n.sim, ncol=3)
      beta.est.ICIRR.filename <- paste0(res.dir, "beta_ests_Z_ICIRR_n_",nn,"_a2_",alpha.20.t.text,"_delta4_",delta4, ".csv")

############storage matrix for Buzkova-Lumley beta 
      beta.BL.mat <- matrix(c(0), nrow=n.sim, ncol=3)
      beta.est.BL.filename <- paste0(res.dir, "beta_ests_Z_BL_n_",nn,"_a2_",alpha.20.t.text,"_delta4_",delta4, ".csv")


###########storage matrix for delta (4 covariates, Z1=X1, Z2=X2, Z3=X3, Z4)
      delta.mat <- matrix(c(0), nrow=n.sim, ncol=4)
      delta.est.filename <- paste0(res.dir, "delta_ests_Z_n_",nn,"_a2_",alpha.20.t.text,"_delta4_",delta4, ".csv")


################################simulation parameters
########X1: P(X1=1), X1 Bernoulli
      p.X1 <- 0.5

########X2: f(X2) ~ N(mu.X2, sd.X2)
      mu.X2 <- 0
      sd.X2 <- 0.5

########X3(t): f(X3(t) | X1, X2) ~ N(mu.X3.t, sd.X3.t)
#mu.X3.t = E[X3(t) | X1, X2, t] = nu0.X3 + nu1.X3*X1 + nu2.X3*X2
      nu0.X3 <- 0
      nu1.X3 <- 0.25
      nu2.X3 <- 0.5
      sd.X3.t <- 0.5

########Z4(t): f(Z4(t) | X1, X2, X3) ~ N(mu.Z4.t, sd.Z4.t)
#mu.Z4.t = E[Z4(t) | X1, X2, X3, t] = nu0.Z4 + nu1.Z4*X1 + nu2.Z4*X2 + nu3.Z4*X3(t)
      nu0.Z4 <- 0
      nu1.Z4 <- 0
      nu2.Z4 <- 0
      nu3.Z4 <- -1
      sd.Z4.t <- 0.75

########V: V = (1-X1)*U1 + X1*U2
#U1 ~ gamma(mean=1, var=0.5) and U2 ~ gamma(mean=1, var=0.8)
      var.U1 <- 0.5
      mean.U1 <- 1
      scale.U1 <- var.U1/mean.U1
      shape.U1 <- mean.U1/scale.U1

      var.U2 <- 0.8
      mean.U2 <- 1
      scale.U2 <- var.U2/mean.U2
      shape.U2 <- mean.U2/scale.U2

########Y(t) ~ N(mu.Y.t, sd.Y)
#mu.Y.t = E[Y(t) | X1, X2, X3(t), Z4(t), V] = alpha.10.t + alpha.20.t*V + beta1.Y*X1 + beta2.Y*X2 + beta3.Y*X3(t) + beta4.Y*Z4(t)
      sd.Y <- sqrt(0.25)
      beta1.Y <- 1
      beta2.Y <- -1
      beta3.Y <- 1.5
      beta4.Y <- 0.5

###############dN(t) ~ Pois(rate.t)
#rate.t = V*lambda0.t*exp(delta1*X1 + delta2*X2 + delta3*X3.t + delta4*Z4.t) 
      lambda0.t <- 0.1 
      delta1 <- -0.5
      delta2 <- -0.25
      delta3 <- 0.25


set.seed(032922)
##########################################################################################
######################################start simulation
##########################################################################################

      for(ii.sim in 1:n.sim){

####simulate baseline covariates X1 and X2
#X1 is bernoulli with probability p.X1
         X1 <- rbinom(n=nn,size=1, prob=p.X1)

#X2 is normal with mean=mu.X2, sd=sd.X2
         X2 <- rnorm(n=nn, mean=mu.X2, sd=sd.X2)

##X3.t is normal and time-varying 
#E[X3(t) | X1, X2, t] = nu0.X3 + nu1.X3*X1 + nu2.X3*X2
         mu.X3.t <- nu0.X3 + nu1.X3*X1 + nu2.X3*X2
#rows are id, columns are time
         X3.t <- t(sapply(mu.X3.t, function(x){rnorm(n=n.times, mean=x, sd=sd.X3.t)}))

##Z4.t is normal and time-varying
#E[Z4(t) | X1, X2, X3(t), t] = nu0.Z4 + nu1.Z4*X1 + nu2.Z4*X2 + nu3.Z4*X3.t
         Z4.sim.f <-function(id){
              mu.Z4.t <- nu0.Z4 + nu1.Z4*X1[id] + nu2.Z4*X2[id] + nu3.Z4*X3.t[id,] 
              Z4 <- rnorm(length(mu.Z4.t), mean=mu.Z4.t, sd=sd.Z4.t)
              Z4
         }

#rows are id, columns are t
         Z4.t <-t(sapply(c(1:nn), Z4.sim.f))

######simulate latent variable V
## V = (1-X1)*U1 + X1*U2
#U1 ~ gamma(mean=1, var=0.5)
         U1 <- rgamma(n=nn,shape=shape.U1,scale= scale.U1)

#U2 ~ gamma(mean=1, var=0.8)
         U2 <- rgamma(n=nn,shape=shape.U2,scale= scale.U2)

         V <- (1-X1)*U1 + X1*U2


###############simulate outcome Y.t
#E[Y(t) | X1, X2, X3(t), Z4(t), V] = alpha.10.t + alpha.20.t*V + beta1.Y*X1 + beta2.Y*X2 + beta3.Y*X3(t) + beta4.Y*Z4(t)

         Y.sim.f <-function(id){
              mu.Y.t <- alpha.10.t + alpha.20.t*V[id] + beta1.Y*X1[id] + beta2.Y*X2[id] + beta3.Y*X3.t[id,] + beta4.Y*Z4.t[id,] 
              Y <- rnorm(length(mu.Y.t), mean=mu.Y.t, sd=sd.Y)
              Y
         }

#rows are id, columns are t
         Y.t <-t(sapply(c(1:nn), Y.sim.f))


###############simulate visit times dN(t)
#rate.t = V*lambda0.t*exp(delta1*X1 + delta2*X2 + delta3*X3.t + delta4*Z4.t) 

         dNi.f <- function(id){
              rate.t <- V[id]*lambda0.t*exp(delta1*X1[id] + delta2*X2[id] + delta3*X3.t[id,] + delta4*Z4.t[id,])
              rate.t <- sapply(rate.t,function(x){min(x,1)}) #to prevent probs >1
              dNi <- rbinom(size=1,n=length(rate.t),prob=rate.t)
              dNi
         }

###rows are ids, columns are t
         dNi.t <- t(sapply(c(1:nn), dNi.f))

##############simulate censoring time
         Ci <-runif(min=1,max=tau,n=nn)

###number of visits per patient
         m.f <- function(id){
              visits <- dNi.t[id,]
              cens <- Ci[id]
              sum(visits[time_discrete<=cens])
         }

         m <- sapply(c(1:nn),m.f)

#############################################complete data
         complete.data <- data.frame(ID=rep(c(1:nn),n.times), t=sort(rep(time_discrete,nn)), m=m, Y=c(Y.t),
                          X1=rep(X1, n.times), X2=rep(X2, n.times), X3=c(X3.t), Z4=c(Z4.t), C=rep(Ci, n.times), dNi=c(dNi.t))

         complete.data <- complete.data[order(complete.data$ID),]

######################################START FUNCTION############################
###############################################################################
         InformTimes.f <- function(complete.data){
#############################################################

###################################plugging in NA for Y when dNi=0 and dropping rows when t>Ci
              obs.data <- complete.data[complete.data$t<=complete.data$C,]
              obs.data$Y[obs.data$dNi==0] <- NA

####################################create observed subset with only observed rows
########################################################################
# Analysis dataset: sim.data
# ID = participant id, t = visit time
# X1 = binary baseline, X2 = continuous baseline, (confounders) 
# X3 = continuous time-dependent (exposure)
# Z4 = continuous time-dependent (auxiliary covariate for observation-time model)
# Y = continuous
# C = censoring time
# m = number of visits (id-level)
########################################################################
              sim.data <- obs.data[obs.data$dNi==1,]
              sim.data <- sim.data[, c("ID", "t", "m", "Y", "X1", "X2", "X3", "Z4", "C")]
              sim.data <- sim.data[with(sim.data, order(ID, t)), ]  #sorted by id, and by time within id
              max.time <- max(sim.data$C)

#########################################################
# Baseline dataset: baseData
# baseData = baseline data for each id
##########################################################
              baseData <- ddply(obs.data, .(ID), function(x) x[1, ])
              baseData <- baseData[, c("ID", "m", "X1", "X2", "C")]


#################################################################
# udt = vector of unique visit times in the dataset
#################################################################
              udt <- unique(sort(sim.data$t[sim.data$t>0]))  #vector of unique times with >=1 visit


##############################################################################
#### expand time-varying covariates to full set of length(udt) x nn rows ####
##############################################################################
              testdata3 <- complete.data
              testdata3$dNi[testdata3$t>testdata3$C] <- 0 
              testdata3$Y[testdata3$dNi==0] <- NA 
              testdata4 <- testdata3[(testdata3$t %in% udt),]  #only keeping times with >0 visits


  ###################################################################
  # delta.hat (observation-model covariates Z)
  ###################################################################

              delta.f <- function(delta){
                   exp_delta <- function(tt){exp(delta[1]*baseData$X1+delta[2]*baseData$X2+delta[3]*testdata4$X3[testdata4$t==tt]+delta[4]*testdata4$Z4[testdata4$t==tt])}
      
                   numer1 <- sapply(sim.data$t, function(u){sum( (baseData$X1*exp_delta(u))[u<=baseData$C], na.rm=T) } )
                   numer2 <- sapply(sim.data$t, function(u){sum( (baseData$X2*exp_delta(u))[u<=baseData$C], na.rm=T) } )
                   numer3 <- sapply(sim.data$t, function(u){sum( (testdata4$X3[testdata4$t==u]*exp_delta(u))[u<=baseData$C], na.rm=T) } )
                   numer4 <- sapply(sim.data$t, function(u){sum( (testdata4$Z4[testdata4$t==u]*exp_delta(u))[u<=baseData$C], na.rm=T) } )
                   denom <- sapply(sim.data$t, function(u){sum( (exp_delta(u))[u<=baseData$C], na.rm=T) } )
                   Vbar <- cbind(numer1/denom, numer2/denom, numer3/denom, numer4/denom)

                   bigV <-cbind(sim.data$X1, sim.data$X2, sim.data$X3, sim.data$Z4)
                   temp <- colSums((bigV-Vbar)/nn, na.rm=T)
                   temp
              }

              delta <- c(delta1, delta2, delta3, delta4)

              delta.hat <- nleqslv(delta, delta.f)$x


#################estimated dLam(t) under Z1+Z2+Z3+Z4 (Z1=X1, Z2=X2, Z3=X3) in observation-time model evaluated at solution for delta 
              exp_delta <- function(u){exp(delta.hat[1]*baseData$X1+
		           delta.hat[2]*baseData$X2+delta.hat[3]*testdata4$X3[testdata4$t==u]+delta.hat[4]*testdata4$Z4[testdata4$t==u])}

#sum over subjects of exp_delta for each time, result length of unique visit times udt
              denom_delta <- sapply(udt, function(u){sum( (exp_delta(u))[u<=baseData$C], na.rm=T) } )

#estimated dLam(t) under X1+X2+X3+Z4
              estlam.t.delta <- sapply(1:length(udt), function(u) sum( ((sim.data$t==udt[u])/denom_delta[u])) )


  ###################################################################
  # gamma.hat (outcome-model covariates X ONLY)
  ###################################################################

              gamma.f <- function(gamma){
                   exp_gamma <- function(tt){exp(gamma[1]*baseData$X1+gamma[2]*baseData$X2+gamma[3]*testdata4$X3[testdata4$t==tt])}
      
                   numer1 <- sapply(sim.data$t, function(u){sum( (baseData$X1*exp_gamma(u))[u<=baseData$C], na.rm=T) } )
                   numer2 <- sapply(sim.data$t, function(u){sum( (baseData$X2*exp_gamma(u))[u<=baseData$C], na.rm=T) } )
                   numer3 <- sapply(sim.data$t, function(u){sum( (testdata4$X3[testdata4$t==u]*exp_gamma(u))[u<=baseData$C], na.rm=T) } )
                   denom <- sapply(sim.data$t, function(u){sum( (exp_gamma(u))[u<=baseData$C], na.rm=T) } )
                   Vbar <- cbind(numer1/denom, numer2/denom, numer3/denom)

                  bigV <-cbind(sim.data$X1, sim.data$X2, sim.data$X3)
                  temp <- colSums((bigV-Vbar)/nn, na.rm=T)
                  temp
              }

              gamma <- c(delta1, delta2, delta3)

              gamma.hat <- nleqslv(gamma, gamma.f)$x


#################estimated dLam(t) under X1+X2+X3 in observation-time model evaluated at solution for gamma 
              exp_gamma <- function(u){exp(gamma.hat[1]*baseData$X1+
	   	         gamma.hat[2]*baseData$X2+gamma.hat[3]*testdata4$X3[testdata4$t==u])}

#sum over subjects of exp_gamma for each time, result length of unique visit times udt
              denom_gamma <- sapply(udt, function(u){sum( (exp_gamma(u))[u<=baseData$C], na.rm=T) } )

#estimated dLam(t) under X1+X2+X3
              estlam.t.gamma <- sapply(1:length(udt), function(u) sum( ((sim.data$t==udt[u])/denom_gamma[u])) )


####estimated Ybar_star (closest neighbor)
####plugging in 0 for ids with 0 visits, b/c will be multiplied by m=0
#################################see 
              Y_star <- function(t){
                   sapply(baseData$ID, function(n){
                   if(n %in% sim.data$ID){
                     tail(sim.data$Y[sim.data$ID==n][(abs(sim.data$t[sim.data$ID==n]-t)==min(abs(sim.data$t[sim.data$ID==n]-t)))],1)
                     }else{
                      0}       
                   })
              }

###########################################
###########################################
# Estimation methods
###########################################
###########################################

              piCi.delta <- sapply(baseData$ID, function(n, t){sum( (exp(delta.hat[1]*baseData$X1[baseData$ID==n]+
		            delta.hat[2]*baseData$X2[baseData$ID==n]+delta.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n]+
                            delta.hat[4]*testdata4$Z4[testdata4$t==t & testdata4$ID==n])*
		            estlam.t.delta)[t<=baseData$C[baseData$ID==n]], na.rm=T) }, t=udt )

              piCi.gamma <- sapply(baseData$ID, function(n, t){sum( (exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		            gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])*
		            estlam.t.gamma)[t<=baseData$C[baseData$ID==n]], na.rm=T) }, t=udt )


              S1 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m/piCi.gamma))*(u<=baseData$C))})
              S2 <- sapply(udt, function(u){ mean((exp_gamma(u)*((baseData$m^2)/piCi.gamma))*(u<=baseData$C))})

              P1 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m*(baseData$m-1)/piCi.gamma^2))*(u<=baseData$C))})
              P2 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m*(baseData$m-1)^2/piCi.gamma^2))*(u<=baseData$C))})

              S1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m/piCi.gamma))*(u<=baseData$C))})
              S1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m/piCi.gamma))*(u<=baseData$C))})
              S1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m/piCi.gamma))*(u<=baseData$C))})

              P1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m*(baseData$m-1)/piCi.gamma^2))*(u<=baseData$C))}) 
              P1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m*(baseData$m-1)/piCi.gamma^2))*(u<=baseData$C))}) 
              P1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi.gamma^2))*(u<=baseData$C))})

              Ybar_star.numer <- sapply(udt, function(u) {mean((Y_star(u)*exp_gamma(u)*(baseData$m/piCi.gamma))*(u<=baseData$C)) } )
              Ybar <- Ybar_star.numer/S1

              Xbar1 <- S1X1/S1
              Xbar2 <- S1X2/S1
              Xbar3 <- S1X3/S1

####B(A)^{-1}
              B.Ainv.f <- function(u){
                   A <- matrix(c(S1[udt==u], P1[udt==u], S2[udt==u], P2[udt==u]), nrow=2,ncol=2, byrow=TRUE)
                   B <- matrix(c(S1X1[udt==u], P1X1[udt==u],S1X2[udt==u], P1X2[udt==u], S1X3[udt==u], P1X3[udt==u]), ncol=2, byrow=TRUE) 

                   B.Ainv <- B%*%ginv(A)
                   B.Ainv
              }


####numerator of ICIRR beta estimate
              beta.ICIRR.num.f <- function(dd){
                   u <- dd['t'] #time
                   m <- dd['m']
                   id <- dd['ID']
                   X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
                   Yval <- dd['Y']
                   inv.rho <- ((exp_gamma(u)/piCi.gamma)/(exp_delta(u)/piCi.delta))[which(baseData$ID==id)]  
  
                   (X.mat - B.Ainv.f(u)%*%matrix(c(1,m), ncol=1))*(Yval-Ybar[which(udt==u)])*inv.rho
              }

####denominator of ICIRR beta estimate
              beta.ICIRR.den.f <- function(dd){
                   u <- dd['t'] #time
                   m <- dd['m']
                   id <- dd['ID']
                   X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
                   Xbar.mat <- matrix(c(Xbar1[which(udt==u)], Xbar2[which(udt==u)], Xbar3[which(udt==u)]), ncol=1)
                   inv.rho <- ((exp_gamma(u)/piCi.gamma)/(exp_delta(u)/piCi.delta))[which(baseData$ID==id)] 

                   b.den <- (X.mat - B.Ainv.f(u)%*%matrix(c(1,m), ncol=1))%*%t(X.mat-Xbar.mat)*inv.rho
                   list(b.den)
              }

              beta.est.ICIRR.num <-  apply(apply(sim.data,1,beta.ICIRR.num.f),1,sum)

              beta.est.ICIRR.denom.L <- apply(sim.data,1,beta.ICIRR.den.f)
              beta.est.ICIRR.denom.A <- array(unlist(beta.est.ICIRR.denom.L), dim = c(nrow(beta.est.ICIRR.denom.L[[1]][[1]]), ncol(beta.est.ICIRR.denom.L[[1]][[1]]), length(beta.est.ICIRR.denom.L)))
              beta.est.ICIRR.denom <- apply(beta.est.ICIRR.denom.A,c(1,2),sum)

####beta ICIRR estimate
              beta.est.ICIRR <- ginv(beta.est.ICIRR.denom)%*%matrix(beta.est.ICIRR.num,ncol=1)



##############################################################################
#######################################################Buzkova & Lumley 2009 method
##############################################################################

#######for BL-weighted terms
              H0 <- sapply(udt, function(u){ sum(exp_gamma(u)*(u<=baseData$C))})/nn
              Hbar.X1 <- sapply(udt, function(u){ sum(exp_gamma(u)*baseData$X1*(u<=baseData$C))/nn})/H0
              Hbar.X2 <- sapply(udt, function(u){ sum(exp_gamma(u)*baseData$X2*(u<=baseData$C))/nn})/H0
              Hbar.X3 <- sapply(udt, function(u){ sum(exp_gamma(u)*testdata4$X3[testdata4$t==u]*(u<=baseData$C))/nn})/H0

              Ybar_star.numer.BL <- sapply(udt, function(u) {mean((Y_star(u)*exp_gamma(u))*(u<=baseData$C)) } )
              Ybar.BL <- Ybar_star.numer.BL/H0

########for BL beta estimation
              beta.num.BL.f <- function(dd){
                   u <- dd['t'] #time
                   m <- dd['m']
                   id <- dd['ID']
                   X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
                   Xbar.BL.mat <- matrix(c(Hbar.X1[which(udt==u)], Hbar.X2[which(udt==u)], Hbar.X3[which(udt==u)]), ncol=1)
                   Yval <- dd['Y']
                   inv.rho <- (exp_gamma(u)/exp_delta(u))[which(baseData$ID==id)]  
                   (X.mat - Xbar.BL.mat)*(Yval-Ybar.BL[which(udt==u)])*inv.rho
              }

              beta.den.BL.f <- function(dd){
                   u <- dd['t'] #time
                   m <- dd['m']
                   id <- dd['ID']
                   X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
                   Xbar.BL.mat <- matrix(c(Hbar.X1[which(udt==u)], Hbar.X2[which(udt==u)], Hbar.X3[which(udt==u)]), ncol=1)
                   inv.rho <- (exp_gamma(u)/exp_delta(u))[which(baseData$ID==id)]
                   b.den.BL <- (X.mat-Xbar.BL.mat)%*%t(X.mat-Xbar.BL.mat)*inv.rho
                   list(b.den.BL)
              }

              beta.est.num.BL <- apply(apply(sim.data,1,beta.num.BL.f),1,sum)

              beta.est.denom.BL.L <- apply(sim.data,1,beta.den.BL.f)
              beta.est.denom.BL.A <- array(unlist(beta.est.denom.BL.L), dim = c(nrow(beta.est.denom.BL.L[[1]][[1]]), ncol(beta.est.denom.BL.L[[1]][[1]]), length(beta.est.denom.BL.L)))
              beta.est.denom.BL <- apply(beta.est.denom.BL.A,c(1,2),sum)

              beta.est.BL <- ginv(beta.est.denom.BL)%*%matrix(beta.est.num.BL,ncol=1)


###estimates
              ests.out <- list(c(delta.hat), c(beta.est.ICIRR), c(beta.est.BL))
              ests.out
         }  #end function


########################run functions and save estimates

       Orig.Ests <- InformTimes.f(complete.data=complete.data)

       delta.mat[ii.sim,] <- Orig.Ests[[1]]
       write.csv(delta.mat, delta.est.filename,row.names=FALSE)
 
       beta.ICIRR.mat[ii.sim,] <- Orig.Ests[[2]]
       write.csv(beta.ICIRR.mat, beta.est.ICIRR.filename,row.names=FALSE)

       beta.BL.mat[ii.sim,] <- Orig.Ests[[3]]
       write.csv(beta.BL.mat, beta.est.BL.filename,row.names=FALSE)


print(ii.sim)

      } ##end simulated dataset iterate
    }  ##end delta4 iterate
  }  ##end alpha.20(t) iterate
}  ## end nn iterate
