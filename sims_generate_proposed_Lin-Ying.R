####################################################################################
#Simulation study for 
#"Time-Varying Latent Effect Models for Repeated Measurements to Address
#Informative Observation Times in the U.S. Medicare Minimum Data Set"  
#By Shardell M, Chen C, and Falvey J
#University of Maryland School of Medicine

#Three covariates: X1 and X2 are time-invariant, X3 is time-varying
#Comparing proposed method to Lin and Ying (2001) (LY)
####################################################################################

##directory of simulation results
res.dir <- "./sims_results/"

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
alpha.20.t.mat <- cbind((1/(1+time_discrete)), sin(time_discrete), rep(0, length(time_discrete)))
alpha.20.t.text.vec <- c("1.over.(1+t)","sin(t)", "0")

##################possible values for gamma3
gamma3.vec <- c(-0.5, 0, 0.5)

##################iterate over different sample size
for(ii.nn in 1:length(nn.vec)){

#sample size per simulation:
  nn <- nn.vec[ii.nn]

################iterate over alpha.20(t) and gamma3)

  for(ii.alpha2 in 1:ncol(alpha.20.t.mat)){

    alpha.20.t <- alpha.20.t.mat[,ii.alpha2]
    alpha.20.t.text <- alpha.20.t.text.vec[ii.alpha2]

    for(ii.gamma3 in 1:length(gamma3.vec)){

      gamma3 <- gamma3.vec[ii.gamma3]

############storage matrices for beta and se (3 covariates)
      beta.mat <- matrix(c(0), nrow=n.sim, ncol=3)
      beta.se.mat <- matrix(c(0), nrow=n.sim, ncol=3)

      beta.est.filename <- paste0(res.dir, "beta_ests_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")
      beta.se.filename <- paste0(res.dir, "beta_ses_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")

###########storage matrices for gamma and se (3 covariates)
      gamma.mat <- matrix(c(0), nrow=n.sim, ncol=3)
      gamma.se.mat <- matrix(c(0), nrow=n.sim, ncol=3)

      gamma.est.filename <- paste0(res.dir, "gamma_ests_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")
      gamma.se.filename <- paste0(res.dir, "gamma_ses_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")


############storage matrices for LY beta and se (3 covariates)
      beta.LY.mat <- matrix(c(0), nrow=n.sim, ncol=3)
      beta.se.LY.mat <- matrix(c(0), nrow=n.sim, ncol=3)

      beta.est.LY.filename <- paste0(res.dir, "beta_ests_LY_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")
      beta.se.LY.filename <- paste0(res.dir, "beta_ses_LY_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")


################storage matrices for scriptA and se 
      dscriptA1.mat <- matrix(0, nrow=n.sim, ncol=n.times)
      dscriptA2.mat <- matrix(0, nrow=n.sim, ncol=n.times)

      scriptA1.se.mat <- matrix(0, nrow=n.sim, ncol=n.times)
      scriptA2.se.mat <- matrix(0, nrow=n.sim, ncol=n.times)

      dscriptA1.est.filename <- paste0(res.dir, "dscriptA1_ests_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")
      dscriptA2.est.filename <- paste0(res.dir, "dscriptA2_ests_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")

      scriptA1.ses.filename <- paste0(res.dir, "scriptA1_ses_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")
      scriptA2.ses.filename <- paste0(res.dir, "scriptA2_ses_n_",nn,"_a2_",alpha.20.t.text,"_gamma3_",gamma3, ".csv")


################################simulation parameters
########X1: P(X1=1), X1 Bernoulli(p.X1)
      p.X1 <- 0.5

########X2: f(X2) ~N(mu.X2, sd.X2) 
      mu.X2 <- 0
      sd.X2 <- 0.5

########X3(t): f(X3(t) | X1, X2) ~ N(mu.X3.t, sd.X3.t)
#mu.X3.t = E[X3(t) | X1, X2, t] = nu0.X3 + nu1.X3*X1 + nu2.X3*X2
      nu0.X3 <- 0
      nu1.X3 <- 0.25
      nu2.X3 <- 0.5
      sd.X3.t <- 0.5

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
#mu.Y.t = E[Y(t) | X1, X2, X3(t), V] = alpha.10.t + alpha.20.t*V + beta1.Y*X1 + beta2.Y*X2 + beta3.Y*X3(t)
      sd.Y <- sqrt(0.25)
      beta1.Y <- 1
      beta2.Y <- -1
      beta3.Y <- 1


###############dN(t) ~ Pois(rate.t)
#rate.t = V*lambda0.t*exp(gamma1*X1 + gamma2*X2 + gamma3*X3.t) 
      lambda0.t <- 0.1 
      gamma1 <- -0.5
      gamma2 <- -0.25


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

######simulate latent variable V
## V = (1-X1)*U1 + X1*U2
#U1 ~ gamma(mean=1, var=0.5)
         U1 <- rgamma(n=nn,shape=shape.U1,scale= scale.U1)

#U2 ~ gamma(mean=1, var=0.8)
         U2 <- rgamma(n=nn,shape=shape.U2,scale= scale.U2)

         V <- (1-X1)*U1 + X1*U2


###############simulate outcome Y.t
#E[Y(t) | X1, X2, X3(t), V] = alpha.10.t + alpha.20.t*V + beta1.Y*X1 + beta2.Y*X2 + beta3.Y*X3(t)

         Y.sim.f <-function(id){
              mu.Y.t <- alpha.10.t + alpha.20.t*V[id] + beta1.Y*X1[id] + beta2.Y*X2[id] + beta3.Y*X3.t[id,]
              Y <- rnorm(length(mu.Y.t), mean=mu.Y.t, sd=sd.Y)
              Y
         }

#rows are id, columns are t
         Y.t <-t(sapply(c(1:nn), Y.sim.f))


###############simulate visit times dN(t)
#rate.t = V*lambda0.t*exp(gamma1*X1 + gamma2*X2 + gamma3*X3.t) 

         dNi.f <- function(id){
              rate.t <- V[id]*lambda0.t*exp(gamma1*X1[id] + gamma2*X2[id] + gamma3*X3.t[id,])
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
                          X1=rep(X1, n.times), X2=rep(X2, n.times), X3=c(X3.t), C=rep(Ci, n.times), dNi=c(dNi.t))

         complete.data <- complete.data[order(complete.data$ID),]

###################################plugging in NA for Y when dNi=0 and dropping rows when t>Ci
         obs.data <- complete.data[complete.data$t<=complete.data$C,]
         obs.data$Y[obs.data$dNi==0] <- NA


#############################################create observed subset with only observed rows
########################################################################
# Analysis dataset: sim.data
# ID = participant id, t = visit time
# X1 = binary baseline, X2 = continuous baseline, X3 = continuous time-dependent
# Y = continuous
# C = censoring time
# m = number of visits (id-level)
########################################################################
         sim.data <- obs.data[obs.data$dNi==1,]
         sim.data <- sim.data[, c("ID", "t", "m", "Y", "X1", "X2", "X3", "C")]
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
# gamma.hat (observation-time model covariates)
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

         gamma <- c(gamma1, gamma2, gamma3)

         gamma.hat <- nleqslv(gamma, gamma.f)$x


#####save gamma.hat
         gamma.mat[ii.sim,] <- gamma.hat
         write.csv(gamma.mat, gamma.est.filename,row.names=FALSE)
 

#################estimated dLam(t) under X1+X2+X3 in observation-time model evaluated at solution for gamma 
         exp_gamma <- function(u){exp(gamma.hat[1]*baseData$X1+
	   	         gamma.hat[2]*baseData$X2+gamma.hat[3]*testdata4$X3[testdata4$t==u])}

#sum over subjects of exp_gamma for each time, result length of unique visit times udt
         denom_gamma <- sapply(udt, function(u){sum( (exp_gamma(u))[u<=baseData$C], na.rm=T) } )

#estimated dLam(t) under X1+X2+X3
         estlam.t.gamma <- sapply(1:length(udt), function(u) sum( ((sim.data$t==udt[u])/denom_gamma[u])) )


####estimated Ybar_star (closest neighbor)
####plugging in 0 for ids with 0 visits, b/c will be multiplied by m=0
############################################ 
         Y_star <- function(t){
              sapply(baseData$ID, function(n){
              if(n %in% sim.data$ID){
                tail(sim.data$Y[sim.data$ID==n][(abs(sim.data$t[sim.data$ID==n]-t)==min(abs(sim.data$t[sim.data$ID==n]-t)))],1)
                }else{
                 0}       
              })
         }

########################################################
########################################################
# Proposed Estimation Method
########################################################
########################################################

         piCi <- sapply(baseData$ID, function(n, t){sum( (exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
	  	         gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])*
		         estlam.t.gamma)[t<=baseData$C[baseData$ID==n]], na.rm=T) }, t=udt )

         S1 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m/piCi))*(u<=baseData$C))})
         S2 <- sapply(udt, function(u){ mean((exp_gamma(u)*((baseData$m^2)/piCi))*(u<=baseData$C))})

         P1 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})
         P2 <- sapply(udt, function(u){ mean((exp_gamma(u)*(baseData$m*(baseData$m-1)^2/piCi^2))*(u<=baseData$C))})

         S1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m/piCi))*(u<=baseData$C))})
         S1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m/piCi))*(u<=baseData$C))})
         S1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m/piCi))*(u<=baseData$C))})

         P1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))}) 
         P1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))}) 
         P1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})


         Ybar_star.numer <- sapply(udt, function(u) {mean((Y_star(u)*exp_gamma(u)*(baseData$m/piCi))*(u<=baseData$C)) } )
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

####numerator of beta estimate
         beta.num.f <- function(dd){
              u <- dd['t'] #time
              m <- dd['m']
              X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
              Yval <- dd['Y']

              (X.mat - B.Ainv.f(u)%*%matrix(c(1,m), ncol=1))*(Yval-Ybar[which(udt==u)])
         }

####denominator of beta estimate
         beta.den.f <- function(dd){
              u <- dd['t'] #time
              m <- dd['m']
              X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
              Xbar.mat <- matrix(c(Xbar1[which(udt==u)], Xbar2[which(udt==u)], Xbar3[which(udt==u)]), ncol=1)

              b.den <- (X.mat - B.Ainv.f(u)%*%matrix(c(1,m), ncol=1))%*%t(X.mat-Xbar.mat)
              list(b.den)
         }

         beta.est.num <-  apply(apply(sim.data,1,beta.num.f),1,sum)

         beta.est.denom.L <- apply(sim.data,1,beta.den.f)
         beta.est.denom.A <- array(unlist(beta.est.denom.L), dim = c(nrow(beta.est.denom.L[[1]][[1]]), ncol(beta.est.denom.L[[1]][[1]]), length(beta.est.denom.L)))
         beta.est.denom <- apply(beta.est.denom.A,c(1,2),sum)

####beta estimate
         beta.est <- ginv(beta.est.denom)%*%matrix(beta.est.num,ncol=1)
         beta.mat[ii.sim,] <- beta.est
         write.csv(beta.mat, beta.est.filename,row.names=FALSE)

################################estimating script A
#######A^{-1}
         Ainv.f <- function(u){
              A <- matrix(c(S1[udt==u], P1[udt==u], S2[udt==u], P2[udt==u]), nrow=2,ncol=2, byrow=TRUE) 
              Ainv <- ginv(A)
              Ainv
         }

         d.scriptA.f <- function(u){
              sim.sub <- sim.data[which(sim.data$t==u),]
              X.mat <- t(cbind(sim.sub$X1, sim.sub$X2, sim.sub$X3))

              val1 <- sum(sim.sub$Y - as.vector(t(beta.est)%*%X.mat))/nn
              val2 <- sum(sim.sub$m*(sim.sub$Y - as.vector(t(beta.est)%*%X.mat)))/nn

              d.scriptA <- Ainv.f(u)%*%matrix(c(val1, val2), ncol=1)
              d.scriptA
         }


#######one value for each unique time
         d.scriptA <- sapply(udt, d.scriptA.f)
         scriptA <- apply(d.scriptA,1,cumsum)

         dscriptA1.mat[ii.sim,which(time_discrete %in% udt)] <- d.scriptA[1,]
         dscriptA2.mat[ii.sim,which(time_discrete %in% udt)] <- d.scriptA[2,]

########write d.scriptA         
         write.csv(dscriptA1.mat, dscriptA1.est.filename,row.names=FALSE)
         write.csv(dscriptA2.mat, dscriptA2.est.filename,row.names=FALSE) 

################################################estimating the variance of estimated beta
################################################ ginv(D) %*% Sigma %*% ginv(t(D))

########################estimate D

         D <- beta.est.denom/nn

################final asymptotic form of Psi 

         S1.X1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X1*(baseData$m/piCi))*(u<=baseData$C))})
         S1.X1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X2*(baseData$m/piCi))*(u<=baseData$C))})
         S1.X1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*testdata4$X3[testdata4$t==u]*(baseData$m/piCi))*(u<=baseData$C))})
         S1.X2X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*baseData$X2*(baseData$m/piCi))*(u<=baseData$C))})
         S1.X2X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*testdata4$X3[testdata4$t==u]*(baseData$m/piCi))*(u<=baseData$C))})
         S1.X3X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*testdata4$X3[testdata4$t==u]*(baseData$m/piCi))*(u<=baseData$C))})

         P1.X1X1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X1*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))}) 
         P1.X1X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X2*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})
         P1.X1X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})
         P1.X2X2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*baseData$X2*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})
         P1.X2X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})
         P1.X3X3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi^2))*(u<=baseData$C))})

         S1.piX1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m/piCi^2))*(u<=baseData$C))})
         S1.piX2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m/piCi^2))*(u<=baseData$C))})
         S1.piX3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m/piCi^2))*(u<=baseData$C))})

         P1.piX1 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*(baseData$m*(baseData$m-1)/piCi^3))*(u<=baseData$C))}) 
         P1.piX2 <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*(baseData$m*(baseData$m-1)/piCi^3))*(u<=baseData$C))}) 
         P1.piX3 <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)/piCi^3))*(u<=baseData$C))})


##################################estimating d.F

         d.F.f <- function(u){

              val1 <- sum((sim.data$t==u))/nrow(baseData)
              val2 <- sum(sim.data$m*(sim.data$t==u))/nrow(baseData)

              d.F <- Ainv.f(u)%*%matrix(c(val1, val2), ncol=1)
              d.F
         }

####one value for each unique time
         d.F <- sapply(udt, d.F.f)
         F <- apply(d.F,1,cumsum)


##################################################psi.i (asymptotic representation of gamma)

         H0 <- sapply(udt, function(u){ sum(exp_gamma(u)*(u<=baseData$C))})/nn
         Hbar.X1 <- sapply(udt, function(u){ sum(exp_gamma(u)*baseData$X1*(u<=baseData$C))/nn})/H0
         Hbar.X2 <- sapply(udt, function(u){ sum(exp_gamma(u)*baseData$X2*(u<=baseData$C))/nn})/H0
         Hbar.X3 <- sapply(udt, function(u){ sum(exp_gamma(u)*testdata4$X3[testdata4$t==u]*(u<=baseData$C))/nn})/H0

##############function for the numerator of estimating equation for gamma
         psi.it.num.f <- function(n,t){

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$ID==n & testdata4$t==t])
              H <- cbind(Hbar.X1[which(udt==t)], Hbar.X2[which(udt==t)], Hbar.X3[which(udt==t)])

              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              dM.it <- dNi - (exp(gamma.hat[1]*baseData$X1[baseData$ID==n] +
	    	    gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$ID==n & testdata4$t==t])*
		    estlam.t.gamma[which(udt==t)])*(t<=baseData$C[baseData$ID==n])

              apply((X.mat - H)*dM.it,2,sum)
         }


         psi.i.num <- sapply(baseData$ID, psi.it.num.f, t=udt)


###############denominator of estimating equation for gamma - this flattens each id's outer(X-H, X-H) matrix
         G.it.f <- function(n,t){

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$t==t & testdata4$ID==n])
              H <- cbind(Hbar.X1[which(udt==t)], Hbar.X2[which(udt==t)], Hbar.X3[which(udt==t)])

              the.outer2 <- t(apply((X.mat-H),1, function(x){(outer(x,x))}))

              the.rate <- (exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
	  	        gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])*
		        estlam.t.gamma[which(udt==t)])*(t<=baseData$C[baseData$ID==n])
              G.it <- the.outer2*the.rate

              G.i <- apply(G.it,2,sum)
              G.i
         }

#####average across ids and "unflatten" id-specific matrix

         G.i <- sapply(baseData$ID, G.it.f,t=udt)
         G <- matrix(apply(G.i,1,mean),nrow=3,ncol=3)


#######################################matrix, length(gamma.hat) by number of participants
         psi.i <- solve(G) %*% psi.i.num

         gamma.var <- var(t(psi.i))/nn
         gamma.se <- sqrt(diag(gamma.var))

         gamma.se.mat[ii.sim,] <- gamma.se

##write standard errors of gamma
         write.csv(gamma.se.mat, gamma.se.filename,row.names=FALSE)

################################phi.it (asymptotic representation of lambda(t)

###################################first need dphi.it 

         dphi.term1.it.f <- function(n,t){
 
              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              dM.it <- dNi -  (exp(gamma.hat[1]*baseData$X1[baseData$ID==n] +
	  	     gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$ID==n & testdata4$t==t])*
		     estlam.t.gamma[which(udt==t)])*(t<=baseData$C[baseData$ID==n])

              dM.it/H0
         }

#dim row=length(udt)*col=length(id) then flattened to make 
#length = length(udt) * length(id) sorted by id, then t
         dphi.term1.it <- c(sapply(baseData$ID, dphi.term1.it.f, t=udt))


#dim length(udt) * p
         H <- cbind(Hbar.X1, Hbar.X2, Hbar.X3)

##dim length(udt) by length(id)
         dphi.term2.it <- apply(H,2, function(x){x*estlam.t.gamma})%*%psi.i

#organized by id, then t within id
         dphi.term2.it <- c(dphi.term2.it)

#################dphi for each id and each time
###length = length(udt) * length(id)
         dphi.it <- dphi.term1.it - dphi.term2.it


#########################################theta.j(xi, Ci), asymptototic representation of pi(Ci, xi)
####################################first need dtheta.j(xi, t) 
         dtheta.it.f <- function(u){

              X.mat <- cbind(baseData$X1,baseData$X2,testdata4$X3[testdata4$t==u])

              exp.gammahat.x <- exp(gamma.hat[1]*X.mat[,1] + gamma.hat[2]*X.mat[,2] + gamma.hat[3]*X.mat[,3])

#Xi, Ci are rows, psi.j are columns:      
              dtheta.ij.u <- outer(exp.gammahat.x, dphi.it[testdata4$t==u])    +  
                            estlam.t.gamma[which(udt==u)]*(diag(exp.gammahat.x)%*%X.mat)%*%psi.i  

              dtheta.ij.u
         }

#########################dtheta.ij.t; 3-dim array, length of ids x length of ids at each timepoint. 
         dtheta.ijt.L <- lapply(udt, dtheta.it.f)
         dtheta.ijt.A <- array(unlist(dtheta.ijt.L), dim = c(nrow(dtheta.ijt.L[[1]]), ncol(dtheta.ijt.L[[1]]), length(dtheta.ijt.L)))


###############################################sum of phi until censoring time for each id
         phi.i.C.f <- function(dd){
              n <- dd['ID']
              C <- dd['C']

              sum(dphi.it[which(testdata4$ID==n & testdata4$t<=C)])
         }

###############phi.i.C
##################length number of ids
         phi.i.C <- apply(baseData,1,phi.i.C.f)

##############################################theta.ij.C (dtheta.ij.t summed up to C)

         theta.ij.C.f <- function(dd){
              n <- dd['ID']
              C <- dd['C']

              dtheta.ijt.A.sub <- dtheta.ijt.A[which(baseData$ID==n), , which(udt<=C)] 
              dtheta.ij.A.sub <- rowSums(dtheta.ijt.A.sub)
              dtheta.ij.A.sub
         }


########################################theta.ij.C: length of ids (xi,ci) x length of ids (phi.j, psi.j)
         theta.ij.C <- t(apply(baseData,1,theta.ij.C.f))

##################################################d.Mstar.it 

         d.Mstar.it.f <- function(n,t){
  
              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$ID==n & testdata4$t==t])
              pi.C <- piCi[baseData$ID==n]
              exp.gamma <- exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		              gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$ID==n & testdata4$t==t])
              d.scriptA.1 <- d.scriptA[1,]
              d.scriptA.2 <- d.scriptA[2,]

              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              Yval <- testdata4$Y[testdata4$ID==n]
              Yval[is.na(Yval)] <- 0  ##will multiply by (dNi=0)              

              d.Mstar.it <- (Yval - c(X.mat%*%beta.est))*dNi - (t<=C)*m*exp.gamma*d.scriptA.1/pi.C - (t<=C)*m*(m-1)*exp.gamma*d.scriptA.2/pi.C^2
              d.Mstar.it
         }

###dim row=length(udt)*col=length(id)
###then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Mstar.it <- c(sapply(baseData$ID, d.Mstar.it.f, t=udt))

####################################d.Mtildestar.it 
         d.Mtildestar.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              pi.C <- piCi[baseData$ID==n]
              exp.gamma <- exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		            gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])
              d.F.1 <- d.F[1,]
              d.F.2 <- d.F[2,]

              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              d.Mtildestar.it <- dNi - (t<=C)*m*exp.gamma*d.F.1/pi.C - (t<=C)*m*(m-1)*exp.gamma*d.F.2/pi.C^2
              d.Mtildestar.it
         }


###dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Mtildestar.it <- c(sapply(baseData$ID, d.Mtildestar.it.f, t=udt))


##################################setting up for d.Qtilde1.it

         S1.pi <- sapply(udt, function(u){ mean(exp_gamma(u)*(baseData$m/piCi^2)*(u<=baseData$C))})
         S2.pi <- sapply(udt, function(u){ mean(exp_gamma(u)*(baseData$m^2/piCi^2)*(u<=baseData$C))})

         P1.pi <- sapply(udt, function(u){ mean(exp_gamma(u)*(baseData$m*(baseData$m-1)/piCi^3)*(u<=baseData$C))})
         P2.pi <- sapply(udt, function(u){ mean(exp_gamma(u)*(baseData$m*(baseData$m-1)^2/piCi^3)*(u<=baseData$C))})


######################### m.i * d.Mstar.it (term that has same expected value)
 
         mi.d.Mstar.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$t==t & testdata4$ID==n])
              pi.C <- piCi[baseData$ID==n]
              exp.gamma <- exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		             gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])
              d.scriptA.1 <- d.scriptA[1,]
              d.scriptA.2 <- d.scriptA[2,]

              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              Yval <- testdata4$Y[testdata4$ID==n]
              Yval[is.na(Yval)] <- 0  ##will multiply by (dNi=0)      

              mi.d.Mstar.it <- m*(Yval - c(X.mat%*%beta.est))*dNi - (t<=C)*m*m*exp.gamma*d.scriptA.1/pi.C - (t<=C)*m*(m-1)^2*exp.gamma*d.scriptA.2/pi.C^2
              mi.d.Mstar.it
         }

######################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         mi.d.Mstar.it <- c(sapply(baseData$ID,mi.d.Mstar.it.f, t=udt))


         S2X1 <- sapply(udt, function(u){ mean(exp_gamma(u)*baseData$X1*(baseData$m^2/piCi)*(u<=baseData$C))})
         S2X2 <- sapply(udt, function(u){ mean(exp_gamma(u)*baseData$X2*(baseData$m^2/piCi)*(u<=baseData$C))})
         S2X3 <- sapply(udt, function(u){ mean(exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m^2/piCi)*(u<=baseData$C))})

         P2X1 <- sapply(udt, function(u){ mean(exp_gamma(u)*baseData$X1*(baseData$m*(baseData$m-1)^2/piCi^2)*(u<=baseData$C))}) 
         P2X2 <- sapply(udt, function(u){ mean(exp_gamma(u)*baseData$X2*(baseData$m*(baseData$m-1)^2/piCi^2)*(u<=baseData$C))}) 
         P2X3 <- sapply(udt, function(u){ mean(exp_gamma(u)*testdata4$X3[testdata4$t==u]*(baseData$m*(baseData$m-1)^2/piCi^2)*(u<=baseData$C))})


#############length udt x id (sorted by id, then t within id) (j column of theta.ij.C, contribution to asymptotic):
#for each column j time t, taking average over i (rows of theta.i.C)
         S1.pi.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$m/piCi^2)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
          }

         S1.pi.theta.j <- apply(testdata4,1, S1.pi.theta.j.f)

#for each column j time t, taking average over i (rows of theta.i.C)
         S2.pi.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$m^2/piCi^2)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         S2.pi.theta.j <- apply(testdata4,1, S2.pi.theta.j.f)

#for each column j time t, take over i (rows of theta.i.C)
         P1.pi.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$m*(baseData$m-1)/piCi^3)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         P1.pi.theta.j <- apply(testdata4,1, P1.pi.theta.j.f)

#for each column j time t, taking average over i (rows of theta.i.C)
         P2.pi.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$m*(baseData$m-1)^2/piCi^3)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         P2.pi.theta.j <- apply(testdata4,1, P2.pi.theta.j.f)


###########################################################d.Qtilde1.it

         d.Qtilde1.it.f <- function(n,t){
  
              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              d.Mstar <- d.Mstar.it[testdata4$ID==n]
              m.d.Mstar <- mi.d.Mstar.it[testdata4$ID==n]

              S1.Xmat <- cbind(S1X1[udt==t], S1X2[udt==t],S1X3[udt==t])
              S2.Xmat <- cbind(S2X1[udt==t], S2X2[udt==t],S2X3[udt==t])

              P1.Xmat <- cbind(P1X1[udt==t], P1X2[udt==t],P1X3[udt==t])
              P2.Xmat <- cbind(P2X1[udt==t], P2X2[udt==t],P2X3[udt==t])

              d.scriptA.1 <- d.scriptA[1,(udt==t)]
              d.scriptA.2 <- d.scriptA[2,(udt==t)]

              S1.pi.theta <- S1.pi.theta.j[testdata4$ID==n]
              S2.pi.theta <- S2.pi.theta.j[testdata4$ID==n]
              P1.pi.theta <- P1.pi.theta.j[testdata4$ID==n]
              P2.pi.theta <- P2.pi.theta.j[testdata4$ID==n]

              d.Q1.it <- P2[udt==t]*d.Mstar - P1[udt==t]*m.d.Mstar -
                         c((P2[udt==t]*S1.Xmat - P1[udt==t]*S2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.scriptA.1 +
                         (P2[udt==t]*S1.pi.theta - P1[udt==t]*S2.pi.theta)*d.scriptA.1 - 
                         c((P2[udt==t]*P1.Xmat - P1[udt==t]*P2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.scriptA.2 +
                         2*(P2[udt==t]*P1.pi.theta - P1[udt==t]*P2.pi.theta)*d.scriptA.2

              d.Qtilde1.it <- d.Q1.it/(P2[udt==t]*S1[udt==t] - P1[udt==t]*S2[udt==t])
              d.Qtilde1.it
         }


###dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Qtilde1.it.mat <- sapply(baseData$ID, d.Qtilde1.it.f, t=udt) #use this for scriptA ses
         d.Qtilde1.it <- c(d.Qtilde1.it.mat)


####################################################################d.Qtilde2.it

         d.Qtilde2.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              d.Mstar <- d.Mstar.it[testdata4$ID==n]
              m.d.Mstar <- mi.d.Mstar.it[testdata4$ID==n]

              S1.Xmat <- cbind(S1X1[udt==t], S1X2[udt==t],S1X3[udt==t])
              S2.Xmat <- cbind(S2X1[udt==t], S2X2[udt==t],S2X3[udt==t])

              P1.Xmat <- cbind(P1X1[udt==t], P1X2[udt==t],P1X3[udt==t])
              P2.Xmat <- cbind(P2X1[udt==t], P2X2[udt==t],P2X3[udt==t])

              d.scriptA.1 <- d.scriptA[1,(udt==t)]
              d.scriptA.2 <- d.scriptA[2,(udt==t)]

              S1.pi.theta <- S1.pi.theta.j[testdata4$ID==n]
              S2.pi.theta <- S2.pi.theta.j[testdata4$ID==n]
              P1.pi.theta <- P1.pi.theta.j[testdata4$ID==n]
              P2.pi.theta <- P2.pi.theta.j[testdata4$ID==n]

              d.Q2.it <- S2[udt==t]*d.Mstar - S1[udt==t]*m.d.Mstar -
                      c((S2[udt==t]*S1.Xmat - S1[udt==t]*S2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.scriptA.1 +
                      (S2[udt==t]*S1.pi.theta - S1[udt==t]*S2.pi.theta)*d.scriptA.1 - 
                      c((S2[udt==t]*P1.Xmat - S1[udt==t]*P2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.scriptA.2 +
                      2*(S2[udt==t]*P1.pi.theta - S1[udt==t]*P2.pi.theta)*d.scriptA.2


              d.Qtilde2.it <- d.Q2.it/(S2[udt==t]*P1[udt==t] - S1[udt==t]*P2[udt==t])
              d.Qtilde2.it
         }


######################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Qtilde2.it.mat <- sapply(baseData$ID, d.Qtilde2.it.f, t=udt)  #for scriptA ses
         d.Qtilde2.it <- c(d.Qtilde2.it.mat)


############################################m.i * d.Mtildestar.it (term that has same expected value)

         mi.d.Mtildestar.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              pi.C <- piCi[baseData$ID==n]
              exp.gamma <- exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		        gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$t==t & testdata4$ID==n])

              d.F.1 <- d.F[1,(udt==t)]
              d.F.2 <- d.F[2,(udt==t)]

              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              d.Mtildestar.it <- m*dNi - (t<=C)*m*m*exp.gamma*d.F.1/pi.C - (t<=C)*m*(m-1)^2*exp.gamma*d.F.2/pi.C^2
              d.Mtildestar.it
         }

######################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         mi.d.Mtildestar.it <- c(sapply(baseData$ID,mi.d.Mtildestar.it.f, t=udt))


######################################################d.Rtilde1.it
         d.Rtilde1.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              d.Mtildestar <- d.Mtildestar.it[testdata4$ID==n]
              m.d.Mtildestar <- mi.d.Mtildestar.it[testdata4$ID==n]

              S1.Xmat <- cbind(S1X1[udt==t], S1X2[udt==t],S1X3[udt==t])
              S2.Xmat <- cbind(S2X1[udt==t], S2X2[udt==t],S2X3[udt==t])

              P1.Xmat <- cbind(P1X1[udt==t], P1X2[udt==t],P1X3[udt==t])
              P2.Xmat <- cbind(P2X1[udt==t], P2X2[udt==t],P2X3[udt==t])

              d.F.1 <- d.F[1,(udt==t)]
              d.F.2 <- d.F[2,(udt==t)]

              S1.pi.theta <- S1.pi.theta.j[testdata4$ID==n]
              S2.pi.theta <- S2.pi.theta.j[testdata4$ID==n]
              P1.pi.theta <- P1.pi.theta.j[testdata4$ID==n]
              P2.pi.theta <- P2.pi.theta.j[testdata4$ID==n]

              d.R1.it <- P2[udt==t]*d.Mtildestar - P1[udt==t]*m.d.Mtildestar -
                      c((P2[udt==t]*S1.Xmat - P1[udt==t]*S2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.F.1 +
                      (P2[udt==t]*S1.pi.theta - P1[udt==t]*S2.pi.theta)*d.F.1 - 
                      c((P2[udt==t]*P1.Xmat - P1[udt==t]*P2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.F.2 +
                      2*(P2[udt==t]*P1.pi.theta - P1[udt==t]*P2.pi.theta)*d.F.2

              d.Rtilde1.it <- d.R1.it/(P2[udt==t]*S1[udt==t] - P1[udt==t]*S2[udt==t])
              d.Rtilde1.it
         }

######################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Rtilde1.it <- c(sapply(baseData$ID, d.Rtilde1.it.f, t=udt))


###################################################d.Rtilde2.it

         d.Rtilde2.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              d.Mtildestar <- d.Mtildestar.it[testdata4$ID==n]
              m.d.Mtildestar <- mi.d.Mtildestar.it[testdata4$ID==n]

              S1.Xmat <- cbind(S1X1[udt==t], S1X2[udt==t],S1X3[udt==t])
              S2.Xmat <- cbind(S2X1[udt==t], S2X2[udt==t],S2X3[udt==t])

              P1.Xmat <- cbind(P1X1[udt==t], P1X2[udt==t],P1X3[udt==t])
              P2.Xmat <- cbind(P2X1[udt==t], P2X2[udt==t],P2X3[udt==t])

              d.F.1 <- d.F[1,(udt==t)]
              d.F.2 <- d.F[2,(udt==t)]

              S1.pi.theta <- S1.pi.theta.j[testdata4$ID==n]
              S2.pi.theta <- S2.pi.theta.j[testdata4$ID==n]
              P1.pi.theta <- P1.pi.theta.j[testdata4$ID==n]
              P2.pi.theta <- P2.pi.theta.j[testdata4$ID==n]

              d.R2.it <- S2[udt==t]*d.Mtildestar - S1[udt==t]*m.d.Mtildestar -
                      c((S2[udt==t]*S1.Xmat - S1[udt==t]*S2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.F.1 +
                      (S2[udt==t]*S1.pi.theta - S1[udt==t]*S2.pi.theta)*d.F.1 - 
                      c((S2[udt==t]*P1.Xmat - S1[udt==t]*P2.Xmat)%*%matrix(psi.i[,baseData$ID==n],ncol=1))*d.F.2 +
                      2*(S2[udt==t]*P1.pi.theta - S1[udt==t]*P2.pi.theta)*d.F.2

              d.Rtilde2.it <- d.R2.it/(S2[udt==t]*P1[udt==t] - S1[udt==t]*P2[udt==t])
              d.Rtilde2.it
         }

######################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.Rtilde2.it <- c(sapply(baseData$ID, d.Rtilde2.it.f, t=udt))


############################S1.piX1.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         S1.piX1.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$X1*baseData$m/piCi^2)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         S1.piX1.theta.j <- apply(testdata4,1, S1.piX1.theta.j.f)

############################S1.piX2.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         S1.piX2.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$X2*baseData$m/piCi^2)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         S1.piX2.theta.j <- apply(testdata4,1, S1.piX2.theta.j.f)

############################S1.piX3.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         S1.piX3.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(testdata4$X3[testdata4$t==u]*baseData$m/piCi^2)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         S1.piX3.theta.j <- apply(testdata4,1, S1.piX3.theta.j.f)

############################P1.piX1.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         P1.piX1.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$X1*baseData$m*(baseData$m-1)/piCi^3)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         P1.piX1.theta.j <- apply(testdata4,1, P1.piX1.theta.j.f)

############################P1.piX2.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         P1.piX2.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(baseData$X2*baseData$m*(baseData$m-1)/piCi^3)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         P1.piX2.theta.j <- apply(testdata4,1, P1.piX2.theta.j.f)

############################P1.piX3.theta.j: for each column j time t, taking average over i (rows of theta.i.C)
         P1.piX3.theta.j.f <- function (dd){
              n <- dd['ID']
              u <- dd['t']

              mean(exp_gamma(u)*(testdata4$X3[testdata4$t==u]*baseData$m*(baseData$m-1)/piCi^3)*(u<=baseData$C)*theta.ij.C[,which(baseData$ID==n)])
         }

         P1.piX3.theta.j <- apply(testdata4,1, P1.piX3.theta.j.f)


##################setting up components for bigPsi
         S1.XXmat.f <- function(t){
              list(matrix(c(S1.X1X1[udt==t],S1.X1X2[udt==t], S1.X1X3[udt==t],
              S1.X1X2[udt==t], S1.X2X2[udt==t], S1.X2X3[udt==t],
              S1.X1X3[udt==t], S1.X2X3[udt==t], S1.X3X3[udt==t]),ncol=3,nrow=3))
         }

         S1.XXmat.L <- sapply(udt,S1.XXmat.f) 
         S1.XXmat <- array(unlist(S1.XXmat.L), dim = c(nrow(S1.XXmat.L[[1]]), ncol(S1.XXmat.L[[1]]), length(S1.XXmat.L)))


         P1.XXmat.f <- function(t){
              list(matrix(c(P1.X1X1[udt==t],P1.X1X2[udt==t], P1.X1X3[udt==t],
              P1.X1X2[udt==t], P1.X2X2[udt==t], P1.X2X3[udt==t],
              P1.X1X3[udt==t], P1.X2X3[udt==t], P1.X3X3[udt==t]),ncol=3,nrow=3))
         }

         P1.XXmat.L <- sapply(udt,P1.XXmat.f) 
         P1.XXmat <- array(unlist(P1.XXmat.L), dim = c(nrow(P1.XXmat.L[[1]]), ncol(P1.XXmat.L[[1]]), length(P1.XXmat.L)))


#####################################bigPsi for variance estimation
         bigPsi.it.f <- function(n,t){

              m <- baseData$m[baseData$ID==n]
              C <- baseData$C[baseData$ID==n]

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$t==t & testdata4$ID==n])

              d.Mstar <- d.Mstar.it[testdata4$ID==n]
              m.d.Mstar <- mi.d.Mstar.it[testdata4$ID==n]

              d.Mtildestar <- d.Mtildestar.it[testdata4$ID==n]
              m.d.Mtildestar <- mi.d.Mtildestar.it[testdata4$ID==n]

              d.Qtilde1 <- d.Qtilde1.it[testdata4$ID==n]
              d.Qtilde2 <- d.Qtilde2.it[testdata4$ID==n]

              d.Rtilde1 <- d.Rtilde1.it[testdata4$ID==n]
              d.Rtilde2 <- d.Rtilde2.it[testdata4$ID==n]

              S1.Xmat <- cbind(S1X1[udt==t], S1X2[udt==t],S1X3[udt==t])
              S2.Xmat <- cbind(S2X1[udt==t], S2X2[udt==t],S2X3[udt==t])

              S1.piXmat <- cbind(S1.piX1[udt==t], S1.piX2[udt==t],S1.piX3[udt==t])
              P1.piXmat <- c(P1.piX1[udt==t], P1.piX2[udt==t],P1.piX3[udt==t])

              S1.piX.theta.mat <- cbind(S1.piX1.theta.j[testdata4$ID==n & testdata4$t==t], 
                                        S1.piX2.theta.j[testdata4$ID==n & testdata4$t==t],
                                        S1.piX3.theta.j[testdata4$ID==n & testdata4$t==t])
              P1.piX.theta.mat <- cbind(P1.piX1.theta.j[testdata4$ID==n & testdata4$t==t], 
                                        P1.piX2.theta.j[testdata4$ID==n & testdata4$t==t],
                                        P1.piX3.theta.j[testdata4$ID==n & testdata4$t==t])

              P1.Xmat <- cbind(P1X1[udt==t], P1X2[udt==t],P1X3[udt==t])
              P2.Xmat <- cbind(P2X1[udt==t], P2X2[udt==t],P2X3[udt==t])

              d.scriptA.1 <- d.scriptA[1,(udt==t)]
              d.scriptA.2 <- d.scriptA[2,(udt==t)]

              d.F.1 <- d.F[1,(udt==t)]
              d.F.2 <- d.F[2,(udt==t)]

              Xbar.mat <- cbind(Xbar1[ which(udt==t)], Xbar2[ which(udt==t)], Xbar3[ which(udt==t)])
              Xbar.beta <- c(Xbar.mat%*%beta.est)
              ybar <- Ybar[which(udt==t)]

              psi <- psi.i[,baseData$ID==n]

              S1.XXmat.psi <- t(apply(S1.XXmat, 3, function(x){matrix(psi, nrow=1)%*%x}))  
              P1.XXmat.psi <- t(apply(P1.XXmat, 3, function(x){matrix(psi, nrow=1)%*%x}))

####term 1 X*d.Mstar
              term1 <- X.mat*d.Mstar

####term 2 -X*(ybar-beta*xbar)*d.Mtildestar
              term2 <- -X.mat*(ybar-Xbar.beta)*d.Mtildestar

##term 3 -S1X*dQtilde1
              term3 <- -S1.Xmat*d.Qtilde1

##term 4 -S1XX*d.scriptA.1*psi
              term4 <- -S1.XXmat.psi*d.scriptA.1    

##term 5 -P1X*dQtilde2
              term5 <- -P1.Xmat*d.Qtilde2

##term 6 -P1XX*d.scriptA.2*psi
              term6 <- -P1.XXmat.psi*d.scriptA.2     

##term 7 S1X*(ybar-beta*xbar)*dRtilde1
              term7 <- S1.Xmat*(ybar-Xbar.beta)*d.Rtilde1

##term 8 S1XX*(ybar-beta*xbar)*d.F.1*psi
              term8 <- S1.XXmat.psi*(ybar-Xbar.beta)*d.F.1    

##term 9 P1X*(ybar-beta*xbar)*dRtilde2 
              term9 <- P1.Xmat*(ybar-Xbar.beta)*d.Rtilde2

##term 10 P1XX*(ybar-beta*xbar)*d.F.2*psi
              term10 <- P1.XXmat.psi*(ybar-Xbar.beta)*d.F.2    

##term 11 S1.piX*d.scriptA.1*theta.XC
              term11 <- S1.piX.theta.mat*d.scriptA.1

##term 12 2*P1.piX*d.scriptA.2*theta.XC
              term12 <- 2*P1.piX.theta.mat*d.scriptA.2

##term 13 -S1.piX*(ybar-beta*xbar)*d.F.1*theta.XC
              term13 <- -S1.piX.theta.mat*(ybar-Xbar.beta)*d.F.1

##term 14 -2*P1.piX*(ybar-beta*xbar)*d.F.2*theta.XC 
              term14 <- -2*P1.piX.theta.mat*(ybar-Xbar.beta)*d.F.2

              bigPsi.it <- term1 + term2 + term3 + term4 + term5 + term6 + term7 +
                        term8 + term9 + term10 + term11 + term12 + term13 + term14  

              apply(bigPsi.it,2,sum)
         }

         bigPsi.i <- t(sapply(baseData$ID,bigPsi.it.f,t=udt))


################################################Beta estimate variance estimation
         Sigma <- var(bigPsi.i)

         beta.var <- ginv(D)%*%Sigma%*%ginv(t(D))/nrow(baseData)

         beta.se <- sqrt(diag(beta.var))

         beta.se.mat[ii.sim,] <- beta.se
         write.csv(beta.se.mat, beta.se.filename,row.names=FALSE)


#################################scriptA variance estimation
         d.J.f <- function(u){
    
              S1X.mat <- c(S1X1[which(udt==u)], S1X2[which(udt==u)], S1X3[which(udt==u)])
              S2X.mat <- c(S2X1[which(udt==u)], S2X2[which(udt==u)], S2X3[which(udt==u)])

              d.J <- -Ainv.f(u)%*%rbind(S1X.mat, S2X.mat)*estlam.t.gamma[which(udt==u)]
              list(d.J)
         }


####one 2 x p matrix for each unique time
         d.J.L <- sapply(udt, d.J.f)
         d.J.A <- array(unlist(d.J.L), dim = c(nrow(d.J.L[[1]]), ncol(d.J.L[[1]]), length(d.J.L)))
         J.t <- apply(d.J.A,c(1:2),cumsum) #length(udt) x 2 x p


         Qtilde1.it <- apply(d.Qtilde1.it.mat,2,cumsum)
         Qtilde2.it <- apply(d.Qtilde2.it.mat,2,cumsum)


         bigOmega.tt.f <- function(u){

              bigPsi.it <- rbind(Qtilde1.it[which(udt==u),], Qtilde2.it[which(udt==u),]) + J.t[which(udt==u),,]%*%t(bigPsi.i)

              list(var(t(bigPsi.it)))
         }

         bigOmega.tt <- sapply(udt,bigOmega.tt.f)


######SE's for scriptA (*not* d.scriptA)
         scriptA.ses <- do.call(rbind,lapply(bigOmega.tt,function(x){sqrt(diag(x))}))/sqrt(nn)
         scriptA1.se.mat[ii.sim,which(time_discrete %in% udt)] <- scriptA.ses[,1]
         scriptA2.se.mat[ii.sim,which(time_discrete %in% udt)] <- scriptA.ses[,2]

         write.csv(scriptA1.se.mat, scriptA1.ses.filename,row.names=FALSE)
         write.csv(scriptA2.se.mat, scriptA2.ses.filename,row.names=FALSE) 


##############################################################################
#######################################################Lin & Ying 2001 method
##############################################################################

##########for LY-weighted terms
         Ybar_star.numer.LY <- sapply(udt, function(u) {mean((Y_star(u)*exp_gamma(u))*(u<=baseData$C)) } )
         Ybar.LY <- Ybar_star.numer.LY/H0


######for LY beta estimation
         beta.num.LY.f <- function(dd){
              u <- dd['t'] #time
              m <- dd['m']
              X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
              Xbar.LY.mat <- matrix(c(Hbar.X1[which(udt==u)], Hbar.X2[which(udt==u)], Hbar.X3[which(udt==u)]), ncol=1)
              Yval <- dd['Y']

              (X.mat - Xbar.LY.mat)*(Yval-Ybar.LY[which(udt==u)])
         }

         beta.den.LY.f <- function(dd){
              u <- dd['t'] #time
              m <- dd['m']
              X.mat <- matrix(c(dd['X1'], dd['X2'], dd['X3']), ncol=1)
              Xbar.LY.mat <- matrix(c(Hbar.X1[which(udt==u)], Hbar.X2[which(udt==u)], Hbar.X3[which(udt==u)]), ncol=1)

              b.den.LY <- (X.mat-Xbar.LY.mat)%*%t(X.mat-Xbar.LY.mat)
              list(b.den.LY)
         }

         beta.est.num.LY <- apply(apply(sim.data,1,beta.num.LY.f),1,sum)

         beta.est.denom.LY.L <- apply(sim.data,1,beta.den.LY.f)
         beta.est.denom.LY.A <- array(unlist(beta.est.denom.LY.L), dim = c(nrow(beta.est.denom.LY.L[[1]][[1]]), ncol(beta.est.denom.LY.L[[1]][[1]]), length(beta.est.denom.LY.L)))
         beta.est.denom.LY <- apply(beta.est.denom.LY.A,c(1,2),sum)

         beta.est.LY <- ginv(beta.est.denom.LY)%*%matrix(beta.est.num.LY,ncol=1)
         beta.LY.mat[ii.sim,] <- beta.est.LY

         write.csv(beta.LY.mat, beta.est.LY.filename,row.names=FALSE)


############for LY dscriptA estimation
         d.scriptA.LY.f <- function(u){
              sim.sub <- sim.data[which(sim.data$t==u),]
              X.mat <- t(cbind(sim.sub$X1, sim.sub$X2, sim.sub$X3))

              val1 <- sum(sim.sub$Y - as.vector(t(beta.est.LY)%*%X.mat))/nn  

              d.scriptA.LY <- val1/H0[which(udt==u)]
              d.scriptA.LY
         }

##########one value for each unique time
         d.scriptA.LY <- sapply(udt, d.scriptA.LY.f)


###############################################d.M.it.LY (martingale for LY joint model) 

         d.M.it.LY.f <- function(n,t){
  
              C <- baseData$C[baseData$ID==n]

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$ID==n & testdata4$t==t])
              exp.gamma <- exp(gamma.hat[1]*baseData$X1[baseData$ID==n]+
		        gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$ID==n & testdata4$t==t])
 
              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              Yval <- testdata4$Y[testdata4$ID==n]
              Yval[is.na(Yval)] <- 0  ##will multiply by (dNi=0)              

              d.M.it.LY <- (Yval - c(X.mat%*%beta.est.LY))*dNi - (t<=C)*exp.gamma*d.scriptA.LY
              d.M.it.LY
         }

#################dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.M.it.LY <- c(sapply(baseData$ID, d.M.it.LY.f, t=udt))


################################################d.scriptM.it.LY (martingale for observation model) 
         d.scriptM.it.LY.f <- function(n,t){
 
              dNi <- as.numeric(t %in% sim.data$t[sim.data$ID==n])

              d.scriptM.it <- dNi -  (exp(gamma.hat[1]*baseData$X1[baseData$ID==n] +
		           gamma.hat[2]*baseData$X2[baseData$ID==n]+gamma.hat[3]*testdata4$X3[testdata4$ID==n & testdata4$t==t])*
		           estlam.t.gamma[which(udt==t)])*(t<=baseData$C[baseData$ID==n])

              d.scriptM.it
         }

###dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.scriptM.it.LY <- c(sapply(baseData$ID, d.scriptM.it.LY.f, t=udt))


#####################################################d.R.it.LY for variance estimation
         d.R.it.LY.f <- function(n,t){

              C <- baseData$C[baseData$ID==n]

              d.M.LY <- d.M.it.LY[testdata4$ID==n]
              d.scriptM.LY <- d.scriptM.it.LY[testdata4$ID==n]

              Xbar.LY.mat <- cbind(Hbar.X1[which(udt==t)], Hbar.X2[which(udt==t)], Hbar.X3[which(udt==t)])
 
              Xbar.beta.LY <- c(Xbar.LY.mat%*%beta.est.LY)
              ybar.LY <- Ybar.LY[which(udt==t)]
 
              d.R.it.LY <- d.M.LY - (ybar.LY-Xbar.beta.LY)*d.scriptM.LY

              d.R.it.LY
         }

###dim row=length(udt)*col=length(id) then flattened to make 
###length = length(udt) * length(id) sorted by id, then t
         d.R.it.LY <- c(sapply(baseData$ID, d.R.it.LY.f, t=udt))


###########estimate D.LY
         D.LY <- beta.est.denom.LY/nn


#######################LY-weighted terms for H.LY and Omega.LY
         H1.X1X1.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X1)*(u<=baseData$C))})/H0
         H1.X1X2.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*baseData$X2)*(u<=baseData$C))})/H0
         H1.X1X3.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X1*testdata4$X3[testdata4$t==u])*(u<=baseData$C))})/H0
         H1.X2X2.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*baseData$X2)*(u<=baseData$C))})/H0
         H1.X2X3.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*baseData$X2*testdata4$X3[testdata4$t==u])*(u<=baseData$C))})/H0
         H1.X3X3.LY <- sapply(udt, function(u){ mean((exp_gamma(u)*testdata4$X3[testdata4$t==u]*testdata4$X3[testdata4$t==u])*(u<=baseData$C))})/H0

         H1.XXmat.diff.LY.f <- function(t){
              list(matrix(c(H1.X1X1.LY[udt==t],H1.X1X2.LY[udt==t], H1.X1X3.LY[udt==t],
                H1.X1X2.LY[udt==t], H1.X2X2.LY[udt==t], H1.X2X3.LY[udt==t],
                H1.X1X3.LY[udt==t], H1.X2X3.LY[udt==t], H1.X3X3.LY[udt==t]),ncol=3,nrow=3) - 
                outer(c(Hbar.X1[udt==t],Hbar.X2[udt==t], Hbar.X3[udt==t]), c(Hbar.X1[udt==t],Hbar.X2[udt==t], Hbar.X3[udt==t])))
         }

         H1.XXmat.diff.LY.L <- sapply(udt,H1.XXmat.diff.LY.f) 
         H1.XXmat.diff.LY <- array(unlist(H1.XXmat.diff.LY.L), 
                             dim = c(nrow(H1.XXmat.diff.LY.L[[1]]), ncol(H1.XXmat.diff.LY.L[[1]]), length(H1.XXmat.diff.LY.L)))


######################################for L-Y H matrix
         H.it.LY.f <- function(dd){
              u <- dd['t'] #time
              ybar.LY <- Ybar.LY[which(udt==u)]
              Yval <- dd['Y']

              H1.XXmat.diff.LY[,,which(udt==u)]*(Yval - ybar.LY)
  
              H.it.LY <- H1.XXmat.diff.LY[,,which(udt==u)]*(Yval - ybar.LY)
              list(H.it.LY)
         }

         H.it.LY.L <- apply(sim.data,1,H.it.LY.f)
         H.it.LY.A <- array(unlist(H.it.LY.L), 
                      dim = c(nrow(H.it.LY.L[[1]][[1]]), ncol(H.it.LY.L[[1]][[1]]), length(H.it.LY.L)))
         H.LY <- apply(H.it.LY.A,c(1,2),sum)/nn


#########################for L-Y Omega matrix
         Omega.it.LY.f <- function(dd){
              u <- dd['t'] #time

              Omega.it.LY <-  H1.XXmat.diff.LY[,,which(udt==u)]
 
              list(Omega.it.LY)
         }

         Omega.it.LY.L <- apply(sim.data,1,Omega.it.LY.f)
         Omega.it.LY.A <- array(unlist(Omega.it.LY.L), 
                             dim = c(nrow(Omega.it.LY.L[[1]][[1]]), ncol(Omega.it.LY.L[[1]][[1]]), length(Omega.it.LY.L)))
         Omega.LY <- apply(Omega.it.LY.A,c(1,2),sum)/nn

         H.Omega.inv.LY <- H.LY%*%ginv(Omega.LY)


#######################for L-Y V matrix
         V.it.LY.f <- function(n,t){

              X.mat <- cbind(baseData$X1[baseData$ID==n],baseData$X2[baseData$ID==n],  testdata4$X3[testdata4$t==t & testdata4$ID==n])
              Xbar.LY.mat <- cbind(Hbar.X1[which(udt==t)], Hbar.X2[which(udt==t)], Hbar.X3[which(udt==t)])
 
              d.R.i.LY <- d.R.it.LY[testdata4$ID==n]
              d.scriptM.i.LY <- d.scriptM.it.LY[testdata4$ID==n]  

              term1 <- (X.mat - Xbar.LY.mat)*d.R.i.LY
              term2 <- -t(H.Omega.inv.LY%*%t(X.mat - Xbar.LY.mat))*d.scriptM.i.LY
 
              V.it.LY <- term1 + term2
              apply(V.it.LY,2,sum)
         }

         V.i.LY <- t(sapply(baseData$ID,V.it.LY.f,t=udt))
         V.LY <- var(V.i.LY)


#############################L-Y SEs of beta estimate
         beta.var.LY <- ginv(D.LY)%*%V.LY%*%ginv(t(D.LY))/nrow(baseData)
         beta.se.LY <- sqrt(diag(beta.var.LY))

         beta.se.LY.mat[ii.sim,] <- beta.se.LY
         write.csv(beta.se.LY.mat, beta.se.LY.filename,row.names=FALSE)

      } ##end simulated dataset iterate
    }  ##end gamma3 iterate
  }  ##end alpha.20(t) iterate
}  ## end nn iterate
