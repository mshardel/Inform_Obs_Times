####################################################################################
#Tables of simulation results for 
#"Time-Varying Latent Eect Models for Repeated Measurements to Address
#Informative Observation Times in the U.S. Medicare Minimum Data Set"  
#By Shardell M, Chen C, and Falvey J
#University of Maryland School of Medicine

#Three outcome model covariates of interest: X1 and X2 are time-invariant, X3 is time-varying
#Four observation model covariates: Z1=X1, Z2=X2, Z3=X3, and Z4 is time-varying
#Comparing proposed inverse conditional intensity rate ratio (ICIRR) weighted
#method to Buzkova and Lumley (2009) (BL)

#Creates Tables S.2 and S.3.
####################################################################################

#read in names of files
sim.dir <- ".\\wt_sims_results"

#creates directory for post-processed simulation results
cur.dir <- ".\\wt_sims_results\\wt_sims_postprocess"
ifelse(!dir.exists(file.path(cur.dir)), dir.create(file.path(cur.dir)), FALSE)

the.files <- sort(list.files(sim.dir))
the.files <- the.files[grepl(".csv", the.files)]

################################# beta ests ##################################
################################# Table S.2 ##################################

beta.true <- c(1, -1, 1)

#############################function to read in file contents, compute means and SEs
######################################################################################
####function to read in contents and post-process
read.beta_ests.f <- function(x){

  dd.ests <- read.csv(paste0(sim.dir,"\\",x[1]))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)     
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)
  
  the.means <- formatC(the.means,format='f',digits=3) 
  the.emp.ses <- formatC(the.emp.ses,format='f',digits=3)
  
  the.names <- rep(c("MeanEst", "ESE"), ncol(dd.ests))
 
  results <- c(rbind(the.means, the.emp.ses))
  names(results) <- the.names 
  return(results)         
}


#######################################################################
###################ICIRR beta ests (proposed)##########################
#######################################################################

###subset of files names with beta_ests with ICIRR

the.files.beta <- the.files[grepl("beta_ests_Z_ICIRR", the.files)]


###get n, a2, delta4
pre.the.files.beta_ests.split <- unlist(strsplit(the.files.beta, ".csv", fixed=T))
the.files.beta_ests.split <- strsplit(pre.the.files.beta_ests.split, "_", fixed=T)
n.vec <-  unlist(lapply(the.files.beta_ests.split,function(x){x[6]}))
a2.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[8]}))
a2.vec[a2.vec=="1.over.(1+t)"] <- "1/(1+t)"
delta4.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[10]}))


#######################################################################

               #results of proposed method ICIRR#  

beta.est.results.proposed.ICIRR <- t(sapply(the.files.beta,  read.beta_ests.f))

params.proposed.ICIRR <- cbind(rep("ICIRR", nrow(beta.est.results.proposed.ICIRR)),n.vec, a2.vec, delta4.vec)
colnames(params.proposed.ICIRR) <- c("Method", "n", "alpha2","delta4")


full.beta.est.results.proposed.ICIRR <- cbind(params.proposed.ICIRR,beta.est.results.proposed.ICIRR) 
#############################################################################


#############################################################################
#######################beta ests (BL)########################################
#############################################################################

###subset of files names with beta_ests with BL

the.files.beta_BL <- the.files[grepl("beta_ests_Z_BL", the.files)]


###get n, a2, delta4
pre.the.files.beta_ests.split <- unlist(strsplit(the.files.beta_BL, ".csv", fixed=T))
the.files.beta_ests.split <- strsplit(pre.the.files.beta_ests.split, "_", fixed=T)
n.vec <-  unlist(lapply(the.files.beta_ests.split,function(x){x[6]}))
a2.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[8]}))
a2.vec[a2.vec=="1.over.(1+t)"] <- "1/(1+t)"
delta4.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[10]}))
#############################################################################

               #results of BL method# 

beta.est.results.BL <- t(sapply(the.files.beta_BL, read.beta_ests.f))

params.BL <- cbind(rep("BL", nrow(beta.est.results.BL)),n.vec, a2.vec, delta4.vec)
colnames(params.BL) <- c("Method", "n", "alpha2","delta4")

full.beta.est.results.BL <- cbind(params.BL,beta.est.results.BL) 

colnames(full.beta.est.results.BL) <- paste0("BL", colnames(full.beta.est.results.BL))
#############################################################################

#######################combine beta sims results and write####################

all.beta.est.results <- cbind(full.beta.est.results.proposed.ICIRR, full.beta.est.results.BL[,-c(1:4)])[,-c(1)] 
all.beta.est.results.df <- data.frame(all.beta.est.results)
all.beta.est.results.df[,c(2)] <- factor(all.beta.est.results.df[,c(2)])
all.beta.est.results.df[,c(1,3)] <- apply(all.beta.est.results.df[,c(1,3)], 2, as.numeric)

oo <- with(all.beta.est.results.df, order(n, alpha2, delta4))

all.beta.est.results <- all.beta.est.results[oo,]

############################# Table S.2 #####################################
write.csv(all.beta.est.results, paste0(cur.dir,"\\TableS.2_beta_ICIRR-weighted-proposed_Buzkova-Lumley.csv"),row.names=F)
#############################################################################


#############################################################################
#####delta ests from Lin et al (2000) (note: does not depend on alpha2)######
########################### Table S.3 #######################################
#############################################################################


###subset of file names with delta
the.files.delta <- the.files[grepl("delta_ests", the.files)]

###subset of files with one constant a2
the.files.delta.sub <- the.files.delta[a2.vec==a2.vec[1]]

delta.true.mat <- cbind(rep(-0.5, length(the.files.delta.sub)),rep(-0.25, length(the.files.delta.sub)), rep(0.25, length(the.files.delta.sub)),
              as.numeric(delta4.vec[a2.vec==a2.vec[1]]))

 
######read in file contents, compute means and SEs

####function to read in contents of scalar parameters and post-process
read.delta_ests.f <- function(x){

  dd.ests <- read.csv(paste0(sim.dir,"\\",x[1]))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)
  the.means <- formatC(the.means,format='f',digits=3) 
  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)
  the.emp.ses <- formatC(the.emp.ses,format='f',digits=3) 

  delta.true <- delta.true.mat[which(the.files.delta.sub[1]==x[1]),]

  the.names <- rep(c("Est", "ESE"), ncol(dd.ests))
  results <- c(rbind(the.means, the.emp.ses))

  names(results) <- the.names 
  return(results)         
}

##########################################################################

               #results of Lin et al (2000) delta estimation#

delta.est.results <- t(sapply(the.files.delta.sub, read.delta_ests.f))

params.sub <- cbind(n.vec[a2.vec==a2.vec[1]], delta4.vec[a2.vec==a2.vec[1]])
colnames(params.sub) <- c("n", "delta4")

param.sub <- apply(params.sub,2,as.numeric)

oo <- with(data.frame(param.sub), order(n, delta4))

full.delta.est.results <- cbind(params.sub,delta.est.results) 
full.delta.est.results <- full.delta.est.results[oo,]

###################################Table S.3 ###################################
write.csv(full.delta.est.results, paste0(cur.dir,"\\TableS.3_delta_ICIRR-weighted-proposed_Buzkova-Lumley.csv"),row.names=F)
##########################################################################

