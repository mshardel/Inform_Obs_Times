####################################################################################
#Tables and Figures of simulation results for 
#"Time-Varying Latent Eect Models for Repeated Measurements to Address
#Informative Observation Times in the U.S. Medicare Minimum Data Set"  
#By Shardell M, Chen C, and Falvey J
#University of Maryland School of Medicine

#Three covariates: X1 and X2 are time-invariant, X3 is time-varying
#Comparing proposed method to Lin and Ying (2001) (LY)
#Observation model coefficients estimated using Lin et al (2000)

#Creates Figures 1, 2, S.1, and S2; and Tables 1 and S.1.
####################################################################################
library(ggplot2)
library(gridExtra)

#read in names of files
sim.dir <- ".\\sims_results"

#creates directory for post-processed simulation results
cur.dir <- ".\\sims_results\\sims_postprocess"
ifelse(!dir.exists(file.path(cur.dir)), dir.create(file.path(cur.dir)), FALSE)

the.files <- sort(list.files(sim.dir))
the.files <- the.files[grepl(".csv", the.files)]

#################################beta ests##################################
################################# Table 1 ##################################

beta.true <- c(1, -1, 1)

############################################################################
############################beta ests (proposed)############################
############################################################################
###subset of file names with beta_ests without LY
the.files.beta_ests <- the.files[grepl("beta_ests", the.files) & !grepl("LY", the.files)]

###subset of file names with beta_ses without LY
the.files.beta_ses <- the.files[grepl("beta_ses", the.files) & !grepl("LY", the.files)]

###combining estimate and SE filenames
the.files.beta <- cbind(the.files.beta_ests, the.files.beta_ses)

###get n, a2, gamma3
pre.the.files.beta_ests.split <- unlist(strsplit(the.files.beta_ests, ".csv", fixed=T))
the.files.beta_ests.split <- strsplit(pre.the.files.beta_ests.split, "_", fixed=T)
n.vec <-  unlist(lapply(the.files.beta_ests.split,function(x){x[4]}))
a2.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[6]}))
a2.vec[a2.vec=="1.over.(1+t)"] <- "1/(1+t)"
gamma3.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[8]}))

#############################read in file contents, compute means and SEs

####function to read in contents of scalar parameters and post-process
read.beta_ests.f <- function(x){

  dd.ests <- read.csv(paste0(sim.dir,"\\",x[1]))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)
  the.means <- formatC(the.means,format='f',digits=3)  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)
  the.emp.ses <- formatC(the.emp.ses,format='f',digits=3)

  dd.ses <- read.csv(paste0(sim.dir,"\\",x[2]))
  the.num.ses <- apply(dd.ses,2, mean, na.rm=TRUE)  
  the.num.ses <- formatC(the.num.ses,format='f',digits=3)  
   
  lcl <- (dd.ests - 2*dd.ses)
  lcl <- lcl[!is.na(lcl[,1]),]

  ucl <- (dd.ests + 2*dd.ses)
  ucl <- ucl[!is.na(ucl[,1]),]

  cover.lcl <- t(apply(lcl,1,function(x){x<=beta.true}))
  cover.ucl <- t(apply(ucl,1,function(x){x>=beta.true}))

  cover <- cover.lcl*cover.ucl
  cover.perc <- apply(cover,2, mean)*100
  cover.perc <- formatC(cover.perc,format='f',digits=1)

  the.names <- rep(c("Est", "SE", "ESE", "Cov %"), ncol(dd.ests))
 
  results <- c(rbind(the.means, the.num.ses, the.emp.ses, cover.perc))
  names(results) <- the.names 
  return(results)         
}
##########################################################################

                #results of proposed method#

beta.est.results.proposed <- t(apply(the.files.beta, 1, read.beta_ests.f))

params.proposed <- cbind(rep("Proposed", nrow(beta.est.results.proposed)),n.vec, a2.vec, gamma3.vec)
colnames(params.proposed) <- c("Method", "n", "alpha2","gamma3")

full.beta.est.results.proposed <- cbind(params.proposed,beta.est.results.proposed) 
##########################################################################


##########################################################################
#############################beta ests (LY)###############################
##########################################################################

###subset of files names with beta_ests with LY

the.files.beta_ests_LY <- the.files[grepl("beta_ests", the.files) & grepl("LY", the.files)]

###subset of file names with beta_ses with LY
the.files.beta_ses_LY <- the.files[grepl("beta_ses", the.files) & grepl("LY", the.files)]

###combining estimate and SE filenames
the.files.beta_LY <- cbind(the.files.beta_ests_LY, the.files.beta_ses_LY)

###get n, a2, gamma3
pre.the.files.beta_ests.split <- unlist(strsplit(the.files.beta_ests_LY, ".csv", fixed=T))
the.files.beta_ests.split <- strsplit(pre.the.files.beta_ests.split, "_", fixed=T)
n.vec <-  unlist(lapply(the.files.beta_ests.split,function(x){x[5]}))
a2.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[7]}))
a2.vec[a2.vec=="1.over.(1+t)"] <- "1/(1+t)"
gamma3.vec <- unlist(lapply(the.files.beta_ests.split,function(x){x[9]}))
##########################################################################

                #results of LY method#

beta.est.results_LY <- t(apply(the.files.beta_LY, 1, read.beta_ests.f))

params.LY <- cbind(rep("LY", nrow(beta.est.results_LY)),n.vec, a2.vec, gamma3.vec)
colnames(params.LY) <- c("Method", "n", "alpha2","gamma3")

full.beta.est.results_LY <- cbind(params.LY,beta.est.results_LY) 
##############################################################################

#######################combine beta sims results and write####################

all.beta.est.results <- rbind(full.beta.est.results.proposed, full.beta.est.results_LY) 

all.beta.est.results.df <- data.frame(all.beta.est.results)
all.beta.est.results.df[,c(1,3)] <- apply(all.beta.est.results.df[,c(1,3)], 2, factor)
all.beta.est.results.df[,c(2,4)] <- apply(all.beta.est.results.df[,c(2,4)], 2, as.numeric)

oo <- with(all.beta.est.results.df, order(rev(Method), n, alpha2, gamma3))

all.beta.est.results <- all.beta.est.results[oo,]

############################# Table 1 #########################################
write.csv(all.beta.est.results, paste0(cur.dir,"\\Table1_beta_proposed_Lin-Ying.csv"),row.names=F)
#############################################################################


#############################################################################
#####gamma ests from Lin et al (2000) (note: does not depend on alpha2)######
########################### Table S.1 #######################################
#############################################################################

###subset of file names with gamma
the.files.gamma_ests <- the.files[grepl("gamma_ests", the.files)]

###subset of file names with gamma
the.files.gamma_ses <- the.files[grepl("gamma_ses", the.files)]

###combining estimate and SE filenames
the.files.gamma <- cbind(the.files.gamma_ests, the.files.gamma_ses)

###subset of files with one constant a2
the.files.gamma.sub <- the.files.gamma[a2.vec==a2.vec[1],]

gamma.true.mat <- cbind(rep(-0.5, nrow(the.files.gamma.sub)),rep(-0.25, nrow(the.files.gamma.sub)), 
              as.numeric(gamma3.vec[a2.vec==a2.vec[1]]))

 
#######################read in file contents, compute means and SEs

####function to read in contents of scalar parameters and post-process
read.gamma_ests.f <- function(x){

  dd.ests <- read.csv(paste0(sim.dir,"\\",x[1]))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)
  the.means <- formatC(the.means,format='f',digits=3) 
  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)
  the.emp.ses <- formatC(the.emp.ses,format='f',digits=3) 

  dd.ses <- read.csv(paste0(sim.dir,"\\",x[2]))
  the.num.ses <- apply(dd.ses,2, mean, na.rm=TRUE)  
  the.num.ses <- formatC(the.num.ses,format='f',digits=3)

   
  lcl <- (dd.ests - 2*dd.ses)
  lcl <- lcl[!is.na(lcl[,1]),]

  ucl <- (dd.ests + 2*dd.ses)
  ucl <- ucl[!is.na(ucl[,1]),]

  gamma.true <- gamma.true.mat[which(the.files.gamma.sub[,1]==x[1]),]

  cover.lcl <- t(apply(lcl,1,function(x){x<=gamma.true}))
  cover.ucl <- t(apply(ucl,1,function(x){x>=gamma.true}))

  cover <- cover.lcl*cover.ucl
  cover.perc <- apply(cover,2, mean)*100
  cover.perc <- formatC(cover.perc,format='f',digits=1)

  #the.names <- rep(c("Est", "SE", "ESE", "Cov %"), ncol(dd.ests))
  #results <- c(rbind(the.means, the.num.ses, the.emp.ses, cover.perc))

  the.names <- rep(c("Est", "SE"), ncol(dd.ests))
  results <- c(rbind(the.means, the.num.ses))

  names(results) <- the.names 
  return(results)         
}
##############################################################################

                #results of Lin et al (2000) gamma estimation#

gamma.est.results <- t(apply(the.files.gamma.sub, 1, read.gamma_ests.f))

params.sub <- cbind(n.vec[a2.vec==a2.vec[1]], gamma3.vec[a2.vec==a2.vec[1]])
colnames(params.sub) <- c("n", "gamma3")

param.sub <- apply(params.sub,2,as.numeric)

oo <- with(data.frame(param.sub), order(n, gamma3))

full.gamma.est.results <- cbind(params.sub,gamma.est.results) 
full.gamma.est.results <- full.gamma.est.results[oo,]

###################################Table S.1 ###################################
write.csv(full.gamma.est.results, paste0(cur.dir,"\\TableS.1_gamma_proposed_Lin-Ying.csv"),row.names=F)
#############################################################################


#############################################################################
############################script A2 ests###################################
############################Figures 2 and S.2 ###############################
#############################################################################

###subset of files names with dscriptA2_ests
##NOTE: take the cumsum to derive scriptA2_ests and empirical SEs
the.files.dscriptA2_ests <- the.files[grepl("dscriptA2_ests", the.files)]

###subset of file names with scriptA2_ses
the.files.scriptA2_ses <- the.files[grepl("scriptA2_ses", the.files)]

###combining estimate and SE filenames
the.files.scriptA2 <- cbind(the.files.dscriptA2_ests, the.files.scriptA2_ses)

lambda0.t <- 0.1
tau <-2 #max follow-up time
time_discrete <- seq(0.02,tau,0.02)
alpha.10.t.true <- 1 + sin(time_discrete)

oo <- with(all.beta.est.results.df, order(rev(Method), n, alpha2, gamma3))
##########################################################################

#####################Deriving SEs for A1 and A2###########################
scriptAx_CI_lims.f <- function(x){

  dd.ests.d <- read.csv(paste0(sim.dir,"\\",x[1]))
  dd.ests <- t(apply(dd.ests.d, 1, cumsum))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)

  dd.ses <- read.csv(paste0(sim.dir,"\\",x[2]))
  dd.ses <- dd.ses[!is.na(dd.ses[,1]),]
  the.num.ses <- apply(dd.ses,2, mean)  
   
  num.lcl <- (the.means - 2*the.num.ses)
  num.ucl <- (the.means + 2*the.num.ses)

  cl.lims <- c(min(num.lcl[c(1:93)]), max(num.ucl[c(1:93)]))
  return(cl.lims)
}

A2.ylims.500.pre <- apply(the.files.scriptA2[n.vec=="500",],1,scriptAx_CI_lims.f)  
A2.ylims.500 <- c(min(A2.ylims.500.pre[1,]), max(A2.ylims.500.pre[2,]))

A2.ylims.250.pre <- apply(the.files.scriptA2[n.vec=="250",],1,scriptAx_CI_lims.f) 
A2.ylims.250 <- c(min(A2.ylims.250.pre[1,]), max(A2.ylims.250.pre[2,]))


read.scriptA2_ests.f <- function(x,ylim=NULL){

  dd.ests.d <- read.csv(paste0(sim.dir,"\\",x[1]))
  dd.ests <- t(apply(dd.ests.d, 1, cumsum))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)

  dd.ses <- read.csv(paste0(sim.dir,"\\",x[2]))
  dd.ses <- dd.ses[!is.na(dd.ses[,1]),]
  the.num.ses <- apply(dd.ses,2, mean)  
   
  num.lcl <- (the.means - 2*the.num.ses)
  num.ucl <- (the.means + 2*the.num.ses)


  the.n <- n.vec[which(the.files.dscriptA2_ests==x[1])]
  the.gamma3 <-  gamma3.vec[which(the.files.dscriptA2_ests==x[1])]
  a2.true <- sapply(time_discrete, function(t){eval(parse(text = a2.vec[which(the.files.dscriptA2_ests==x[1])]))})
  the.a2 <- a2.vec[which(the.files.dscriptA2_ests==x[1])]  

  scriptA2.true <- cumsum(a2.true*lambda0.t)

  ##create a figure   
  result.data <- data.frame(time_discrete, scriptA2.true, the.means, num.lcl, num.ucl)

  myYlab = bquote( A[2](t))
  mytitle = bquote(gamma["03"]~"="~.(the.gamma3)*","~alpha[20](t)~"="~.(the.a2))

  pp <- ggplot(result.data[1:93,], aes(time_discrete, scriptA2.true)) +
      geom_line(aes(time_discrete,scriptA2.true), linetype="solid") +
      geom_line(aes(time_discrete,the.means), linetype="dashed") +
      geom_ribbon(aes(ymin=num.lcl, ymax=num.ucl), alpha=0.5, fill = "darkgray")+
      #geom_line(aes(time_discrete,num.lcl), linetype="dotdash") +
      #geom_line(aes(time_discrete,num.ucl), linetype="dotdash") +
      coord_cartesian(xlim = c(min(time_discrete), max(time_discrete)),ylim=ylim) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         plot.title = element_text(size=10)) +
      xlab("t") + ylab(myYlab)+ 
      ggtitle(mytitle)
      
      return(pp)
}

########################### Figure S.2 ###################################################
A2.plots.250 <- apply(the.files.scriptA2[n.vec=="250",][oo[1:9],],1,read.scriptA2_ests.f,ylim=A2.ylims.250)
A2.plots.250.grob <- marrangeGrob(A2.plots.250, nrow=3,ncol=3, top=NULL)
ggsave(file=paste0(cur.dir,"\\","FigureS.2_A2_plots_250_proposed.pdf"), A2.plots.250.grob, width=6.6, height=8.0, units="in")
ggsave(file=paste0(cur.dir,"\\","FigureS.2_A2_plots_250_proposed.png"), A2.plots.250.grob, width=6.6, height=8.0, units="in")

########################### Figure 2 ####################################################
A2.plots.500 <- apply(the.files.scriptA2[n.vec=="500",][oo[1:9],],1,read.scriptA2_ests.f, ylim=A2.ylims.500)
A2.plots.500.grob <- marrangeGrob(A2.plots.500, nrow=3,ncol=3, top=NULL)
ggsave(file=paste0(cur.dir,"\\","Figure2_A2_plots_500_proposed.pdf"), A2.plots.500.grob, width=6.6, height=8.0, units="in")
ggsave(file=paste0(cur.dir,"\\","Figure2_A2_plots_500_proposed.png"), A2.plots.500.grob, width=6.6, height=8.0, units="in")
##########################################################################

#######################################################################
##########################script A1 ests ##############################
##########################Figures 1 and S.1############################
#######################################################################

###subset of files names with dscriptA1_ests
##NOTE: Take the cumsum to derive scriptA1_ests and empirical SEs
the.files.dscriptA1_ests <- the.files[grepl("dscriptA1_ests", the.files)]

###subset of file names with scriptA1_ses
the.files.scriptA1_ses <- the.files[grepl("scriptA1_ses", the.files)]

###combining estimate and SE filenames
the.files.scriptA1 <- cbind(the.files.dscriptA1_ests, the.files.scriptA1_ses)

lambda0.t <- 0.1
tau <-2 #max follow-up time
time_discrete <- seq(0.02,tau,0.02)
alpha.10.t.true <- 1 + sin(time_discrete)

A1.ylims.500.pre <- apply(the.files.scriptA1[n.vec=="500",],1,scriptAx_CI_lims.f)  
A1.ylims.500 <- c(min(A1.ylims.500.pre[1,]), max(A1.ylims.500.pre[2,]))

A1.ylims.250.pre <- apply(the.files.scriptA1[n.vec=="250",],1,scriptAx_CI_lims.f) 
A1.ylims.250 <- c(min(A1.ylims.250.pre[1,]), max(A1.ylims.250.pre[2,]))


read.scriptA1_ests.f <- function(x,ylim=NULL){

  dd.ests.d <- read.csv(paste0(sim.dir,"\\",x[1]))
  dd.ests <- t(apply(dd.ests.d, 1, cumsum))

  the.means <- apply(dd.ests,2, mean, na.rm=TRUE)  
  the.emp.ses <-   apply(dd.ests,2, sd, na.rm=TRUE)

  dd.ses <- read.csv(paste0(sim.dir,"\\",x[2]))
  dd.ses <- dd.ses[!is.na(dd.ses[,1]),]
  the.num.ses <- apply(dd.ses,2, mean)  
   
  num.lcl <- (the.means - 2*the.num.ses) 
  num.ucl <- (the.means + 2*the.num.ses)
 

  the.n <- n.vec[which(the.files.dscriptA1_ests==x[1])]
  the.gamma3 <-  gamma3.vec[which(the.files.dscriptA1_ests==x[1])]
  the.a2 <- a2.vec[which(the.files.dscriptA1_ests==x[1])]  

  scriptA1.true <- cumsum(alpha.10.t.true*lambda0.t)

  ##create a figure   
  result.data <- data.frame(time_discrete, scriptA1.true, the.means, num.lcl, num.ucl)

  myYlab = bquote( A[1](t))
  mytitle = bquote(gamma["03"]~"="~.(the.gamma3)*","~alpha[20](t)~"="~.(the.a2))

  pp <- ggplot(result.data[1:93,], aes(time_discrete, scriptA1.true)) +
      geom_line(aes(time_discrete,scriptA1.true), linetype="solid") +
      geom_line(aes(time_discrete,the.means), linetype="dashed") +
      geom_ribbon(aes(ymin=num.lcl, ymax=num.ucl), alpha=0.5, fill = "darkgray")+
      #geom_line(aes(time_discrete,num.lcl), linetype="dotdash") +
      #geom_line(aes(time_discrete,num.ucl), linetype="dotdash") +
      coord_cartesian(xlim = c(min(time_discrete), max(time_discrete)),ylim=ylim) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         plot.title = element_text(size=10)) +
      xlab("t") + ylab(myYlab)+ 
      ggtitle(mytitle)
 
  return(pp)       
}


########################### Figure S.1 ####################################################
A1.plots.250 <- apply(the.files.scriptA1[n.vec=="250",][oo[1:9],],1,read.scriptA1_ests.f,ylim=A1.ylims.250)
A1.plots.250.grob <- marrangeGrob(A1.plots.250, nrow=3,ncol=3, top=NULL)
ggsave(file=paste0(cur.dir,"\\","FigureS.1_A1_plots_250_proposed.pdf"), A1.plots.250.grob, width=6.6, height=8.0, units="in")
ggsave(file=paste0(cur.dir,"\\","FigureS.1_A1_plots_250_proposed.png"), A1.plots.250.grob, width=6.6, height=8.0, units="in")


########################### Figure 1 #####################################################
A1.plots.500 <- apply(the.files.scriptA1[n.vec=="500",][oo[1:9],],1,read.scriptA1_ests.f,ylim=A1.ylims.500)
A1.plots.500.grob <- marrangeGrob(A1.plots.500, nrow=3,ncol=3, top=NULL)
ggsave(file=paste0(cur.dir,"\\","Figure1_A1_plots_500_proposed.pdf"), A1.plots.500.grob, width=6.6, height=8.0, units="in")
ggsave(file=paste0(cur.dir,"\\","Figure1_A1_plots_500_proposed.png"), A1.plots.500.grob, width=6.6, height=8.0, units="in")
##########################################################################

