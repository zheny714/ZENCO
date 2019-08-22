# ***README***
# The following script is used to fit Zero-inflated count data to 
# the ZENCO model proposed in the submitted manuscript titled 
# "Modeling Dynamic Correlation in Zero-inflated Bivariate Count 
# Data with Applications to Single-cell RNA Sequencing Data"


library(rjags)
library(purrr)
library(MASS)


###########################################################
# Zero Inflated Bivariate Negative Binomial Data Simulation 
# tau0 = t0, tau1 = t1
###########################################################
sim.Z <- function(t0, t1){
  N <- 1000
  phi <- c(4,4) 
  mu <- c(15,15)
  lambda0 <- 0.1 # mean of poisson
  b0 = 0.1401 # parameter of dropout rate
  b1 = -0.02393 # can be fitted from real data 
  tau0 = t0
  tau1 = t1
  # dropout rate
  p <- c(exp(b0+b1*mu[1])/(1+exp(b0+b1*mu[1])), exp(b0+b1*mu[2])/(1+exp(b0+b1*mu[2]))) 
  cp <- mu
  # x3
  mu3 <- 15
  phi3 <- 4 
  x3 <- rnbinom(N, size=1/phi3, mu=mu3)
  p3 <- rbernoulli(N, p=exp(b0+b1*mu3)/(1+exp(b0+b1*mu3))) 
  z3 <- rpois(N, lambda0) # poisson
  x3 <- (1-p3)*x3+p3*z3  
  # if x3 dropout, x3 is mu3 
  rho.mu3 <- (exp(tau0+tau1*mu3)-1)/(exp(tau0+tau1*mu3)+1)
  rho <- (exp(tau0+tau1*x3)-1)/(exp(tau0+tau1*x3)+1)
  rho <- (1-p3)*rho+p3*rho.mu3  
  # quadratic
  cp[1] <- 1/phi[1]
  cp[2] <- 1/phi[2] 
  z <- matrix(rep(0,2*N), ncol=2)
  y <- matrix(rep(0,2*N), ncol=2)
  for(i in 1:N){
    Sigma <- matrix(c(1,rho[i],rho[i],1),2,2)
    z[i, ] <- mvrnorm(n = 1, rep(0, 2), Sigma)
    for(j in 1:2){
      y[i,j] <- rpois(1,qgamma(pnorm(z[i,j]),cp[j],cp[j])*mu[j])
    }
  }
  
  # Add "Zero"
  z1 <- rpois(N, lambda0) # poisson
  z2 <- rpois(N, lambda0)
  p1 <- rbernoulli(N, p=p[1])
  p2 <- rbernoulli(N, p=p[2])
  y1 <- (1-p1)*y[,1]+p1*z1
  y2 <- (1-p2)*y[,2]+p2*z2
  y <- cbind(y1, y2)
  return(cbind(y,x3))
}

#####################################################
# Use the simulated data to model dynamic correlation
# with four different methods
#####################################################



##############################################
# Method 1: Without considering Zero-inflation
##############################################
fun1 <- function(dat){
  ### JAGS ###
  Ind.string <-"
  model{
  for (i in 1:N){
  rho[i] <- (exp(tau0+tau1*x[i])-1)/(exp(tau0+tau1*x[i])+1)
  sigma[1,1,i] <- 1/(1-rho[i]*rho[i])
  sigma[2,2,i] <- 1/(1-rho[i]*rho[i])
  sigma[1,2,i] <- -rho[i]/(1-rho[i]*rho[i])
  sigma[2,1,i] <- -rho[i]/(1-rho[i]*rho[i])
  Z[i,1:2] ~ dmnorm(rep(0,2),sigma[,,i])
  
  for (j in 1:2){
  Y[i,j] ~ dpois(lambda[i,j])
  lambda[i,j] <- qgamma(pnorm(Z[i,j],0,1), inverphi[j], inverphi[j])*mu[j]
  }
  }
  tau0 ~ dnorm(0, 250)
  tau1 ~ dnorm(0, 250)
  for (i in 1:2){
  mu[i] ~ dlnorm(0, 1)
  inverphi[i] ~ dgamma(1, 0.01)
  }
  }
  "
  y.ind <- cbind(dat[,1], dat[,2])
  x3.ind <- dat[,3]
  Ind.spec <- textConnection(Ind.string)
  jags_data = list(N=length(dat[,1]),Y=y.ind, x=x3.ind)
  phi.ind <- c((var(y.ind[,1])/mean(y.ind[,1])-1)/mean(y.ind[,1]), (var(y.ind[,2])/mean(y.ind[,2])-1)/mean(y.ind[,2]))
  jags_inits = list(mu=c(mean(y.ind[,1]), mean(y.ind[,2])), inverphi=1/phi.ind, tau0=0, tau1=0)
  jags_model_Ind = jags.model(Ind.spec, data=jags_data, n.adapt=10000, inits=jags_inits, n.chain=3)
  update(jags_model_Ind, 5000)
  samps.coda.Ind <- coda.samples(jags_model_Ind, c('mu', 'inverphi', 'tau0', 'tau1'), n.iter = 50000, thin=10, n.burnin=10000)
  
  # summary(samps.coda.Ind)
  samps.Ind<- summary(samps.coda.Ind)
  return(samps.Ind)
}





#######################################
# Method 2: With Zero-inflation (ZENCO)
#######################################
fun2 <- function(dat){
  ### JAGS ###
  IndZ.string <-"
  model{
  for (i in 1:N){
  p[i, 1] ~ dbern(exp(0.1401-0.02393*mu[1])/(1+exp(0.1401-0.02393*mu[1])))
  p[i, 2] ~ dbern(exp(0.1401-0.02393*mu[2])/(1+exp(0.1401-0.02393*mu[2])))
  p3[i] ~ dbern(exp(0.1401-0.02393*mu3)/(1+exp(0.1401-0.02393*mu3)))
  
  x[i] ~ dpois(lmd3[i])
  h[i] ~ dgamma(inverphi3, inverphi3)
  lmd3[i] <- ifelse(p3[i]==0, h[i]*mu3, 0.1)
  
  
  rho_tmp[i] <- ifelse(p3[i]==0, ((exp(tau0+tau1*x[i])-1)/(exp(tau0+tau1*x[i])+1)), (exp(tau0+tau1*mu3)-1)/(exp(tau0+tau1*mu3)+1))
  rho[i] <- ifelse(abs(rho_tmp[i])==1, rho_tmp[i]/(1+10^-9), rho_tmp[i])
  sigma[1,1,i] <- 1/(1-rho[i]*rho[i])
  sigma[2,2,i] <- 1/(1-rho[i]*rho[i])
  sigma[1,2,i] <- -rho[i]/(1-rho[i]*rho[i])
  sigma[2,1,i] <- -rho[i]/(1-rho[i]*rho[i])
  Z[i,1:2] ~ dmnorm(rep(0,2),sigma[,,i])
  
  for (j in 1:2){
  Y[i,j] ~ dpois(lambda[i,j])
  lambda[i,j] <- ifelse(p[i,j]==0, qgamma(pnorm(Z[i,j],0,1),inverphi[j],inverphi[j])*mu[j], 0.1)
  }
  }
  tau0 ~ dnorm(0, 250)
  tau1 ~ dnorm(0, 250)
  mu3 ~ dlnorm(0, 1)
  inverphi3 ~ dgamma(1, 0.01)
  for (i in 1:2){
  mu[i] ~ dlnorm(0, 1)
  inverphi[i] ~ dgamma(1, 0.01)
  }
  }
  "
  y.indz <- cbind(dat[,1], dat[,2])
  x3.indz <- dat[,3]
  IndZ.spec <- textConnection(IndZ.string)
  jags_data = list(N=length(dat[,1]),Y=y.indz, x=x3.indz)
  phi.indz <- c((var(y.indz[,1])/mean(y.indz[,1])-1)/mean(y.indz[,1]), (var(y.indz[,2])/mean(y.indz[,2])-1)/mean(y.indz[,2]))
  jags_inits = list(mu=c(mean(y.indz[,1]), mean(y.indz[,2])), inverphi=1/phi.indz, mu3=mean(x3.indz), inverphi3=1/((var(x3.indz)/mean(x3.indz)-1)/mean(x3.indz)),tau0=0, tau1=0)
  jags_model_IndZ = jags.model(IndZ.spec, data=jags_data, n.adapt=10000, inits=jags_inits, n.chain=3)
  update(jags_model_IndZ, 5000)
  samps.coda.IndZ <- coda.samples(jags_model_IndZ, c('mu', 'inverphi', 'tau0', 'tau1', 'mu3', 'inverphi3'), n.iter = 50000, thin=10, n.burnin=10000)
  
  # summary(samps.coda.IndZ)
  samps.IndZ<- summary(samps.coda.IndZ)
  return(samps.IndZ)
}



######################################
# Method 3: TLA
######################################

# Calculate TLA using bootstrap standard error
TLA <- function(data123, B){
  est.boot <- rep(0, B)
  for (i in 1:B){
    idx <- sample(1:length(data123), length(data123), replace=T)
    tmp <- data123[idx]
    est.boot[i] <- mean(tmp) 
  }
  # bootstrap standard error
  se.boot <- sqrt(sum((est.boot-mean(est.boot))^2/(B-1)))
  # TLA
  return(mean(data123)/se.boot)
}

# p-value from permutation
TLA.pval <- function(x123){
  # standardize data
  x1.std <- scale(x123[,1])
  x2.std <- scale(x123[,2])
  x3.std <- scale(x123[,3])
  
  # three product moment
  x.product <- x1.std*x2.std*x3.std
  
  
  # library(boot)
  # boot(data=cbind(y,x3), statistic=betfun,R=100)
  # betfun = function(data, ind){  
  #  return(mean(data[ind,1]*data[ind,2]*data[ind,3])) 
  # }
  
  # Calcualte TLA from real data 
  TLA.real <- TLA(x.product, 500)
  
  # Use permutation test to calcualte p-value
  # H0: tau1=0
  TLA.perm <- rep(0, 1000)
  for (j in 1:1000){
    x3.perm <- sample(x3.std)
    x.perm <- x1.std*x2.std*x3.perm
    TLA.perm[j] <- TLA(x.perm, 500)
  }
  # hist(TLA.perm)
  # abline(v=TLA.real, col="blue")
  return(2*min(mean(TLA.perm<TLA.real), mean(TLA.perm>TLA.real)))
}

fun3 <- function(dat){
  p.TLA<- TLA.pval(dat)
  return(p.TLA)
}


######################################
# Method 4: CNM-Full
######################################
fun4 <- function(dat){
  library("LiquidAssociation")
  # add 10^-3 to zero count
  dat[dat==0] <- 10^-3
  # log transformation
  dat <- log(dat)
  # standardize
  x1.cnm <- scale(dat[,1])
  x2.cnm <- scale(dat[,2])
  x3.cnm <- scale(dat[,3])
  x.cnm <- cbind(x1.cnm, x2.cnm, x3.cnm)
  colnames(x.cnm)<-c("Gene1", "Gene2", "Gene3")
  # p-value of b5
  p.b5 <- CNM.full(x.cnm)@output[8,4]
  return(p.b5)
}



######################################
# Wrapper Function
######################################
wrapper <- function(t00, t11){
  dt <- sim.Z(t00,t11)
  re1 <- fun1(dt)
  re2 <- fun2(dt)
  re3 <- fun3(dt)
  re4 <- fun4(dt)
  # Store Results from Four Methods
  f1<-paste("~/ZENCO.RData", sep="")
  save(re1, re2, re3, re4, file=f1)
}


# wrapper(tau0, tau1)
# example: wrapper(0, 0.05)
