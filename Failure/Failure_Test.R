# measurement alignment test with Fisher & Karl (2019)

##### 0. Package importing and dataset loading

# package load
library(lavaan)
library(sirt)
library(foreach)
library(parallel)
library(doParallel)
library(psych)
library(MASS)

# set some constants
fits <- c('rmsea.scaled','srmr','cfi.scaled')
var.help <- c('vaccine_attitudes_2_midneutral',
              'vaccine_attitudes_9_midneutral',
              'antiexpert_1','conspirational_think_1')

# data load
data <- read.csv('data_filtered.csv')
DATA <- data

# how many group?
table(data$UserLanguage)

# sort by group
data <- data[order(data$UserLanguage),]

# extract group names
countries <- labels(table(data$UserLanguage))[[1]]
# number of  groups to be examined
n.include <- length (countries)

##### 1. MG-CFA
#####
# first, let's start with MGCFA
# let's focus on help for now for alignment
cfa_model.help<- '
              help =~ vaccine_attitudes_9_midneutral + 
 vaccine_attitudes_2_midneutral + 
antiexpert_1+conspirational_think_1'

#  invariance test
# configural invariance
# use WLSMV (ordinal scale)
fit.help.configural <- cfa (cfa_model.help, data, group = 'UserLanguage',
                            estimator='WLSMV')
fitMeasures(fit.help.configural)[fits]
#rmsea.scaled         srmr   cfi.scaled 
#   0.12913445   0.04880355   0.88495049 

# metric invariance
fit.help.metric <- cfa (cfa_model.help, data, group = 'UserLanguage',
                            estimator='WLSMV', group.equal='loadings')
fitMeasures(fit.help.metric)[fits]
#rmsea.scaled         srmr   cfi.scaled 
#0.10157052   0.06747898   0.82587208   
# changes
fitMeasures(fit.help.metric)[fits]-fitMeasures(fit.help.configural)[fits]
# -0.02756393   0.01867543  -0.05907842 
# unacceptable

# scalar invariance
fit.help.scalar <- cfa (cfa_model.help, data, group = 'UserLanguage',
                        estimator='WLSMV', group.equal=c('loadings','intercepts'))
fitMeasures(fit.help.scalar)[fits]
#rmsea.scaled         srmr   cfi.scaled 
#0.2022793    0.1345629    0.0000000  
# changes
fitMeasures(fit.help.scalar)[fits]-fitMeasures(fit.help.metric)[fits]
# 0.1007088    0.0670839   -0.8258721
# unacceptable

# residual invariance 
fit.help.residual <- cfa (cfa_model.help, data, group = 'UserLanguage',
                        estimator='WLSMV', group.equal=c('loadings','intercepts',
                                                         'residuals'))
fitMeasures(fit.help.residual)[fits]
#rmsea.scaled         srmr   cfi.scaled 
#0.1816056    0.2168882    0.0000000 
fitMeasures(fit.help.residual)[fits]-fitMeasures(fit.help.scalar)[fits]
#rmsea.scaled         srmr   cfi.scaled 
#-0.02067372   0.08232531   0.00000000


##### 2. Measurement Alignment
#####
# Then, let's perform measurement alignment

# help
# extract cfa parameters
par.help <- invariance_alignment_cfa_config(dat = data[,var.help], 
                                          group = data$UserLanguage)
# do alignment
# following the suggested threshold values in Fisher & Karl (2019)
mod.help <- invariance.alignment(lambda = par.help$lambda, nu =
                                   par.help$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))

# test performance
mod.help$es.invariance['R2',]
#  loadings intercepts 
# 0.9344503  0.9888906 


# item-level test
cmod <- invariance_alignment_constraints(mod.help, lambda_parm_tol = .4, nu_parm_tol = .22)
summary(cmod)
# lambda noninvariance item = 55.4%
# nu noninvariance item = 42.9%
# acceptable

##### 3. Monte Carlo Simulations
#####
# Monte Carlo simulation

# define function for simulation/iteration
# simulation function
# times: how many time to repeat the same test?
# n: n for each test (e.g., 100, 200, 500)
# data: data to be tested
# model: CFA model
# lv: name of the latent variable e.g., help
# n.include: number of groups. in this case, 8 (countries)
# gruops: names of groups to be tested. e.g., data$country
# items: items (variable names) to be tested. an array
# par: parameters for cfa (created by invariance_alignment_cfa_config)
# seed: random seed
simulation <- function(n,data,model,lv,
                           n.include,groups,items,par,seed=1){

  # create a matrix to return results
  # correlation and R2s
  cor.mean <- 0
  cor.var <- 0
  R2.loading <- 0
  R2.intercept <- 0
  
  # begin simulation
  set.seed(seed)
  G <- n.include # number of groups
  I <- length(items) # number of items
    
  # lambda, nu, and error_var to be created for simulation
  err_var.cle <- matrix(1, nrow=G,ncol=I)
  
  # Create stimulated data
  # enter group mu and sigma
  data$Y <- rowMeans(data[,items])
  mu<-scale(aggregate(x=data$Y,
                        by = list(groups),
                        FUN=mean, na.rm=T)[,2])[,1]
  sigma <- (aggregate(x=data$Y,
                        by = list(groups),
                        FUN=sd, na.rm=T)[,2])
  N <- rep(n,G)
  
  # do simulation for data creation
  dat <- invariance_alignment_simulate(
      par$nu,par$lambda,err_var.cle,mu,sigma,N
    )
  
  # do CFA. parameter extraction for alignment
  par.simul <- invariance_alignment_cfa_config(dat = dat[,items], 
                                               group = dat$group,
                                               estimator = 'WLSMV')
    

  #cfa.test <-cfa(model,dat,estimator='WLSMV',group='group')
  #ipars <- parameterEstimates(cfa.test)
  
  # then, do alignment
  mod1.simul <- invariance.alignment(lambda = par.simul$lambda, nu =
                                       par.simul$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25),
                                     optimizer='nlminb')
    
  # do CFA. scalar invariance model for further calculation
  cfa.simul <- cfa(model,dat,estimator='WLSMV',group='group',
                     group.equal=c('loadings','intercepts'),meanstructure=T)
    
  # get group mean
  params.simul <- parameterEstimates(cfa.simul)
  alpha.simul <- params.simul[(params.simul$op=='~1')&(params.simul$lhs==lv),'est']
    
  # group mean correlation (Muthen 2018)
  correlation <- corr.test(alpha.simul,mod1.simul$pars$alpha0,method='spearman')$r
  
  # get group intercept
  psi.simul <- params.simul[(params.simul$op=='~~')&(params.simul$lhs==lv),'est']
  correlation.psi <- corr.test(psi.simul,mod1.simul$pars$psi0,method='spearman')$r
    
  cor.mean <- correlation
  cor.var <- correlation.psi
    
  # R2. The extent to which non-invariances in loadings and intercepts were absorbed?
  # ideally >= 75-80%
  R2.loading <- mod1.simul$es.invariance['R2',1]
  R2.intercept <- mod1.simul$es.invariance['R2',2]
  
  # make matrix
  to.return <-cbind(cor.mean,cor.var,R2.loading,R2.intercept)
  to.return <- data.matrix(to.return)
  
  return(to.return)
  
}



# use five cores
# and repeat each test 500 times
cores <- 5
times <- 10

# create threads to distribute tasks
cl <- parallel::makeCluster(cores,type='FORK')
doParallel::registerDoParallel(cl)

# start simulation with n = 100
# time measure as well

#start_100 <-Sys.time()
n <- 500 
now <- foreach (i = seq(1,times)) %dopar%{
  # do simulation
  #simulation <- function(n,data,model,lv,
  #                       n.include,groups,items,par,seed=1)
  simulation(n ,data,cfa_model.help,'help',
             n.include, data$UserLanguage,var.help,par.help,i)
  #  message(sprintf('%d',i))
}
#end_100<-Sys.time()
#elapsed_100 <- end_100 - start_100
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_100 <- now[[1]]
  }else{
    simulate_100 <- rbind(simulate_100,now[[i]])
  }
}
# save n = 100
write.csv(data.frame(simulate_100),file='simulate_100.csv',row.names = FALSE)

# results for n = 100
print(describe(simulate_100),digits=4)

# do the sam ewith n = 200
start_200 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # do simulation
  #simulation <- function(n,data,model,lv,
  #                       n.include,groups,items,par,seed=1)
  simulation(200,data,cfa_model.help,'help',
             n.include, data$UserLanguage,var.help,par.help,i)
  #  message(sprintf('%d',i))
}
end_200<-Sys.time()
elapsed_200 <- end_200 - start_200
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_200 <- now[[1]]
  }else{
    simulate_200 <- rbind(simulate_200,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_200),file='simulate_200.csv',row.names = FALSE)
print(describe(simulate_200),digits=4)

# then n= 500
start_500 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # do simulation
  #simulation <- function(n,data,model,lv,
  #                       n.include,groups,items,par,seed=1)
  simulation(500,data,cfa_model.help,'help',
             n.include, data$UserLanguage,var.help,par.help,i)
  #  message(sprintf('%d',i))
}
end_500<-Sys.time()
elapsed_500 <- end_500 - start_500
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_500 <- now[[1]]
  }else{
    simulate_500 <- rbind(simulate_500,now[[i]])
  }
}
# save n = 500
write.csv(data.frame(simulate_500),file='simulate_500.csv',row.names = FALSE)
print(describe(simulate_500),digits=4)

# terminate multiprocessing
parallel::stopCluster(cl)

# in all cases, cor ≥ 95%. Good
