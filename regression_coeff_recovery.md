#coeff recovery in regression analysis using brms

Here we show an example for a simulation where we generate data based on known 'true' coeff and recover thoese using brms.
This simulation allows us to explore how sample-size and effect-size influence our uncertinty estimates (aka CI).
Its also a first step before going to a power analysis.
For each study:
sample individual parameters based on known coeff and their sds.
simulate data based on the individual parameters
estimate the parameters based on the "observed" data

We can repeat this for a few studies, and then see what size of CI/pd we observe, and whether our point-estimate is close to the real values used to sample individual parameters.

```
#Aim: get the results of N studies trying to estimate a true effect (intercept and slope).
#the fixed and scale parameters are used to generate individual parameters for a specific study
#then this study estimates the fixed effect while also including CI and pd estimates.


rm(list=ls())
library(brms)
library(cmdstanr)
library(bayestestR)
library(dplyr)
library(lme4)
library(lmerTest)

####some preparations with individual parameters and data simulation----------------
#simulation configuration
Nsubjects    = 50
Ntrials      = 100
Nstudies     = 10


#function for data generation
generate_df = function(Nsubjects, Ntrials, eta) {
  df=data.frame()
  
  for (i in 1:Nsubjects){
    subject        = rep(i, each = Ntrials)
    trial          = seq(1, Ntrials)
    x0             = rep(1, each = Ntrials) #for the intercept
    x1             = sample.int(2,Ntrials,replace=T)-1 #condition coded 0 or 1
    beta_0         = rep(eta[i,1],Ntrials)
    beta_1         = rep(eta[i,2],Ntrials)
    y              = cbind(x0,x1)%*%eta[i,]+rnorm(Ntrials,0,1)
    
    df = rbind(df,data.frame(subject,trial,beta_0,beta_1,x1,y))
  }
  return(df)
}


#that's the model we'll be using
fit_prior = c(
  set_prior(
    prior = "normal(0,0.2)",
    class = "b",
    coef = "Intercept"
  ),
  set_prior(
    prior = "normal(0,0.2)",
    class = "b",
    coef = "x1")
)

model_empty <- brm(y ~ 0+Intercept+x1 +(1+x1 |subject),
                   data   = data.frame(subject=1,
                                       x1=1,
                                       y=1),
                   prior  = fit_prior,
                   backend= "cmdstanr",
                   chains = 1,
                   silent = 0)

#initialize data frames for simulation output
df_bayes=df_freq=list()
df_bayes[['x0']]=df_bayes[['x1']]=data.frame()
df_freq [['x0']]=df_freq [['x1']]=data.frame()



####simulate studies-----------------------

for (study in 1: Nstudies){
  #a current study begins here
  
  #we sample some participant from the population. the means are the 'fixed effects' and the sds are the 'random effect'.
  #note that this also has a covariance parameter for the association between the intercept and slope in the population
  mysigma = matrix(c(1.0,0.2,
                     0.2,1.0),
                   2,2) #varcov matrix where diag is the variance for each param, and off-diag is the covar
  
  eta    = MASS::mvrnorm(n     = Nsubjects,
                         mu    = c(0,.20),  #these are the fixed effects for the intercept and slope,
                         Sigma = mysigma       #variance-covariance matrix
                         )
  
  
  #then we 'run' the study and observe the data
  df=generate_df(Nsubjects,Ntrials,eta)
  
  #Bayes estimation
  model <- update(model_empty,
                  newdata   = df,
                     prior  = fit_prior,
                     backend= "cmdstanr",
                     iter   =1100,
                     warmup =1000,
                     chains = 10,
                     silent = 0,
                     refresh = 0,
                     cores=10)
  
  df_bayes[['x0']]=rbind(df_bayes[['x0']],describe_posterior(model,parameters='b_Intercept'))
  df_bayes[['x1']]=rbind(df_bayes[['x1']],describe_posterior(model,parameters= 'b_x1'))

  
  #freq estimation
  model = lmer(y~1+x1+(1+x1|subject),
               data=df) 
  df_freq[['x0']]=rbind(df_freq[['x0']],summary(model)$coefficients[1,])
  df_freq[['x1']]   =rbind(df_freq[['x1']],summary(model)$coefficients[2,])
}

```
