# bayesian_analysis_in_r
This includes some typical R code that we use in the lab.

Find here guidelines for how to describe your analysis in a paper:

<a href="https://www.nature.com/articles/s41562-021-01177-7/tables/1">Bayesian Analysis Reporting Guidelines</a>


Find here guidelines for how to describe your effect:

<a href="https://www.frontiersin.org/articles/10.3389/fpsyg.2019.02767/full">Indices of Effect Existence and Significance in the Bayesian Framework</a>

## BRMS package

brms is a package which wraps stan code for regression analysis in the common base R format of y~x like in lm and glm functions.

### Installation
To use brms you need to install both brms and rstan (stan for R). We start by deleting any disturbing files:
```
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
```
We then install rstan
```
install.packages("rstan")
```
We verify the installation was succesfull by running an example model using data which is within the package:
```
example(stan_model, package = "rstan", run.dontrun = TRUE)
```
We now install brms package:
```
install.packages("brms")
```
After installing rstan and brms you are good to go, but may still get extra speed by using "cmdstanr" package. 
This package enables your model to run without some unneccessary waste of time as is described <a href="http://mc-stan.org/cmdstanr/articles/cmdstanr.html#introduction-1">here</a> 

**Steps to install cmdstanr:**

1)Restart R session

2)Install cmdstan R package (check the link above for updates from the time of this being written)
```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```
3)Call the library
```
library(cmdstanr)
```
4)Check your computer is ready for this!
```
check_cmdstan_toolchain()
```
5)Install cmdstan 
```
install_cmdstan(cores = 2)
```
That should be enough but check the outputs in the console.
You can now use the backend = "cmdstanr" argument when calling a brms or a stan fitting function.
Bare in mind that you won't have access to marginal likelihoods if you want to calculate Bayesfactors.
### logistic regression
When having two options for the outcome variables we use family = bernoulli(link = "logit").

When more than two options are possible for the y variable, we use the family = binomial(link = "logit").

The model in equations may look like that:

![image](https://user-images.githubusercontent.com/51457131/160368497-95b8db60-a7fe-4741-8e8f-268b93a9b905.png)

Always define weakly informative priors! a uniform prior on the log-odds scale will be a bi-modal distribution focused on 0 and 1 on the probability scale.
![image](https://user-images.githubusercontent.com/51457131/165359009-ec0c43f4-4930-4e9b-aeb9-637af3825fe1.png)

```
library(brms)
mypriors=c(set_prior(prior = "normal(0,0.2)", class = "b",coef = "Intercept"),
           set_prior(prior = "normal(0,0.2)", class = "b",coef = "reward_oneback"))

 
model= brm(stay ~ 0 + Intercept+reward_oneback+(1+reward_oneback| subject), 
           data = df%>%filter(reoffer_ch==T), 
           family = bernoulli(link = "logit"),
           warmup = 1000,
           iter = 2000,    
           cores =4,
           chains=4,
           prior=mypriors)
save(model,file='./data/empirical_data/brms_weakly_informative_priors.rdata')
#use the update function to run a similar model without having it to recompile the model. You may use it when giving newdata or different priors as can be seen here <a href="https://rdrr.io/cran/brms/man/update.brmsfit.html">here</a>
fit_only_intercept = update(model,formula. = -reward_oneback)
```
### exmine  mcmc
```
plot(model)
pairs(model)
```

### examine posterior distrbution - statistics
```
library(bayestestR)
# posterior estimates
describe_posterior(model)
#define a range for the rope (region of practical equivalence) and the function will return the % of the posterior inside that region.
#read this article to learn more about the ROPE and the use of bayesian regressions:
#Indices of Effect Existence and Significance in the Bayesian Framework (Makowski et al., 2019)
describe_posterior(model,rope_range=c(-0.1,0.1)) 
#bayesfactor vs null point
mybayes=bayesfactor_parameters(model)
mybayes
plot(mybayes)

#bayesfactor vs ROPE
mybayes=bayesfactor_parameters(model,null=c(-0.0132, 0.0132))
plot(mybayes)
```

### examine posterior distrbution - visualization 
```
#plot fixed effect (e.g., reward_oneback)
plot(conditional_effects(model,effect='reward_oneback'), plot = FALSE)[[1]] + 
  ggtitle("(A) empirical data") + 
  xlab("previous-outcome") +
  ylab("pStay") + theme_bw()+
  scale_x_discrete(labels=c("Unrewarded","Rewarded"))

#plot historgram for posterior samples
library(RLR)
posterior_samples = insight::get_parameters(model)
my_posteriorplot(x       = posterior_samples$b_reward_oneback,
                 myxlim  = c(-.75,+.25),
                 my_vline= 0, 
                 myxlab  = expression(beta['previous-outcome']),
                 mycolor = "blanchedalmond")+ylab('density')
```

### prior and posterior predictive checks

```
#set empirical y estimates (should be numeric)
y   =as.numeric(df$stay)

#set panel factor (e.g., unrewarded vs rewarded - should be a factor var)
group_vec=df$reward_oneback

#posterior_predictive_check
yrep=posterior_predict(model)
ppc_stat_grouped(y, yrep, group_vec$reward_oneback)
```
### model comparision
you might want to comapre which model fits better to your data. 
e.g., you might have several brms models, and would like to find which set of predictors is most sutiable. 
check this video to understand what loo does: https://youtu.be/Re-2yVd0Mqk?t=315
using loo:
```
library(loo)
#run loo
model2 <- loo(model2)
model1 <- loo(model1)
model0 <- loo(model0)

#compare with loo
loo_compare(model2,model1,model0)

#stacking weights
stack_models <- cbind(model2$pointwise[,"elpd_loo"],
                      model1$pointwise[,"elpd_loo"],
                      model0$pointwise[,"elpd_loo"])
stacking_weights(stack_models)
```

using bridge_sampler:
```
#compare with bayesfactor_models
BF_model_comparison <- bayesfactor_models(model1,model2, denominator =model0)
BF_model_comparison
BF_model_inclusion =bayesfactor_inclusion(BF_model_comparison_unch)

```

### contrasts
after fitting you model you might want to run some specific contrasts, which can be easily computed using emmeans.
e.g., say we predict 'stay' with 'reward_oneback' and 'condition'. let's also assume that 'condition' is a factor variable with three levels A,B and C.
then brms will produce two dummy variabels for 'condition', but you might want to be more specific. You might for example want to check whether the reward effect in C is different from both A and B, or check seperately whether the reward effect is different for A vs B, A vs C and B vs C with out refitting the model.

```
#brms fit
model= brm(stay ~ reward_oneback*condition+(reward_oneback*condition| subject), 
           data = df, 
           family = bernoulli(link = "logit"))

#create and plot an emmeans object
library(emmeans)
em=emmeans::emmeans(model,~reward_n1back*cue_n1back)
plot(em)

#produce your own comparisions by given weights to each of the six levels (2 for reward_oneback X 3 for condition)
cont=
emmeans::contrast(em, list('A   vs. B'=c(-1, 1, 1,-1, 0, 0),
                           'A   vs. C'=c(-1, 1, 0, 0, 1,-1),
                           'B   vs. C'=c( 0, 0,-1, 1, 1,-1),
                           'A+B vs. C'=c(-1, 1,-1, 1, 2,-2)))
```
### power analysis
Unlike in NHST (null-hypothesis-significance-testing), Bayesian analysis focuses on parameter estimation. Thus, we are not focusing on the chance of obtaining a specific point-estimate (e.g point-null hypothesis) but are satisfied with **precisly** estimating our parameter of interest. Our "criterion of interest" will therefore be whether or not our Credible Interval (CI) is narrow enough for us to be certain about our results. 

The following code simulates a power analysis for a simple linear regression model having only an intercept (it takes a while to run).
```
library(brms)
library(cmdstanr)
library(bayestestR)
library(dplyr)
#Power analysis
possible_Nsubjects = seq(25, 100, by = 25)
Ntrials = 50
Nsimulations = 50
# Define prior model
df    = data.frame(y=0) # data has just the minimal number of rows to let brm know what the various
# predictors look like (levels, nesting...)

#priors are based on preliminary findings or on "ideal" results
prior = c(set_prior(
  prior = "normal(0,0.5)",class = "b",coef = "Intercept"),
  set_prior(
  prior = "cauchy(0,0.5)",class = "sigma"))

#data generating model using the priors
model = brm(y~0 + Intercept,prior=prior,sample_prior = 'only',data=df,backend = "cmdstanr")
#define df for results of criterion
power_df_all=data.frame()
#sanity check for hdi
hdi(model)
hdi(rnorm(100000,0,0.5))

# Function that makes x data
generate_x = function(Nsubjects, Ntrials) {
    df = data.frame(
      subject = rep(seq(1, Nsubjects), each = Ntrials),
      trial = rep(seq(1, Ntrials), Nsubjects)
    )
    return(df)
  }

for (Nsubjects in possible_Nsubjects) {
  power_df=data.frame(sample_size =rep(Nsubjects,Nsimulations),CI_width=rep(NA,Nsimulations))
  #we create x data for each sample size using this function
  x_data = generate_x(Nsubjects = Nsubjects, Ntrials=Ntrials)
  
  #every draw from the posterior distribution is like a new study
  #the posterior_predictive_distribution (PPD) will have the structure of Ndraws X Nobservations.
  #thus, the number of rows refers to the number of studies on which we will loop.
  PPD = posterior_predict(model, newdata = x_data, ndraw = Nsimulations) 
  
  for (i in 1:nrow(PPD)) {
    # make full data by combining the x data and y data
    temp_data = x_data
    
    #adding each iteration the whole row as the y variable to the x data
    temp_data$y <- PPD[i, ]
    
    # estimate sub model for the current draw from the posterior
    # use priors that you put in the prereg - NOT the ones for the gen-model.
    fit_prior = c(set_prior(
      prior = "normal(0,0.2)",class = "b",coef = "Intercept"),
      set_prior(
        prior = "cauchy(0,0.2)",class = "sigma"))
    
    temp_model <- brm(y~0 + Intercept, data = temp_data, 
                      sample_prior = FALSE,
                      prior = fit_prior, 
                      backend = "cmdstanr")
    
    # get HDI / estimate / ROPE according to your criterion
    CI_width=hdi(temp_model,parameters="Intercept")[[4]] - hdi(temp_model,parameters="Intercept")[[3]]
    power_df[i,"CI_width"] = CI_width
  }
  power_df_all=rbind(power_df_all,power_df)
}

#check power criterion
precision_criterion = 0.15 #the minimal CI width we want
power_df_all= power_df_all%>%mutate(criterion=(CI_width<precision_criterion)*1)
power = power_df_all%>%group_by(sample_size)%>%summarise(power=mean(criterion))
power_df_all%>%group_by(sample_size)%>%summarise(mean(CI_width))
```
### logistic regression power analysis
More commonly, you will run a logistic regression which means your y variable is dichotomous and thus the link function between the predictor and the y variable is a binomial or a bernoulli function and not a gaussian. The output variable will correspond to the logit or log-odds and can then be transformed to probability.

Here is an example for power analysis with a logistic regression:
```
library(brms)
library(cmdstanr)
library(bayestestR)
library(dplyr)
library(RLR)
#Power analysis
possible_Nsubjects = seq(25, 100, by = 25)
Ntrials = 125
Nsimulations = 50
# Define prior model
# data has just the minimal number of rows to let brm know what the various
# predictors look like (levels, nesting...)
df    = data.frame(subject=rep(1,2),generalize=rep(1,2),reward_oneback=as.factor(c(0,1)))
#priors are based on preliminary findings or on "ideal" results
prior = c(set_prior(
  prior = "normal(0.1,0.15)",class = "b",coef = "Intercept"),
  set_prior(
    prior = "normal(0.3,0.2)",class = "b",coef = "reward_oneback"))
#data generating model using the priors
model = brm(generalize~0 + Intercept+reward_oneback,
            prior=prior,
            sample_prior = 'only',
            data=df,
            family = bernoulli(link="logit"),
            backend = "cmdstanr")

#sanity check for hdi
hdi(model)
hdi(rnorm(100000,0.3,0.2)) #put here the mean and sd of your prior
# Function that makes x data
generate_x = function(Nsubjects, Ntrials) {
  df = data.frame(
    subject = rep(seq(1, Nsubjects), each = Ntrials),
    trial = rep(seq(1, Ntrials), Nsubjects),
    reward_oneback = as.factor(sample(c(0,1),Nsubjects*Ntrials,replace = TRUE))
  )
  return(df)
}
#define df for results of criterion
power_df_all=data.frame()
for (Nsubjects in possible_Nsubjects) {
  power_df=data.frame(sample_size =rep(Nsubjects,Nsimulations),Median=rep(NA,Nsimulations),CI_width=rep(NA,Nsimulations),zero_criterion=rep(NA,Nsimulations))
  #we create x data for each sample size using this function
  x_data = generate_x(Nsubjects = Nsubjects, Ntrials=Ntrials)
  #every draw from the posterior distribution is like a new study
  #the posterior_predictive_distribution (PPD) will have the structure of Ndraws X Nobservations.
  #thus, the number of rows refers to the number of studies on which we will loop.
  PPD = posterior_predict(model, newdata = x_data, ndraw = Nsimulations)
  for (i in 1:nrow(PPD)) {
    # make full data by combining the x data and y data
    temp_data = x_data
    #adding each iteration the whole row as the y variable to the x data
    temp_data$generalize <- PPD[i, ]
    # estimate sub model for the current draw from the posterior
    # use priors that you put in the prereg - NOT the ones for the gen-model.
    fit_prior = c(
      set_prior(
        prior = "normal(0,0.2)",
        class = "b",
        coef = "Intercept"
      ),
      set_prior(
        prior = "normal(0,0.2)",
        class = "b",
        coef = "reward_oneback1"),
      )
    temp_model <- brm(generalize ~ 0+Intercept+reward_oneback,
                      data = temp_data,
                      sample_prior = FALSE,
                      prior = fit_prior,
                      family = bernoulli(link="logit"),
                      backend = "cmdstanr")
    # get HDI / estimate / ROPE according to your criterion
    Median=fixef(temp_model,pars="reward_oneback1")[[1]]
    power_df[i,"Median"] = Median
    CI_width=hdi(temp_model,parameters="reward_oneback1")[[4]] - hdi(temp_model,parameters="reward_oneback1")[[3]]
    power_df[i,"CI_width"] = CI_width
    #check if zero is in the hdi by multiplying the two edges of the hdi
    zero_criterion=(hdi(temp_model,parameters="reward_oneback1")[[4]] * hdi(temp_model,parameters="reward_oneback1")[[3]])>0
    power_df[i,"zero_criterion"] = zero_criterion
  }
  power_df_all=rbind(power_df_all,power_df)
}
#check power criterion
precision_criterion = 0.2 #the minimal CI width we want
power_df_all= power_df_all%>%mutate(criterion=(CI_width<precision_criterion)*1)
save(power_df_all,file="power_df_all.Rdata")
power = power_df_all%>%group_by(sample_size)%>%summarise(power=mean(criterion))
CI_width_by_sample_size = power_df_all%>%group_by(sample_size)%>%summarise(power=mean(CI_width))
```

### reliability
Check <a href="https://discourse.mc-stan.org/t/computing-icc-like-reliability-within-unit-variance-for-hierarchical-models/20729">this</a> thread for an explanation

There are two versions for the same analysis, one where we manually divide the data and one where we add session as a factor so that brms automatically divides it.

Notice also that when we create the random factors in two separate parantheses we assume that there is no correlation between the two.

Example:

```
library(dplyr)
library(brms) 
rm(list=ls())

#get data
df=read.csv('./data/empirical_data/df.csv')
#add variables

df=df%>%mutate(coherence=if_else(coherence==1,1,-1),reveal=if_else(reveal==1,1,-1),session=as.factor(session),
               intercept1=if_else(session=="session1",1,0),intercept2=if_else(session=="session2",1,0),
               reveal1=if_else(session=="session1",reveal,0),reveal2=if_else(session=="session2",reveal,0))

# df=df%>%mutate(coherence=if_else(coherence==1,1,-1),reveal=if_else(reveal==1,0.5,-0.5),session=as.factor(session))

#formulas
formula_reliability =  coherence ~ 0+intercept1+intercept2+reveal1+reveal2+(0+intercept1+intercept2+reveal1+reveal2|subject)

#model
model_reliability =
  brm(
    formula = formula_reliability,
    data = df,
    family = bernoulli(link = "logit"),
    warmup = 1000,
    iter = 1200,
    chains = 20,
    cores = 20,
    seed = 123,
    backend ="cmdstanr"
  )
save(model_reliability, file = 'data/brms/model_reliability.Rdata')
```
