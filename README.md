# bayesian_analysis_in_r
This includes some typical R code that we use in the lab

## BRMS package

### logistic regression
```
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
#use the update function to run a similar model.
fit_only_intercept = update(model,formula. = -reward_oneback)
```
### exmine  mcmc
```
mcmc_trace(model)
mcmc_pairs(model)
```

### examine posterior distrbution - statistics
```
# posterior estimates
describe_posterior(model)

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
```
library(brms)
library(cmdstanr)
library(bayestestR)
library(dplyr)
#Power analysis
possible_Nsubjects = seq(10, 100, by = 10)
Ntrials = 100
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
power_df=data.frame(simul=seq(1:Nsimulations),CI_width=rep(NA,Nsimulations))
#sanity check for hdi
hdi(model)
hdi(rnorm(100000,0,0.5))

# Function that makes x data
generate_x = function(Nsubjects, Ntrials) {
    df = data.frame(
      subject = rep(seq(1, Nsubjects), each = Ntrials),
      trial = rep(seq(1, Ntrials), Nsubjects),
      RT = rnorm(
        Nsubjects * Ntrials,
        400,
        50
      )
    )
    return(df)
  }

for (Nsubjects in possible_Nsubjects) {
  #we create x data for each sample size using this function
  x_data = generate_x(Nsubjects = Nsubjects, Ntrials=Ntrials)
  
  #every draw from the posterior distribution is like a new study
  #the posterior_predictive_distribution (PPD) will have the structure of Ndraws X Nobservations.
  #thus, the number of rows refers to the number of studies on which we will loop.
  PPD = posterior_predict(model, newdata = x_data, ndraw = Nsimulations, replace = TRUE) 
  
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
}
```
