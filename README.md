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

