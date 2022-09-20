## examine priors (intercept only)

we can start by initilazing some brms model. this is just to compile a model we can later change. not data is only a placeholder required by brms. we will only sample from the priors here.

```
library(brms)
library(bayestestR)
rm(list=ls())


####initial preperation and model compiliation
#load data but note this is technical - we are not going to actually use the data here since brms will work under "prior=only"
load('./df.rdata')
df=df_n

#brms regressions argument
myformula=brmsformula(
  rt    ~ 0+Intercept,
  sigma ~ 0+Intercept,
  beta  ~ 0+Intercept
)

#brms family argument
myfamily=exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

#sample
model<-
  brm(myformula,
      family      = myfamily,
      data        = df, 
      sample_prior='only',
      iter        = 1,  
      warmup      = 1,
      cores       = 1, 
      chains      = 1, 
      backend='cmdstan')

#we can look into stan if we need to
stancode(model)
```

now we can visullay examine some prior. this is important since we are dealing with a 'log' link meaning we need to set the prior on the log scale

```
# now we can change the prior and see what happens
# here we will set meanrt prior around 500ms, sigma around 50ms, tau around 150ms
# first lets just visuallt and manully examine the priors for sigma and tau
# dealing with the fact that we need to set the priors after the log transformation

# examining our priors for sigma:
# we will set a tau mean and sd (putting care into the log transformation)
# and then visually examine the prior and use common sense. 
# usually sigma estimates should get around 50ms values as a prior for an intercept

priormean = log(50/1000) #this will be our tau prior mean for 50ms
priorsd   = 1.5 #this will be our tau prior sd - change it to see if you get what you need

x=seq(-4,4,0.001) #just to generate a prior plot in the next line
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) #this is the prior in log estimates
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 #lets see what will be the tau in ms at the 20,40,60 and quantile of the prior.


#examining our priors for tau (tahe same way as sigma)
# usually tau estimates should get around 150ms values as a prior for an intercept

priormean = log(150/1000) #this will be our tau prior mean for 50ms
priorsd   = 1 #this will be our tau prior sd - change it to see if you get what you need

x=seq(-4,4,0.001) #just to generate a prior plot in the next line
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) #this is the prior in log estimates
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 #lets see what will be the tau in ms at the 20,40,60 and quantile of the prior.
```

now we can test these priors by actually sampling rt in brms from these priors

```
###sandbox------
#now lets set some priors and examine our simulated rts
myprior  = c(
  #intercept
  set_prior(prior="normal(0.5,       0.1)",       class="b", coef="Intercept", dpar=""),
  set_prior(prior="normal(-2.995732, 1.5)",  class="b", coef="Intercept", dpar="sigma"),
  set_prior(prior="normal(-1.89712,  1)",  class="b", coef="Intercept", dpar="beta")
)

model=update(model,prior=myprior,newdata=df,iter=100000)

#check you got estimates around what we asked for (meanrt ~0.5sec, sigma~0.05sec, tau~0.15sec). Remember bo_tau=exp(b0_tau) and same for sigma.
describe_posterior(model,ci=.80) #(try to back-transform your coeff to see they make sense in ms. you can also do that for the 95%CI)

#check prior predictive checks
brms::pp_check(model,
               ndraws=20, #define how many samples to use. note that each sample generates a distrbution,
               prefix='ppd' #so brms wont plot the empirical data
               )

```
