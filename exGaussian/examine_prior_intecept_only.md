# examine priors (intercept only)
This code will walk you through prior testing for exGaussian sampling. The idea here is to find good priors to use when sampling exgaussian posteriors with rt data.

## manullay examine some priors
Lets start by visually and manullay examining some possible distrbutions for our priors.
This is important since we are dealing with a 'log' link meaning we need to set the prior on the log scale.
What we want to get here, is exgaussian distrbutions that will get us resnoable estimates for rts.
It would be nice to see the meanrt set aroun 500ms, sigma around 50ms, and tau around 150

#### examining our priors for sigma:
we will set a tau mean and sd (putting care into the log transformation)
and then visually examine the prior and use common sense. 
usually sigma estimates should get around 50ms values as a prior for an intercept
the last line of the next code will get you the 20th,40th,60th, and 8th percntile of the prior distrbution in ms
which should be very informative for us

```
priormean = log(50/1000) #this will be our prior mean for 50ms
priorsd   = 1.5 #this will be our tau prior sd - change it to see if you get what you need

x=seq(-4,4,0.001) #just to generate a prior plot in the next line
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) #this is the prior in log estimates
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 #lets see what will be the estimate in ms at the 20,40,60 and quantile of the prior.
```

#### examining our priors for tau (the same way as sigma)
usually tau estimates should get around 150ms values as a prior for an intercept
```
priormean = log(150/1000) #this will be our tau prior mean for 50ms
priorsd   = 1 #this will be our tau prior sd - change it to see if you get what you need

x=seq(-4,4,0.001) #just to generate a prior plot in the next line
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) #this is the prior in log estimates
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 #lets see what will be the tau in ms at the 20,40,60 and quantile of the prior.
```

## sampling data in brms usingthe defined priors

```
library(brms)
library(bayestestR)

#brms data argument (note this is a place holder since we will sample from the prior)
df=data.frame(rt=1)

#brms regressions argument
myformula=brmsformula(
  rt    ~ 0+Intercept,
  sigma ~ 0+Intercept,
  beta  ~ 0+Intercept
)

#brms family argument
myfamily=exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

myprior  = c(
  #intercept
  set_prior(prior="normal(0.5,       0.1)",       class="b", coef="Intercept", dpar=""),
  set_prior(prior="normal(-2.995732, 1.5)",  class="b", coef="Intercept", dpar="sigma"),
  set_prior(prior="normal(-1.89712,  1)",  class="b", coef="Intercept", dpar="beta")
)

#sample
model<-
  brm(myformula,
      myprior,
      family      = myfamily,
      data        = df, 
      sample_prior='only',
      iter        = 1,  
      warmup      = 1,
      cores       = 1, 
      chains      = 1, 
      backend='cmdstan')
```
now we can check we got estimates around what we asked for (meanrt ~0.5sec, sigma~0.05sec, tau~0.15sec). Remember bo_tau=exp(b0_tau) and same for sigma.

```
describe_posterior(model,ci=.80) #(try to back-transform your coeff to see they make sense in ms. you can also do that for the 95%CI)
```

we can plot some simulated RTs. each line in the next plot will be an exgaussian dist from a single prior draw.
so its ok if this 'dance' a bit, but should be centered overall around the ms position we want
```
#check prior predictive checks
brms::pp_check(model,
               ndraws=20, #define how many samples to use. note that each sample generates a distrbution,
               prefix='ppd' #so brms wont plot the empirical data
               )
```

# sandbox
you can try different priors without having to compile your model (just quicker)
###sandbox------
#now lets set some priors and examine our simulated rts
myprior  = c(
  #intercept
  set_prior(prior="normal(0.5,       0.1)",       class="b", coef="Intercept", dpar=""),
  set_prior(prior="normal(-2.995732, 1.5)",  class="b", coef="Intercept", dpar="sigma"),
  set_prior(prior="normal(-1.89712,  1)",  class="b", coef="Intercept", dpar="beta")
)

model=update(model,prior=myprior,newdata=df,iter=100000)
describe_posterior(model,ci=.80)
brms::pp_check(model,
               ndraws=20, #define how many samples to use. note that each sample generates a distrbution,
               prefix='ppd' #so brms wont plot the empirical data
               )

```
