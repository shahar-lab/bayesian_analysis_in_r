# examine priors (intercept only)
This code will walk you through prior testing for exGaussian sampling. The idea here is to find good priors to use when sampling exgaussian posteriors with rt data.

## manually examine some priors
Lets start by visually and manullay examining some possible distrbutions for our priors.
This is important since we are dealing with a 'log' link meaning we need to set the prior on the log scale.
What we want to get here, is exgaussian distrbutions that will get us resnoable estimates for rts.
It would be nice to see the meanrt set aroun 500ms, sigma around 50ms, and tau around 300ms.

#### examining our priors for sigma:
We will set a  mean and sd (putting care into the log transformation)
and then visually examine the prior and use common sense. 
usually sigma estimates should get around 50ms values as a prior for an intercept
the last line of the next code here will get you the 20th,40th,60th, and 80th percntile of the prior distrbution in ms
which should be very informative for us

```
priormean = log(50/1000) 
priorsd   = 0.5

#plot in log scale
x=seq(-4,4,0.001) 
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) 

#get some quantiles in ms scale
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 
```

#### examining our priors for tau (the same way as sigma)
usually tau estimates should get around 300ms values as a prior for an intercept
```
priormean = log(300/1000) 
priorsd   = 0.5

#plot in log scale
x=seq(-4,4,0.001) 
curve(dnorm(x,mean=priormean,sd=priorsd), from=-10, to=10) 

#get some quantiles in ms scale
exp(qnorm(c(.20,.40,.60,.80),mean=priormean,sd=priorsd))*1000 
```

## sampling data in brms using the defined priors

```
#### sample rts using brms and only priors
library(brms)
library(bayestestR)
rm(list=ls())
#brms data argument (note this is a place holder since we will sample from the prior)
df=data.frame(rt=rep(1,10000))

#brms regressions argument
myformula=brmsformula(
  rt    ~ 0+Intercept,
  sigma ~ 0+Intercept,
  beta  ~ 0+Intercept
)

#brms family argument
myfamily=exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

#brms prior argument
myprior  = c(
  #intercept
  set_prior(prior="normal(0.5,       0.1)",       class="b", coef="Intercept", dpar=""),
  set_prior(prior="normal(-2.995732, 0.5)",  class="b", coef="Intercept", dpar="sigma"),
  set_prior(prior="normal(-1.203973,  0.5)",  class="b", coef="Intercept", dpar="beta")
)

#compile
model<-
  brm(myformula,
      prior=myprior,
      sample_prior='only', 
      family      = myfamily,
      data        = df, 
      sample_prior='only',
      iter        = 1,  
      warmup      = 1,
      cores       = 1, 
      chains      = 1, 
      backend='cmdstan')

#sample
model=update(model,prior=myprior,newdata=df,iter=100000)


```
now we can check we got estimates around what we asked for (meanrt about 0.5sec, sigma about 0.05sec, tau about 0.15sec). Remember bo_tau=exp(b0_tau) and same for sigma.

```
describe_posterior(model,ci=.80) #(try to back-transform your coeff to see they make sense in ms. you can also do that for the 95%CI)
```

we can plot some simulated RTs. each line in the next plot will be an exgaussian dist from a single prior draw.
so its ok if this 'dance' a bit, but should be centered overall around the ms position we want
```
#check prior predictive checks
library(ggplot2)
brms::pp_check(model,
               ndraws=20, #define how many samples to use. note that each sample generates a distrbution,
               prefix='ppd' #so brms wont plot the empirical data
               )+coord_cartesian(xlim = c(-5, 5))


# another why to check the simulated data
preds <- posterior_predict(model, newdata = data.frame(rt=1))

library(ggdist)
ggplot() + 
  stat_slab(aes(x = preds[,1], y = "")) +
  xlim(-1,2)

```
