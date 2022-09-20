go over this only after you read the intercept prior sampling. this adds some linear effects


```
library(brms)
library(bayestestR)
rm(list=ls())

###visually test some priors
#we want the effects on the mean rt to be around about 50-100ms and 25-50 from tau
# examine tau prior for the effect
priormean = log((0+150)/1000) - log(150/1000) 
priorsd   = 1

#get some quantiles in ms scale with the effect (b0+b1*1)
x=exp(rnorm(1000,-1.89712,1)+rnorm(1000,0,1))
quantile(x,c(.20,.40,.60,.80))*1000 


#compare to some quantiles in ms scale without the effect (b0+b1*0)
x=exp(rnorm(1000,-1.89712,1))
quantile(x,c(.20,.40,.60,.80))*1000 




####prior sampling

#load data 
#note that only the number of observations and the indep are going to have a role here. 
#the actual empirical rt will not have an influence since we will use prior sampling).
load('./df.rdata')
df=df_n

df$stage=df$stage - 1
df$group <- relevel(df$group, ref = 'control') # make sure the baseline group is control



#brms regressions argument
myformula=brmsformula(
  rt    ~ 0+Intercept+stage*group,
  sigma ~ 0+Intercept,
  beta  ~ 0+Intercept+stage*group
)

#brms family argument
myfamily=exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

#brms prior argument (you can examine brms defualts using get_prior(myformula, family      = myfamily,data=df))
myprior  = c(
  
  #meanrt 
  set_prior(prior="normal(0.5,0.1)", class="b", coef="Intercept"      ,dpar=""),
  set_prior(prior="normal(0  ,0.1)", class="b", coef="groupADHD"      ,dpar=""),
  set_prior(prior="normal(0  ,0.1)", class="b", coef="stage"          ,dpar=""),
  set_prior(prior="normal(0  ,0.1)", class="b", coef="stage:groupADHD",dpar=""),
  
  #sigma
  set_prior(prior="normal(-2.995732, 1.5)", class="b", coef="Intercept", dpar="sigma"),
  
  #tau
  set_prior(prior="normal(-1.89712,  1)", class="b", coef="Intercept"      ,dpar="beta"),
  set_prior(prior="normal(0       ,  1)", class="b"  , coef="groupADHD"      ,dpar="beta"),
  set_prior(prior="normal(0       ,  1)", class="b"  , coef="stage"          ,dpar="beta"),
  set_prior(prior="normal(0       ,  1)", class="b"  , coef="stage:groupADHD",dpar="beta")
  
)

#compile
model<-
  brm(myformula,
      family      = myfamily,
      prior       = myprior,
      sample_prior='only', 
      data        = df, 
      iter        = 1,  
      warmup      = 1,
      cores       = 1, 
      chains      = 1, 
      backend='cmdstan')

#sample
model=update(model,prior=myprior,newdata=df,iter=100000)

#check estimates
describe_posterior(model,ci=.80) #(try to back-transform your coeff to see they make sense in ms. you can also do that for the 95%CI)

#check prior predictive checks
library(ggplot2)
brms::pp_check(model,
               ndraws=20, #define how many samples to use. note that each sample generates a distrbution,
               prefix='ppd' #so brms wont plot the empirical data
)+coord_cartesian(xlim = c(-10, 10))

```
