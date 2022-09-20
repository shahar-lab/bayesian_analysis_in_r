go over this only after you read the intercept prior sampling. this adds some linear effects


```
library(brms)
library(bayestestR)
rm(list=ls())

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
  set_prior(prior="normal(-2.995732, 0.5)", class="b", coef="Intercept", dpar="sigma"),
  
  #tau
  set_prior(prior="normal(-1.89712,  0.5)", class="b", coef="Intercept"      ,dpar="beta"),
  set_prior(prior="normal(0       ,  0.5)", class="b"  , coef="groupADHD"      ,dpar="beta"),
  set_prior(prior="normal(0       ,  0.5)", class="b"  , coef="stage"          ,dpar="beta"),
  set_prior(prior="normal(0       ,  0.5)", class="b"  , coef="stage:groupADHD",dpar="beta")
  
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
