## Example 1: Intercept only model
Parameter recovery to a single set of exGaussian data. 
This includes three examples: 
(1) estimate parameters using ML in gamlss
(2) estimate posterior using brms 
(3) estimate posterior using brms with a 'log' link for tau and sigma (this can help mcmc sampling but requires using exp() to get the natural scale of milliseconds)

```
####simulate data from true mu=400,sigma=50,tau=150
N   =5000
rt  =gamlss::rexGAUS(N, mu = 400, sigma = 50, nu = 150)
df  =data.frame(rt)

####fit with ML
model<-gamlss(rt~1, data=df,
              family=exGAUS(mu.link = "identity", 
                            sigma.link = "identity", 
                            nu.link = "identity"))
summary(model)




####fit with brms 
model<-brm( 
  
  brmsformula(
    rt    ~ 1,
    sigma ~ 1,
    beta  ~ 1
  ), 
  
  data = df,
  warmup = 500,
  iter = 1000,  
  cores =1, 
  chains=1,
  backend='cmdstan',
  
  family = exgaussian(link = "identity", 
                      link_sigma = "identity", 
                      link_beta = "identity"),
  
  prior=c(prior(student_t(3, 523.1, 125.1),class=Intercept, dpar=""),#you can plot this using plot(dstudent_t(seq(0,1000,1),df=3,mu=523.1,sigma=125.1))
          prior(student_t(3, 100, 50),     class=Intercept, dpar="beta"),
          prior(student_t(3, 0, 2.5),      class=Intercept, dpar="sigma"))
)
  

library(bayestestR)
describe_posterior(model)




####fit with brms using log link (and then you need to convert the parameters back using exp())
model<-brm( 
  brmsformula(
    rt    ~ 1,
    sigma ~ 1,
    beta  ~ 1
  ), 
  data = df,warmup = 500,iter = 1000,  cores =1, chains=1,
  family = exgaussian(link = "identity", 
                      link_sigma = "log", 
                      link_beta = "log"),
  backend='cmdstan')

samples = insight::get_parameters(model)



#examine your recovered parameters (note meanRT=mu+tau)

samples = samples%>%mutate(RTmean=b_Intercept,
                           mu  =b_Intercept-exp(b_beta_Intercept),
                           tau =exp(b_beta_Intercept))


print(paste('RTmean_pred =', median(samples$RTmean))) #you can use hdi(samples$b_Intercept) to get CI
print(paste('mu_pred ='    , median(samples$mu)))
print(paste('tau_pred ='   , median(samples$tau)))
```
library(brms)
library(bayestestR)
library(dplyr)
library(insight)

#recover exGaussian parameters mu and tau for a single data set (no indepndent variables - only intercept)
N   =5000
rt  =rexGAUS(N, mu = 400, sigma = 50, nu = 150)
df  =data.frame(rt)

#fit with brms and obtain the posterior samples
model<-brm( 
  brmsformula(
  rt    ~ 1,
  sigma ~ 1,
  beta  ~ 1
), 
data = df,warmup = 500,iter = 1000,  cores =1, chains=1,
family = exgaussian(),
backend='cmdstan')

samples = insight::get_parameters(model)

#examine your recovered paramters 
describe_posterior(model)

samples = samples%>%mutate(RTmean=b_Intercept,
                           mu  =b_Intercept-exp(b_beta_Intercept),
                           tau =exp(b_beta_Intercept))


print(paste('RTmean_pred =', median(samples$RTmean))) #you can use hdi(samples$b_Intercept) to get CI
print(paste('mu_pred ='    , median(samples$mu)))
print(paste('tau_pred ='   , median(samples$tau)))
```


## Exmaple 2: One indepndent variable
```
#### single independent categorical variable -----

#simulate data from true mu {b0=400,b1=100},sigma=50,tau={b0=150,b1=50}
N   =5000
x1  =rbinom(N,size=1,prob=.5) #or if you want this to be continuous use something likernorm(N,3,1)
rt  =rexGAUS(N, mu = 400+100*x1, sigma = 50, nu = 150+50*x1)
df  =data.frame(x1,rt)



#fit with brms and obtain the posterior samples
model<-brm( 
  brmsformula(
    rt    ~ 1+x1,
    sigma ~ 1,
    beta  ~ 1+x1
  ), 
  data = df,warmup = 500,iter = 1000,  cores =1, chains=1,
  family = exgaussian(),
  backend='cmdstan')

samples = insight::get_parameters(model)



#recover regression parameters

samples = samples%>%mutate(#first calculate the intercept parameter b0 for RTmean, mu and tau
                           RTmean_b0=b_Intercept,
                           mu_b0    =b_Intercept-exp(b_beta_Intercept),
                           tau_b0   =exp(b_beta_Intercept),
                           
                           #now calculate the slope parameter b1 for RTmean, mu and tau
                           RTmean_b1=b_Intercept+b_x1 - RTmean_b0,
                           mu_b1    =b_Intercept+b_x1-exp(b_beta_Intercept+b_beta_x1)-mu_b0,
                           tau_b1   =exp(b_beta_Intercept+b_beta_x1)-tau_b0
                           )


print(paste('mu_b0 (intercept)=',median(samples$mu_b0))) #you can use hdi(samples$...) to get CI
print(paste('mu_b1 (slope)='    ,median(samples$mu_b1)))

print(paste('tau_b0 (intercept)=',median(samples$tau_b0)))
print(paste('tau_b1 (slope)='    ,median(samples$tau_b1)))



#examine center estimates for mu, and tau (e.g., in case you want to plot this) 
print(paste('RTmean (x1 is 0) =', median(samples$RTmean_b0),
            'RTmean (x1 is 1) =', median(samples$RTmean_b1+samples$RTmean_b0))) 
print(paste('mu     (x1 is 0) =', median(samples$mu_b0),
            'mu     (x1 is 1) =', median(samples$mu_b1+samples$mu_b0))) 
print(paste('tau    (x1 is 0) =', median(samples$tau_b0),
            'tau    (x1 is 1) =', median(samples$tau_b1+samples$tau_b0))) 
```


OLD:
פרשנות הפרמטרים



שלב אחד הוא מחשב את הפרמטרים לפי הרגרסיה
mu =inter+beta_1*x
beta=inter+beta_1*x

שלב שני הוא מחשב את הmu לפי
mu = mu - exp(beta)

שלב שלישי הוא מחשב את הטאו לפי
tau = exp(beta)

הערה - פונקציית הלייקליהוד מקבלת למדא ולא טאו, ולכן לה יש גם המרה של Inverse שאותנו לא מעניינת


קובץ לשחזור פרמטרים באקסגאוסיין

```
library(gamlss)
library(brms)
library(bayestestR)
library(see)
library(dplyr)
theme_set(theme_modern())
N=5000

#simulate exGaussian data
x1=rexGAUS(N/2, mu = 400, sigma = 50, nu = 150)
x2=rexGAUS(N/2, mu = 500, sigma = 50, nu = 175)
df=data.frame(
condition=c(rep(0,N/2),rep(1,N/2)),
rt       =c(x1,x2)
)


df%>%group_by(condition)%>%summarise(mean(rt))
#fit with brms
model<-brm( 
  brmsformula(
  rt    ~ 1+condition,
  sigma ~ 1,
  beta  ~ 1+condition
), 
data = df,
warmup = 500,
iter = 1000,    
cores =1,
chains=1,
family = exgaussian(),
backend='cmdstan')


describe_posterior(model)
posterior_samples = insight::get_parameters(model)
names(posterior_samples)

posterior_samples = posterior_samples%>%mutate(mu_1  =b_Intercept-exp(b_beta_Intercept),
                                               tau_1 =exp(b_beta_Intercept),
                                               mu_2  =b_Intercept+b_condition-exp(b_beta_Intercept+b_beta_condition),
                                               tau_2 =exp(b_beta_Intercept+b_beta_condition)
                                               )
hist(posterior_samples$mu_1)
hist(posterior_samples$tau_1)
hist(posterior_samples$mu_2)
hist(posterior_samples$tau_2)

```
