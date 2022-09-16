## Example 1: Parameter recovery to a single set of exGaussian data. Note that for clarity this is generated without noise.
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
