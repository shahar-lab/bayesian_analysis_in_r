```

#### intercept only model -----

#simulate data from true mu=400,sigma=50,tau=150
N   =5000
rt  =gamlss::rexGAUS(N, mu = 400, sigma = 50, nu = 150)
df  =data.frame(rt)

#fit with ML
model<-gamlss(rt~1, data=df,
              family=exGAUS(mu.link = "identity", 
                            sigma.link = "identity", 
                            nu.link = "identity"))
summary(model)


#fit with brms 
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

#fit with brms using log link (and then you need to convert the parameters back using exp())

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
