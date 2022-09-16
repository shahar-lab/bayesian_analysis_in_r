## Example 1: 
We will simulate data with three indepndent varabels (two continuous  and one categorial). 
Then we will generate a depndent using some fixed regression parameters.
Finaly we will use both frequentist and bayesian modeling to recover thoese regression parameters from the data

```
#Aim:Generate and recover linear regression parameters using freq and bayes models

####data with two continuous  and one categorical variables
N=100 
x1   =rnorm(N,mean=50, sd=9)
x2   =rnorm(N,mean=200,sd=64)
x3   =rbinom(N,size=1, prob=0.7)
error=rnorm(100,0,16)

#Generate the dependent variable (b0=150, b1=-4, b2=2.5, b3=5)
y    =150-(4*x1)+(2.5*x2)+(5*x3)+error
df   =data.frame(x1,x2,x3,y)

#Frequentist model
model=lm(y~x1+x2+x3,data=df)
summary(model)
autoplot(m2)

#Bayesian model
model=brm(y~x1+x2+x3,
          data=df,
          backend='cmdstan')

bayestestR::describe_posterior(model)

```
## Example 2
Same as before, only now we are goinf to use to con indepndent varabels and their interaction.
So now the third regrssion coeff is the interaction parameter.

```
#Aim:Generate and recover linear regression parameters using freq and bayes models

####data with two continuous and interaction
N=100 
x1   =rnorm(N,mean=50, sd=9)
x2   =rnorm(N,mean=200,sd=64)
x3   =x1*x2
error=rnorm(100,0,16)

#Generate the dependent variable (b0=150, b1=-4, b2=2.5, b3=5)
y    =150-(4*x1)+(2.5*x2)+(5*x3)+error
df   =data.frame(x1,x2,x3,y)

#Frequentist model
model=lm(y~x1*x2,data=df)
summary(model)
autoplot(model)

#Bayesian model
model=brm(y~x1*x2,
          data=df,
          backend='cmdstan')

bayestestR::describe_posterior(model)
```
