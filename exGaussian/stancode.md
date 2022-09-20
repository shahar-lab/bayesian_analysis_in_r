# This is the model block stan uses when using an exGaussian family with a 'log' ling for both beta and sigma

So we can see that a few important thing are happening here:
(1) mu is the RT mean (not the exguassian 'mu' which is [mu - beta] here)
(2) sigma has an exp() transformation. So if stan sampled sigma = 0 for one itreation, exp_mod_normal_lpdf will get exp(0). <br>
Also important here is that for some reason the sigma transfomation applys only when sampling with data, but not when sampling only priors, when prior='only'
(3) Beta is the tau parameter in the exgaussian.Like sigma, it has an exp() transformation. This is true fro both when sampling with empirical data and when using prior='only'.
(4) The likelihood function actually gets lambda which is inv(beta) but this doesn't affect us since we are usually intrested in beta (tau), not lambda

When interpeting the coeff:
for the mean, just use the regular coeff you get. b0_meanrt=b0

for tau, you need to use exp(). Note this has to be on the whole estimate. E.g., if we have beta=b0+b1*x1 
you can have b0_tau=exp(bo_tau), and the b1_tau=exp(b0_tau+b1_tau*1)-exp(b0_tau). (exp(b1_tau) won't work here - important!).

for mu (the exgaussian mu) you need to substract tau from the mean.
so b0_mu=bo_meanrt-exp(b0_tau) or b1_mu=b1_meanrt-[exp(b0_tau+b1_tau*1)-exp(b0_tau)]

```
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = X * b;
    // initialize linear predictor term
    vector[N] sigma = X_sigma * b_sigma;
    // initialize linear predictor term
    vector[N] beta = X_beta * b_beta;
    for (n in 1 : N) {
      // apply the inverse link function
      sigma[n] = exp(sigma[n]);
    }
    for (n in 1 : N) {
      // apply the inverse link function
      beta[n] = exp(beta[n]);
    }
    target += exp_mod_normal_lpdf(Y | mu - beta, sigma, inv(beta));
  }
  // priors including constants
  target += lprior;
}
```
