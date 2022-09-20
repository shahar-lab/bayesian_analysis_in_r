# This is the model block stan uses when using an exGaussian family with a 'log' link for both beta and sigma

So we can see that a few important thing are happening here:<br/>
(1) mu is the RT mean (not the exguassian 'mu' which is [mu - beta] here)<br/>
(2) sigma has an exp() transformation. So if stan sampled sigma = 0 for one itreation, exp_mod_normal_lpdf will get exp(0). <br/>
(3) Beta is the tau parameter in the exgaussian.Like sigma, it has an exp() transformation. This is true fro both when sampling with empirical data and when using prior='only'.<br/>
(4) The likelihood function actually gets lambda which is inv(beta) but this doesn't affect us since we are usually intrested in beta (tau), not lambda<br/>
<br/>
When interpeting the coeff:<br/>
for the mean, just use the regular coeff you get. b0_meanrt=b0<br/>
<br/>
for tau, you need to use exp(). Note this has to be on the whole estimate. E.g., if we have beta=b0+b1*x1 <br/>
you can have b0_tau=exp(bo_tau), and the b1_tau=exp(b0_tau+b1_tau*1)-exp(b0_tau). (exp(b1_tau) won't work here - important!).<br/>
<br/>
for mu (the exgaussian mu) you need to substract tau from the mean.<br/>
so b0_mu=bo_meanrt-exp(b0_tau) or b1_mu=b1_meanrt-[exp(b0_tau+b1_tau*1)-exp(b0_tau)]<br/>
<br/>

```
parameters {
  vector[K] b; // population-level effects
  vector[K_sigma] b_sigma; // population-level effects
  vector[K_beta] b_beta; // population-level effects
}
transformed parameters {
  real lprior = 0; // prior contributions to the log posterior
  lprior += student_t_lpdf(b[1] | 3, 0.5, 0.1);
  lprior += student_t_lpdf(b_sigma[1] | 3, 1.051271, 1.01005);
  lprior += student_t_lpdf(b_beta[1] | 3, 1.161834, 1.077884);
}
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
