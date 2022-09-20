# This is the model block stan uses when using an exGaussian family with a 'log' ling for both beta and sigma

So we can see that a few important thing are happening here:

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
