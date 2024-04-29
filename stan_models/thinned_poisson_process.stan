
    data{
      int<lower = 1> N;
      vector[N] X;
      vector[N] Z;
      // Only arrays can store integer values
      // https://mc-stan.org/docs/stan-users-guide/indexing-efficiency.html
      array[N] int y;
    }
    parameters{
      real alpha;
      real beta;
      real gamma;
      real delta;
    }
    model{
      //priors
      target += normal_lpdf(alpha | 0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(gamma | 0,10);
      target += normal_lpdf(delta | 0,10);

      // likelihood
      // log parameterization so we dont have to exponentiate coefs before passing into poisson distribution
      // https://mc-stan.org/docs/functions-reference/poisson-distribution-log-parameterization.html
      target += poisson_log_lpmf(y | alpha + beta * X + gamma + delta * Z);
    }
    generated quantities{
      vector[N] lambda_rep;
      vector[N] lambda_rep_log;
      vector[N] b_rep;
      vector[N] lambda_bias;
      array[N] int y_rep;
      lambda_rep = exp(alpha + beta * X);
      lambda_rep_log = alpha + beta * X;
      //.*abs(alpha + beta * X);
      b_rep = exp(gamma + delta * Z);
      //.*abs(gamma + delta * Z);
      // Thin the process
      //lambda_bias = lambda_rep .* b_rep;
      lambda_bias = exp(alpha + beta * X + gamma + delta * Z);
      // User guide has some examples of log poisson generation and poisson PPCs
      y_rep = poisson_log_rng(alpha + beta * X + gamma + delta * Z);
    }
