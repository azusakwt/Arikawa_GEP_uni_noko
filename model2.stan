// generated with brms 2.14.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_GEP;  // number of observations
  vector[N_GEP] Y_GEP;  // response variable
  int<lower=1> K_GEP;  // number of population-level effects
  matrix[N_GEP, K_GEP] X_GEP;  // population-level design matrix
  int<lower=1> N_PPFD;  // number of observations
  vector[N_PPFD] Y_PPFD;  // response variable
  // data for spline s(month,k=4,bs="cc")
  int nb_PPFD_1;  // number of bases
  int knots_PPFD_1[nb_PPFD_1];  // number of knots
  // basis function matrices
  matrix[N_PPFD, knots_PPFD_1[1]] Zs_PPFD_1_1;
  int<lower=1> K_sigma_PPFD;  // number of population-level effects
  matrix[N_PPFD, K_sigma_PPFD] X_sigma_PPFD;  // population-level design matrix
  int<lower=1> N_TEMP;  // number of observations
  vector[N_TEMP] Y_TEMP;  // response variable
  int<lower=1> K_TEMP;  // number of population-level effects
  matrix[N_TEMP, K_TEMP] X_TEMP;  // population-level design matrix
  // data for spline s(month,k=6,by=state,bs="cc",id=0)Desertified
  int nb_TEMP_1;  // number of bases
  int knots_TEMP_1[nb_TEMP_1];  // number of knots
  // basis function matrices
  matrix[N_TEMP, knots_TEMP_1[1]] Zs_TEMP_1_1;
  // data for spline s(month,k=6,by=state,bs="cc",id=0)Vegetated
  int nb_TEMP_2;  // number of bases
  int knots_TEMP_2[nb_TEMP_2];  // number of knots
  // basis function matrices
  matrix[N_TEMP, knots_TEMP_2[1]] Zs_TEMP_2_1;
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_GEP = K_GEP - 1;
  matrix[N_GEP, Kc_GEP] Xc_GEP;  // centered version of X_GEP without an intercept
  vector[Kc_GEP] means_X_GEP;  // column means of X_GEP before centering
  int Kc_sigma_PPFD = K_sigma_PPFD - 1;
  matrix[N_PPFD, Kc_sigma_PPFD] Xc_sigma_PPFD;  // centered version of X_sigma_PPFD without an intercept
  vector[Kc_sigma_PPFD] means_X_sigma_PPFD;  // column means of X_sigma_PPFD before centering
  int Kc_TEMP = K_TEMP - 1;
  matrix[N_TEMP, Kc_TEMP] Xc_TEMP;  // centered version of X_TEMP without an intercept
  vector[Kc_TEMP] means_X_TEMP;  // column means of X_TEMP before centering
  vector[nresp] Y[N];  // response array
  for (i in 2:K_GEP) {
    means_X_GEP[i - 1] = mean(X_GEP[, i]);
    Xc_GEP[, i - 1] = X_GEP[, i] - means_X_GEP[i - 1];
  }
  for (i in 2:K_sigma_PPFD) {
    means_X_sigma_PPFD[i - 1] = mean(X_sigma_PPFD[, i]);
    Xc_sigma_PPFD[, i - 1] = X_sigma_PPFD[, i] - means_X_sigma_PPFD[i - 1];
  }
  for (i in 2:K_TEMP) {
    means_X_TEMP[i - 1] = mean(X_TEMP[, i]);
    Xc_TEMP[, i - 1] = X_TEMP[, i] - means_X_TEMP[i - 1];
  }
  for (n in 1:N) {
    Y[n] = transpose([Y_GEP[n], Y_PPFD[n], Y_TEMP[n]]);
  }
}
parameters {
  vector[Kc_GEP] b_GEP;  // population-level effects
  real Intercept_GEP;  // temporary intercept for centered predictors
  real<lower=0> sigma_GEP;  // residual SD
  real Intercept_PPFD;  // temporary intercept for centered predictors
  // parameters for spline s(month,k=4,bs="cc")
  // standarized spline coefficients
  vector[knots_PPFD_1[1]] zs_PPFD_1_1;
  real<lower=0> sds_PPFD_1_1;  // standard deviations of spline coefficients
  vector[Kc_sigma_PPFD] b_sigma_PPFD;  // population-level effects
  real Intercept_sigma_PPFD;  // temporary intercept for centered predictors
  vector[Kc_TEMP] b_TEMP;  // population-level effects
  real Intercept_TEMP;  // temporary intercept for centered predictors
  // parameters for spline s(month,k=6,by=state,bs="cc",id=0)Desertified
  // standarized spline coefficients
  vector[knots_TEMP_1[1]] zs_TEMP_1_1;
  real<lower=0> sds_TEMP_1_1;  // standard deviations of spline coefficients
  // parameters for spline s(month,k=6,by=state,bs="cc",id=0)Vegetated
  // standarized spline coefficients
  vector[knots_TEMP_2[1]] zs_TEMP_2_1;
  real<lower=0> sds_TEMP_2_1;  // standard deviations of spline coefficients
  real<lower=0> sigma_TEMP;  // residual SD
  cholesky_factor_corr[nresp] Lrescor;  // parameters for multivariate linear models
}
transformed parameters {
  // actual spline coefficients
  vector[knots_PPFD_1[1]] s_PPFD_1_1;
  // actual spline coefficients
  vector[knots_TEMP_1[1]] s_TEMP_1_1;
  // actual spline coefficients
  vector[knots_TEMP_2[1]] s_TEMP_2_1;
  // compute actual spline coefficients
  s_PPFD_1_1 = sds_PPFD_1_1 * zs_PPFD_1_1;
  // compute actual spline coefficients
  s_TEMP_1_1 = sds_TEMP_1_1 * zs_TEMP_1_1;
  // compute actual spline coefficients
  s_TEMP_2_1 = sds_TEMP_2_1 * zs_TEMP_2_1;
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_GEP] mu_GEP = Intercept_GEP + Xc_GEP * b_GEP;
    // initialize linear predictor term
    vector[N_PPFD] mu_PPFD = Intercept_PPFD + rep_vector(0.0, N_PPFD) + Zs_PPFD_1_1 * s_PPFD_1_1;
    // initialize linear predictor term
    vector[N_PPFD] sigma_PPFD = Intercept_sigma_PPFD + Xc_sigma_PPFD * b_sigma_PPFD;
    // initialize linear predictor term
    vector[N_TEMP] mu_TEMP = Intercept_TEMP + Xc_TEMP * b_TEMP + Zs_TEMP_1_1 * s_TEMP_1_1 + Zs_TEMP_2_1 * s_TEMP_2_1;
    // multivariate predictor array
    vector[nresp] Mu[N];
    vector[nresp] sigma[N];
    // cholesky factor of residual covariance matrix
    matrix[nresp, nresp] LSigma[N];
    for (n in 1:N_PPFD) {
      // apply the inverse link function
      sigma_PPFD[n] = exp(sigma_PPFD[n]);
    }
    // combine univariate parameters
    for (n in 1:N) {
      Mu[n] = transpose([mu_GEP[n], mu_PPFD[n], mu_TEMP[n]]);
      sigma[n] = transpose([sigma_GEP, sigma_PPFD[n], sigma_TEMP]);
      LSigma[n] = diag_pre_multiply(sigma[n], Lrescor);
    }
    for (n in 1:N) {
      target += multi_normal_cholesky_lpdf(Y[n] | Mu[n], LSigma[n]);
    }
  }
  // priors including all constants
  target += student_t_lpdf(Intercept_GEP | 3, 12.1, 7);
  target += student_t_lpdf(sigma_GEP | 3, 0, 7)
    - 1 * student_t_lccdf(0 | 3, 0, 7);
  target += student_t_lpdf(Intercept_PPFD | 3, 0.4, 2.5);
  target += student_t_lpdf(sds_PPFD_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_PPFD_1_1);
  target += student_t_lpdf(Intercept_sigma_PPFD | 3, 0, 2.5);
  target += student_t_lpdf(Intercept_TEMP | 3, 0.7, 2.5);
  target += student_t_lpdf(sds_TEMP_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_TEMP_1_1);
  target += student_t_lpdf(sds_TEMP_2_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_TEMP_2_1);
  target += student_t_lpdf(sigma_TEMP | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += lkj_corr_cholesky_lpdf(Lrescor | 1);
}
generated quantities {
  // actual population-level intercept
  real b_GEP_Intercept = Intercept_GEP - dot_product(means_X_GEP, b_GEP);
  // actual population-level intercept
  real b_PPFD_Intercept = Intercept_PPFD;
  // actual population-level intercept
  real b_sigma_PPFD_Intercept = Intercept_sigma_PPFD - dot_product(means_X_sigma_PPFD, b_sigma_PPFD);
  // actual population-level intercept
  real b_TEMP_Intercept = Intercept_TEMP - dot_product(means_X_TEMP, b_TEMP);
  // residual correlations
  corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  vector<lower=-1,upper=1>[nrescor] rescor;
  // extract upper diagonal of correlation matrix
  for (k in 1:nresp) {
    for (j in 1:(k - 1)) {
      rescor[choose(k - 1, 2) + j] = Rescor[j, k];
    }
  }
}
