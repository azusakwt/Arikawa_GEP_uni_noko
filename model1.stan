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
  int<lower=1> K_shape_PPFD;  // number of population-level effects
  matrix[N_PPFD, K_shape_PPFD] X_shape_PPFD;  // population-level design matrix
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
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_GEP = K_GEP - 1;
  matrix[N_GEP, Kc_GEP] Xc_GEP;  // centered version of X_GEP without an intercept
  vector[Kc_GEP] means_X_GEP;  // column means of X_GEP before centering
  int Kc_shape_PPFD = K_shape_PPFD - 1;
  matrix[N_PPFD, Kc_shape_PPFD] Xc_shape_PPFD;  // centered version of X_shape_PPFD without an intercept
  vector[Kc_shape_PPFD] means_X_shape_PPFD;  // column means of X_shape_PPFD before centering
  int Kc_TEMP = K_TEMP - 1;
  matrix[N_TEMP, Kc_TEMP] Xc_TEMP;  // centered version of X_TEMP without an intercept
  vector[Kc_TEMP] means_X_TEMP;  // column means of X_TEMP before centering
  for (i in 2:K_GEP) {
    means_X_GEP[i - 1] = mean(X_GEP[, i]);
    Xc_GEP[, i - 1] = X_GEP[, i] - means_X_GEP[i - 1];
  }
  for (i in 2:K_shape_PPFD) {
    means_X_shape_PPFD[i - 1] = mean(X_shape_PPFD[, i]);
    Xc_shape_PPFD[, i - 1] = X_shape_PPFD[, i] - means_X_shape_PPFD[i - 1];
  }
  for (i in 2:K_TEMP) {
    means_X_TEMP[i - 1] = mean(X_TEMP[, i]);
    Xc_TEMP[, i - 1] = X_TEMP[, i] - means_X_TEMP[i - 1];
  }
}
parameters {
  vector[Kc_GEP] b_GEP;  // population-level effects
  real Intercept_GEP;  // temporary intercept for centered predictors
  real<lower=0> shape_GEP;  // shape parameter
  real Intercept_PPFD;  // temporary intercept for centered predictors
  // parameters for spline s(month,k=4,bs="cc")
  // standarized spline coefficients
  vector[knots_PPFD_1[1]] zs_PPFD_1_1;
  real<lower=0> sds_PPFD_1_1;  // standard deviations of spline coefficients
  vector[Kc_shape_PPFD] b_shape_PPFD;  // population-level effects
  real Intercept_shape_PPFD;  // temporary intercept for centered predictors
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
    vector[N_PPFD] shape_PPFD = Intercept_shape_PPFD + Xc_shape_PPFD * b_shape_PPFD;
    // initialize linear predictor term
    vector[N_TEMP] mu_TEMP = Intercept_TEMP + rep_vector(0.0, N_TEMP) + Zs_TEMP_1_1 * s_TEMP_1_1 + Zs_TEMP_2_1 * s_TEMP_2_1;
    for (n in 1:N_PPFD) {
      // apply the inverse link function
      shape_PPFD[n] = exp(shape_PPFD[n]);
    }
    for (n in 1:N_GEP) {
      // apply the inverse link function
      mu_GEP[n] = shape_GEP * exp(-(mu_GEP[n]));
    }
    for (n in 1:N_PPFD) {
      // apply the inverse link function
      mu_PPFD[n] = shape_PPFD[n] * exp(-(mu_PPFD[n]));
    }
    target += gamma_lpdf(Y_GEP | shape_GEP, mu_GEP);
    target += gamma_lpdf(Y_PPFD | shape_PPFD, mu_PPFD);
    target += normal_id_glm_lpdf(Y_TEMP | Xc_TEMP, mu_TEMP, b_TEMP, sigma_TEMP);
  }
  // priors including all constants
  target += student_t_lpdf(b_GEP[1] | 3, 0, 1);
  target += student_t_lpdf(b_GEP[2] | 3, 0, 1);
  target += student_t_lpdf(b_GEP[3] | 3, 0, 1);
  target += student_t_lpdf(b_GEP[4] | 3, 0, 1);
  target += student_t_lpdf(Intercept_GEP | 3, 3, 2);
  target += normal_lpdf(shape_GEP | 0, 1)
    - 1 * normal_lccdf(0 | 0, 1);
  target += normal_lpdf(Intercept_PPFD | 0, 1);
  target += student_t_lpdf(sds_PPFD_1_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_PPFD_1_1);
  target += student_t_lpdf(b_shape_PPFD | 3, 0, 1);
  target += normal_lpdf(Intercept_shape_PPFD | 0, 2);
  target += student_t_lpdf(b_TEMP | 3, 0, 1);
  target += normal_lpdf(Intercept_TEMP | 0, 2);
  target += student_t_lpdf(sds_TEMP_1_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_TEMP_1_1);
  target += student_t_lpdf(sds_TEMP_2_1 | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(zs_TEMP_2_1);
  target += student_t_lpdf(sigma_TEMP | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
}
generated quantities {
  // actual population-level intercept
  real b_GEP_Intercept = Intercept_GEP - dot_product(means_X_GEP, b_GEP);
  // actual population-level intercept
  real b_PPFD_Intercept = Intercept_PPFD;
  // actual population-level intercept
  real b_shape_PPFD_Intercept = Intercept_shape_PPFD - dot_product(means_X_shape_PPFD, b_shape_PPFD);
  // actual population-level intercept
  real b_TEMP_Intercept = Intercept_TEMP - dot_product(means_X_TEMP, b_TEMP);
}
