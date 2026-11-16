// Witch hat in peak-centered local time, with the per-person spreads carried
// on the LOG scale instead of the natural scale.
//
// Identical to witch_toy.stan in data, curve, likelihood and priors on the
// means. The ONLY substantive difference is the three sigma_* parameters and
// how they enter `transformed parameters`. witch_toy.stan writes
//
//     V = mu_V * exp((sigma_V / mu_V) * V_raw)
//
// so the log-scale multiplier is a RATIO whose denominator is itself sampled.
// As mu_V approaches zero that ratio diverges and the posterior funnels: the
// adapted step size collapses and a chain that wanders in cannot climb back
// out. Here the multiplier is a parameter in its own right, so no division by
// a sampled quantity appears anywhere in the model.
//
// ---------------------------------------------------------------------------
// PSI. This file is otherwise a copy of witch_toy_reparam.stan, and the only
// difference is the prior on psi. That file puts a FLAT prior on [0, span[j]],
// the whole record. This one puts a Beta(a, a) on [psi_lo[j], psi_hi[j]],
// where psi_lo is the last negative test before the first positive and psi_hi
// the first negative after the last positive, both supplied as data. Nothing
// else changes -- same parameters, same likelihood, same priors on the means
// and spreads -- so the two arms differ in the psi prior alone.
//
// WHY THE BRACKET. psi lies strictly inside it for every curve in the family,
// as a matter of geometry rather than of assumption: the cone is positive on
// (psi - delt_p, psi + delt_c), so a test that read negative sits at or beyond
// a foot, and the peak is between the feet. witch_toy_reparam.stan's support
// instead leans on the sampling design -- the driver's W_DAYS >= L_0 -- to
// keep the peak inside the record. The bracket is also 2-3x narrower than the
// record: measured over the shipped generator it runs 11-15 days against a
// span of 20-29.
//
// WHY BETA RATHER THAN UNIFORM. At psi = psi_lo the peak would coincide with a
// test that read negative, which is impossible for V > 0; likewise at psi_hi.
// A flat prior gives those endpoints the same density as the middle. Beta(a, a)
// with a > 1 sends the density to zero at both ends, which is correct, and
// puts the mode at the midpoint.
//
// ON THE MIDPOINT, HONESTLY. The true peak does NOT sit at the midpoint of
// this bracket. It sits at delt_p / (delt_p + delt_c) of the way across --
// 0.6 for the generator's 6 and 4 -- and measured over the shipped grid the
// true relative position has a median of 0.57-0.60, with EVERY person above
// 0.5 at l_0 = 1. A symmetric prior therefore pulls psi early, which shortens
// the fitted rise and lengthens the fitted fall, biasing mu_delt_p down and
// mu_delt_c up. It is used anyway because the correct centre is a function of
// delt_p and delt_c -- the very quantities being estimated -- so centering
// there would smuggle the answer into the prior. The bias is kept small by
// keeping `a` small instead: the Beta(a, a) density ratio between 0.6 and 0.5
// is 0.96^(a-1), so 0.96 at a = 2 and 0.85 at a = 5. The prior's job here is
// to supply curvature and a narrow support, not information.
//
// `a` IS DATA, and a = 1 is exactly the uniform. Run with psi_shape = 1,
// psi_lo = 0 and psi_hi = span, this file reduces to witch_toy_reparam.stan
// exactly -- beta(1, 1) contributes a constant to the log density -- which is
// the control arm for anything the bracket is credited with.

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1,upper=J> id[N];
  vector[N] x;
  vector[N] y;
  vector[J] psi_lo;               // last negative before the first positive
  vector[J] psi_hi;               // first negative after the last positive
  real<lower=1> psi_shape;        // Beta(a, a) on psi_unit; 1 is the uniform
  real<lower=0> prior_sd;
}

parameters {
  real<lower=0> mu_V;
  real<lower=0> mu_delt_p;
  real<lower=0> mu_delt_c;

  // Log-scale spreads: dimensionless multipliers, NOT the natural-scale
  // (load units / days) spreads that witch_toy.stan declares here.
  real<lower=0> sigma_log_V;
  real<lower=0> sigma_log_delt_p;
  real<lower=0> sigma_log_delt_c;

  vector[J] V_raw;
  vector[J] delt_p_raw;
  vector[J] delt_c_raw;

  vector<lower=0,upper=1>[J] psi_unit;

  real<lower=0> sigma_eps;
}

transformed parameters {
  // Lognormal non-centering, as before -- but the multiplier no longer
  // divides by a sampled mean, which is the whole point of this file.
  vector[J] V      = mu_V      * exp(sigma_log_V      * V_raw);
  vector[J] delt_p = mu_delt_p * exp(sigma_log_delt_p * delt_p_raw);
  vector[J] delt_c = mu_delt_c * exp(sigma_log_delt_c * delt_c_raw);

  // Beta on [psi_lo[j], psi_hi[j]]. The map from psi_unit is affine, so its
  // Jacobian is constant and needs no adjustment -- the same argument that let
  // witch_toy_reparam.stan rescale by span.
  vector[J] psi = psi_lo + psi_unit .* (psi_hi - psi_lo);
}

model {
  mu_V      ~ normal(6, prior_sd);
  mu_delt_p ~ normal(6, prior_sd);
  mu_delt_c ~ normal(4, prior_sd);

  // prior_sd divided by the prior centre of the matching mean. witch_toy.stan
  // puts normal(0, prior_sd) on a NATURAL-scale spread, which at the prior
  // centre implies a log-scale spread of prior_sd / centre. Dividing here
  // reproduces that implied width, so this arm differs from the baseline in
  // GEOMETRY ONLY. Without the divisor the reparameterised arm would also
  // carry a ~6x wider prior and the experiment would confound the two.
  sigma_log_V      ~ normal(0, prior_sd / 6);
  sigma_log_delt_p ~ normal(0, prior_sd / 6);
  sigma_log_delt_c ~ normal(0, prior_sd / 4);

  V_raw      ~ normal(0, 1);
  delt_p_raw ~ normal(0, 1);
  delt_c_raw ~ normal(0, 1);

  sigma_eps ~ normal(0, 0.5);

  // The one line that separates this file from witch_toy_reparam.stan, where
  // psi_unit is bounded but never given a prior and so picks up Stan's
  // implicit uniform. At psi_shape = 1 this IS that uniform.
  psi_unit ~ beta(psi_shape, psi_shape);

  for (n in 1:N){
    int j = id[n];
    real t = x[n] - psi[j];
    real mu_n;

    if (t > -delt_p[j] && t <= 0){
      mu_n = V[j] * (1 + t / delt_p[j]);
    }
    else if (t > 0 && t < delt_c[j]){
      mu_n = V[j] * (1 - t / delt_c[j]);
    }
    else{
      mu_n = 0;
    }

    y[n] ~ normal(mu_n, sigma_eps);
  }
}

generated quantities {
  real mu_duration = mu_delt_p + mu_delt_c;

  // Natural-scale spreads, recovered by inverting this file's mapping, so the
  // two arms' sigma_V / sigma_delt_p / sigma_delt_c can be read against each
  // other and against the generating SIGMA_IND directly.
  real sigma_V      = sigma_log_V      * mu_V;
  real sigma_delt_p = sigma_log_delt_p * mu_delt_p;
  real sigma_delt_c = sigma_log_delt_c * mu_delt_c;
}
