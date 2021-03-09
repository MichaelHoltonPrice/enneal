# clear the workspace
rm(list=ls())

# Define some functions used in testing

# Calculate the un-normalized negative log of the multivariate Gaussian
# density.
unnorm_neg_log_gauss_dens <- function(x,mu=0,covMat=1,invCovMat=NA) {
  # To speed computation, invCovMat (the inverse covariance matrix) can be input
  # instead of the covariance matrix (covMat). covMat is ignored if invCovMat
  # is not NA
  if(all(is.na(invCovMat))) {
    invCovMat <- solve(covMat)
  }
  v <- as.matrix(x-mu)
  return(0.5*t(v) %*% invCovMat %*% v)
}

# Parameters for testing scalar function
mu_scalar <- -1.5
covMat_scalar <- 0.2
invCovMat_scalar <- 1/covMat_scalar
x0_scalar <- 0

# (1) Scalar input to do_mh_sampling_at_temp [new chain]
# Make sure we can sample a new chain and that an object of class mh_chain is
# returned; also check the values and shapes of the elements of the returned
# object, chain.
expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  NA
)

expect_equal(
  class(chain),
  "mh_chain"
)

expect_equal(
  names(chain),
  c("theta0",
    "eta0",
    "temp",
    "prop_scale",
    "eta_best",
    "theta_best",
    "accept_vect",
    "eta_vect",
    "num_samp",
    "eta",
    "theta",
    "num_samp_vect",
    "neg_log_obj_func",
    "temp",
    "save_theta",
    "theta_mat")
)

expect_equal(
  chain$theta0,
  x0_scalar
)

expect_equal(
  chain$eta0,
  unnorm_neg_log_gauss_dens(x0_scalar,mu=mu_scalar,covMat_scalar)
)

expect_equal(
  chain$prop_scale,
  0.1
)

expect_equal(
  length(chain$eta_best),
  1
)

expect_equal(
  length(chain$theta_best),
  1
)

expect_equal(
  length(chain$accept_vect),
  chain$num_samp_vect
)

expect_equal(
  length(chain$eta_vect),
  chain$num_samp_vect
)

expect_equal(
  length(chain$eta),
  1
)

expect_equal(
  length(chain$theta),
  1
)

expect_equal(
  chain$num_samp_vect,
  400
)

expect_equal(
  chain$neg_log_obj_func,
  unnorm_neg_log_gauss_dens
)

expect_equal(
  chain$temp,
  2
)

expect_equal(
  chain$prop_scale,
  0.1
)

expect_equal(
  chain$save_theta,
  TRUE
)

expect_equal(
  dim(chain$theta_mat),
  c(1,chain$num_samp_vect)
)

# Check that, as applicable, each error that can be thrown for new chains is
# in fact thrown (and also that the first error on the following line, which
# involves an incorrect first input, is also thrown).
expect_error(
  chain <- do_mh_sampling_at_temp("cannot_use_strings",
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  paste0("init should be a starting parameter vector or continued chain of ",
         "class mh_mcmc")
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  "num_samp must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  num_samp=400,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  "neg_log_obj_func must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  prop_scale=0.1,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  "temp must be input for new chains"
)
expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  "prop_scale must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_scalar,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=c(0.1,0.2),
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar,
                                  save_theta=TRUE),
  "prop_scale must be a scalar or have the same length as theta"
)

# A dummy function that returns infinity
infinity_func <- function(x){Inf}
expect_error(
  chain <- do_mh_sampling_at_temp(0,
                                  num_samp=400,
                                  neg_log_obj_func=infinity_func,
                                  temp=2,
                                  prop_scale=0.1,
                                  save_theta=TRUE),
  paste0("The negative log objective function is not finite for the input ",
         "initialization vector theta0")
)


# (2) Scalar input to do_mh_sampling_at_temp [continued chain]
# Make sure we can sample when initializing with a continuing chain and that an
# object of class mh_chain is returned; also check the values and shapes of the
# elements of the returned object, chain2.
expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_scalar),
  NA
)

expect_equal(
  class(chain2),
  "mh_chain"
)

expect_equal(
  names(chain2),
  c("theta0",
    "eta0",
    "temp",
    "prop_scale",
    "eta_best",
    "theta_best",
    "accept_vect",
    "eta_vect",
    "num_samp",
    "eta",
    "theta",
    "num_samp_vect",
    "neg_log_obj_func",
    "temp",
    "save_theta",
    "theta_mat")
)

expect_equal(
  chain2$theta0,
  x0_scalar
)

expect_equal(
  chain2$eta0,
  unnorm_neg_log_gauss_dens(x0_scalar,mu=mu_scalar,covMat_scalar)
)

expect_equal(
  chain2$prop_scale,
  0.1
)

expect_equal(
  length(chain2$eta_best),
  1
)

expect_equal(
  length(chain2$theta_best),
  1
)

expect_equal(
  length(chain2$accept_vect),
  sum(chain2$num_samp_vect)
)

expect_equal(
  length(chain2$eta_vect),
  sum(chain2$num_samp_vect)
)

expect_equal(
  length(chain2$eta),
  1
)

expect_equal(
  length(chain2$theta),
  1
)

expect_equal(
  chain2$num_samp_vect,
  c(400,400)
)

expect_equal(
  chain2$neg_log_obj_func,
  unnorm_neg_log_gauss_dens
)

expect_equal(
  chain2$temp,
  2
)

expect_equal(
  chain2$save_theta,
  TRUE
)

expect_equal(
  dim(chain2$theta_mat),
  c(1,sum(chain$num_samp_vect))
)

# Check that, as applicable, each error that can be thrown for continuing chains
# is in fact thrown.
expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                   neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                   mu=mu_scalar,
                                   invCovMat=invCovMat_scalar),

  "neg_log_log_obj_func should not be input for continuing chains."
)

expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                   save_theta=TRUE,
                                   mu=mu_scalar,
                                   invCovMat=invCovMat_scalar),
  "save_theta should not be input for continuing chains."
)

# (3) Vector input to do_mh_sampling_at_temp [new chain]
# Parameters for testing vector function
mu_vector <- c(-1.5,0.25)
covMat_vector <- diag(c(1,2))
invCovMat_vector <- solve(covMat_vector)
x0_vector <- c(0,0)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  NA
)

expect_equal(
  class(chain),
  "mh_chain"
)
expect_equal(
  names(chain),
  c("theta0",
    "eta0",
    "temp",
    "prop_scale",
    "eta_best",
    "theta_best",
    "accept_vect",
    "eta_vect",
    "num_samp",
    "eta",
    "theta",
    "num_samp_vect",
    "neg_log_obj_func",
    "temp",
    "save_theta",
    "theta_mat")
)

expect_equal(
  chain$theta0,
  x0_vector
)

expect_equal(
  chain$eta0,
  unnorm_neg_log_gauss_dens(x0_vector,mu=mu_vector,covMat_vector)
)

expect_equal(
  chain$prop_scale,
  c(0.1,0.1)
)

expect_equal(
  length(chain$eta_best),
  1
)

expect_equal(
  length(chain$theta_best),
  2
)

expect_equal(
  length(chain$accept_vect),
  chain$num_samp_vect
)

expect_equal(
  length(chain$eta_vect),
  chain$num_samp_vect
)

expect_equal(
  length(chain$eta),
  1
)

expect_equal(
  length(chain$theta),
  2
)

expect_equal(
  chain$num_samp_vect,
  400
)

expect_equal(
  chain$neg_log_obj_func,
  unnorm_neg_log_gauss_dens
)

expect_equal(
  chain$temp,
  2
)

expect_equal(
  chain$save_theta,
  TRUE
)

expect_equal(
  dim(chain$theta_mat),
  c(2,chain$num_samp_vect)
)

# Check that, as applicable, each error that can be thrown for new chains is
# in fact thrown (and also that the first error on the following line, which
# involves an incorrect first input, is also thrown).
expect_error(
  chain <- do_mh_sampling_at_temp("cannot_use_strings",
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  paste0("init should be a starting parameter vector or continued chain of ",
         "class mh_mcmc")
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  "num_samp must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  num_samp=400,
                                  temp=2,
                                  prop_scale=0.1,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  "neg_log_obj_func must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  prop_scale=0.1,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  "temp must be input for new chains"
)
expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  "prop_scale must be input for new chains"
)

expect_error(
  chain <- do_mh_sampling_at_temp(x0_vector,
                                  num_samp=400,
                                  neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                  temp=2,
                                  prop_scale=c(0.1,0.2,0.3),
                                  mu=mu_vector,
                                  invCovMat=invCovMat_vector,
                                  save_theta=TRUE),
  "prop_scale must be a scalar or have the same length as theta"
)

# (4) Vector input to do_mh_sampling_at_temp [continued chain]
# Make sure we can sample when initializing with a continuing chain and that an
# object of class mh_chain is returned; also check the values and shapes of the
# elements of the returned object, chain2.
expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                  mu=mu_scalar,
                                  invCovMat=invCovMat_vector),
  NA
)

expect_equal(
  class(chain2),
  "mh_chain"
)

expect_equal(
  names(chain2),
  c("theta0",
    "eta0",
    "temp",
    "prop_scale",
    "eta_best",
    "theta_best",
    "accept_vect",
    "eta_vect",
    "num_samp",
    "eta",
    "theta",
    "num_samp_vect",
    "neg_log_obj_func",
    "temp",
    "save_theta",
    "theta_mat")
)

expect_equal(
  chain2$theta0,
  x0_vector
)

expect_equal(
  chain2$eta0,
  unnorm_neg_log_gauss_dens(x0_vector,mu=mu_vector,covMat_vector)
)

expect_equal(
  chain2$prop_scale,
  c(0.1,0.1)
)

expect_equal(
  length(chain2$eta_best),
  1
)

expect_equal(
  length(chain2$theta_best),
  2
)

expect_equal(
  length(chain2$accept_vect),
  sum(chain2$num_samp_vect)
)

expect_equal(
  length(chain2$eta_vect),
  sum(chain2$num_samp_vect)
)

expect_equal(
  length(chain2$eta),
  1
)

expect_equal(
  length(chain2$theta),
  2
)

expect_equal(
  chain2$num_samp_vect,
  c(400,400)
)

expect_equal(
  chain2$neg_log_obj_func,
  unnorm_neg_log_gauss_dens
)

expect_equal(
  chain2$temp,
  2
)

expect_equal(
  chain2$save_theta,
  TRUE
)

expect_equal(
  dim(chain2$theta_mat),
  c(2,sum(chain$num_samp_vect))
)

# Check that, as applicable, each error that can be thrown for continuing chains
# is in fact thrown.
expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                   neg_log_obj_func=unnorm_neg_log_gauss_dens,
                                   mu=mu_vector,
                                   invCovMat=invCovMat_vector),

  "neg_log_log_obj_func should not be input for continuing chains."
)

expect_error(
  chain2 <- do_mh_sampling_at_temp(chain,
                                   save_theta=TRUE,
                                   mu=mu_vector,
                                   invCovMat=invCovMat_vector),
  "save_theta should not be input for continuing chains."
)
