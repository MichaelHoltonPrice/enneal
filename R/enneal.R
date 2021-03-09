#' Do Metropolis-Hastings sampling at a given temperature
#'
#' Do Metropolis-Hastings sampling given an input negative log of some target
#' distribution (probably a negative log-likelihood) at a given temperature. The
#' chain can be initialized either with a starting parameter vector or with a
#' chain created by a previous call to do_mh_sampling_at_temp. For continuing
#' chains, some (but not all) of the inputs from the previously used values can
#' be input.
#'
#' @param init The initializing object (either a parameter vector or a chain
#'   created by a previous call to do_mh_sampling_at_temp).
#' @param num_samp The number of samples to make.
#' @param neg_log_obj_func The function that calculates the negative log of the
#'   function being sampled.
#' @param temp The temperature
#' @param prop_scale The standard deviation(s) to use for the proposal
#'   distribution (either a vector that is the same length as the parameter
#'   vector, theta, or a scalar that is used for all parameters).
#' @param save_theta Whether or not to save and return the accepted parameter
#'   values (assumed False if not input). Cannot be over-ridden for continuing
#'   chains.
#' @param ... Variables required by neg_log_obj_func.
#' @return A list-like object of class "mh_chain" with sampling information
#' @export
do_mh_sampling_at_temp <- function(init,
                           num_samp=NA,
                           neg_log_obj_func=NA,
                           temp=NA,
                           prop_scale=NA,
                           save_theta=NA,
                           ...) {

  # TODO: explore saving the ... parameters for use in continuing chains.

  # Is this a new chain, or are we continuining an existing chain?
  if ("mh_chain" %in% class(init)) {
    new_chain <- FALSE
  } else if (is.numeric(init)) {
    new_chain <- TRUE
  } else {
    stop(paste0("init should be a starting parameter vector or continued ",
                "chain of class mh_mcmc"))
  }

  if(new_chain) {
    # Do error checking on the inputs
    if (is.na(num_samp)) {
      stop("num_samp must be input for new chains")
    }

    if(class(neg_log_obj_func) != "function") {
      stop("neg_log_obj_func must be input for new chains")
    }

    if(is.na(temp)) {
      stop("temp must be input for new chains")
    }

    # prop_scale must always be input
    if (any(is.na(prop_scale))) {
      stop("prop_scale must be input for new chains")
    }

    if (length(prop_scale) != 1 && length(prop_scale) != length(init)) {
      stop("prop_scale must be a scalar or have the same length as theta")
    }

    if (length(prop_scale) == 1) {
      if (length(init) != 1) {
        prop_scale <- rep(prop_scale,length(init))
      }
    }

    # If not input, save_theta is set to False. The default in the function call
    # is NA so we can ensure that the user has not over-ridden the value for
    # continuing chains.
    if(is.na(save_theta)) {
      save_theta <- False
    }

    # For a new chain, calculate the starting negative log-likelihood and set
    # theta0 equal to the input parameter vector.
    theta0 <- init
    eta0 <- neg_log_obj_func(theta0,...)
    if (!is.finite(eta0)) {
      stop(paste0("The negative log objective function is not finite for the ",
                  "input initialization vector theta0"))
    }

    # Set eta_best and theta_best
    theta_best <- theta0
    eta_best   <- eta0

    # Set theta and eta, the currently accepted variables, to theta0 and eta0
    eta <- eta0
    theta <- theta0
  } else {
    # Currently, resetting the sampling function is not supported inside this
    # function. It can be rest by directly editing the function inside the
    # mh_chain object outside the function. Allowing it to be reset seems likely
    # to cause mis-uses more often than it is needed.
    if (class(neg_log_obj_func) != "function") {
       neg_log_obj_func <- init$neg_log_obj_func
    } else {
      stop("neg_log_log_obj_func should not be input for continuing chains.")
    }

    # For an existing chain, get eta0, theta0, eta_best, and theta_best from
    # init
    eta0       <- init$eta0
    theta0     <- init$theta0
    eta_best   <- init$eta_best
    theta_best <- init$theta_best

    # Set theta and eta, the currently accepted variables, using the last values
    # from the input chain.
    #num_prev_samples <- init$total_samples
    eta   <- init$eta
    theta <- init$theta

    # If necessary, set variables using the values stored in init for continuing
    # chains or, for new chains, throw errors if necessary variables are not
    # input.
    if (is.na(num_samp)) {
      num_samp <- init$num_samp
    }

    if (is.na(num_samp)) {
      num_samp <- init$num_samp
    }


    # The temperature can always be over-ridden
    if(is.na(temp)) {
      temp <- init$temp
    }

    # The proposal scale can always be over-ridden
    if(any(is.na(temp))) {
      prop_scale <- init$prop_scale
    }

    # The number of samples can always be over-ridden
    if(is.na(num_samp)) {
      num_samp <- init$num_samp_vect[length(init$num_samp_vect)]
    }

    # Currently, at least, the value of save_theta cannot be overwritten.
    if(is.na(save_theta)) {
      save_theta <- init$save_theta
    } else {
      stop("save_theta should not be input for continuing chains.")
    }

  }

  # accept_vect is a vector of length num_samp that records whether each sample
  # was accepted or rejected
  accept_vect <- rep(NA,num_samp)

  # eta_vect is a vector of length num_samp that records the value of the
  # negative log-likelihood (or, if direct is TRUE, the sampled function).
  eta_vect <- rep(NA,num_samp)

  # If necessary, initialize the matrix of save parameter values, theta_mat,
  # which has dimensions Number of Parameters x Number of Samples.
  if(save_theta) {
    theta_mat <- matrix(NA,length(theta0),num_samp)
  }

  for(n in 1:num_samp) {
    # Create a new proposal parameter
    theta_prop <- theta + rnorm(length(theta))*prop_scale
    eta_prop <- neg_log_obj_func(theta_prop,...)

    if(!is.finite(eta_prop)) {
      accept <- F
    } else {
      alpha <- min(1,exp(-(eta_prop-eta)/temp))
      accept <- runif(1) < alpha
    }

    accept_vect[n] <- accept
    if(!accept) {
      theta_prop <- theta
      eta_prop   <- eta
    } else {
      if(eta_prop < eta_best) {
        theta_best <- theta_prop
        eta_best   <- eta_prop
      }
    }
    eta_vect[n] <- eta_prop

    # Get ready for next sample
    theta <- theta_prop
    eta   <- eta_prop

    if(save_theta) {
      theta_mat[,n] <- theta
    }
  }

  # If this is a new chain, create and return a new list-like object. If this is
  # a continuing chain, update the input chain and return it.
  if(new_chain) {
    output <- list(theta0=theta0,
                   eta0=eta0,
                   temp=temp,
                   prop_scale=prop_scale,
                   eta_best=eta_best,
                   theta_best=theta_best,
                   accept_vect=accept_vect,
                   eta_vect=eta_vect,
                   num_samp=num_samp,
                   eta=eta,
                   theta=theta,
                   num_samp_vect=num_samp,
                   neg_log_obj_func=neg_log_obj_func,
                   temp=temp,
                   save_theta=save_theta)
    if(save_theta) {
      output$theta_mat <- theta_mat
    }
    class(output) <- "mh_chain"
    return(output)
  } else {
    # Not a new chain

    # Update modified values
    init$eta_best <- eta_best
    init$theta_best <- theta_best
    init$accept_vect <- c(init$accept_vect,accept_vect)
    init$eta_vect <- c(init$eta_vect,accept_vect)
    init$num_samp_vect <- c(init$num_samp_vect,num_samp)
    return(init)
  }

  stop("This point should never be reached")
}