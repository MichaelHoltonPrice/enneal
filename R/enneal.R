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
#' @param neg_log_cost_func The function that calculates the negative log of the
#'   function being sampled.
#' @param num_samp The number of samples to make.
#' @param temp The temperature
#' @param prop_scale The standard deviation(s) to use for the proposal
#'   distribution (either a vector that is the same length as the parameter
#'   vector, theta, or a scalar that is used for all parameters).
#' @param save_theta Whether or not to save and return the accepted parameter
#'   values (assumed False if not input). Cannot be over-ridden for continuing
#'   chains.
#' @param ... Variables required by neg_log_cost_func.
#' @return A list-like object of class "mh_chain" with sampling information
#' @export
do_mh_sampling_at_temp <- function(init,
                           neg_log_cost_func=NA,
                           num_samp=NA,
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

    if(class(neg_log_cost_func) != "function") {
      stop("neg_log_cost_func must be input for new chains")
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
    eta0 <- neg_log_cost_func(theta0,...)
    if (!is.finite(eta0)) {
      stop(paste0("The negative log cost function is not finite for the ",
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
    if (class(neg_log_cost_func) != "function") {
       neg_log_cost_func <- init$neg_log_cost_func
    } else {
      stop("neg_log_cost_func should not be input for continuing chains.")
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

    # The temperature can always be over-ridden
    if(is.na(temp)) {
      temp <- init$temp
    }

    # The proposal scale can always be over-ridden
    if(any(is.na(prop_scale))) {
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
    eta_prop <- neg_log_cost_func(theta_prop,...)

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
                   neg_log_cost_func=neg_log_cost_func,
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
    init$eta <- eta
    init$theta <- theta
    init$eta_best <- eta_best
    init$theta_best <- theta_best
    init$accept_vect <- c(init$accept_vect,accept_vect)
    init$eta_vect <- c(init$eta_vect,eta_vect)
    init$num_samp_vect <- c(init$num_samp_vect,num_samp)
    return(init)
  }

  stop("This point should never be reached")
}

#' Parallel tempering for function minimization
#'
#' Minimize an input cost function using parallel tempering. The cost function
#' is often a negative log-likelihood, and is exponentiated at the accept/reject
#' steps so that each chain is sampling the likelihood function raised to the
#' temperature of the chain. The tempering consists of two steps:
#'
#' (1) Do Metropolis-Hastings sampling of each chain samps_per_cyc (samples per
#'     cycle) times.
#' (2) Randomly choose two adjacent temperatures to attempt a swap.
#'
#' These two steps are performed num_cyc (number of cycles) times. A precise
#' description of these steps follows.
#'
#' Let k = 1...K index temperatures in the vector temp_vect. Each cycle of the
#' algorithm consists of a set of within-in chain Metropolis samples for each
#' temperature, followed by a single attempt to randomly swap samples across
#' chains for a single, randomly chosen pair of adjacent temperatures (k and
#' k+1). M within chain samples are made by calling the function
#' enneal::do_mh_sampling_at_temp. Following these within-chain samples, an
#' attempt is made to swap samples between one pair of adjacent temperatures,
#' also using a Metropolis criterion. The pair for which this is attempted is
#' randomly chosen. Let k be the lower temperature. The weighting accorded the
#' non-swapped configuration is exp(-eta_k/T_k)*exp(-eta_kp1/T_kp1), where kp1
#' stands for k+1 (k plus one) and eta_k (etc.) is the cost function evaluated
#' for the current value of the chain, theta_k. The weighting accorded the
#' swapped configuration is exp(-eta_kp1/T_k)*exp(-eta_k/T_kp1). The swap is
#' always accepted if the swapped configuration has a higher weighting;
#' otherwise, it is accepted with a probability equal to the proportion of these
#' weightings. To be precise, the Metropolis acceptance ratio is
#'
#' a = min(1,exp(-(eta_kp1-eta_k)*(1/T_k-1/T_kp1)))
#'
#' The swap is accepted with probability a, and otherwise rejected. The cycle of
#' within-chain sampling followed by a single swap attempt is repeated num_cyc
#' times. The minimum value of this procedure can be extracted from the output
#' by calling the function par_temper_best.
#' TODO: implement par_temper_best
#'
#' @param theta0 The starting point for sampling (used to initialize all chains)
#' @param neg_log_cost_func The negative log of the cost function
#' @param samps_per_cyc Number of within chain Metropolis samples per cycle
#' @param temp_vect Vector of temperatures for each chain (ordered highest to
#'   lowest)
#' @param prop_scale A speficiation of the standard deviations to use for
#'   proposing new parameter values. prop_scale should be a scalar, a vector of
#'   length(theta), or a matrix with dimensions length(theta) by
#'   length(temp_vect).
#' @param num_cyc Number of cycles (a swap is attempted after each cycle).
#' @param ... Variables required by neg_log_cost_func
#'
#' @return An object of class \code{par_temper} that consists of (a) chains (the
#'   sampled chains) and (b) swap_mat, a matrix summarizing the results of the
#'   swap attempts.

#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export
par_temper <- function(theta0,
                       neg_log_cost_func,
                       samps_per_cyc=200,
                       temp_vect=10^(rev(seq(-1,1,by=.1))),
                       prop_scale=1,
                       num_cyc=100,
                       ...) {

  chains = list()
  # prop_scale should be a scalar, a vector of length theta, or a matrix with
  # dimensions length theta x length temp_vect.
  if (is.vector(prop_scale)) {
    if(length(prop_scale) == 1) {
      # A scalar
      prop_scale <- matrix(prop_scale,length(theta0),length(temp_vect))
    } else if (length(prop_scale) == length(theta0)) {
      prop_scale <- replicate(length(temp_vect),prop_scale)
    } else {
      stop(paste0("If prop_scale is a vector, it should be length 1 or the ",
                  "same length as theta0"))
    }
  } else if (is.matrix(prop_scale)) {
    if(dim(prop_scale) != c(length(theta0),length(temp_vect))) {
      stop(paste0("If prop_scale is a matrix, it should have dimensions ",
                  "length(theta0) by length(temp_vect)"))
    }
  } else {
    stop("prop_scale should be a scalar, vector, or matrix")
  }

  swap_mat <- matrix(NA,3,num_cyc)
  # Iterate over number of cycles
  for (cc in 1:num_cyc) {
    for (k in 1:length(temp_vect)) {
      if (cc == 1) {
        # Start new chains
        chains[[k]] <-
        do_mh_sampling_at_temp(theta0,
                               num_samp=samps_per_cyc,
                               neg_log_cost_func=neg_log_cost_func,
                               temp=temp_vect[k],
                               prop_scale=prop_scale[,k],
                               save_theta=T,
                               ...)
      } else {
        # Extend chains for this cycle
        chains[[k]] <- do_mh_sampling_at_temp(chains[[k]],...)
      }
    }

    # Randomly choose an adjacent pair of temperatures to attempt a swap
    k <- sample(1:(length(temp_vect)-1),1)

    # Extract the current negative log cost from the two adjacent chains
    eta_k     <- chains[[k  ]]$eta
    eta_kp1   <- chains[[k+1]]$eta

    # Extract the temperatures of the two adjacent chains
    T_k     <- temp_vect[k]
    T_kp1   <- temp_vect[k+1]

    # Calculate the swap probability
    a_swap <- exp(-(eta_kp1-eta_k)*(1/T_k-1/T_kp1))
    if(a_swap >= 1) {
      accept <- T
    } else {
      accept <- runif(1) <= a_swap
    }

    if(accept) {
      theta_k     <- chains[[k  ]]$theta
      theta_kp1   <- chains[[k+1]]$theta

      # Swap eta and theta (there is no need to update eta_best and theta_best
      # here).
      chains[[k  ]]$eta <- eta_kp1
      chains[[k+1]]$eta <- eta_k
      chains[[k  ]]$theta <- theta_kp1
      chains[[k+1]]$theta <- theta_k
    }

    swap_mat[1,cc] <- samps_per_cyc*cc
    swap_mat[2,cc] <- k
    if(accept) {
      swap_mat[3,cc] <- 1
    } else {
      swap_mat[3,cc] <- 0
    }
  }
  output <- list(chains=chains,swap_mat=swap_mat)
  return(output)
}