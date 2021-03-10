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
                       samps_per_cyc=NA,
                       temp_vect=NA,
                       prop_scale=NA,
                       num_cyc=NA,
                       ...) {

  chains = list()
  # prop_scale should be a scalar, a vector of length theta, or a matrix with
  # dimensions length theta x length temp_vect.
  if (any(is.na(prop_scale))) {
    stop("prop_scale must be input for new chains")
  }

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
    print("----")
    print(cc)
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

#    plot(chains[[1]]$eta_vect)
#    readline(prompt="Press [enter] to continue")
#    for (k in 1:length(temp_vect)) {
#      print(chains[[k]]$eta_best)
#      print(chains[[k]]$eta_best)
#    }
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