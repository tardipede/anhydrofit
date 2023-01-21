#' Fit Exponential model to anhydrobiosis recovery data
#'
#' @param data a data frame table with observations in rows and variables in columns
#' @param column_IDs a vector with the numeric indexes of the column containing in order: population, time, moving individuals, total individuals. Default: c(1,2,3,4).
#' @param prop numeric value between 0 and 1 indicating the quantile for the calculation of Q. Default: 0.90.
#' @param n.iter number of iterations for the JAGS model. Default: 1000000.
#' @param ... additional arguments for \code{R2Jags::jags}.
#' @return a JAGS model

fit_model_exponential <- function(data,
                                  column_IDs = c(1, 2, 3, 4),
                                  prop = 0.9,
                                  n.iter = 1000000,
                                  ...) {

    # Rename columns to match the script below
  colnames(data)[column_IDs] <- c("population", "time", "moving", "total")

  # Set the population factors in alphabetical order
  data$population <- factor(data$population, levels = sort(unique(data$population)))
  # Save the populations order
  pops_order <- levels(data$population)

  # Prepare data for JAGS
  data.jags <- list(
    time = data$time,
    moving = data$moving,
    total = data$total,
    species_num = as.numeric(as.factor(data$population)),
    Nsp = length(unique(data$population)),
    prop = prop
  )

  f <- tempfile("model.txt")
  sink(f)
  cat("

  # Prepare JAGS model
  model {
    # Priors
    for (sp in 1:Nsp) {                # weak estimate Q and p separately for each species
      p[sp] ~ dunif(0,1)               # low informative prior bound to be >0
      Q[sp] ~ dnorm(0,1.0E-2)T(0,)     # low informative prior
    }
    # Likelihood
    for (i in 1:length(time)){
      lambda[i] <- -log(1-prop)/Q[species_num[i]]
      mu[i] <- (1-(exp(-1*lambda[i]*time[i]))) * p[species_num[i]]
      mu_regularized[i]  <- ifelse(mu[i] == 0, 0.00001, ifelse(mu[i] == 1, 0.99999, mu[i])) # This is to avoid the model to crash when p is exactly 0 or 1
      moving[i]~dbinom(mu_regularized[i],total[i])
    }} ")

  closeAllConnections()
  # Run JAGS model


  results.jags <- R2jags::jags(
    data = data.jags,
    parameters.to.save = c("p", "Q"),
    model.file = f,
    n.iter = n.iter,...
  )

  attr(results.jags, "data") <- data
  attr(results.jags, "params") <- list(prop = prop)
  attr(results.jags, "num_to_pop") <- data.frame(population = pops_order, number = 1:length(pops_order))
  attr(results.jags, "model") <- "exponential"
  attr(results.jags, "type") <- "simple"

  return(results.jags)
}


#' Fit Weibull model to anhydrobiosis recovery data
#'
#' @param data a data frame table with observations in rows and variables in columns
#' @param column_IDs a vector with the numeric indexes of the column containing in order: population, time, moving individuals, total individuals. Default: c(1,2,3,4).
#' @param prop numeric value between 0 and 1 indicating the quantile for the calculation of Q. Default: 0.90.
#' @param n.iter number of iterations for the JAGS model. Default: 1000000.
#' @param ... additional arguments for \code{R2Jags::jags}.
#' @return a JAGS model

fit_model_weibull <- function(data,
                              column_IDs = c(1, 2, 3, 4),
                              prop = 0.9,
                              n.iter = 1000000,
                               ...) {

  # Rename columns to match the script below
  colnames(data)[column_IDs] <- c("population", "time", "moving", "total")

  # Set the population factors in alphabetical order
  data$population <- factor(data$population, levels = sort(unique(data$population)))
  # Save the populations order
  pops_order <- levels(data$population)

  # Prepare data for JAGS
  data.jags <- list(
    time = data$time,
    moving = data$moving,
    total = data$total,
    species_num = as.numeric(as.factor(data$population)),
    Nsp = length(unique(data$population)),
    prop = prop
  )

  f <- tempfile("model.txt")
  sink(f)
  cat("

  # Prepare JAGS modelz
  model {
    # Priors
    for (sp in 1:Nsp) {                # weak estimate Q and p separately for each species
      p[sp] ~ dunif(0,1)               # low informative prior bound to be >0
      Q[sp] ~ dnorm(0,1.0E-2)T(0,)     # low informative prior
      sh_log[sp] ~ dnorm(0,1)        # weakly informative prior centered at 0
      sh[sp] <- exp(sh_log[sp])

    }
    # Likelihood


    for (i in 1:length(time)){


      upper_lambda[i] <- pow((-log(1-prop)),1/sh[species_num[i]])
      lambda[i] <- upper_lambda[i]/Q[species_num[i]]

      mu[i] <- (1-(exp(-1*pow(lambda[i]*time[i],sh[species_num[i]])))) * p[species_num[i]]

      mu_regularized[i]  <- ifelse(mu[i] == 0, 0.00001, ifelse(mu[i] == 1, 0.99999, mu[i])) # This is to avoid the model to crash when p is exactly 0 or 1
      moving[i]~dbinom(mu_regularized[i],total[i])
    }} ")

  closeAllConnections()
  # Run JAGS model


  results.jags <- R2jags::jags(
    data = data.jags,
    parameters.to.save = c("p", "Q", "sh_log"),
    model.file = f,
    n.iter = n.iter, ...
  )

  attr(results.jags, "data") <- data
  attr(results.jags, "params") <- list(prop = prop)
  attr(results.jags, "num_to_pop") <- data.frame(population = pops_order, number = 1:length(pops_order))
  attr(results.jags, "model") <- "weibull"
  attr(results.jags, "type") <- "simple"

  return(results.jags)
}



#' Wrapper function to fit both Exponential and Weibull models
#'
#' @param model string indicating the model to fit ("exponential" or "weibull")
#' @param data a data frame table with observations in rows and variables in columns
#' @param column_IDs a vector with the numeric indexes of the column containing in order: population, time, moving individuals, total individuals. Default: c(1,2,3,4).
#' @param prop numeric value between 0 and 1 indicating the quantile for the calculation of Q. Default: 0.90.
#' @param n.iter number of iterations for the JAGS model. Default: 1000000.
#' @param ... additional arguments for \code{R2Jags::jags}.
#' @return a JAGS model
#' @export
#'
#' @examples
#' \dontrun{
#' data = data.frame(group = c("A","A","A","B","B","B"),
#'        time = c(1,2,24,1,2,24),
#'        moving = c(0,1,4,2,5,16),
#'        tot = rep(20,6))
#' mod = fit_model(data, model = "weibull")}

fit_model <- function(data,
                      column_IDs = c(1, 2, 3, 4),
                      prop = 0.9,
                      n.iter = 1000000,
                      model = c("weibull","exponential"),
                      ...) {

  if(model == "weibull"){results.jags = fit_model_weibull(data=data, column_IDs, prop=prop, n.iter=n.iter,...)}
  if(model == "exponential"){results.jags = fit_model_exponential(data=data, column_IDs, prop=prop, n.iter=n.iter,...)}

  return(results.jags)
}
