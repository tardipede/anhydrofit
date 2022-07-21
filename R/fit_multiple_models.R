#' Fit both Exponential and Weibull models to each experimental group
#'
#' @importFrom crayon blue green red magenta
#'
#' @param data a data frame table with observations in rows and variables in columns
#' @param column_IDs a vector with the numeric indexes of the column containing in order: population, time, moving individuals, total individuals. Default: c(1,2,3,4).
#' @param prop numeric value between 0 and 1 indicating the quantile for the calculation of Q. Default: 0.90.
#' @param n.iter number of iterations for the JAGS model. Default: 1000000.
#' @param ... additional arguments for \code{R2Jags::jags}.
#' @return a list of models
#' @export
#' @examples
#' \dontrun{
#' data = data.frame(group = c("A","A","A","B","B","B"),
#'        time = c(1,2,24,1,2,24),
#'        moving = c(0,1,4,2,5,16),
#'        tot = rep(20,6))
#' mod = fit_multiple_models(data)}

fit_multiple_models = function(data,
                               column_IDs = c(1, 2, 3, 4),
                               prop = 0.9,
                               n.iter = 1000000,
                               ...){

  population = NA

  # Rename columns to match the script below
  colnames(data)[column_IDs] <- c("population", "time", "moving", "total")
  # Set the population factors in alphabetical order
  data$population <- factor(data$population, levels = sort(unique(data$population)))
  # Save the populations order
  pops_order <- levels(data$population)


  results = list()
  for(i in 1:length(pops_order)){
    cat(paste0(magenta("\r Running models for: "), blue(paste0(pops_order[i], "\n"))))
    data_temp = subset(data, population == pops_order[i])

    cat(green("\r  Exponential model \n"))
    exp_temp = fit_model(data_temp,n.iter=n.iter,column_IDs=column_IDs,model="exponential", quiet = T, ...)
    Rhats = exp_temp$BUGSoutput$summary[,"Rhat"]
    ESS = exp_temp$BUGSoutput$summary[,"n.eff"]
    n_points = exp_temp$BUGSoutput$n.keep * exp_temp$BUGSoutput$n.chains
    check_condition = all(Rhats<1.1, (ESS/n_points)>0.5)
    if(check_condition == FALSE){cat(red("\r  The model Rhat and ESS seems to be too low, check the model output and if needed increase n.iter \n"))}

    cat(green("\r  Weibull model \n"))
    wei_temp = fit_model(data_temp,n.iter=n.iter,column_IDs=column_IDs,model="weibull", quiet = T, ...)
    Rhats = wei_temp$BUGSoutput$summary[,"Rhat"]
    ESS = wei_temp$BUGSoutput$summary[,"n.eff"]
    n_points = wei_temp$BUGSoutput$n.keep * wei_temp$BUGSoutput$n.chains
    check_condition = all(Rhats<1.02, (ESS/n_points)>0.5)
    if(check_condition == FALSE){cat(red("\r  The model Rhat and ESS seems to be too low, check the model output and if needed increase n.iter \n"))}


    models = list(exponential = exp_temp, weibull = wei_temp)
    results[[i]] = models
    names(results)[i] = pops_order[i]}

  attr(results, "data") <- data
  attr(results, "type") <- "multiple"
  attr(results, "params") <- list(prop = prop)


  return(results)

}
