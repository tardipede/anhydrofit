#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`



# If the model has log_sh among its parameters, extract it
# Otherwise return a sequence of zeros as long as the posterior sample
.extract_log_sh = function(model){

  if("sh_log" %in% model$parameters.to.save){
    sh_chains = data.frame(do.call(rbind,(as.mcmc(model))))$sh_log}

  if(!("sh_log" %in% model$parameters.to.save)){
    n_points = model$BUGSoutput$n.keep * model$BUGSoutput$n.chains
    sh_chains = data.frame(sh_log = rep(0,n_points))}

  return(sh_chains)

}

#' @importFrom coda as.mcmc
#' @importFrom stringr str_match_all
.extract_chains_model_multiple = function(model, DIC_treshold = 1, show.table = TRUE, interval_ahi = c(0,24)){

  # Extract DICs from the models
  DIC_list = lapply(model, FUN = function(x){c(x[[1]]$BUGSoutput$DIC,x[[2]]$BUGSoutput$DIC)})
  DIC_tab = do.call(rbind,DIC_list)
  colnames(DIC_tab) = c("exponential","weibull")
  if(show.table == TRUE){.print_nice_DIC_table(DIC_tab , DIC_treshold)}

  # Choose the models to keep
  models_to_keep = (((DIC_tab[,1]-DIC_tab[,2])>DIC_treshold)*1)+1

  # Extract best models
  models_extracted = list()
  for (i in 1:length(model)){
    models_extracted[[i]] = model[[i]][[models_to_keep[i]]]}
  names(models_extracted) = names(model)


  # Extract Q and p chains
  Qs = do.call(cbind,lapply(models_extracted, FUN = function(x){data.frame(do.call(rbind,(as.mcmc(x))))$Q}))
  colnames(Qs) = names(models_extracted)

  ps = do.call(cbind,lapply(models_extracted, FUN = function(x){data.frame(do.call(rbind,(as.mcmc(x))))$p}))
  colnames(ps) = names(models_extracted)

  logshs = do.call(cbind,lapply(models_extracted, FUN = function(x){.extract_log_sh(x)}))
  colnames(logshs) = names(models_extracted)

  AhIs = do.call(cbind,lapply(models_extracted, FUN = function(x){.area_UC_model_inmultiple(x, limits = interval_ahi)}))
  colnames(AhIs) = names(models_extracted)



  results = list(Qs = Qs, ps = ps, log_shs = logshs, AhIs = AhIs)
  return(results)

}

#' @importFrom coda as.mcmc
#' @importFrom stringr str_match_all
.extract_chains_model_simple = function(model, interval_ahi = c(0,24)){

  # Extract chains
  chains = data.frame(do.call(rbind,as.mcmc(model)))

  # Extract Q and p posteriors and name them with the population names
  Qs = chains[, grep("Q", colnames(chains))]
  if(!is.null(ncol(Qs))){
    Qs = Qs[, order(as.numeric(unlist(str_match_all(colnames(Qs), "[0-9]+"))))]
    colnames(Qs) = attributes(model)$num_to_pop$population }
  if(is.null(ncol(Qs))){Qs = data.frame(Q = Qs)}

  ps = chains[, grep("p", colnames(chains))]
  if(!is.null(ncol(ps))){
    ps = ps[, order(as.numeric(unlist(str_match_all(colnames(ps), "[0-9]+"))))]
    colnames(ps) = attributes(model)$num_to_pop$population }
  if(is.null(ncol(ps))){ps = data.frame(p = ps)}

  if(attributes(model)$model == "weibull"){
    logshs = chains[, grep("sh_log", colnames(chains))]
    if(!is.null(ncol(logshs))){
      logshs = logshs[, order(as.numeric(unlist(str_match_all(colnames(logshs), "[0-9]+"))))]
      colnames(logshs) = attributes(model)$num_to_pop$population }
    if(is.null(ncol(logshs))){logshs = data.frame(p = logshs)}}

  if(attributes(model)$model == "exponential"){
    columns_n = ncol(ps)
    rows_n = nrow(ps)
    logshs = data.frame(matrix(rep(0,columns_n*rows_n), nrow = rows_n))
    colnames(logshs) = colnames(ps)}

  # Calculate AhI
  prop <- attributes(model)$params$prop
  if(length(attributes(model)$num_to_pop$population)==1){
                    temp_AhI = list(.area_UC_list(Q = Qs[,1], p = ps[,1], log_sh = logshs[,1], prop = prop, limits = interval_ahi))}

  if(length(attributes(model)$num_to_pop$population)>1){
    temp_AhI <- list()
    for (i in 1:length(attributes(model)$num_to_pop$population)) {
      temp_AhI[[i]] <- .area_UC_list(
        Q = Qs[, attributes(model)$num_to_pop$population[i]],
        p = ps[, attributes(model)$num_to_pop$population[i]],
        log_sh = logshs[, attributes(model)$num_to_pop$population[i]],
        prop = prop, limits = interval_ahi
      )
    }}
  AhIs <- data.frame(do.call(cbind, temp_AhI))
  colnames(AhIs) <- attributes(model)$num_to_pop$population


  results = list(Qs = Qs, ps = ps, log_shs = logshs, AhIs = AhIs)
  return(results)

}



## Other general helpers

#' @importFrom stats quantile sd
.my_summary_function <- function(x) {
  results <- c(mean(x, na.rm = T), sd(x, na.rm = T), quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  names(results) <- c("mean", "sd", "2.5%", "25%", "50%", "75%", "9.75%")
  return(results)
}

#' @importFrom bayestestR pd_to_p p_direction
.bayes_testing <- function(x, y) {
  pval = as.numeric(pd_to_p(p_direction(x - y)))
  return(pval)
}



### Function for 2d p-direction and effect size

## Helper functions ##
.get_angle = function(dx,dy){
  return(atan2(dy,dx) * (180 / pi)+180)}

.correct_angle = function(x){
  if(x<0){x = 360+x}
  if(x>360){x = x-360}
  return(x)}

.get_angle_in_interval = function(x,m, get.logical = FALSE){
  min_a = min(.correct_angle((m-90)),.correct_angle((m+90)))
  max_a = max(.correct_angle((m-90)),.correct_angle((m+90)))

  if(get.logical == FALSE){
    prop = mean((x > min_a) & (x < max_a))
    if(prop < 0.5){prop = 1-prop}
    return(prop)}

  if(get.logical == TRUE){
    prop = ((x > min_a) & (x < max_a))
    if(mean(prop) < 0.5){prop = !prop}
    return(prop)}
}


.es_2d = function(group1, group2, Qnorm){

  # Calculate distances between groups centroids
  diff_c1 = mean(group1[,1],na.rm=T) - mean(group2[,1],na.rm=T)
  diff_c2 = mean((group1[,2]/Qnorm),na.rm=T) - mean((group2[,2]/Qnorm), na.rm=T)
  dist_c = sqrt(diff_c1^2 + diff_c2^2)
  return(dist_c)
}


#' @importFrom bayestestR area_under_curve
.area_UC <- function(Q, p, log_sh, limits = c(0, 24), step = 0.1, prop = 0.9) {

  x <- seq(limits[1], limits[2], by = step)
  lambda <- ((-log(1 - prop))^(1/exp(log_sh))) / Q
  y <- (1 - exp(-(lambda * x)^exp(log_sh))) * p
  AUC_abs <- area_under_curve(x, y, method = c("trapezoid"))
  AUC_max <- diff(limits)
  AUC_rel <- AUC_abs / AUC_max

  return(AUC_rel)
}

.area_UC_list <- function(Q, p, log_sh, prop, limits = c(0, 24)) {
  data <- data.frame(Q = Q, p = p, log_sh = log_sh)
  results <- as.numeric(unlist(apply(data, MARGIN = 1, FUN = function(x) {
    .area_UC(Q = x[1], p = x[2], log_sh = x[3], prop = prop, limits = limits)
  })))
  return(results)
}


.area_UC_model_inmultiple = function(model, limits){
  prop = attributes(model)$params$prop
  Qs = data.frame(do.call(rbind,(as.mcmc(model))))$Q
  ps = data.frame(do.call(rbind,(as.mcmc(model))))$p
  sh_log = .extract_log_sh(model)

  UC = .area_UC_list(Q = Qs, p = ps, log_sh = sh_log, prop = prop, limits = limits)
  return(UC)

}



######

#Print a coloured DIC summary table
#' @importFrom knitr kable
#' @importFrom crayon blue
.print_nice_DIC_table = function(DIC_tab, DIC_treshold){
  DIC_tab = data.frame(DIC_tab)
  DIC_tab$delta_DIC = DIC_tab$exponential-DIC_tab$weibull
  DIC_tab$model = c("Exponential", "Weibull")[(((DIC_tab[,1]-DIC_tab[,2])>DIC_treshold)*1)+1]
  colnames(DIC_tab) = c("DIC Exponential", "DIC Weibull", "delta DIC", "Model selected")
  cat(c("",blue(paste(gsub(":", "-", kable(DIC_tab)), "\n"))))     }


# originally from aplpack package, plotting functions removed
#' @importFrom grDevices chull
#' @importFrom stats aggregate median setNames
.plothulls_ <- function(x, y, fraction, n.hull = 1,
                        col.hull, lty.hull, lwd.hull, density=0, ...){
  # function for data peeling:
  # x,y : data
  # fraction.in.inner.hull : max percentage of points within the hull to be drawn
  # n.hull : number of hulls to be plotted (if there is no fractiion argument)
  # col.hull, lty.hull, lwd.hull : style of hull line
  # plotting bits have been removed, BM 160321
  # pw 130524
  if(ncol(x) == 2){ y <- x[,2]; x <- x[,1] }
  n <- length(x)
  if(!missing(fraction)) { # find special hull
    n.hull <- 1
    if(missing(col.hull)) col.hull <- 1
    if(missing(lty.hull)) lty.hull <- 1
    if(missing(lwd.hull)) lwd.hull <- 1
    x.old <- x; y.old <- y
    idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
    for( i in 1:(length(x)/3)){
      x <- x[-idx]; y <- y[-idx]
      if( (length(x)/n) < fraction ){
        return(cbind(x.hull,y.hull))
      }
      idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx];
    }
  }
  if(missing(col.hull)) col.hull <- 1:n.hull
  if(length(col.hull)) col.hull <- rep(col.hull,n.hull)
  if(missing(lty.hull)) lty.hull <- 1:n.hull
  if(length(lty.hull)) lty.hull <- rep(lty.hull,n.hull)
  if(missing(lwd.hull)) lwd.hull <- 1
  if(length(lwd.hull)) lwd.hull <- rep(lwd.hull,n.hull)
  result <- NULL
  for( i in 1:n.hull){
    idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
    result <- c(result, list( cbind(x.hull,y.hull) ))
    x <- x[-idx]; y <- y[-idx]
    if(0 == length(x)) return(result)
  }
  result
} # end of definition of plothulls
#################################


# Make hulls dataframe for plotting
.get_hulls = function(data, fraction){

  group = NULL

  results = list()
  for(i in 1:length(unique(data$group))){

    data_temp = subset(data, group == unique(data$group)[i])
    the_matrix = matrix(data = c(data_temp$p, data_temp$Q), ncol = 2)
    hull_temp = setNames(data.frame(.plothulls_(the_matrix, fraction = fraction)), nm = c("p", "Q"))
    hull_temp$group = rep(unique(data$group)[i], nrow(hull_temp))
    results[[i]] = hull_temp

  }

  hull_data = do.call(rbind,results)
  return(hull_data)

}

