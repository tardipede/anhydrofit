#' Analyse fitted models from \code{fit_model} and \code{fit_multiple_models}
#'
#' @param model the output model from \code{fit_model} or \code{fit_multiple_models}
#' @param pairwise.test if test the pairwise differences between each group
#' @param Q.norm.factor normalization factor for calculating effect magnitide between treatment, if NA is calculated automatically ad 99 percentile of all Qs. Default: NA.
#' @param interval_ahi vector of 2 indicating the time interval to calculate AhI. Default: c(0,24)
#' @param DIC_treshold minimum difference needed in DIC to select the more complex model (Weibull). Default: 1.
#' @param get.chains get as output the parameter chains extracted from model. Defauls: FALSE.
#' @export
#'
#' @examples
#' \dontrun{
#' data = data.frame(group = c("A","A","A","B","B","B"),
#'        time = c(1,2,24,1,2,24),
#'        moving = c(0,1,4,2,5,16),
#'        tot = rep(20,6))
#' mod1 = fit_model(data, model = "weibull")
#' analyse_model(mod1)
#'
#' mod2 = fit_multiple_models(data)
#' analyse_model(mod2)}


analyse_model = function(model,
                         pairwise.test = TRUE,
                         Q.norm.factor = NA,
                         interval_ahi = c(0,24),
                         DIC_treshold = 1,
                         get.chains = FALSE){
  requireNamespace("R2jags")
  requireNamespace("bayestestR")

  # Extract chains
  if(attributes(model)$type == "simple"){chains = .extract_chains_model_simple(model, interval_ahi = interval_ahi)}
  if(attributes(model)$type == "multiple"){chains = .extract_chains_model_multiple(model,DIC_treshold = DIC_treshold,
                                                                                   interval_ahi = interval_ahi)}

  ps = chains$ps
  Qs = chains$Qs
  log_shs = chains$log_shs
  AhIs = chains$AhIs


  # Calculate AhI
  prop = attributes(model)$params$prop

  # Calculate summary tables
  p_summary = t(apply(ps, MARGIN = 2, FUN = .my_summary_function))
  colnames(p_summary) = paste0("p ", colnames(p_summary))

  Q_summary = t(apply(Qs, MARGIN = 2, FUN = .my_summary_function))
  colnames(Q_summary) = paste0("Q ", colnames(Q_summary))

  AhI_summary = t(apply(AhIs, MARGIN = 2, FUN = .my_summary_function))
  colnames(AhI_summary) = paste0("AhI ", colnames(AhI_summary))



  summary_table = cbind(p_summary,Q_summary, AhI_summary)
  results = list(summary = summary_table)


  if(pairwise.test == TRUE){

       # Create a table with the possible pairwise combinations
       populations = colnames(ps)
       pairwise.table = expand.grid(pop1 = populations, pop2 = populations)
       pairwise.table = pairwise.table[pairwise.table$pop1 != pairwise.table$pop2, ]
       pairwise.table = data.frame(t(apply(pairwise.table, MARGIN = 1, FUN = sort)))
       pairwise.table = pairwise.table[duplicated(pairwise.table), ]

       pairwise.table = data.frame( group1 = pairwise.table[,2],
                                    group2 = pairwise.table[,1],
                                    p_pval = rep(NA, nrow(pairwise.table)),
                                    Q_pval = rep(NA, nrow(pairwise.table)),
                                    AhI_pval = rep(NA, nrow(pairwise.table)),
                                    biv_pval = rep(NA, nrow(pairwise.table)),
                                    biv_es = rep(NA, nrow(pairwise.table)))

       # Get Qnorm
       if(is.na(Q.norm.factor)){Q.norm.factor = as.numeric(quantile(unlist(Qs),prob=c(0.99)))}

       # Fill the table
       for (i in 1:nrow(pairwise.table)) {
         pairwise.table[i, 3] = .bayes_testing(ps[, pairwise.table[i, 1]], ps[, pairwise.table[i, 2]])
         pairwise.table[i, 4] = .bayes_testing(Qs[, pairwise.table[i, 1]], Qs[, pairwise.table[i, 2]])
         pairwise.table[i, 5] = .bayes_testing(AhIs[, pairwise.table[i, 1]], AhIs[, pairwise.table[i, 2]])
         pairwise.table[i, 6] <- bayestestR::pd_to_p(p_direction_2d(cbind(ps[, pairwise.table[i, 1]],Qs[, pairwise.table[i, 1]]),
                                                        cbind(ps[, pairwise.table[i, 2]],Qs[, pairwise.table[i, 2]])))
         pairwise.table[i, 7] <- .es_2d(cbind(ps[, pairwise.table[i, 1]],Qs[, pairwise.table[i, 1]]),
                                        cbind(ps[, pairwise.table[i, 2]],Qs[, pairwise.table[i, 2]]),
                                   Qnorm = Q.norm.factor)}
     results[[2]] = pairwise.table
     names(results)[2] = "pairwise.table"
     attr(results, "Qnorm") = Q.norm.factor}

  if(get.chains == FALSE){return(results)}
  if(get.chains == TRUE){return(chains)}

}


#' Plot p and Q from fitted models from \code{fit_model} and \code{fit_multiple_models}
#'
#' @importFrom ggplot2 ggplot aes theme_bw  scale_y_reverse theme xlab ylab scale_fill_distiller geom_polygon element_blank geom_point element_blank median_hilow
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by summarise
#'
#' @param model the output model from \code{fit_model} or \code{fit_multiple_models}
#' @param interval_ahi vector of 2 indicating the time interval to calculate AhI. Default: c(0,24)
#' @param DIC_treshold minimum difference needed in DIC to select the more complex model (Weibull). Default: 1.
#' @export

plot_model = function(model, interval_ahi = c(0,24), DIC_treshold = 1){

  group = median_AhI = p = Q = x = y = NA

  # Extract chains
  if(attributes(model)$type == "simple"){chains = .extract_chains_model_simple(model, interval_ahi = interval_ahi)}
  if(attributes(model)$type == "multiple"){chains = .extract_chains_model_multiple(model,
                                                                                   interval_ahi = interval_ahi,
                                                                                   DIC_treshold = DIC_treshold,
                                                                                   show.table = FALSE)}
  ps = data.frame(chains$ps)
  Qs = data.frame(chains$Qs)
  AhI = data.frame(chains$AhIs)

  # Prepare the data for plotting
  plot_data_p = gather(ps,"group","p")
  plot_data_Q = gather(Qs,"group","Q")
  plot_data_AhI = gather(AhI,"group","AhI")
  plot_data = cbind(plot_data_p,plot_data_Q[,2],plot_data_AhI[,2])
  colnames(plot_data)[3:4] = c("Q","AhI")


  # Aggregate treatments by their median AhI
  data_ggregated = aggregate(plot_data$AhI,by = list(plot_data$group), FUN = median)
  colnames(data_ggregated) = c("group", "median_AhI")

  # Generate hull shapes
  hull_data = .get_hulls(plot_data, fraction = 0.8)


  # Add median AhI column for fill
  hull_data = merge(hull_data,data_ggregated )

  # Get summary to plot mean points
  summary_table = plot_data %>% group_by(group) %>%
    summarise(x = mean(p), y = median(Q))




 # Make the plot

  ggplot(hull_data, aes(x=p, y=Q))+
    theme_bw()+
    geom_polygon(aes(fill=median_AhI, group = group),alpha = 0.20, col=NA, show.legend = FALSE)+
    scale_y_reverse()+
    geom_point(data = summary_table,
               aes(x=x, y=y))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_fill_distiller(palette = "Spectral", direction = 1 )+
    geom_text_repel(data = summary_table,
                    aes(x = x, y = y, label = group),
                    min.segment.length = 0, size = 4,
                    segment.colour = "black", force_pull = 1, box.padding = 1, max.overlaps =100)+
    xlab("Survival (p)") + ylab("Recovery speed (Q)")



}

#' Plot AhI violin plots from fitted models from \code{fit_model} and \code{fit_multiple_models}
#'
#' @importFrom ggplot2 ggplot aes theme_bw  scale_y_continuous theme xlab ylab scale_fill_distiller geom_violin stat_summary
#' @importFrom tidyr gather
#'
#' @param model the output model from \code{fit_model} or \code{fit_multiple_models}
#' @param interval_ahi vector of 2 indicating the time interval to calculate AhI. Default: c(0,24)
#' @param DIC_treshold minimum difference needed in DIC to select the more complex model (Weibull). Default: 1.
#' @export

plot_ahi = function(model, interval_ahi = c(0,24), DIC_treshold = 1){
  group = median_AhI = NA

  requireNamespace("tidyverse")

  # Extract chains
  if(attributes(model)$type == "simple"){chains = .extract_chains_model_simple(model)}
  if(attributes(model)$type == "multiple"){chains = .extract_chains_model_multiple(model, show.table = FALSE)}
  AhI = chains$AhIs

  plot_data = gather(data.frame(AhI),"group","AhI")

  # Aggregate treatments by their median AhI
  data_ggregated = aggregate(plot_data$AhI,by = list(plot_data$group), FUN = median)
  colnames(data_ggregated) = c("group", "median_AhI")

  # relevel the group factor
  plot_data$group = factor(plot_data$group, levels = data_ggregated[order(data_ggregated$median_AhI, decreasing = T),1])

  # Add median AhI column for fill
  plot_data = merge(plot_data,data_ggregated )

  # Make the plot
  ggplot(plot_data)+
    theme_bw()+
    geom_violin(aes(x=group, y=AhI, fill=median_AhI), scale = "width", alpha=0.5, show.legend = FALSE, col=NA)+
    scale_fill_distiller(palette = "Spectral", direction = 1 )+
    stat_summary(aes(x=group, y=AhI),fun=median, colour="black", geom="point", size = 3)+
    stat_summary(aes(x=group, y=AhI),fun.data=median_hilow, colour="black", geom="linerange")+
    scale_y_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))+
    theme(panel.grid.minor = element_blank())+
    ylab("Anhydrobiosis Index") + xlab("")
}


