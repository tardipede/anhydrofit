#' Make a plot showing predicted values along with the input data
#'
#' @importFrom ggplot2 ggplot theme_bw geom_ribbon geom_point aes facet_wrap theme element_blank
#' @importFrom scales alpha
#' @importFrom coda as.mcmc
#'
#' @param model the output model from \code{fit_model} or \code{fit_multiple_models}
#' @return a saved pdf file named model_fits.pdf
#' @export


plot_fits = function(model){
  if(attributes(model)$type == "simple"){chains = .extract_chains_model_simple(model)}
  if(attributes(model)$type == "multiple"){chains = .extract_chains_model_multiple(model, show.table = FALSE)}
  ps = data.frame(chains$ps)
  Qs = data.frame(chains$Qs)
  shs = exp(data.frame(chains$log_shs))

  prop = attributes(model)$params$prop

  data = attributes(model)$data
  data = data[,c("population","time","moving","total")]

  max_time = max(data$time)*1.1

  # Format the population names to avoid issues
  data$population = gsub(" ", "_", data$population)
  data$population = gsub("[.]", "_", data$population)
  data$population = gsub("[-]", "_", data$population)
  data$prob = data$moving/data$total

  data$low_ci = .wilson_score_interval(s = data$moving, n = data$total, type="low")
  data$high_ci = .wilson_score_interval(s = data$moving, n = data$total, type="high")

  colnames(ps) = gsub(" ", "_", colnames(ps))
  colnames(ps) = gsub("[.]", "_", colnames(ps))
  colnames(ps) = gsub("[-]", "_", colnames(ps))

  colnames(Qs) = gsub(" ", "_", colnames(Qs))
  colnames(Qs) = gsub("[.]", "_", colnames(Qs))
  colnames(Qs) = gsub("[-]", "_", colnames(Qs))

  colnames(shs) = gsub(" ", "_", colnames(shs))
  colnames(shs) = gsub("[.]", "_", colnames(shs))
  colnames(shs) = gsub("[-]", "_", colnames(shs))

  results = list()
  for(i in 1:ncol(ps)){

    ps_temp = ps[,i]
    Qs_temp = Qs[,i]
    sh_temp = shs[,i]

    x = seq(0,max_time,length.out = 500)
    predicted_temp = matrix(, nrow = length(ps_temp), ncol = length(x))
    for(j in 1:length(ps_temp)){
      lambda = ((-log(1 - prop))^(1/sh_temp[j])) / Qs_temp[j]
      y = (1 - exp(-(lambda * x)^sh_temp[j])) * ps_temp[j]
      predicted_temp[j,] = y
    }

    df_temp = data.frame(t(apply(predicted_temp, MARGIN = 2, FUN = function(x){quantile(x, probs = c(0.027, 0.975), na.rm=T)})))
    df_temp$time = x
    df_temp$population = rep(colnames(ps)[i], nrow(df_temp))
    results[[i]] = df_temp}

    plot_data = do.call(rbind,results)
    colnames(plot_data)[1:2] = c("low","high")

    ggplot(plot_data)+
      theme_bw()+
      geom_ribbon(aes(x=time, ymin=low, ymax=high), fill = "#60effc", col="#60effc",alpha=0.5)+
      facet_wrap(.~population)+
      geom_linerange(data = data, aes(x = time, ymin = low_ci, ymax = high_ci), col="#09919e")+
      geom_point(data=data, aes(x=time, y=prob), shape = 21, size=3, fill="#60effc",alpha=0.5)+
      theme(panel.grid.minor = element_blank())
}
