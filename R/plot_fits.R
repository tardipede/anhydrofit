#' Save as pdf file in the working directory with plots showing the model fits
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics layout lines par points
#' @importFrom scales alpha
#'
#' @param model the output model from \code{fit_model} or \code{fit_multiple_models}
#' @return a saved pdf file mamed model_fits.pdf
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
  data$population = gsub(" ", ".", data$population)
  data$prob = data$moving/data$total

  pdf("model_fits.pdf", width = 8.27, height = 11.69)
  group <- "Title"
  layout(matrix(c(1:12), 6, 2) )
  par(mar = c(2.5, 4.1,2.5, 2.1),oma = c(2, 2, 2, 2))

  for(i in 1:ncol(ps)){
    data_temp = data[data$population==colnames(ps)[i],]

    sample_post = sample(1:nrow(ps), size = 100)

    ps_temp = ps[sample_post,i]
    Qs_temp = Qs[sample_post,i]
    sh_temp = shs[sample_post,i]

    plot(x = c(0,50), y=c(0,1), type="n", xlab = "Time", ylab = "p", main = colnames(ps)[i])

    for(j in 1:length(sample_post)){

      x = seq(0,50,0.1)
      lambda = ((-log(1 - prop))^(1/sh_temp[j])) / Qs_temp[j]
      y = (1 - exp(-(lambda * x)^sh_temp[j])) * ps_temp[j]
      lines(x=x, y=y, type="l", col= alpha("grey",0.35))
  }

    points(x = data_temp$time, y = data_temp$prob, col = alpha("red", 0.70))



  }

  dev.off()
}
