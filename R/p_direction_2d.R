#' Calculated p_direction of differences in angular direction between two groups
#'
#' @importFrom circular circular plot.circular arrows.circular
#'
#' @param group1 a matrix with 2 columns with x and y coordinates
#' @param group2 a matrix with 2 columns with x and y coordinates
#' @param make.plot make a direction plot of differences between groups. Defauls: FALSE.

#' @return a p-direction value
#' @export


p_direction_2d = function(group1, group2, make.plot = FALSE){

  # Calcuate difference between groups
  diff_var1 = group1[,1] - group2[,1]
  diff_var2 = group1[,2] - group2[,2]

  # Standardize the differences
  diff_var1_s = diff_var1/sd(diff_var1)
  diff_var2_s = diff_var2/sd(diff_var2)

  # Get angles
  diff_angles = .get_angle(diff_var2_s, diff_var1_s)

  # Get median angle
  diff_circ = circular(diff_angles, units="degrees", template="none")
  median_angle = median(diff_circ, na.rm = T)

  # Calculate p-direction
  p_dir = .get_angle_in_interval(diff_circ,median_angle)

  # Plot angles
  if(make.plot == TRUE){
    plot.circular(diff_circ, stack = T, pch = 20, sep = 0.1, shrink = 1.8, bins = 360,
                  main = paste0("p-direction = ", round(p_dir, 3)))

    arrows.circular(median_angle)
    arrows.circular(median_angle-90, col="red", length = 0)
    arrows.circular(median_angle+90, col="red", length = 0)}

  return(p_dir)
}
