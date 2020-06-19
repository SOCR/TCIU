#' @title plotPhaseDistributions
#' @description  Creating plot result for K-space transformation
#'
#' @details This function provides a visualization method to show results of K-space transformation
#' This function can choose to generate all plots related to the K-space transformation.
#' Or we could choose features code to select features we like to see to generate K-space transformation result.
#'
#' @param dataFT data inherited from function "kSpaceTransform", K-space transformed dataset.
#' @param dataColnames feature names for the transformed dataset
#' @param size size of the text inside the plot
#' @param option two options avaliable: "All": to plot for all features, "Select": to select several features you wish to show
#' @param select_feature if "option" is setted as "Select", then put inside a sequence to indicate the features you wish to show with this function.
#'
#' @return This function returns a plot based on ggplot2 showing results of K-space transformation
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' require(ggplot2)
#' require(dplyr)
#' require(tidyr)
#' 
#' # Showing all results based on K-space transformation
#' ## Make sure to change the height and width of your plot screen so that the plot can fit your monitor
#' plotPhaseDistributions(FT_Belgium, dataColnames)
#' 
#' # Only show a few features from the results. Here we choose to generate the first 6 results
#' plotPhaseDistributions(FT_Belgium, dataColnames,option = "Select",size = 10,select_feature = c(1,2,3,4,5,6))
#' 
plotPhaseDistributions <- function (dataFT, dataColnames, size=10, option="ALL",select_feature=NULL,...) {
  df.phase <- as.data.frame(Re(dataFT$phases))
  df.phase %>% gather() %>% head()
  colnames(df.phase) <- dataColnames
  phaseDistributions <- gather(df.phase)
  colnames(phaseDistributions) <- c("Feature", "Phase")
  if (is.null(size)) size=10

  # map the value as our x variable, and use facet_wrap to separate by the key column:
  if(option=="ALL"){
    ggplot(phaseDistributions, aes(Phase)) +
      # geom_histogram(bins = 10) +
      geom_histogram(aes(y=..density..), bins = 10) +
      facet_wrap( ~Feature, scales = 'free_x') +
      xlim(-pi, pi) +
      theme(strip.text.x = element_text(size = size, colour = "black", angle = 0))}
  else if(option=="Select"){
    choosefeat<-dataColnames[select_feature]
    NphaseDistributions<-phaseDistributions[which(phaseDistributions$Feature %in% choosefeat),]
    ggplot(NphaseDistributions, aes(Phase)) +
      # geom_histogram(bins = 10) +
      geom_histogram(aes(y=..density..), bins = 10) +
      facet_wrap( ~Feature, scales = 'free_x') +
      xlim(-pi, pi) +
      theme(strip.text.x = element_text(size = size, colour = "black", angle = 0))}
  else{
    return("Wrong input on option, Please try again")
  }
}
