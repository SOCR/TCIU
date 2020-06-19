#' @title changetitle
#' @description  Reshape the title for plots
#' @details This funciton reshapes a vector of characters inside. This vector often serves as a list of titles
#' for different plots. The aim for this function is to narrow down the length of each element in this vector
#' so that the name will fit into plots. This function will put a line break(i.e. "\\n") each time after counting
#' 25 characters.
#'
#' @param titl vector contains titles that need to be reshaped
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @examples
#' head(dataColnames)
#' newtitles<-changetitle(dataColnames)
#' head(newtitles)
#'
#' @export
#'
changetitle <- function(titl) {
  newtitle<-rep(NA,length(titl))
  for (i in 1:length(titl)) {
    tempchar<-substring(titl[i],1:nchar(titl[i]),1:nchar(titl[i]))
    for (j in 1:8) {
      k=25*j
      tempchar[which(tempchar==" ")[which(tempchar==" ")>k][1]]<-"\n"
    }
    tempchar<-paste(tempchar,collapse = "")
    newtitle[i]<-tempchar
  }
  return(newtitle)
}
