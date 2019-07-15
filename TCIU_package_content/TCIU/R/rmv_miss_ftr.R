#' @title rmv_miss_ftr
#' @description Removes features with all missing values
#'
#' @details This function deletes features that are missing for all time points
#' associated with a certain country.
#' 
#' @param countryName give a country name that is shown on the original dataset
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' require(dplyr)
#' 
#' #take "Sweden" as an example
#' head(rmv_miss_ftr("Sweden"))[,1:10]
rmv_miss_ftr = function(countryName = "Bulgaria"){
  DataSuperSample =
    time_series %>%
    filter(country == countryName) %>%
    select_if( function(col) sum(is.na(col)) != length(col) ) %>%
    select(-time, -country) %>%
    as.matrix() %>%
    cleardata() %>%
    as.data.frame()%>%
    # remove feature that only has one value
    select_if(function(col) sum(is.na(col)) == 0 ) %>%
    # remove feature that all the values are 0
    select_if(function(col) sum(col) != 0 ) %>%
    splinecreate() %>%
    as.data.frame()
  return(DataSuperSample)
}

##Concern: The name of country "German" has changed inside the dataset once.
