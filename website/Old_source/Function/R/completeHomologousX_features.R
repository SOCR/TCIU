#' @title completeHomologousX_features
#' @description Ensure the XReg predictors are homologous
#' @details This function ensures the XReg predictors for ALL 31 EU countries are homologous
#' \strong{Notice:} This function mainly serves as a functional part for other functions. So we
#' seldomly apply this function by itself.
#'
#' @param list_of_dfs a list of datasets of different countries
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#'
#' @export
#'
completeHomologousX_features <- function (list_of_dfs) {
  # delete features that are missing at all time points
  for (j in 1:length(list_of_dfs)) {
    print(paste0("Pre-processing Country: ...", countryNames[j], "... "))
    data = list_of_dfs[[j]]
    data = data[ , colSums(is.na(data)) != nrow(data)]%>%
      select(-time, -country)%>%
      as.matrix()%>%
      cleardata()->DataMatrix
    DataMatrix = DataMatrix[ , colSums(is.na(DataMatrix)) == 0] # remove features with only 1 value
    DataMatrix = DataMatrix[ , colSums(DataMatrix) != 0] %>% # remove features with all values=0
      # Supersample 72 --*5--> 360 timepoints
      splinecreate()%>%
      as.data.frame()->DataSuperSample # super-Sample the data
    # remove some of features
    DataSuperSample = DataSuperSample[, -c(50:80)]; dim(X)  # 360 167
    # ensure full-rank design matrix, DataSuperSample
    DataSuperSample <-
      DataSuperSample[ , qr(DataSuperSample)$pivot[seq_len(qr(DataSuperSample)$rank)]]
    print(paste0("dim()=(", dim(DataSuperSample)[1], ",", dim(DataSuperSample)[2], ") ..."))
    # update the current DF/Country
    list_of_dfs_CommonFeatures[[j]] <- DataSuperSample
  }

  # Identify All Xreg features that are homologous (same feature columns) across All 31 countries
  # Identify Common Columns (freatures)
  comCol <- Reduce(intersect, lapply(list_of_dfs_CommonFeatures, colnames))
  list_of_dfs_CommonFeatures <- lapply(list_of_dfs_CommonFeatures, function(x) x[comCol])

  for (j in 1:length(list_of_dfs_CommonFeatures)) {
    list_of_dfs_CommonFeatures[[j]] <- subset(list_of_dfs_CommonFeatures[[j]], select = comCol)
    print(paste0("dim(", countryNames[j], ")=(", dim(list_of_dfs_CommonFeatures[[j]])[1],
                 ",", dim(list_of_dfs_CommonFeatures[[j]])[2], ")!"))  # 72 * 197
  }
  return(list_of_dfs_CommonFeatures)
}
