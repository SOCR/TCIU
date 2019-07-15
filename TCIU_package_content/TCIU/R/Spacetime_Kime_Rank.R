#' @title Spacetime_Kime_Rank
#' @description Generate numerous LASSO regression analysis results based on 1. "Weak" or "Strong" signals 2. different
#' Space Kime transformation methods
#' 
#' @details This function \code{Spacetime_Kime_Rank} is a function that takes in a data frame X and a vector Y. X is the 
#' data fram for all predicor variables and Y is the vector for the response variable. Then the function will base on two
#' other parameters "signal" and "phase" to determine which features will be used on creating ARIMAX models and also 
#' whether SpaceKime transformation methods will be applied on X and Y. It will generate numerous results. And the result
#' that is returned eventually is deterimined on parameter "result". More details can be viewed on the parameters part.
#' 
#' @param X a data frame X containing all predictors. X needs to contain features both from ARIMA time series model and
#' also features for the external regressors(i.e. the "xreg" part).
#' @param Y a vector Y containing information of the response result. 
#' @param signal two options for this parameter: "Weak" and "Strong". Set to "weak" will analyse the model only based on 
#' ARIMA features (i.e. features with AR and MA parts). Set to "strong" will analyse model based on all features
#' @param phase three options: "Normal", "Nil" and "Swap". Using "Normal" will not conduct SpaceKime transformation.
#' "Nil" will conduct "Nil-Phase" transformation and "Swap" will conduct "Swapped-Phase" transformation on the data.
#' @param result There are generally 5 choices to show 5 different results. If result=1, then the MSE plot of LASSO
#' cross-validation will be shown. If result=2 then the coefficient of the top 9 sailent features under LASSO penalty
#' will be shown. If result=3 then the ranking comparison with the observed value against the predicted value (under
#' different SpaceKime transformation results) will be shown. Also the correlation will also be shown. If result=4 a plot
#' similiar to the result of 3 will be shown. But this time the plot will also include the information about whether a 
#' country is top 30 country or not. And the name of different countries will also be shown on this plot to give readers
#' a much clearer view. If result=5, then a list will be shown, this list contains dimension of predictors used in our model,
#' information about first 9 highest features, correlation between predicted and true results, all features importance
#' comparison, and other results based on the choice of another parameter "lam".
#' Moreover, there are another two options, result=6 and result=7. These two options will only work when we choose phase
#' parameter as "Nil". Result=6 will show a plot based on phase distribution.
#'  
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
Spacetime_Kime_Rank <- function(X,Y,signal="Weak",phase="Normal",result=1,Seed=4321,lam="min") {
  X = X[ , colSums(is.na(X)) == 0]
  if(signal=="Weak"){
    X<-X[,1:378]
  }
  if(signal=="Strong"){
    X<-X[,1:386]
  }
  if(phase=="Nil"){
    FT_aggregate_arima_vector_country_ranking_df <-
      kSpaceTransform(aggregate_arima_vector_country_ranking_df, inverse = FALSE, reconPhases = NULL)
    if(result==6){
      return(plotPhaseDistributions(FT_aggregate_arima_vector_country_ranking_df,
                                    colnames(aggregate_arima_vector_country_ranking_df), size=4, cex=0.1))
    }
    IFT_FT_aggregate_arima_vector_country_ranking_df <-
      kSpaceTransform(FT_aggregate_arima_vector_country_ranking_df$magnitudes,
                      TRUE, FT_aggregate_arima_vector_country_ranking_df$phases)
    if(result==7){
      return(PS_ARIMAX_remodel(aggregate_arima_vector_country_ranking_df[,1:ncol(aggregate_arima_vector_country_ranking_df)-1],
                               Y,NULL,result = "reconstruct",rename_Y = "Yvec")[1:4])
    }
    IFT_NilPhase_FT_aggregate_arima_vector<-PS_ARIMAX_remodel(aggregate_arima_vector_country_ranking_df[,
                                                                                                        1:ncol(aggregate_arima_vector_country_ranking_df)-1],Y,NULL,result = "reconstruct",rename_Y = "Yvec")[[6]]
    rownames(IFT_NilPhase_FT_aggregate_arima_vector)<-rownames(X)
  }
  if(phase=="Swap"){
    FT_aggregate_arima_vector_country_ranking_df <-
      kSpaceTransform(aggregate_arima_vector_country_ranking_df, inverse = FALSE, reconPhases = NULL)
    swappedPhase_FT_aggregate_arima_vector <- FT_aggregate_arima_vector_country_ranking_df$phases
    IFT_SwappedPhase_FT_aggregate_arima_vector <- array(complex(), c(dim(temp_Data)[1], dim(temp_Data)[2]))
    set.seed(Seed)   # sample randomly Phase-columns for each of the 131 covariates (X)
    swappedPhase_FT_aggregate_arima_vector1 <- as.data.frame(cbind(
      swappedPhase_FT_aggregate_arima_vector[ ,
                                              sample(ncol(swappedPhase_FT_aggregate_arima_vector[ , 1:378]))],  # mix ARIMA signature phases
      swappedPhase_FT_aggregate_arima_vector[ ,
                                              sample(ncol(swappedPhase_FT_aggregate_arima_vector[ , 379:386]))],# mix the meta-data phases
      swappedPhase_FT_aggregate_arima_vector[ , 387]))                          # add correct Outcome phase
    swappedPhase_FT_aggregate_arima_vector <- swappedPhase_FT_aggregate_arima_vector1

    colnames(swappedPhase_FT_aggregate_arima_vector) <- colnames(temp_Data)
  }

  text0<-dim(X)
  set.seed(Seed)
  if(phase=="Normal"){
    Xmat<-X
    CVlassotext<-c("CV LASSO (using only Timeseries data): Number of Nonzero (Active) Coefficients")
    result2text<-c("Top 9 salient features (LASSO penalty)")
  }
  if(phase=="Nil"){
    Xmat<-IFT_NilPhase_FT_aggregate_arima_vector
    CVlassotext<-c("(Spacekime, Nil-phase) CV LASSO: Number of Nonzero (Active) Coefficients")
    result2text<-c("(Spacekime) Top 9 salient features (LASSO penalty)")
  }
  if(phase=="Swap"){
    Xmat<-swappedPhase_FT_aggregate_arima_vector
    CVlassotext<-c("(Spacekime, Swapped-Phases) CV LASSO: Number of Nonzero (Active) Coefficients")
    result2text<-c("(Spacekime, Swapped-Phases) Top 9 salient features (LASSO penalty)")
  }
  CVlasso<-cv.glmnet(data.matrix(Xmat[ , 1:ncol(Xmat)]), Y, alpha = 1, parallel=TRUE)
  if(result==1){
    plot(CVlasso)
    mtext(CVlassotext, side=3, line=2.5)
  }
  PREDlasso<-predict(CVlasso,s=ifelse(lam=="min",CVlasso$lambda.min,lam),newx=data.matrix(Xmat[,1:ncol(Xmat)]))
  Coeflist<-coef(CVlasso,s=ifelse(lam=="min",CVlasso$lambda.min,lam))
  Coeflist <- data.frame(Coeflist@Dimnames[[1]][Coeflist@i+1],Coeflist@x)
  names(Coeflist) <- c('Feature','EffectSize')
  text1<-arrange(Coeflist, -abs(EffectSize))[2:10, ] ##This place Normal and Nil phase are not the same. See through it.
  text2<-cor(Y, predLASSO_lim[, 1])
  text3<-varImp(CVlasso, lambda = ifelse(lam=="min",CVlasso$lambda.min,lam))
  if(result==2){
    COEFplot<-coef(CVlasso, s = ifelse(lam=="min",CVlasso$lambda.min,lam)) %>%  # "lambda.min"
      broom::tidy() %>%
      filter(row != "(Intercept)") %>%
      top_n(100, wt = abs(value)) %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE, aes(size = abs(value))) +
      ggtitle(result2text) +
      xlab("Effect-size") +
      ylab(NULL)
    return(COEFplot)
  }
  if(phase=="Normal"){
    validation_lim <- data.frame(matrix(NA, nrow = dim(PREDlasso)[1], ncol=2), row.names=countryNames)
    validation_lim [ , 1] <- Y; validation_lim[ , 2] <- PREDlasso[, 1]
    colnames(validation_lim) <- c("Orig_Y", "LASSO")
    NN<-3
    LMylab<-c("LASSO (42*9 +8) param model")
    LMmain<-c("Observed (X) vs. LASSO-Predicted (Y) Overall Country Ranking, cor=%.02f")
    Xin_result4<-validation_lim[ , 1]
    Yin_result4<-validation_lim[ , 2]
    result4title<-c("Spacetime LASSO Predicted (y) vs. Observed (x) Overall Country Ranking, cor=%.02f")
    result4Y_T<-c( "Spacetime LASSO Predicted")
  }
  if(phase=="Nil"){
    validation_lim <- cbind(PREDlasso[, 1],
                            Xmat[ , 387], Y)
    colnames(validation_lim) <- c("predLASSO_kime", "IFT_NilPhase_FT_Y", "Orig_Y")
    rownames(validation_lim)[11] <- "Germany"
    NN<-4
    LMylab<-c("IFT_NilPhase predLASSO_kime")
    LMmain<-c("Observed (x) vs. IFT_NilPhase Predicted (y) Overall Country Ranking, cor=%.02f")
    Xin_result4<-Y
    Yin_result4<-PREDlasso
    result4title<-c("Spacekime LASSO Predicted, Nil-Phases (y) vs. Observed (x) Overall Country Ranking, cor=%.02f")
    result4Y_T<-c("Spacekime LASSO Predicted, using Nil-Phases")
  }
  if(phase=="Swap"){
    set.seed(Seed)
    cvLASSO_lim = cv.glmnet(data.matrix(X[ , 1:ncol(X)]), Y, alpha = 1, parallel=TRUE)
    predLASSO_lim <-  predict(cvLASSO_lim, s = ifelse(lam=="min",cvLASSO_lim$lambda.min,lam), # cvLASSO_lim$lambda.min,
                              newx = data.matrix(X[ , 1:ncol(X)]))
    IFT_NilPhase_FT_aggregate_arima_vector<-PS_ARIMAX_remodel(aggregate_arima_vector_country_ranking_df[,
                                                                                                        1:ncol(aggregate_arima_vector_country_ranking_df)-1],Y,NULL,result = "reconstruct",rename_Y = "Yvec")[[6]]
    rownames(IFT_NilPhase_FT_aggregate_arima_vector)<-rownames(X)
    set.seed(Seed)
    cvLASSO_nil_kime = cv.glmnet(data.matrix(IFT_NilPhase_FT_aggregate_arima_vector[ , 1:ncol(IFT_NilPhase_FT_aggregate_arima_vector)]),
                                 Y, alpha = 1, parallel=TRUE)
    predLASSO_nil_kime <-  predict(cvLASSO_nil_kime, s = ifelse(lam=="min",cvLASSO_nil_kime$lambda.min,lam),
                                   newx = data.matrix(IFT_NilPhase_FT_aggregate_arima_vector[ , 1:ncol(IFT_NilPhase_FT_aggregate_arima_vector)]))
    validation_lim <- cbind(PREDlasso[ , 1],predLASSO_lim[, 1], predLASSO_nil_kime[, 1], Y)
    colnames(validation_kime_swapped) <- c("predLASSO_IFT_SwappedPhase","predLASSO (spacetime)", "predLASSO_IFT_NilPhase", "Orig_Y")
    NN<-5
    LMylab<-c("predLASSO_spacekime_swapped Country Overall Ranking")
    LMmain<-c("Observed (x) vs. Kime IFT_SwappedPhase_FT_Y (y) Overall Country Ranking, cor=%.02f")
    Xin_result4<-Y
    Yin_result4<-PREDlasso
    result4title<-c("Spacekime LASSO Predicted, Swapped-Phases (y) vs. Observed (x) Overall Country Ranking, cor=%.02f")
    result4Y_T<-c("Spacekime LASSO Predicted, using Swapped-Phases")
  }
  head(validation_nil_kime)
  text4<-dim(validation_lim)
  text5<-head(validation_lim)
  validation_lim <- as.data.frame(cbind(validation_lim, ifelse (validation_lim[, 1]<=30, 1, 0)))
  colnames(validation_lim)[NN] <- "Top30Rank"
  text6<-head(validation_lim)
  text7<-cor(validation_lim[ , 1], validation_lim[, NN-1])
  linFit_lim <- lm(validation_lim[ , 1] ~ validation_lim[, 2])
  if(result==3){
    plot(validation_lim[ , 1] ~ validation_lim[, NN-1],
         col="blue", xaxt='n', yaxt='n', pch = 16, cex=3,
         xlab="Observed Country Overall Ranking", ylab=LMylab,
         main = sprintf(LMmain,
                        cor(validation_lim[ , 1], validation_lim[, NN-1])))
    abline(linFit_lim, lwd=3, col="red")
  }
  if(result==4){
    myPlot <- ggplot(as.data.frame(validation_lim), aes(x=Xin_result4,
                                                        y=Yin_result4, label=rownames(validation_lim))) +
      geom_smooth(method='lm') +
      geom_point() +
      # Color by groups
      # geom_text(aes(color=factor(rownames(validation_lim)))) +
      geom_label_repel(aes(label = rownames(validation_lim),
                           fill = factor(Top30Rank)), color = 'black', size = 5,
                       point.padding = unit(0.3, "lines")) +
      # theme(legend.position = "bottom") +
      theme(legend.position = c(0.1, 0.9),
            legend.text = element_text(colour="black", size=12, face="bold"),
            legend.title = element_text(colour="black", size=14, face="bold")) +
      scale_fill_discrete(name = "Country Overall Ranking",
                          labels = c("Below 30 Rank", "Top 30 Rank")) +
      labs(title=sprintf(result4title,  cor(validation_lim[ , 1], validation_lim[, NN-1])),
           x ="Observed Overall Country Ranking (1 is 'best')",
           y = "Spacetime LASSO Predicted")
    return(myPlot)
  }
  if(result==5){
    RESlist<-list(dimX=text0,Coef_first9=text1,COR_pred_true=text2,LASSO_Varimp=text3,dimLASSO=text4,LASSO_part_feature=text5,LASSO_rename_feature=text6,pred_COR=text7)
    return(RESlist)
  }
}
