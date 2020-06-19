#' @title PS_ARIMAX_remodel
#' @description ARIMAX Remodeling
#' @details Conduct three different phase synthesis methods and create ARIMAX model based on different method.
#' This function can reconstruct the dataset based on three phase synthesis methods:
#' Nil-Phase Synthesis, Swapped-Phase Synthesis and Random-Phase Synthesis.
#' Also it can construct different ARIMAX model based on the phase synthesis method chosen.
#'
#' @param X inherited from function "homologousX_features", we can get parameter X by
#' "homologousX_features(country1,country2)$$X_Country1"
#' @param Y inherited from function "preprocess_ARIMA", we can get parameter Y by "preprocess_ARIMA$Y"
#' @param ts_Y_test time series vector, must be set from function "ts" and
#' must be created based on testing set
#' @param option choose the phase synthesis method used in the reconstruction or ARIMAX model.
#' Three options are valid: "Nil-Phase", "Swapped-Phase", "Random-Phase"
#' @param result choose the result shown by this function. Two options are valid:
#' "ARIMAX": this will show the ARIMAX model created by the chosen phase synthesis method.
#' "reconstruct": this will show the reconstruction procedure based on the chosen phase synthesis method,
#' but the ARIMAX model won't be built.
#' @param rename_Y rename the return result of Y
#' @return If the "result" parameter is setted as "reconstruct", then the reconstruction result will be shown.
#' If setted as "ARIMAX", then the detailed ARIMAX model will be built.
#'
#' @author SOCR team <http://socr.umich.edu/people/>
#' @export
#' 
#' @examples 
#' require(DT)
#' # Show the reconstruction results based on "Nil-Phase" transformation
#' Nil_reconstruct<-PS_ARIMAX_remodel(X_Belgium,Y_Belgium,ts_Y_Belgium_test,result = "reconstruct")
#' ## All different features under the "Nil-Phase"
#' datatable(data.frame(feature_name=Nil_reconstruct[[4]]),fillContainer = TRUE)
#' ## All reconstruct results
#' datatable(data.frame(Nil_reconstruct[[6]]),fillContainer = TRUE)
#' 
#' # Example for showing the ARIMA model based on "Swapped-Phase" transformation
#' Swapped_ARIMAX<-PS_ARIMAX_remodel(X_Belgium,Y_Belgium,ts_Y_Belgium_test,option = "Swapped-Phase",result = "ARIMAX")
#' ## ARIMA information (AR & MA)
#' Swapped_ARIMAX$AR_MA
#' ## Model & time series information
#' Swapped_ARIMAX$ARIMAX_model
#' ## Information for most important features
#' Swapped_ARIMAX$feature_effect
PS_ARIMAX_remodel <- function(X,Y,ts_Y_test,option="Nil-Phase",result="ARIMAX",rename_Y="Y_GDP_Belgium") {
  temp_data<-cbind(X,Y)
  FT<-kSpaceTransform(temp_data, FALSE, NULL)
  if(option=="Nil-Phase"){
    nilPhase_FT_data <- array(complex(real=0, imaginary=0), c(dim(temp_data)[1], dim(temp_data)[2]))
    ##IFT_NilPhase_FT <- array(complex(), c(dim(temp_data)[1], dim(temp_data)[2]))
    #       Invert back to spacetime the FT_Belgium$magnitudes[ , i] signal with nil-phase
    IFT_NilPhase_FT<- Re(kSpaceTransform(FT$magnitudes, TRUE, nilPhase_FT_data))
    colnames(IFT_NilPhase_FT) <- c(colnames(X), rename_Y)

    #result
    if(result=="reconstruct"){
      return(list(
        dim(nilPhase_FT_data),dim(IFT_NilPhase_FT), dim(FT$magnitudes),
        colnames(IFT_NilPhase_FT), head(IFT_NilPhase_FT),IFT_NilPhase_FT)) # head(temp_Data)
    }
    #rename
    IFT_FT<-IFT_NilPhase_FT
  }
  else if(option=="Swapped-Phase"){

    swapped_phase_FT_data <- FT$phases
    colnames(swapped_phase_FT_data) <- c(colnames(X), rename_Y)
    swapped_phase_FT_data1 <- swapped_phase_FT_data
    IFT_SwappedPhase_FT <- array(complex(), c(dim(temp_data)[1], dim(temp_data)[2]))
    set.seed(12345)   # sample randomly Phase-columns for each of the 131 covariates (X)
    #swap_phase_FT_Belgium_indices <- sample(ncol(swapped_phase_FT_Belgium_data)-1)
    # for (j in 1:131) {  # for all coluns of the design Xreg matrix, excluding Y, randomly swap columns phases
    #  swapped_phase_FT_Belgium_data1[ , j] <- swapped_phase_FT_Belgium_data[, swap_phase_FT_Belgium_indices[j]]
    #}
    swapped_phase_FT_data1 <- as.data.frame(cbind(
      swapped_phase_FT_data[ , sample(ncol(swapped_phase_FT_data[ , 1:(ncol(swapped_phase_FT_data)-1)]))],
      swapped_phase_FT_data[ , ncol(swapped_phase_FT_data)]))
    swapped_phase_FT_data <- swapped_phase_FT_data1
    colnames(swapped_phase_FT_data)[ncol(swapped_phase_FT_data)] <- rename_Y
    colnames(swapped_phase_FT_data)

    # Invert back to spacetime the FT_Belgium$magnitudes[ , i] signal using the feature swapped phases
    IFT_SwappedPhase_FT <- Re(kSpaceTransform(FT$magnitudes, TRUE, swapped_phase_FT_data))

    colnames(IFT_SwappedPhase_FT) <- c(colnames(X), rename_Y)

    #result
    if(result=="reconstruct"){
      return(list(
        dim(swapped_phase_FT_data) ,
        # [1] 360 132
        dim(swapped_phase_FT_data), dim(FT$phases),
        dim(IFT_SwappedPhase_FT), dim(FT$magnitudes),
        colnames(IFT_SwappedPhase_FT), tail(IFT_SwappedPhase_FT),IFT_SwappedPhase_FT)) # tail(temp_Data)
    }
    #rename
    IFT_FT<-IFT_SwappedPhase_FT
  }
  else if(option=="Random-Phase"){

    #randPhase_FT_data <- array(complex(), c(dim(temp_data)[1], dim(temp_data)[2]))
    IFT_RandPhase_FT <- array(complex(), c(dim(temp_data)[1], dim(temp_data)[2]))
    randPhase_FT_data <- FT$phases
    for (i in 1:(dim(randPhase_FT_data)[2] -1)) {
      if (i < dim(randPhase_FT_data)[2]) {
        set.seed(12345)   # sample randomly Phases for each of the 131 predictors covariates (X)
        randPhase_FT_data[ , i] <- FT$phases[sample(nrow(FT$phases)), i]
      } else {   } # for the Y outcome (Last Column) - do not change the phases of the Y
    }
    #       Invert back to spacetime the FT_Belgium$magnitudes[ , i] signal with avg-phase
    IFT_RandPhase_FT <- Re(kSpaceTransform(FT$magnitudes, TRUE, randPhase_FT_data))
    colnames(IFT_RandPhase_FT) <- c(colnames(X), rename_Y)

    #result
    if(result=="reconstruct"){
      return(  list(
        dim(randPhase_FT_data),    # ;  head(randPhase_FT_data)
        # [1] 360 132
        dim(IFT_RandPhase_FT), dim(FT$magnitudes) ,
        colnames(IFT_RandPhase_FT), tail(IFT_RandPhase_FT), # tail(temp_Data)
        dim(IFT_RandPhase_FT), head(Re(IFT_RandPhase_FT)), tail(Re(IFT_RandPhase_FT))))
    }
    #rename
    IFT_FT<-IFT_RandPhase_FT
  }
  if(result=="ARIMAX"){
    #construction of time-series analysis
    # Perform ARIMAX modeling on IFT_NilPhase_FT_Belgium; report (p,d,q) params and quality metrics AIC/BIC
    # library(forecast)
    IFT_FT_Y_train <- IFT_FT[1:300, ncol(IFT_FT)]
    IFT_FT_Y_test <- IFT_FT[301:360]

    # Training and Testing Data Covariates explaining the longitudinal outcome (Y)
    IFT_FT_X_train <- as.data.frame(IFT_FT)[1:300, 1:(ncol(IFT_FT)-1)]; dim(IFT_FT_X_train)
    IFT_FT_X_test <- as.data.frame(IFT_FT)[301:360, 1:(ncol(IFT_FT)-1)]; dim(IFT_FT_X_test)

    # Outcome Variable to be ARIMAX-modelled, as a timeseries
    ts_IFT_FT_Y_train <-
      ts(IFT_FT_Y_train, start=c(2000,1), end=c(2014, 20), frequency = 20)

    set.seed(1234)
    modArima_IFT_FT_Y_train <-
      auto.arima(ts_IFT_FT_Y_train, xreg=as.matrix(IFT_FT_X_train))

    pred_arimax <- forecast(modArima_IFT_FT_Y_train, xreg = as.matrix(IFT_FT_X_test))
    pred_arimax_2015_2017 <-
      ts(pred_arimax$mean, frequency=20, start=c(2015,1), end=c(2017,20))
    # alternatively:
    # pred_arimax_1_0_1_2015_2017 <- predict(modArima_IFT_NilPhase_FT_Belgium_Y_train,
    #                                              n.ahead = 3*20, newxreg = IFT_NilPhase_FT_Belgium_X_test)$pred
    sort(modArima_IFT_FT_Y_train$coef)[1:10]
    # Labor cost for LCI (compensation of employees plus taxes minus subsidies), effect=-1.5972295
    #                                                                      ar1, effect=-1.4231617
    #                                 Labor cost other than wages and salaries, effect=-1.2213214
    #                                                                      ma1, effect=-0.9869571
    #                                                                      ar2, effect=-0.8591937
    #           Unemployment , Females, From 15-64 years, From 12 to 17 months, effect=-0.7075454
    #             Unemployment , Total, From 15-64 years, From 18 to 23 months, effect=-0.5797656
    #               Unemployment , Males, From 15-64 years, from 3 to 5 months, effect=-0.5026139
    #                                                                     sar2, effect=-0.3450866
    #             Unemployment , Males, From 15-64 years, from 24 to 47 months, effect=-0.2965540
    return(list(
      training_set_length=length(IFT_FT_Y_train),
      testing_set_length=length(IFT_FT_Y_test),
      AR_MA=modArima_IFT_FT_Y_train$arma,
      ARIMAX_model=pred_arimax_2015_2017,
      feature_effect=sort(modArima_IFT_FT_Y_train$coef)[1:10],
      Correlation=cor(pred_arimax$mean, ts_Y_test),
      Mean=mean(pred_arimax_2015_2017)))}
  else{
    print(paste0("Wrong input for option or result, please try again.") )
    return(NULL)
  }
}
