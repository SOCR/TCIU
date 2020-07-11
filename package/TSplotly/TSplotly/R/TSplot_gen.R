#' @title TSplot_gen
#' @description \code{plot_ly} method working on time series analysis. This is a more general function that can work on any
#' fitted ARIMA(X) model with any time format. Also, a list of time series data can be applied into this function
#' so that new time lines can be directly drawn in the plot without referring to the result of \code{ADDline} function.
#'
#' @param origin_t Number of periods of original time series data you wish to include in the plot write all if all
#' periods should be included
#' @param ARIMAmodel ARIMA model created by function \code{auto.arima()}
#' @param XREG if using ARIMAX model, put in the regularized X matrix
#' @param NEWtitle title for this plot
#' @param Ylab label of Y axis
#' @param Xlab label of X axis
#' @param plot_labels To include a specific labels for each time point of all training time series data and
#' predicted ARIMA(X) result. A vector should be applied to this parameter that has the same length of chosen
#' periods of training time series data (i.e. parameter \code{origin_t}) along with predicted time series periods.
#' @param ts_original label for original time series line
#' @param ts_forecast label for forecasted time series line
#' @param title_size size of the title
#' @param ts_list applying this function can help you draw more time lines into the original plot. A list should
#' be applied to this parameter which contains all extra time series data that you wish to draw on the original
#' plot. Each element on this list should be created by function \code{ts()}
#' @param ts_labels when drawing extra time lines with parameter \code{ts_list}. You could create specific labels
#' for each time points. A list with the same shape of list in \code{ts_list} should be applied. Each element in
#' this list should contain time labels corresponding with the list in \code{ts_list}
#' @param ts_names Creating labels for each extra time lines you draw. Labels will appear on the legend of the plot.
#' @param COLO Specifying colors for each new line that is drawn.
#'
#' @details
#' The function \code{TSplot_gen} is based on package \code{plotly}. It applies \code{plot_ly} function to create
#' interactive plot for time-series analysis result. It requires a fitted model by function \code{auto.arima}.
#' If you are fitting an ARIMA model with external regressors (i.e. \code{Xreg}), then you must put inside the
#' external regressors again.
#'
#' @return a plot result created by \code{plot_ly} function
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#' @export
#'
#' @import forecast zoo plotly
#' @importFrom stats time
#'
#' @examples
#' # Create time labels
#' require(zoo)
#' require(forecast)
#' require(plotly)
#'
#' tl1<-as.yearmon(time(modArima_train$x))[(length(modArima_train$x)-48+1):length(modArima_train$x)]
#' tl2<-as.yearmon(time(forecast(modArima_train,xreg = as.matrix(X_test))$mean))
#' tl<-as.character(c(tl1,tl2))
#' 
#' # Create list and other information for new lines
#' TSlist<-list(MCSI_Data_monthAvg_ts_Y_test)
#' TSlabel<-list(as.character(as.yearmon(time(TSlist[[1]]))))
#' TSname<-c("Original result")
#' 
#' # Put them into related parameters
#' TSplot_gen(48,modArima_train,as.matrix(X_test),title_size = 8,ts_original = "Original time series",
#'         ts_forecast = "Predicted time series",plot_labels = tl, #labels of original plot
#'         ts_list = TSlist,ts_names = TSname,ts_labels = TSlabel,COLO = "black")
#'
TSplot_gen<- function(origin_t,ARIMAmodel,XREG=NULL,NEWtitle="Result",Ylab="Value",Xlab="Time",plot_labels=NULL
                  ,ts_original="original time series",ts_forecast="forecasted time series",title_size=10,
                  ts_list="empty",ts_labels=NULL,ts_names=NULL,COLO=NULL) {
  tsmodel<-forecast(ARIMAmodel, xreg = XREG)
  if(origin_t=="all"){
    TIME=1
  }
  else{
    TIME=(length(tsmodel$x)-origin_t+1)
  }
  includetime<-c(tsmodel$x[TIME:length(tsmodel$x)],rep(NA,length(tsmodel$mean)))
  includetime2<-c(rep(NA,length((time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$mean)
  includetime3<-c(rep(NA,length((time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$lower[,1])
  includetime4<-c(rep(NA,length((time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$upper[,1])
  includetime5<-c(rep(NA,length((time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$lower[,2])
  includetime6<-c(rep(NA,length((time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$upper[,2])
  alltime<-c((time(tsmodel$x)[TIME:length(tsmodel$x)]),(time(tsmodel$mean)))
  TSP<-plot_ly(type="scatter",mode="lines")%>%
    layout(title= list(text=NEWtitle,font=list(family = "Times New Roman",size = title_size,color = "black" )),
           paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
           xaxis = list(title = Xlab,
                        gridcolor = 'rgb(255,255,255)',
                        showgrid = TRUE,
                        showline = FALSE,
                        showticklabels = TRUE,
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE),
           yaxis = list(title = Ylab,
                        gridcolor = 'rgb(255,255,255)',
                        showgrid = TRUE,
                        showline = FALSE,
                        showticklabels = TRUE,
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE))%>%
    add_lines(x=alltime,text=plot_labels,y=includetime,name=ts_original,line=list(color="green"))%>%
    add_lines(x=alltime,text=plot_labels,y=includetime5,name="95% lower bound",line=list(color="powderblue"))%>%
    add_trace(x=alltime,text=plot_labels,y=includetime6,type="scatter",mode="lines",line=list(color="powderblue"),fill = 'tonexty',fillcolor="powderblue",name="95% upper bound")%>%
    add_lines(x=alltime,text=plot_labels,y=includetime3,name="80% lower bound",line=list(color="lightpink"))%>%
    add_trace(x=alltime,text=plot_labels,y=includetime4,type="scatter",mode="lines",line=list(color="lightpink"),fill = 'tonexty',fillcolor="lightpink",name="80% upper bound")%>%
    add_lines(x=alltime,text=plot_labels,y=includetime2,name=ts_forecast,line=list(color="red"))
  #Add new lines to TSP:
  if(ts_list!="empty"){
    for (i in 1:length(ts_list)) {
      tsd<-ts_list[[i]]
      tsl<-ts_labels[[i]]
      tsn<-ts_names[i]
      Color<-COLO[i]
      TSP<-add_trace(TSP,x=time(tsd),text=tsl,type="scatter",mode="lines",
                     y=tsd,name=tsn,line=list(Color))
    }
  }
  return(TSP)
}
