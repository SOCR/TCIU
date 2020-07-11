#' @title TSplot
#' @description \code{plot_ly} method working on time-series analysis. Work only for month & year dataset (i.e. dataset that can
#' satisfy the format required by function \code{as.yearmon} from package \code{zoo})
#'
#' @param origin_t Number of periods of original time-series data you wish to include in the plot
#' write all if all periods should be included
#' @param ARIMAmodel ARIMA model created by function \code{auto.arima()}
#' @param XREG if using ARIMAX model, put in the regularized X matrix
#' @param NEWtitle title for this plot
#' @param Ylab label of Y axis
#' @param Xlab label of X axis
#' @param ts_original label for original time series line
#' @param ts_forecast label for forecasted time series line
#' @param title_size size of the title
#'
#' @details
#' The function \code{TSplot} is based on package \code{plotly}. It applies \code{plot_ly} function to create
#' interactive plot for time-series analysis result. It requires a fitted model by function \code{auto.arima}.
#' If you are fitting an ARIMA model with external regressors (i.e. \code{Xreg}), then you must put inside the
#' external regressors again.
#'
#' @return a plot result created by plot_ly() function
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#' @export
#'
#' @import forecast zoo plotly
#' @importFrom stats time
#'
#' @examples
#' require(forecast)
#' require(zoo)
#' require(plotly)
#'
#' # Create time-series data
#' MCSI_Data_monthAvg_ts_Y <- ts(Y, start=c(1978,1), end=c(2018, 12), frequency = 12)
#'
#' # Apply ARIMAX model
#' modArima <- auto.arima(MCSI_Data_monthAvg_ts_Y, xreg=X)
#'
#' # Create plot_ly results
#' # 48 means that there will be 48 periods from the original
#' # time-series dataset that is included in the plot result.
#' # You could also change this to "all" to see all original dataset in a single plot.
#' TSplot(48,modArima,X_new,title_size = 8,ts_original = "Original time series",
#' ts_forecast = "Predicted time series")
#'
TSplot<- function(origin_t,ARIMAmodel,XREG=NULL,NEWtitle="Result",Ylab="Value",Xlab="Time(Month/Year)"
                  ,ts_original="original time series",ts_forecast="forecasted time series",title_size=10) {
  tsmodel<-forecast(ARIMAmodel, xreg = XREG)
  if(origin_t=="all"){
    TIME=1
  }
  else{
    TIME=(length(tsmodel$x)-origin_t+1)
  }
  includetime<-c(tsmodel$x[TIME:length(tsmodel$x)],rep(NA,length(tsmodel$mean)))
  includetime2<-c(rep(NA,length(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$mean)
  includetime3<-c(rep(NA,length(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$lower[,1])
  includetime4<-c(rep(NA,length(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$upper[,1])
  includetime5<-c(rep(NA,length(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$lower[,2])
  includetime6<-c(rep(NA,length(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]))),tsmodel$upper[,2])
  alltime<-c(as.yearmon(time(tsmodel$x)[TIME:length(tsmodel$x)]),as.yearmon(time(tsmodel$mean)))
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
    add_lines(x=alltime,text=as.character(alltime),y=includetime,name=ts_original,line=list(color="green"))%>%
    add_lines(x=alltime,text=as.character(alltime),y=includetime5,name="95% lower bound",line=list(color="powderblue"))%>%
    add_trace(x=alltime,text=as.character(alltime),y=includetime6,type="scatter",mode="lines",line=list(color="powderblue"),fill = 'tonexty',fillcolor="powderblue",name="95% upper bound")%>%
    add_lines(x=alltime,text=as.character(alltime),y=includetime3,name="80% lower bound",line=list(color="lightpink"))%>%
    add_trace(x=alltime,text=as.character(alltime),y=includetime4,type="scatter",mode="lines",line=list(color="lightpink"),fill = 'tonexty',fillcolor="lightpink",name="80% upper bound")%>%
    add_lines(x=alltime,text=as.character(alltime),y=includetime2,name=ts_forecast,line=list(color="red"))
  return(TSP)
}
