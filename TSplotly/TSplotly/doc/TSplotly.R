## ----Installation of the CBDA package, eval = FALSE----------------------
#  # Installation from the Windows binary (recommended for Windows systems)
#  ## Set the directory to the working space contains "TSplotly_1.1.0.tar.gz"
#  system("R CMD INSTALL TSplotly")
#  # Installation from the source (recommended for Macs and Linux systems)
#  install.packages("~/TSplotly_1.1.0.tar.gz", repos = NULL, type = "source")

## ----Installation of the CBDA package from CRAN, eval = FALSE------------
#  install.packages("TSplotly")

## ----warning=FALSE,message=FALSE-----------------------------------------
require(TSplotly)
require(zoo)
require(ggplot2)
require(plotly)
require(forecast)

# Creating time series data
MCSI_Data_monthAvg_ts_Y <- ts(Y, start=c(1978,1), end=c(2018, 12), frequency = 12)

# Applying ARIMAX model
modArima <- auto.arima(MCSI_Data_monthAvg_ts_Y, xreg=X)

# Creating plot_ly results
## 48 means that there will be 48 periods from the original
## time series dataset that is included in the plot result.
## You could also change this to "all" to see all original dataset in a single plot.
TSplot(48,modArima,X_new,title_size = 8,ts_original = "Original time series",
ts_forecast = "Predicted time series")

## ----warning=FALSE,message=FALSE-----------------------------------------
#Create labels for training time series data and ARIMAX result (48 periods of training data included)
require(zoo)
#Time labels for training data
time_label1<-as.yearmon(time(MCSI_Data_monthAvg_ts_Y))[(length(MCSI_Data_monthAvg_ts_Y)-48+1):length(MCSI_Data_monthAvg_ts_Y)]
#Time labels for ARIMAX model(need to fit model first)
time_pred<-forecast(modArima,xreg = X_new)
time_label2<-as.yearmon(time(time_pred$mean))

time_label<-as.character(c(time_label1,time_label2))

TSplot_gen(48,modArima,X_new,title_size = 8,ts_original = "Original time series",
ts_forecast = "Predicted time series", #inculde labels inside
plot_labels = time_label)

## ----warning=FALSE,message=FALSE-----------------------------------------
# Step 1: creating the base plot
## Creating time labels
tl1<-as.yearmon(time(modArima_train$x))[(length(modArima_train$x)-48+1):length(modArima_train$x)]
tl2<-as.yearmon(time(forecast(modArima_train,xreg = as.matrix(X_test))$mean))
tl<-as.character(c(tl1,tl2))
  
Tempplot<-TSplot_gen(48,modArima_train,as.matrix(X_test),title_size = 8,ts_original = "Original time series",
ts_forecast = "Predicted time series")

# Show base plot if no other elements(labels, new time lines, etc)is included
Tempplot

## ----warning=FALSE,message=FALSE-----------------------------------------
# Step 2: including new lines and labels
## Creating list and other information for new lines
TSlist<-list(MCSI_Data_monthAvg_ts_Y_test)
TSlabel<-list(as.character(as.yearmon(time(TSlist[[1]]))))
TSname<-c("Original result")

## Put them into related parameters
TSplot_gen(48,modArima_train,as.matrix(X_test),title_size = 8,ts_original = "Original time series",
           ts_forecast = "Predicted time series",plot_labels = tl, #labels of original plot
            ts_list = TSlist,ts_names = TSname,ts_labels = TSlabel,COLO = "black")

## ----warning=FALSE,message=FALSE-----------------------------------------
require(forecast)

#Firstly create a base plotly plot
Tempplot<-TSplot(48,modArima_train,as.matrix(X_test),title_size = 8,ts_original = "Original time series",
ts_forecast = "Predicted time series")

# Generate a new line with ADDline function
newline<-ADDline(TS = MCSI_Data_monthAvg_ts_Y_test,linetype = "TS",Name = "Original Result")

## Put the new line into our plot
Tempplot%>%
  add_lines(x=newline$X,text=newline$TEXT,y=newline$Y,name=newline$NAME,line=list(color="grey"))

## ----warning=FALSE,fig.width=14,fig.height=10,out.width=1920,out.height=1080,message=FALSE----
ggplot(MCSI_Data_monthAvg_melt[MCSI_Data_monthAvg_melt$series!="INCOME", ],
       aes(YYYYMM, value)) +
  geom_line(aes(linetype=series, colour = series), size=2) +
  geom_point(aes(shape=series, colour = series), size=0.3) +
  geom_smooth(aes(colour = series), se = TRUE) +
  coord_trans(y="log10") +
  xlab("Time (monthly)") + ylab("Index Values (log-scale)") +
  scale_x_date(date_breaks = "12 month", date_labels =  "%m-%Y")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=20))+ theme(legend.position="top")

## ----warning=FALSE,message=FALSE-----------------------------------------
PYdf<-GtoP_trans(MCSI_Data_monthAvg_melt[MCSI_Data_monthAvg_melt$series!="INCOME", ],NAME="series",X="YYYYMM",Y="value")
PYdf$INCOME<-NULL
#Log10 transformation
PYdf<-log10(PYdf)

## ----warning=FALSE,message=FALSE-----------------------------------------
#Create an interactive list
updatemenus <- list(
  list(
    active = -1,
    type= 'buttons',
    buttons = list(
      list(
        label = "ALL",
        method = "update",
        args = list(list(visible = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)),
                    list(title = "All Indexes"))),
      list(
        label = "ICS",
        method = "update",
        args = list(list(visible = c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)),
                    list(title = "ICS"))),
      list(
        label = "ICC",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)),
                    list(title = "ICC"))),
      list(
        label = "GOVT",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE)),
                    list(title = "GOVT"))),
      list(
        label = "DUR",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)),
                    list(title = "DUR"))),
      list(
        label = "HOM",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE)),
                    list(title = "HOM"))),
      list(
        label = "CAR",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)),
                    list(title = "CAR"))),
      list(
        label = "AGE",
        method = "update",
        args = list(list(visible = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)),
                    list(title = "AGE"))),
      list(
        label = "EDUC",
        method = "update",
        args = list(list(visible = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)),
                    list(title = "EDUC")))
      )
  )
)
# Apply plot_ly to finish generating result
plot_ly(type="scatter",mode="lines")%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$ICS,name="ICS",line=list(color="powderblue"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$ICC,name="ICC",line=list(color="red"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$GOVT,name="GOVT",line=list(color="green"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$DUR,name="DUR",line=list(color="orange"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$HOM,name="HOM",line=list(color="purple"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$CAR,name="CAR",line=list(color="pink"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$AGE,name="AGE",line=list(color="brown"))%>%
  add_lines(x=as.yearmon(rownames(PYdf)),text=rownames(PYdf),y=PYdf$EDUC,name="EDUC",line=list(color="powderblue"))%>%
  layout(title= list(text="Time series for 8 indexes",font=list(family = "Times New Roman",size = 16,color = "black" )),
           paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
           xaxis = list(title ="Time (monthly)",
                        gridcolor = 'rgb(255,255,255)',
                        showgrid = TRUE,
                        showline = FALSE,
                        showticklabels = TRUE,
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE),
           yaxis = list(title = "Index Values (log-scale)",
                        gridcolor = 'rgb(255,255,255)',
                        showgrid = TRUE,
                        showline = FALSE,
                        showticklabels = TRUE,
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE),
         updatemenus=updatemenus)

