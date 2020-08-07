# See more details in TSplotly package

# Helper function for fmri_time_series
GTSplot <- function(tsdata, NEWtitle = "Result", Ylab = "Value", Xlab = "Time", Unit = NULL, ts_name = NULL, title_size = 10, COLO = NULL) {
    TSP <- plot_ly(type = "scatter", mode = "lines")
    for (i in 1:ncol(tsdata)) {
        tsd <- tsdata[, i]
        tsn <- ts_name[i]
        Col <- COLO[i]
        TSP <- add_trace(TSP, x = time(tsd), text = paste(time(tsd), Unit), type = "scatter", mode = "lines", opacity = 0.75, y = tsd, 
            name = tsn, line = list(color = c(Col)))
    }
    TSP <- TSP %>% layout(title = list(text = NEWtitle, font = list(family = "Times New Roman", size = title_size, color = "black")), 
        paper_bgcolor = "rgb(255,255,255)", plot_bgcolor = "rgb(229,229,229)", xaxis = list(title = Xlab, gridcolor = "rgb(255,255,255)", 
            showgrid = TRUE, showline = FALSE, showticklabels = TRUE, tickcolor = "rgb(127,127,127)", ticks = "outside", zeroline = FALSE), 
        yaxis = list(title = Ylab, gridcolor = "rgb(255,255,255)", showgrid = TRUE, showline = FALSE, showticklabels = TRUE, tickcolor = "rgb(127,127,127)", 
            ticks = "outside", zeroline = FALSE))
    return(TSP)
}

# Helper function for fmri_ts_forecast
TSplot_gen <- function(origin_t, ARIMAmodel, XREG = NULL, periods = NULL, NEWtitle = "Result", Ylab = "Value", Xlab = "Time", plot_labels = NULL, 
    ts_original = "original time series", ts_forecast = "forecasted time series", title_size = 10, ts_list = "empty", ts_labels = NULL, 
    ts_names = NULL, COLO = NULL) {
    tsmodel <- forecast(ARIMAmodel, xreg = XREG, h = periods)
    if (origin_t == "all") {
        TIME = 1
    } else {
        TIME = (length(tsmodel$x) - origin_t + 1)
    }
    includetime <- c(tsmodel$x[TIME:length(tsmodel$x)], rep(NA, length(tsmodel$mean)))
    includetime2 <- c(rep(NA, length((time(tsmodel$x)[TIME:length(tsmodel$x)]))), tsmodel$mean)
    includetime3 <- c(rep(NA, length((time(tsmodel$x)[TIME:length(tsmodel$x)]))), tsmodel$lower[, 1])
    includetime4 <- c(rep(NA, length((time(tsmodel$x)[TIME:length(tsmodel$x)]))), tsmodel$upper[, 1])
    includetime5 <- c(rep(NA, length((time(tsmodel$x)[TIME:length(tsmodel$x)]))), tsmodel$lower[, 2])
    includetime6 <- c(rep(NA, length((time(tsmodel$x)[TIME:length(tsmodel$x)]))), tsmodel$upper[, 2])
    alltime <- c((time(tsmodel$x)[TIME:length(tsmodel$x)]), (time(tsmodel$mean)))
    TSP <- plot_ly(type = "scatter", mode = "lines") %>% layout(title = list(text = NEWtitle, font = list(family = "Times New Roman", 
        size = title_size, color = "black")), paper_bgcolor = "rgb(255,255,255)", plot_bgcolor = "rgb(229,229,229)", xaxis = list(title = Xlab, 
        gridcolor = "rgb(255,255,255)", showgrid = TRUE, showline = FALSE, showticklabels = TRUE, tickcolor = "rgb(127,127,127)", 
        ticks = "outside", zeroline = FALSE), yaxis = list(title = Ylab, gridcolor = "rgb(255,255,255)", showgrid = TRUE, showline = FALSE, 
        showticklabels = TRUE, tickcolor = "rgb(127,127,127)", ticks = "outside", zeroline = FALSE)) %>% add_lines(x = alltime, 
        text = plot_labels, y = includetime, name = ts_original, line = list(color = "green")) %>% add_lines(x = alltime, text = plot_labels, 
        y = includetime5, name = "95% lower bound", line = list(color = "powderblue")) %>% add_trace(x = alltime, text = plot_labels, 
        y = includetime6, type = "scatter", mode = "lines", line = list(color = "powderblue"), fill = "tonexty", fillcolor = "powderblue", 
        name = "95% upper bound") %>% add_lines(x = alltime, text = plot_labels, y = includetime3, name = "80% lower bound", line = list(color = "lightpink")) %>% 
        add_trace(x = alltime, text = plot_labels, y = includetime4, type = "scatter", mode = "lines", line = list(color = "lightpink"), 
            fill = "tonexty", fillcolor = "lightpink", name = "80% upper bound") %>% add_lines(x = alltime, text = plot_labels, 
        y = includetime2, name = ts_forecast, line = list(color = "red"))
    # Add new lines to TSP:
    if (ts_list != "empty") {
        for (i in 1:length(ts_list)) {
            tsd <- ts_list[[i]]
            tsl <- ts_labels[[i]]
            tsn <- ts_names[i]
            Color <- COLO[i]
            TSP <- add_trace(TSP, x = time(tsd), text = tsl, type = "scatter", mode = "lines", y = tsd, name = tsn, line = list(Color))
        }
    }
    return(TSP)
}


# Helper function for fmri_kimesurface
# create continuous hexadecimal color vector for data
fmri_split_ab_bl <- function(vect, option = "vector") {
    if (option == "list") {
        overalllen <- length(which(vect != 0)) + 1
        ab_len <- length(which(vect > 0)) + 1
        bl_len <- length(which(vect < 0)) + 1
        nulllist <- as.list(1:overalllen)
        s <- seq_gradient_pal("#FFFF00", "#FFFFFF")(seq(0, 1, length.out = bl_len))
        for (i in 1:bl_len) {
            nulllist[[i]] <- c((i - 1)/overalllen, s[i])
        }
        s <- seq_gradient_pal("#FFFFFF", "#0000FF")(seq(0, 1, length.out = ab_len))
        for (i in 1:ab_len) {
            nulllist[[bl_len - 1 + i]] <- c((bl_len - 1 + i - 1)/overalllen, s[i])
        }
    } else if (option == "vector") {
        overalllen <- length(which(vect != 0)) + 1
        ab_len <- length(which(vect > 0)) + 1
        bl_len <- length(which(vect < 0)) + 1
        nulllist <- rep(NA, overalllen)
        s <- seq_gradient_pal("#FFFF00", "#FFFFFF")(seq(0, 1, length.out = bl_len))
        for (i in 1:bl_len) {
            nulllist[i] <- s[i]
        }
        s <- seq_gradient_pal("#FFFFFF", "#0000FF")(seq(0, 1, length.out = ab_len))
        for (i in 1:ab_len) {
            nulllist[[bl_len - 1 + i]] <- s[i]
        }
    }
    return(nulllist)
}

