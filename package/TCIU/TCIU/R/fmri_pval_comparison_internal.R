# @title fmri_pval_comparison_internal 
# 
# @description This function servers as an internal function to compare two p-values.
# @param p_val1 the first p value 
# @param p_val2 the second p value 
# 
# @details The function \code{pval_comparison.internal} is used to compare two p-values.  
# @author SOCR team <\url{http://socr.umich.edu/people/}> 
# @examples

fmri_pval_comparison_internal = function(p_val1,
                                         p_val2,
                                         names = c("p_val1", "p_val2"),
                                         breaks = 10) {
    p_val_df = data.frame(p_val1 = as.vector(p_val1), p_val2 = as.vector(p_val2))
    names(p_val_df) = names
    
    testData = p_val_df %>% tidyr::gather(key = "p_value_name", value = "p_value")
    
    
    toPlot <- testData %>% mutate(bin = cut(p_value, pretty(p_value, breaks))) %$% table(bin, p_value_name) %>% as.data.frame() %>% 
        mutate(plotVal = ifelse(p_value_name == names[1], -1 * Freq, Freq))
    
    
    p_val_comparison = toPlot %>% ggplot(aes(x = bin, y = plotVal, fill = p_value_name)) + geom_col() + ggtitle(paste0("Comparison between ", 
        names[1], " and ", names[2])) + theme(plot.title = element_text(lineheight = 0.8, face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 90, 
        hjust = 0.5))
    
    return(p_val_comparison)
    
}
