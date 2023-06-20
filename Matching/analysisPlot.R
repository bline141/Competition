# Generates a density plot for numeric variables and a barplot for factors/logicals
# Distribution is shown for each "therapy"
analysisPlot <- function(df, column, therapy){
  # df: data frame with containing atleast two variables (column, therapy)
  # column: variable for which plot should be made
  # therapy: treatment variable
  
  if (is(pull(df, column))[1] %in% c("numeric", "integer")){
    p <- ggplot(df, aes(x = .data[[column]])) + 
      geom_density(aes(fill = .data[[therapy]], colour = .data[[therapy]]), alpha = 0.6) +
      xlab(column) + 
      scale_fill_brewer("Cohort", palette = "Set1") +
      scale_colour_brewer("Cohort", palette = "Set1") +
      guides(color = guide_legend(override.aes = list(linetype = 0))) +
      theme_bw()
  } else if (is(pull(df, column))[1] %in% c("factor", "logical", "character")){
    p <- ggplot(df, aes(x = .data[[column]])) + 
      geom_bar(aes(fill = .data[[therapy]], colour = .data[[therapy]]), alpha = 0.6, position = "dodge") +
      xlab(column) + 
      scale_fill_brewer("Cohort", palette = "Set1") +
      guides(colour = "none") +
      theme_bw()
  }
  
  return(p)
}
