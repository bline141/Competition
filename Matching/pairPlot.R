# For a given group or pair column: 
#   Calculate the difference for numeric variables and the common pairs (table) 
#     for factors/logicals
pairDiffCalc <- function(df, paircol = "subclass"){
  # df: data frame containing at least a grouping column and one more column
  # paircol: grouping variable as string
  
  df %>% 
    group_by_at(paircol) %>%
    summarise(across(where(~ (is.factor(.x) | is.logical(.x) | is.character(.x))), 
                     function(x){paste0(sort(x), collapse = "_")}),
              across(where(~ (is.numeric(.x) | is.integer(.x))), diff))
}

# Generates a histogram for numeric variables and a barplot for factors/logicals
pairPlot <- function(df, column){
  # df: data frame containing at least one variables (column)
  # column: variable for which plot should be made
  if (is(pull(df, column))[1] %in% c("numeric", "integer")){
    p <- ggplot(df, aes(x = .data[[column]])) + 
      geom_histogram(bins = round(sqrt(nrow(df)))) +
      geom_vline(xintercept = 0, col = "red") +
      xlab(column) + 
      theme_bw()
  } else if (is(pull(df, column))[1] %in% c("factor", "logical", "character")){
    df[["same_group"]] <- sapply(strsplit(df[[column]], "_"), function(x){length(unique(x)) == 1})
    p <- ggplot(df, aes(x = .data[[column]])) + 
      geom_bar(aes(fill = same_group)) +
      scale_fill_manual(values = c("red", "green4")) +
      xlab(column) + 
      theme_bw()
  }
}
