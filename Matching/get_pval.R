

#' P-value computation
#' A.D. & J.B.
#' 
#' Compute a statistical test for two groups on one variable and return its p-value. 
#' The test specified depends on the variable type (categorical/continuous) and on
#' whether the values are considered paired between the groups or not.
#' **** IMPORTANT: for paired tests, variable must be in order of the pair IDs within each
#' group. *****
#'         Paired    No                  Yes
#' Values           
#'    Continuous      Wilcoxon rank-sum   Wilcoxon signed-rank
#'    Categorical     Chi-squared         McNemar
#'    
#' @param group Vector containing grouping variable. Either numeric or factor so that the order of groups is respected.
#' @param variable Vector containing the values of the variable to compute the test on. Same length as group. 
#' Either numeric, categorical, or factor.
#' @param pair_id Character vector specifying IDs for pairs of observations. Default is NULL, in which case
#' the data are not assumed to be paired and tests for unpaired data are used; if not NULL, this should be
#' a character vector of the same length as group and variable and tests for paired data will be applied.
#' 
#' @return Numeric value if there are two distinct groups, otherwise NA_real_. 
#' 
#' @example
#' df_example <- data.frame(therapy = factor(rep(c("A", "B"), each = 20)),
#'                          age     = c(rnorm(n = 20, mean = 40, sd = 4),
#'                                      rnorm(n = 20, mean = 30, sd = 4)),
#'                          eye_colour = c(rep(c("blue", "brown", "green"), times = c(6, 6, 8)),
#'                                         rep(c("green", "blue", "brown"), times = c(12, 4, 4))),
#'                          stringsAsFactors = FALSE)
#' 
#' get_pval(group = df_example$therapy, variable = df_example$age, test = "wilcox.ranksum")
#' get_pval(group = df_example$therapy, variable = df_example$eye_colour, test = "chisq")
#'
#' @export
get_pval <- function(group, variable, pair_id = NULL, ...) {
  
  ## Prepare testing data -----------------------------------------
  ## Get group IDs
  group_id <- as.numeric(factor(group))
  
  ## Gather all testing data
  df_test <- data.frame(group_id,
                        variable,
                        stringsAsFactors = FALSE)
  
  ## If data are paired, add ID and sort by pairs!!! Very important
  if (!is.null(pair_id)) {
    
    df_test <- df_test %>%
      mutate(pair_id = pair_id) %>%
      arrange(pair_id, group_id)
  }
  
  if (max(group_id) != 2) {
    return(NA_real_)
  }
  
  ## Test for categorical variable ----------------------------------
  
  if (class(variable) %in% c("character", "factor")) {
    
    df_test <- df_test %>%
      mutate(variable_char   = as.character(variable),
             variable_factor = factor(variable_char,
                                      levels = unique(variable_char)))
    levs <- levels(df_test$variable_factor)
    
    if (is.null(pair_id)) {
      
      ## Chi-squared test
      
      tbl <- with(df_test, table(group_id, variable_factor))
      
      p <-  chisq.test(x = tbl,
                       correct = TRUE,
                       ...)$p.value
    } else {
      
      ## McNemar test
      
      tbl <- with(df_test, 
                  table(variable_factor[group_id == 1],
                        variable_factor[group_id == 2]))
      
      p <- mcnemar.test(tbl)$p.value
      
    }
    
  
  ## Test for continuous variable --------------------------------------
    
  } else if (class(variable) %in% c("numeric", "integer")) {
    
    if (is.null(pair_id)) {
      
      ## Wilcoxon rank-sum test
      
      p <- with(df_test,
                wilcox.test(x = variable[group_id == 1],
                            y = variable[group_id == 2],
                       paired = FALSE,
                       ...))$p.value
    } else {
      
      ## Wilcoxon signed-rank test
      
      p <- with(df_test,
                wilcox.test(x = variable[group_id == 1],
                            y = variable[group_id == 2],
                       paired = TRUE,
                       ...))$p.value
      
      
    }
    
    
  } else { stop("For computing p-value, expecting variable of type character, factor, integer, or numeric.")}
  
  return(p)
}
