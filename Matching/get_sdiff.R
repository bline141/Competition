

#' Standardized difference of a categorical or numerical variable between two groups.
#' 
#' Function to compute the standardized difference of a variable
#' between two groups. 
#'
#' @param group Vector containing grouping variable. Either numeric or factor so that the order of groups is respected.
#' @param variable Vector containing the values of the variable to compute the standardized difference on. Same length as group. 
#' Either numeric, categorical, or factor.
#' @param type Optional string specifying whether the variable is categorical ("cat") or numeric ("num"). 
#' Default is categorical if class(variable) %in% c("character", "factor") and numeric if class(variable) %in% c("integer", "numeric").
#' 
#' @return Numeric value if there are two distinct groups, otherwise NA_real_. 
#' Additionally, if the variable is categorical, both groups must have all categories present, otherwise the result is NA_real_.
#' 
#' @example
#' df_example <- data.frame(therapy = factor(rep(c("A", "B"), each = 20)),
#'                          age     = c(rnorm(n = 20, mean = 40, sd = 4),
#'                                      rnorm(n = 20, mean = 30, sd = 4)),
#'                          eye_colour = c(rep(c("blue", "brown", "green"), times = c(6, 6, 8)),
#'                                         rep(c("green", "blue", "brown"), times = c(12, 4, 4))),
#'                          stringsAsFactors = FALSE)
#' 
#' get_sdiff(group = df_example$therapy, variable = df_example$age)
#' get_sdiff(group = df_example$therapy, variable = df_example$eye_colour)
#' 
#' # NA because not all categories appear for each group
#' df_example$eye_colour[1:6] <- "brown"
#' get_sdiff(group = df_example$therapy, variable = df_example$eye_colour)


get_sdiff <- function(group, variable, type = NULL, na.rm = FALSE) {
  
  
    
    if (is.null(type)) {
    type <- ifelse(class(variable) %in% c("character", "factor"),
                   "cat",
                   ifelse(class(variable) %in% c("integer", "numeric"),
                          "num",
                          stop("For get_sdiff() need variable to be of class character, factor, integer, or numeric. Specify what to do with new class.")))
  } else if (!type %in% c("cat", "num")) {
    stop("For get_sdiff() you need to specify variable of type 'cat' or 'num'.")
  }
  
  group_id <- as.numeric(factor(group))
  
  # if there's only one group_id, can't compute a standardized difference
  if (length(unique(group_id)) == 1) {
    return(NA_real_)
  }
  
  ## Categorical -------------------------------------------------
  if (type == "cat") {
    
    ## (A) create proportions of variable per group_id
    ## (B) compute standardized difference
    
    ## (A)
    variable_char <- as.character(variable)
    levs     <- unique(variable_char)
    variable_factor <- factor(variable_char, levels = levs)
    
    tbl   <- table(group_id, variable_factor)
    props <- tbl / rowSums(tbl)
    
    p1 <- props[1, ]
    p2 <- props[2, ]
    
    
    # ## (B)
     if (any(c(p1, p2) == 0) | any(is.nan(c(p1, p2)))) {
       return(NA_real_)
     }
    
    cov_p1 <- diag(p1) - outer(p1, p1)
    cov_p2 <- diag(p2) - outer(p2, p2)
    
    
    S <- (cov_p1 + cov_p2)/2
    
    # don't need all categories, drop one (dropped one = 1 - sum(all others))
    S  <- S[-1, -1]
    p1 <- p1[-1]
    p2 <- p2[-1]

    sdiff <- sqrt( t(p1 - p2) %*% solve(S) %*% (p1 - p2) )
    
  # Continuous ----------------------------------------------------
  } else {
    
    mu    <- c(mean(variable[group_id == 1], na.rm = na.rm),
               mean(variable[group_id == 2], na.rm = na.rm))
    sigma <- c(sd(variable[group_id == 1], na.rm = na.rm),
               sd(variable[group_id == 2], na.rm = na.rm))
  
    sdiff <- (mu[1] - mu[2])/sqrt(sum(sigma^2)/2)
    
  }
  
  return(as.numeric(sdiff))
  
}

