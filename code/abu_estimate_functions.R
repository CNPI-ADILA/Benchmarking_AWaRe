


# Code for analysis functions for estimating AWaRe antibiotic use 





#---------------------------------#
# Function to estimate ratio of total DID draws : infection draws for benchmark countries ####
#----------------------------------#

# do for all 1000 draws for both benchmark countries 
# do this element wise (e.g. x1 to x1)

# Function to do this 
divide_draws_f <- function(df1, df2, id_col = "iso3_code") {
  # Basic checks
  if (nrow(df1) != 1 || nrow(df2) != 1) {
    stop("Both data frames must have exactly one row.")
  }
  if (df1[[id_col]] != df2[[id_col]]) {
    stop("Identifier values do not match.")
  }
  
  # Get draw columns (all except id_col)
  draw_cols <- setdiff(names(df1), id_col)
  
  # Element-wise division using base R
  divided <- df1[draw_cols] / df2[draw_cols]
  
  # Add back the identifier column
  result <- cbind(df1[id_col], divided)
  
  # Ensure it's a data frame with appropriate column names
  names(result)[1] <- id_col
  result <- as.data.frame(result)
  
  return(result)
}





#--------------------------------#
# Function multiply relevant ratio matrix by total infection incidence matrix for each country by cluster assignment #####
#---------------------------------#

multiply_draws_matrix_f <- function(main_df, divided_df, id_col = "iso3_code") {
  # Extract the id column separately
  iso3_codes <- main_df[[id_col]]
  
  # Extract numeric draw columns as matrices
  draw_cols <- setdiff(names(main_df), id_col)
  main_mat <- as.matrix(main_df[draw_cols])
  divided_vec <- as.numeric(divided_df[1, draw_cols])
  
  # Element-wise multiplication (broadcasted across rows)
  scaled_mat <- main_mat * matrix(divided_vec, nrow = nrow(main_mat), ncol = length(divided_vec), byrow = TRUE)
  
  # Combine back into a data frame
  result_df <- data.frame(iso3_code = iso3_codes, scaled_mat)
  names(result_df) <- c(id_col, draw_cols)
  
  return(result_df)
}


#-------------------------------------------#
# Function to loop over datasets, filter to analysis countries & order in the right order to prep for summing ####
#-------------------------------------------#

summing_data_prep_f <- function(datasets, analysis_countries, dataset_names = NULL) {
  output_list <- vector("list", length(datasets))
  
  for (i in seq_along(datasets)) {
    data <- datasets[[i]]
    
    # Filter and reorder
    data_filtered <- data[data$iso3_code %in% analysis_countries, ]
    data_filtered <- data_filtered %>%
      arrange(match(iso3_code, analysis_countries))
    
    # Check structure
    if (ncol(data_filtered) == 1001 && "iso3_code" %in% names(data_filtered)) {
      data_filtered <- data_filtered[, !names(data_filtered) %in% "iso3_code"]
      output_list[[i]] <- data_filtered
    } else {
      name <- if (!is.null(dataset_names)) dataset_names[i] else paste("dataset", i)
      warning(sprintf("Dataset %s does not have the correct format. Skipping.", name))
    }
  }
  
  return(output_list)
}




#------------------------------#
# Function to summarise global estimates in DDD ####
#-------------------------------#

summarise_global_stat_f <- function(data, stat_fn, sum_type_label) {
  data %>%
    select(-iso3_code) %>%
    summarise(across(everything(), ~ stat_fn(.x, na.rm = TRUE))) %>%
    mutate(sum_type = sum_type_label) %>%
    relocate(sum_type)
}



# --------------------------------#
# Function to convert to DID ####
#----------------------------------#

did_convert_f <- function(data, covariate_data, pop_column = "pop_totl") {
  data %>%
    left_join(covariate_data %>% select(iso3_code, {{ pop_column }}), by = "iso3_code") %>%
    mutate(across(starts_with("X"), ~ (.x * 1000) / (.data[[pop_column]] * 365)))
}




#----------------------------------#
# Function for summarising estimates (median, 95%CI) #####
#-----------------------------------#


summarise_draws_f <- function(data, estimate_label = "Total DID", scenario_label = "Scenario 1") {
  data %>%
    rowwise() %>%
    mutate(
      draw_median = median(c_across(-iso3_code), na.rm = TRUE),
      draw_lower  = quantile(c_across(-iso3_code), probs = 0.025, na.rm = TRUE),
      draw_upper  = quantile(c_across(-iso3_code), probs = 0.975, na.rm = TRUE),
      estimate = estimate_label,
      scenario = scenario_label
    ) %>%
    ungroup() %>%
    select(iso3_code, scenario, estimate, draw_median, draw_lower, draw_upper)
}




