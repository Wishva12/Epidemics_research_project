library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)
library(ISOweek)
library(nasapower)
library(xgboost)
library(Metrics)
library(patchwork)

epi <- read_csv("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/Epidemiology_full_dataset.csv")

df_clean <- epi %>%
  mutate(RDHS = case_when(
    RDHS == "NuwaraEliya" ~ "Nuwara Eliya",
    RDHS == "Kalmune" ~ "Kalmunai",
    RDHS == "paha" ~ "Gampaha",
    RDHS == "Gamapaha" ~ "Gampaha",
    RDHS == "Dr Sudath SamaraweeraWERColombo" ~ "Colombo",
    RDHS == "Anuradhapur" ~ "Anuradhapura",
    RDHS == "Mulativu" ~ "Mullaitivu",
    RDHS == "Baticaloa" ~ "Batticaloa",
    RDHS == "Puttlam" ~ "Puttalam",
    TRUE ~ RDHS
  ))

district_shp <- read_sf("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/lka_admbnda_adm2_slsd_20220816.shp")

district_centroids <- district_shp %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  ) %>%
  select(RDHS = ADM2_EN, lon, lat)

df_with_coords <- df_clean %>%
  left_join(district_centroids, by = "RDHS") %>%
  filter(!is.na(lon), !is.na(lat), !is.na(Year), !is.na(Week), Week > 0, Week <= 53) %>%
  mutate(
    week_str = paste0(Year, "-W", sprintf("%02d", Week), "-1"),
    date = ISOweek::ISOweek2date(week_str),
    month = month(date),
    month_name = month(date, label = TRUE)
  )

df_monthly_all <- df_with_coords %>%
  group_by(RDHS, Year, month, month_name, lon, lat) %>%
  summarise(Dengue_A_monthly = sum(Dengue_A, na.rm = TRUE), .groups = "drop")

get_monthly_rainfall_temp <- function(lon, lat, start_date, end_date) {
  tryCatch({
    get_power(
      community = "RE",
      lonlat = c(lon, lat),
      dates = c(start_date, end_date),
      temporal_api = "MONTHLY",
      pars = c("PRECTOTCORR", "T2M")
    )
  }, error = function(e) {
    message("Error fetching data for lon: ", lon, ", lat: ", lat, " - ", e)
    return(NULL)
  })
}

start_date <- "2018-01-01"
end_date <- "2024-05-01"
unique_coords <- df_with_coords %>% distinct(RDHS, lon, lat)

get_all_weather_data <- function(coords_df, start_date, end_date) {
  all_weather <- list()
  
  for (i in 1:nrow(coords_df)) {
    rdhs <- coords_df$RDHS[i]
    lon <- coords_df$lon[i]
    lat <- coords_df$lat[i]
    
    message("Fetching weather data for ", rdhs, " (", i, "/", nrow(coords_df), ")")
    
    weather_data <- get_monthly_rainfall_temp(lon, lat, start_date, end_date)
    
    if (!is.null(weather_data)) {
      weather_data$RDHS <- rdhs
      all_weather[[rdhs]] <- weather_data
    }
    
    Sys.sleep(1)
  }
  
  return(all_weather)
}

all_weather_data <- get_all_weather_data(unique_coords, start_date, end_date)

process_weather_data <- function(weather_list) {
  all_processed <- list()
  
  for (rdhs in names(weather_list)) {
    weather_data <- weather_list[[rdhs]]
    
    long_data <- weather_data %>%
      pivot_longer(cols = JAN:DEC, names_to = "month_name", values_to = "value") %>%
      mutate(month = match(month_name, toupper(month.abb)))
    
    wide_data <- long_data %>%
      select(LON, LAT, YEAR, month, PARAMETER, value) %>%
      pivot_wider(names_from = PARAMETER, values_from = value) %>%
      mutate(
        month_name = factor(month.abb[month], levels = month.abb, ordered = TRUE),
        RDHS = rdhs
      )
    
    all_processed[[rdhs]] <- wide_data
  }
  
  return(bind_rows(all_processed))
}

processed_weather <- process_weather_data(all_weather_data)

df_monthly_all <- df_monthly_all %>%
  rename(YEAR = Year, LON = lon, LAT = lat) %>%
  mutate(month_name = factor(month_name, levels = month.abb, ordered = TRUE))

df_merged_all <- df_monthly_all %>%
  left_join(processed_weather, by = c("RDHS", "YEAR", "month", "month_name", "LON", "LAT"))

# Enhanced outlier capping function
cap_outliers <- function(x, lower_percentile = 0.05, upper_percentile = 0.95) {
  if (all(is.na(x))) return(x)
  
  lower_bound <- quantile(x, lower_percentile, na.rm = TRUE)
  upper_bound <- quantile(x, upper_percentile, upper_percentile, na.rm = TRUE)
  
  x[x < lower_bound] <- lower_bound
  x[x > upper_bound] <- upper_bound
  
  return(x)
}

# Enhanced lag features function with extended lags
add_enhanced_lag_features <- function(df, target_vars, n_lags = 6) {
  df_sorted <- df %>% arrange(RDHS, YEAR, month)
  
  for (target_var in target_vars) {
    for (i in 1:n_lags) {
      lag_name <- paste0(target_var, "_lag", i)
      df_sorted <- df_sorted %>%
        group_by(RDHS) %>%
        mutate(!!lag_name := lag(.data[[target_var]], i)) %>%
        ungroup()
    }
  }
  return(df_sorted)
}

# Enhanced temporal features function
add_enhanced_temporal_features <- function(df) {
  df %>%
    arrange(RDHS, YEAR, month) %>%
    group_by(RDHS) %>%
    mutate(
      # Basic temporal features
      month_sin = sin(2 * pi * month / 12),
      month_cos = cos(2 * pi * month / 12),
      quarter = ceiling(month / 3),
      
      # Enhanced temporal features
      month_sin2 = sin(4 * pi * month / 12),  # 6-month cycle
      month_cos2 = cos(4 * pi * month / 12),
      
      # Trend features
      time_index = row_number(),
      time_squared = time_index^2,
      
      # Moving averages for trend
      dengue_ma3 = zoo::rollmean(Dengue_A_monthly, k = 3, fill = NA, align = "right"),
      dengue_ma6 = zoo::rollmean(Dengue_A_monthly, k = 6, fill = NA, align = "right"),
      
      # Weather moving averages
      rainfall_ma3 = zoo::rollmean(PRECTOTCORR, k = 3, fill = NA, align = "right"),
      temp_ma3 = zoo::rollmean(T2M, k = 3, fill = NA, align = "right"),
      
      # Seasonal indicators
      is_monsoon = ifelse(month %in% c(5, 6, 7, 8, 9, 10), 1, 0),
      is_dry_season = ifelse(month %in% c(1, 2, 3), 1, 0),
      
      # Year effects
      year_normalized = (YEAR - min(YEAR, na.rm = TRUE)) / (max(YEAR, na.rm = TRUE) - min(YEAR, na.rm = TRUE))
    ) %>%
    ungroup()
}

# Function to add interaction terms
add_interaction_features <- function(df) {
  df %>%
    mutate(
      # Weather interactions
      temp_rainfall_interaction = T2M * PRECTOTCORR,
      temp_rainfall_lag1 = T2M_lag1 * PRECTOTCORR_lag1,
      temp_rainfall_lag2 = T2M_lag2 * PRECTOTCORR_lag2,
      
      # Seasonal-weather interactions
      monsoon_rainfall = is_monsoon * PRECTOTCORR,
      monsoon_temp = is_monsoon * T2M,
      dry_season_temp = is_dry_season * T2M,
      
      # Lag interactions
      dengue_temp_lag1 = Dengue_A_monthly_lag1 * T2M_lag1,
      dengue_rainfall_lag1 = Dengue_A_monthly_lag1 * PRECTOTCORR_lag1,
      
      # Temporal interactions
      month_sin_temp = month_sin * T2M,
      month_cos_rainfall = month_cos * PRECTOTCORR,
      
      # Higher order terms
      temp_squared = T2M^2,
      rainfall_squared = PRECTOTCORR^2,
      temp_cubed = T2M^3,
      
      # Ratio features
      temp_rainfall_ratio = ifelse(PRECTOTCORR != 0, T2M / PRECTOTCORR, 0)
    )
}

# Process data with all enhancements
df_all_enhanced <- df_merged_all %>%
  filter(!is.na(Dengue_A_monthly), !is.na(PRECTOTCORR), !is.na(T2M)) %>%
  group_by(RDHS) %>%
  mutate(
    # Apply outlier capping
    Dengue_A_monthly = cap_outliers(Dengue_A_monthly),
    PRECTOTCORR = cap_outliers(PRECTOTCORR),
    T2M = cap_outliers(T2M)
  ) %>%
  ungroup() %>%
  add_enhanced_lag_features(target_vars = c("Dengue_A_monthly", "PRECTOTCORR", "T2M"), n_lags = 6) %>%
  add_enhanced_temporal_features() %>%
  add_interaction_features() %>%
  drop_na()

# Enhanced XGBoost training function
train_enhanced_xgboost_model <- function(data, min_data_points = 30) {
  if (nrow(data) < min_data_points) return(NULL)
  
  # Enhanced feature set
  feature_cols <- c(
    # Basic weather features
    "PRECTOTCORR", "T2M",
    
    # Extended lag features (6 lags instead of 3)
    paste0("PRECTOTCORR_lag", 1:6),
    paste0("T2M_lag", 1:6),
    paste0("Dengue_A_monthly_lag", 1:6),
    
    # Enhanced temporal features
    "month_sin", "month_cos", "month_sin2", "month_cos2",
    "quarter", "time_index", "time_squared", "year_normalized",
    "is_monsoon", "is_dry_season",
    
    # Moving averages (excluding current month to prevent leakage)
    "dengue_ma3", "dengue_ma6", "rainfall_ma3", "temp_ma3",
    
    # Interaction features
    "temp_rainfall_interaction", "temp_rainfall_lag1", "temp_rainfall_lag2",
    "monsoon_rainfall", "monsoon_temp", "dry_season_temp",
    "dengue_temp_lag1", "dengue_rainfall_lag1",
    "month_sin_temp", "month_cos_rainfall",
    "temp_squared", "rainfall_squared", "temp_cubed",
    "temp_rainfall_ratio"
  )
  
  # Filter features that exist in the data
  available_features <- intersect(feature_cols, names(data))
  
  data <- data %>% arrange(YEAR, month)
  n_total <- nrow(data)
  n_train <- floor(0.8 * n_total)
  if (n_total - n_train < 5) n_train <- n_total - 5
  
  train_data <- data[1:n_train, ]
  test_data <- data[(n_train + 1):n_total, ]
  
  # Remove any remaining NA values in feature columns
  train_clean <- train_data[complete.cases(train_data[available_features]), ]
  test_clean <- test_data[complete.cases(test_data[available_features]), ]
  
  if (nrow(train_clean) < 10 || nrow(test_clean) < 3) return(NULL)
  
  dtrain <- xgb.DMatrix(data = as.matrix(train_clean[available_features]), 
                        label = train_clean$Dengue_A_monthly)
  dtest <- xgb.DMatrix(data = as.matrix(test_clean[available_features]), 
                       label = test_clean$Dengue_A_monthly)
  
  # Enhanced parameters with lower learning rate and more regularization
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = 5,               # Slightly reduced to prevent overfitting
    eta = 0.05,                  # Lower learning rate
    subsample = 0.7,             # More aggressive subsampling
    colsample_bytree = 0.7,      # More feature sampling
    alpha = 1,                   # L1 regularization
    lambda = 2,                  # L2 regularization
    min_child_weight = 3,        # Minimum sum of instance weight in child
    gamma = 1                    # Minimum loss reduction for split
  )
  
  set.seed(123)
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 500,               # More rounds with lower learning rate
    watchlist = list(train = dtrain, test = dtest),
    early_stopping_rounds = 25,  # More patience
    verbose = 0
  )
  
  train_pred <- predict(xgb_model, dtrain)
  test_pred <- predict(xgb_model, dtest)
  
  # Calculate metrics
  train_r2 <- 1 - sum((train_clean$Dengue_A_monthly - train_pred)^2) /
    sum((train_clean$Dengue_A_monthly - mean(train_clean$Dengue_A_monthly))^2)
  
  test_r2 <- 1 - sum((test_clean$Dengue_A_monthly - test_pred)^2) /
    sum((test_clean$Dengue_A_monthly - mean(test_clean$Dengue_A_monthly))^2)
  
  return(list(
    model = xgb_model,
    train_data = train_clean,
    test_data = test_clean,
    train_pred = train_pred,
    test_pred = test_pred,
    feature_cols = available_features,
    metrics = data.frame(
      District = unique(data$RDHS),
      Train_RMSE = rmse(train_clean$Dengue_A_monthly, train_pred),
      Test_RMSE = rmse(test_clean$Dengue_A_monthly, test_pred),
      Test_MAE = mae(test_clean$Dengue_A_monthly, test_pred),
      Test_MAPE = mape(test_clean$Dengue_A_monthly, test_pred) * 100,
      Train_R2 = train_r2,
      Test_R2 = test_r2,
      N_Train = nrow(train_clean),
      N_Test = nrow(test_clean),
      N_Features = length(available_features)
    )
  ))
}

train_all_enhanced_models <- function(df) {
  rdhs_units <- unique(df$RDHS)
  all_results <- list()
  all_metrics <- list()
  
  for (rdhs in rdhs_units) {
    message("Training enhanced model for ", rdhs)
    
    rdhs_data <- df %>% filter(RDHS == rdhs)
    result <- train_enhanced_xgboost_model(rdhs_data)
    
    if (!is.null(result)) {
      all_results[[rdhs]] <- result
      all_metrics[[rdhs]] <- result$metrics
    } else {
      message("Insufficient data for ", rdhs)
    }
  }
  
  return(list(
    models = all_results,
    metrics = bind_rows(all_metrics)
  ))
}

# Train enhanced models
cat("Training enhanced models with improved features...\n")
all_enhanced_models <- train_all_enhanced_models(df_all_enhanced)

if (nrow(all_enhanced_models$metrics) > 0) {
  cat("\n=== ENHANCED MODEL PERFORMANCE METRICS ===\n")
  print(all_enhanced_models$metrics %>% arrange(Test_RMSE))
} else {
  print("No enhanced models were successfully trained. Check data availability.")
}

# Enhanced plotting function
create_enhanced_prediction_plot <- function(result, rdhs_name) {
  train_plot_df <- result$train_data %>%
    mutate(
      Date = as.Date(paste(YEAR, month, "01", sep = "-")),
      Actual = Dengue_A_monthly,
      Predicted = result$train_pred,
      Set = "Train",
      Residual = Actual - Predicted
    )
  
  test_plot_df <- result$test_data %>%
    mutate(
      Date = as.Date(paste(YEAR, month, "01", sep = "-")),
      Actual = Dengue_A_monthly,
      Predicted = result$test_pred,
      Set = "Test",
      Residual = Actual - Predicted
    )
  
  plot_df <- bind_rows(train_plot_df, test_plot_df)
  split_date <- min(test_plot_df$Date)
  
  ggplot(plot_df, aes(x = Date)) +
    geom_line(aes(y = Actual, color = "Actual"), size = 0.8) +
    geom_line(aes(y = Predicted, color = "Predicted"), size = 0.8, linetype = "dashed") +
    geom_vline(xintercept = as.numeric(split_date), linetype = "dotdash", 
               color = "blue", size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c("Actual" = "black", "Predicted" = "red")) +
    labs(
      title = paste(rdhs_name),
      subtitle = paste("RMSE:", round(result$metrics$Test_RMSE, 2)),
      y = "Cases", x = "Date", color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = -5),
      plot.title = element_text(size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# Create enhanced plots
if (nrow(all_enhanced_models$metrics) > 0) {
  all_districts <- all_enhanced_models$metrics %>%
    arrange(Test_RMSE) %>%  # Sort by RMSE (lower is better)
    pull(District)
  
  plot_list <- list()
  for (i in 1:length(all_districts)) {
    district <- all_districts[i]
    if (district %in% names(all_enhanced_models$models)) {
      plot_list[[i]] <- create_enhanced_prediction_plot(all_enhanced_models$models[[district]], district)
    }
  }
  
  n_plots <- length(plot_list)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)
  
  if (n_plots > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows)
    
    cat("Displaying enhanced panel chart for all", n_plots, "RDHS units\n")
    cat("Grid layout:", n_rows, "rows x", n_cols, "columns\n\n")
    
    print(combined_plot)
  }
} else {
  print("No enhanced models available for plotting.")
}

# Enhanced feature importance plot
if (nrow(all_enhanced_models$metrics) > 0) {
  best_district <- all_enhanced_models$metrics %>%
    arrange(Test_RMSE) %>%  # Best = lowest RMSE
    head(1) %>%
    pull(District)
  
  if (best_district %in% names(all_enhanced_models$models)) {
    importance_data <- xgb.importance(model = all_enhanced_models$models[[best_district]]$model)
    
    top_features <- as.data.frame(importance_data) %>%
      arrange(desc(Gain)) %>%
      head(15) %>%  # Show top 15 features
      mutate(
        Feature = factor(Feature, levels = rev(Feature)),
        Feature_Type = case_when(
          grepl("lag", Feature) ~ "Lag Features",
          grepl("interaction|ratio|squared|cubed", Feature) ~ "Interaction/Non-linear",
          grepl("sin|cos|month|quarter|season|monsoon", Feature) ~ "Temporal/Seasonal",
          grepl("ma", Feature) ~ "Moving Average",
          grepl("PRECTOTCORR|T2M", Feature) ~ "Weather",
          TRUE ~ "Other"
        )
      )
    
    importance_plot <- ggplot(top_features, aes(x = Feature, y = Gain, fill = Feature_Type)) +
      geom_col() +
      coord_flip() +
      scale_fill_viridis_d(name = "Feature Type") +
      labs(
        title = paste("Enhanced Feature Importance -", best_district, "District"),
        subtitle = paste("Model RMSE:", round(all_enhanced_models$models[[best_district]]$metrics$Test_RMSE, 2)),
        x = "Feature",
        y = "Gain"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    
    print(importance_plot)
  }
} else {
  print("No enhanced models available for feature importance analysis.")
}

# Performance comparison if original models exist
if (exists("all_models") && nrow(all_models$metrics) > 0 && nrow(all_enhanced_models$metrics) > 0) {
  cat("\n=== PERFORMANCE COMPARISON ===\n")
  
  comparison <- all_models$metrics %>%
    select(District, Original_Test_RMSE = Test_RMSE, Original_Test_R2 = Test_R2) %>%
    inner_join(
      all_enhanced_models$metrics %>%
        select(District, Enhanced_Test_RMSE = Test_RMSE, Enhanced_Test_R2 = Test_R2),
      by = "District"
    ) %>%
    mutate(
      RMSE_Improvement = Original_Test_RMSE - Enhanced_Test_RMSE,
      RMSE_Improvement_Pct = (RMSE_Improvement / Original_Test_RMSE) * 100,
      R2_Improvement = Enhanced_Test_R2 - Original_Test_R2
    ) %>%
    arrange(desc(RMSE_Improvement_Pct))
  
  cat("Districts with largest RMSE improvements:\n")
  print(comparison)
  
  cat("\nOverall improvements:\n")
  cat("Average RMSE improvement:", round(mean(comparison$RMSE_Improvement, na.rm = TRUE), 2), "\n")
  cat("Average RMSE improvement (%):", round(mean(comparison$RMSE_Improvement_Pct, na.rm = TRUE), 1), "%\n")
  cat("Average R² improvement:", round(mean(comparison$R2_Improvement, na.rm = TRUE), 3), "\n")
}

cat("\n=== ENHANCED MODEL SUMMARY ===\n")
if (nrow(all_enhanced_models$metrics) > 0) {
  best_district <- all_enhanced_models$metrics %>%
    arrange(Test_RMSE) %>%
    head(1) %>%
    pull(District)
  
  cat("Total RDHS units with enhanced models:", nrow(all_enhanced_models$metrics), "\n")
  cat("Best performing district (lowest RMSE):", best_district, "\n")
  cat("Best Test RMSE:", round(min(all_enhanced_models$metrics$Test_RMSE, na.rm = TRUE), 2), "\n")
  cat("Average Test R²:", round(mean(all_enhanced_models$metrics$Test_R2, na.rm = TRUE), 3), "\n")
  cat("Average Test RMSE:", round(mean(all_enhanced_models$metrics$Test_RMSE, na.rm = TRUE), 2), "\n")
  cat("Average number of features used:", round(mean(all_enhanced_models$metrics$N_Features, na.rm = TRUE), 0), "\n")
} else {
  cat("No enhanced models were successfully trained.\n")
  cat("Please check if there is sufficient data for each RDHS unit.\n")
}


# --- START: CODE TO CALCULATE OVERALL MODEL ACCURACY ---

# Calculate Overall Model Performance Metrics
if (exists("all_enhanced_models") && length(all_enhanced_models$models) > 0) {
  
  # Combine all test set predictions and actuals from every model
  all_test_actuals <- do.call(c, lapply(all_enhanced_models$models, function(res) res$test_data$Dengue_A_monthly))
  all_test_predictions <- do.call(c, lapply(all_enhanced_models$models, function(res) res$test_pred))
  
  # Ensure there are values to calculate
  if (length(all_test_actuals) > 0 && length(all_test_predictions) > 0) {
    cat("\n=== OVERALL MODEL ACCURACY ACROSS ALL DISTRICTS ===\n")
    
    # Calculate overall regression metrics
    overall_rmse <- rmse(all_test_actuals, all_test_predictions)
    overall_mae <- mae(all_test_actuals, all_test_predictions)
    overall_mape <- mape(all_test_actuals, all_test_predictions) * 100
    overall_r2 <- 1 - sum((all_test_actuals - all_test_predictions)^2) / sum((all_test_actuals - mean(all_test_actuals))^2)
    
    # Print the overall metrics
    cat("Overall Test RMSE:", round(overall_rmse, 2), "\n")
    cat("Overall Test MAE:", round(overall_mae, 2), "\n")
    cat("Overall Test MAPE:", round(overall_mape, 2), "%\n")
    cat("Overall Test R-squared (R²):", round(overall_r2, 3), "\n")
    cat("---------------------------------------------------\n")
    cat("Note: These metrics are calculated by combining the test sets from all individual district models to provide a single performance overview.\n")
    
  } else {
    cat("\nCould not calculate overall accuracy. No test data found in the model results.\n")
  }
}

# --- END: CODE TO CALCULATE OVERALL MODEL ACCURACY ---
