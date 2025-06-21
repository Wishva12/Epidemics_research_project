library(sf)
library(tidyverse)
library(osmdata)
library(ggrepel)
library(spatstat)
library(ggmap)
library(mapview)
library(maps)
library(leaflet)
library(maptools)
library(httr)
library(jsonlite)
library(devtools)
#install_github('Chrisjb/basemapR')
library(basemapR)
library(sp)
library(tigris)
# ------------------------------------Nasa Precipitation Monthly----------------
library(nasapower)
library(dplyr)
library(purrr)
library(sf)
library(DT)  # for datatable display
library(readr)


# Step 1: Load your dataset
epi <- read_csv("H:\\Tajith Research\\Tajith Report\\Epidemiology_updated_dataset.csv")

# Step 4: Clean mismatched district names in your dataset
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

# Step 2: Load your district shapefile
district_shp <- read_sf("H:/Tajith Research/lka_adm_20220816_shp/lka_admbnda_adm2_slsd_20220816.shp")



# Step 3: Get centroids of each district
district_centroids <- district_shp %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  ) %>%
  select(RDHS = ADM2_EN, lon, lat)  # Rename ADM2_EN to match your RDHS column

# Step 5: Join centroid coordinates to the cleaned dataset
df_with_coords <- df_clean %>%
  left_join(district_centroids, by = "RDHS")
  
# Create a function to get rainfall data for each district
get_monthly_rainfall_temp <- function(lon, lat, start_date, end_date) {
  data <- tryCatch({
    get_power(
      community = "RE",
      lonlat = c(lon, lat),
      dates = c(start_date, end_date),
      temporal_api = "MONTHLY",
      pars = c("PRECTOTCORR", "T2M")
    )
  }, error = function(e) {
    message("Error fetching data for coordinates: ", lon, ", ", lat)
    return(NULL)
  })
  
  return(data)
}

# Set the date range for your analysis
start_date <- "2018-01-01"
end_date <- "2024-05-01"

# Use purrr::pmap_dfr to loop over each RDHS district
rainfall_temp_all_districts <- purrr::pmap_dfr(
  list(district_centroids$lon, district_centroids$lat),
  function(lon, lat) {
    get_monthly_rainfall_temp(lon, lat, start_date, end_date)
  }
)

# Check the result
head(rainfall_temp_all_districts)

# Optional: Save the result as CSV
write.csv(rainfall_temp_all_districts, "H://Tajith Research//monthly_rainfall_temp_all_districts_2018_2024.csv", row.names = FALSE)

#save.image(file = "H:/Tajith Research/SARIMAX.RData")

#------------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(lubridate)
library(ISOweek)  # for converting year-week to date

# Step 1: Create a date column (e.g., first day of the week)
df_with_coords <- df_with_coords %>%
  mutate(
    week_str = paste0(Year, "-W", sprintf("%02d", Week), "-1"),
    date = ISOweek::ISOweek2date(week_str),
    month = month(date),
    month_name = month(date, label = TRUE)
  )

print(df_with_coords, width = Inf)

# Step 2: Aggregate to monthly level (example: summing Dengue_A cases)
df_monthly_DengA <- df_with_coords %>%
  group_by(RDHS, Year, month, month_name, lon, lat) %>%
  summarise(
    Dengue_A_monthly = sum(Dengue_A, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: View the result
head(df_monthly_DengA)


#------------- Convert rainfall_temp_all_districts from wide to long format

library(tidyverse)  
library(dplyr)
library(tidyr)

long_data <- rainfall_temp_all_districts %>%
  pivot_longer(
    cols = JAN:DEC,
    names_to = "month_name",
    values_to = "value"
  ) %>%
  mutate(
    month = match(month_name, toupper(month.abb))  # Convert JAN to 1, FEB to 2, etc.
  )


wide_data <- long_data %>%
  select(LON, LAT, YEAR, month, PARAMETER, value) %>%
  pivot_wider(
    names_from = PARAMETER,
    values_from = value
  )


#------------------Standardize column names and types
# Clean df_monthly_DengueA
df_monthly_DengA <- df_monthly_DengA %>%
  rename(
    YEAR = Year,
    LON = lon,
    LAT = lat
  )

# Clean wide_data
wide_data <- wide_data %>%
  rename(
    YEAR = YEAR,
    LON = LON,
    LAT = LAT
  )
  
wide_data <- wide_data %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb, ordered = TRUE))


#----------------- ensure the month_name columns match in type and levels
df_monthly_DengA <- df_monthly_DengA %>%
  mutate(month_name = factor(month_name, levels = month.abb, ordered = TRUE))


#----------- Merge with df_monthly_DengueA using left_join()

df_merged <- df_monthly_DengA %>%
  left_join(wide_data, by = c("YEAR", "month", "month_name", "LON", "LAT"))


# Optional: Save the result as CSV
write.csv(df_merged, "H://Tajith Research//data_monthly_2018_2024.csv", row.names = FALSE)



#---------------- for All RDHS
library(dplyr)
library(forecast)
library(purrr)

# Get unique RDHS names
rdhs_list <- unique(df_merged$RDHS)

# Initialize a list to store models
sarimax_models <- list()

# Loop through each RDHS
# Loop through each RDHS
for (rdhs in rdhs_list) {
  
  # Filter data for the current RDHS
  df_sub <- df_merged %>%
    filter(RDHS == rdhs) %>%
    filter(!is.na(Dengue_A_monthly), !is.na(PRECTOTCORR), !is.na(T2M)) %>%
    arrange(YEAR, month)
  
  # Check if there is enough data
  if (nrow(df_sub) >= 24) {
    
    n_total <- nrow(df_sub)
    n_train <- floor(0.8 * n_total)
    
    # Split into training and testing sets
    train_data <- df_sub[1:n_train, ]
    test_data <- df_sub[(n_train + 1):n_total, ]
    
    # Create time series object for training
    y_train <- ts(train_data$Dengue_A_monthly, frequency = 12,
                  start = c(min(train_data$YEAR), min(train_data$month)))
    
    # Exogenous variables for training and testing
    xreg_train <- as.matrix(train_data %>% select(PRECTOTCORR, T2M))
    xreg_test <- as.matrix(test_data %>% select(PRECTOTCORR, T2M))
    
    # Fit SARIMAX model
    fit <- auto.arima(y_train, xreg = xreg_train, seasonal = TRUE)
    
    # Forecast on test data
    forecast_result <- forecast(fit, xreg = xreg_test, h = nrow(test_data))
    
    # Store model and results
    sarimax_models[[rdhs]] <- list(
      model = fit,
      forecast = forecast_result,
      test_actual = test_data$Dengue_A_monthly,
      test_dates = test_data %>% transmute(date = as.Date(paste(YEAR, month, "01", sep = "-")))
    )
  } else {
    warning(paste("Skipping", rdhs, "- not enough data"))
  }
}

# Optional: list names of models fitted
names(sarimax_models)

model_summaries <- lapply(sarimax_models, summary)



#---------------------------------------------------------------------------
#-------------- View a brief summary for each model in a loop
for (rdhs_name in names(sarimax_models)) {
  cat("\n--- Model summary for RDHS:", rdhs_name, "---\n")
  print(summary(sarimax_models[[rdhs_name]]))
}


library(dplyr)
library(tidyr)
library(ggplot2)

# Prepare plot data for all districts (example for one district here)
plot_data_list <- lapply(names(sarimax_models), function(rdhs_name) {
  model_info <- sarimax_models[[rdhs_name]]
  
  # Extract training fitted values and dates
  fit <- model_info$model
  train_ts <- model_info$model$x
  fitted_vals <- as.numeric(fitted(fit))
  
  ts_start <- start(train_ts)
  n_train <- length(train_ts)
  
  train_dates <- seq(as.Date(paste0(ts_start[1], "-", ts_start[2], "-01")), by = "month", length.out = n_train)
  
  # Extract forecast values and test dates
  forecast_vals <- as.numeric(model_info$forecast$mean)
  test_dates <- model_info$test_dates$date
  
  # Actual values (train + test)
  actual_vals <- c(as.numeric(train_ts), model_info$test_actual)
  actual_dates <- c(train_dates, test_dates)
  
  # Combine into a dataframe
  df_plot <- tibble(
    RDHS = rdhs_name,
    date = c(actual_dates, train_dates, test_dates),
    Cases = c(actual_vals, fitted_vals, forecast_vals),
    Type = c(
      rep("Actual", length(actual_vals)),
      rep("Fitted", length(fitted_vals)),
      rep("Forecast", length(forecast_vals))
    )
  )
  
  return(df_plot)
})

plot_data <- bind_rows(plot_data_list)

# Make Type a factor to control colors and legend order
plot_data$Type <- factor(plot_data$Type, levels = c("Actual", "Fitted", "Forecast"))

# Plot
ggplot(plot_data, aes(x = date, y = Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ RDHS, scales = "free_y") +
  scale_color_manual(values = c("Actual" = "blue", "Fitted" = "green", "Forecast" = "red")) +
  labs(
    title = "Dengue Cases: Actual, Fitted, and Forecast by District",
    x = "Date",
    y = "Monthly Dengue Cases",
    color = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")



#------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

plot_data <- lapply(names(sarimax_models), function(rdhs_name) {
  model_info <- sarimax_models[[rdhs_name]]
  
  if (!is.null(model_info$model)) {
    # Training data (Actual & Fitted)
    actual_train <- as.numeric(model_info$model$x)
    fitted_train <- as.numeric(fitted(model_info$model))
    start_train <- start(model_info$model$x)
    n_train <- length(actual_train)
    dates_train <- seq(as.Date(paste0(start_train[1], "-", start_train[2], "-01")), 
                       by = "month", length.out = n_train)
    
    # Test data (Actual & Forecast)
    forecasted <- as.numeric(model_info$forecast$mean)
    actual_test <- as.numeric(model_info$test_actual)
    test_dates <- model_info$test_dates$date
    n_test <- length(forecasted)
    
    # Combine fitted + forecast as a continuous predicted line
    predicted <- c(fitted_train, forecasted)
    predicted_dates <- c(dates_train, test_dates)
    
    # Combine actual values for train + test
    actual <- c(actual_train, actual_test)
    actual_dates <- c(dates_train, test_dates)
    
    tibble(
      RDHS = rdhs_name,
      date = c(actual_dates, predicted_dates),
      Cases = c(actual, predicted),
      Type = rep(c("Actual", "Predicted"), times = c(length(actual_dates), length(predicted_dates)))
    )
  } else {
    NULL
  }
}) %>% bind_rows()

# Plot
ggplot(plot_data, aes(x = date, y = Cases, color = Type)) +
  geom_line(size = 1) +
  facet_wrap(~ RDHS, scales = "free_y") +
  scale_color_manual(values = c("Actual" = "blue", "Predicted" = "red")) +
  labs(
    title = "Dengue Cases: Actual vs Predicted (Fitted + Forecast) by District",
    x = "Date",
    y = "Monthly Dengue Cases",
    color = "Legend"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")




#------------------------------------------------------------------------------

library(dplyr)
library(tibble)

model_summary_table <- tibble(
  RDHS = character(),
  p = integer(),
  d = integer(),
  q = integer(),
  P = integer(),
  D = integer(),
  Q = integer(),
  s = integer(),
  AIC = numeric(),
  BIC = numeric()
)

for (rdhs_name in names(sarimax_models)) {
  model_fit <- sarimax_models[[rdhs_name]]$model
  
  if (!is.null(model_fit)) {
    order <- model_fit$arma[c(1,6,2)]       # p,d,q
    seasonal <- model_fit$arma[c(3,7,4)]    # P,D,Q
    period <- model_fit$arma[5]              # s
    
    model_summary_table <- model_summary_table %>% 
      add_row(
        RDHS = rdhs_name,
        p = order[1],
        d = order[2],
        q = order[3],
        P = seasonal[1],
        D = seasonal[2],
        Q = seasonal[3],
        s = period,
        AIC = AIC(model_fit),
        BIC = BIC(model_fit)
      )
  }
}

print(model_summary_table)


library(dplyr)
library(tibble)
library(broom)

# Define empty tibble with columns
model_summary_with_coefs <- tibble(
  RDHS = character(),
  p = integer(),
  d = integer(),
  q = integer(),
  P = integer(),
  D = integer(),
  Q = integer(),
  s = integer(),
  AIC = numeric(),
  BIC = numeric(),
  Exogenous_Coefficients = character()
)

for (rdhs_name in names(sarimax_models)) {
  model_info <- sarimax_models[[rdhs_name]]
  model_fit <- model_info$model
  
  if (!is.null(model_fit)) {
    order <- model_fit$arma[c(1,6,2)]       # p,d,q
    seasonal <- model_fit$arma[c(3,7,4)]    # P,D,Q
    period <- model_fit$arma[5]             # s
    
    coefs <- broom::tidy(model_fit)
    
    xreg_vars <- colnames(model_fit$xreg)
    
    exog_coefs <- coefs %>%
      filter(term %in% xreg_vars) %>%
      select(term, estimate)
    
    if (nrow(exog_coefs) == 0) {
      exog_coefs <- tibble(term = NA_character_, estimate = NA_real_)
    }
    
    exog_coefs_str <- paste0(exog_coefs$term, "=", round(exog_coefs$estimate, 4), collapse = ", ")
    
    model_summary_with_coefs <- model_summary_with_coefs %>% 
      add_row(
        RDHS = rdhs_name,
        p = order[1],
        d = order[2],
        q = order[3],
        P = seasonal[1],
        D = seasonal[2],
        Q = seasonal[3],
        s = period,
        AIC = AIC(model_fit),
        BIC = BIC(model_fit),
        Exogenous_Coefficients = exog_coefs_str
      )
  }
}

print(model_summary_with_coefs)


library(Metrics)  # For rmse, mae, mape functions

# Create an empty results table
model_accuracy <- tibble(
  RDHS = character(),
  RMSE = numeric(),
  MAE = numeric(),
  MAPE = numeric()
)

# Loop through each model and calculate metrics
for (rdhs_name in names(sarimax_models)) {
  model_info <- sarimax_models[[rdhs_name]]
  
  if (!is.null(model_info)) {
    actual <- as.numeric(model_info$test_actual)
    predicted <- as.numeric(model_info$forecast$mean)
    
    if (length(actual) == length(predicted) && length(actual) > 0) {
      rmse_val <- rmse(actual, predicted)
      mae_val <- mae(actual, predicted)
      mape_val <- mape(actual, predicted) * 100  # Convert to percentage

      model_accuracy <- model_accuracy %>% 
        add_row(
          RDHS = rdhs_name,
          RMSE = rmse_val,
          MAE = mae_val,
          MAPE = mape_val
        )
    }
  }
}

# View accuracy metrics
print(model_accuracy)


