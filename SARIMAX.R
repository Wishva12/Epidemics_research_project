
library(dplyr)
library(readr)
library(sf)
library(tidyr)
library(forecast)
library(ggplot2)
library(purrr)
library(lubridate)
library(ISOweek)

# Load cleaned dataset with dengue, temperature, rainfall
epi <- read_csv("G:/Academic/4000 Level/Research/Coding from sir/Epidemiology_full_dataset.csv")
epi

# Clean mismatched district names
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

# Load district shapefile
shp <- read_sf("G:/Academic/4000 Level/Research/Coding from sir/lka_admbnda_adm2_slsd_20220816.shp")

# Compute centroids for joining spatially
centroids <- shp %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(lon = st_coordinates(centroid)[, 1],
         lat = st_coordinates(centroid)[, 2]) %>%
  select(RDHS = ADM2_EN, lon, lat)

# Merge centroid coordinates to cleaned data
df_with_coords <- df_clean %>%
  left_join(centroids, by = "RDHS") %>%
  mutate(week_str = paste0(Year, "-W", sprintf("%02d", Week), "-1"),
         date = ISOweek::ISOweek2date(week_str),
         month = month(date),
         month_name = month(date, label = TRUE))

# Aggregate to monthly dengue data
df_monthly_DengA <- df_with_coords %>%
  group_by(RDHS, Year, month, month_name, lon, lat) %>%
  summarise(Dengue_A_monthly = sum(Dengue_A, na.rm = TRUE), .groups = "drop")

# Load NASA POWER monthly climate data
rainfall_temp_all_districts <- read_csv("G:/Academic/4000 Level/Research/Coding from sir/New/monthly_rainfall_temp_all_districts_2018_2024.csv")

# Convert to long format and reshape
long_data <- rainfall_temp_all_districts %>%
  pivot_longer(cols = JAN:DEC, names_to = "month_name", values_to = "value") %>%
  mutate(month = match(month_name, toupper(month.abb)))

wide_data <- long_data %>%
  select(LON, LAT, YEAR, month, PARAMETER, value) %>%
  pivot_wider(names_from = PARAMETER, values_from = value)

# Standardize and merge
df_monthly_DengA <- df_monthly_DengA %>%
  rename(YEAR = Year, LON = lon, LAT = lat) %>%
  mutate(month_name = factor(month_name, levels = month.abb, ordered = TRUE))

wide_data <- wide_data %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb, ordered = TRUE))

df_merged <- df_monthly_DengA %>%
  left_join(wide_data, by = c("YEAR", "month", "month_name", "LON", "LAT"))

# Fit SARIMAX models for all RDHS districts with lagged predictors
rdhs_list <- unique(df_merged$RDHS)
sarimax_models <- list()

for (rdhs in rdhs_list) {
  df_sub <- df_merged %>%
    filter(RDHS == rdhs) %>%
    arrange(YEAR, month) %>%
    mutate(PRECTOTCORR_lag1 = lag(PRECTOTCORR, 1),
           T2M_lag1 = lag(T2M, 1)) %>%
    filter(!is.na(Dengue_A_monthly), !is.na(PRECTOTCORR_lag1), !is.na(T2M_lag1))
  
  if (nrow(df_sub) >= 24) {
    y_ts <- ts(df_sub$Dengue_A_monthly, frequency = 12,
               start = c(min(df_sub$YEAR), min(df_sub$month)))
    xreg_matrix <- as.matrix(df_sub %>% select(PRECTOTCORR_lag1, T2M_lag1))
    fit <- auto.arima(y_ts, xreg = xreg_matrix, seasonal = TRUE)
    sarimax_models[[rdhs]] <- fit
  }
}

# Extract model summaries and significance
test_results <- map_dfr(names(sarimax_models), function(rdhs_name) {
  model <- sarimax_models[[rdhs_name]]
  coefs <- coef(summary(model))
  tibble(
    RDHS = rdhs_name,
    AIC = AIC(model),
    BIC = BIC(model),
    MAPE = mean(abs((model$x - fitted(model)) / model$x), na.rm = TRUE) * 100,
    PREC_pval = coefs["xregPRECTOTCORR_lag1", "Pr(>|t|)"],
    T2M_pval = coefs["xregT2M_lag1", "Pr(>|t|)"]
  )
})

colnames(test_results)

test_results <- test_results %>%
  mutate(Urban = ifelse(rdhs %in% urban_districts, "Urban", "Rural"))


# Add urban classification
urban_districts <- c("Colombo", "Gampaha", "Kandy", "Kurunegala", "Jaffna")




test_results <- map_dfr(names(sarimax_models), function(rdhs_name) {
  model <- sarimax_models[[rdhs_name]]
  coefs <- tryCatch(coef(summary(model)), error = function(e) NULL)
  
  if (!is.null(coefs)) {
    tibble(
      RDHS = rdhs_name,
      AIC = AIC(model),
      BIC = BIC(model),
      MAPE = mean(abs((model$x - fitted(model)) / model$x), na.rm = TRUE) * 100,
      PREC_pval = coefs[grep("PRECTOTCORR", rownames(coefs)), "Pr(>|t|)", drop = TRUE],
      T2M_pval = coefs[grep("T2M", rownames(coefs)), "Pr(>|t|)", drop = TRUE]
    )
  } else {
    tibble(
      RDHS = rdhs_name,
      AIC = NA,
      BIC = NA,
      MAPE = NA,
      PREC_pval = NA,
      T2M_pval = NA
    )
  }
})


test_results <- map_dfr(names(sarimax_models), function(rdhs_name) {
  model <- sarimax_models[[rdhs_name]]
  
  # Try extracting coefficient table safely
  coefs <- tryCatch({
    as.data.frame(coef(summary(model)))
  }, error = function(e) NULL)
  
  # Compute MAPE
  mape_val <- tryCatch({
    mean(abs((model$x - fitted(model)) / model$x), na.rm = TRUE) * 100
  }, error = function(e) NA)
  
  # Extract p-values safely
  prec_p <- tryCatch({
    coefs[grep("PRECTOTCORR", rownames(coefs)), "Pr(>|t|)", drop = TRUE]
  }, error = function(e) NA)
  
  t2m_p <- tryCatch({
    coefs[grep("T2M", rownames(coefs)), "Pr(>|t|)", drop = TRUE]
  }, error = function(e) NA)
  
  tibble(
    RDHS = rdhs_name,
    AIC = AIC(model),
    BIC = BIC(model),
    MAPE = mape_val,
    PREC_pval = prec_p,
    T2M_pval = t2m_p
  )
})


colnames(test_results)

# Save results
write.csv(test_results, "G:/Academic/4000 Level/Research/Coding from sir/New/sarimax_model_summary.csv", row.names = FALSE)





library(forecast)
library(dplyr)
library(lubridate)

# 1) Define a helper to forecast & compute MAPE for one district
forecast_district <- function(df, fit_model, train_frac = 0.8) {
  # df: one‐district data.frame sorted by YEAR, month
  n <- nrow(df)
  # split index
  idx_train <- seq_len(floor(n * train_frac))
  idx_test  <- (max(idx_train) + 1):n
  
  # build ts object on full training period
  y_train <- ts(df$Dengue_A_monthly[idx_train],
                freq = 12,
                start = c(df$YEAR[1], df$month[1]))
  x_train <- as.matrix(df[idx_train, c("PRECTOTCORR","T2M")])
  x_test  <- as.matrix(df[idx_test,  c("PRECTOTCORR","T2M")])
  y_test  <- df$Dengue_A_monthly[idx_test]
  
  # if you haven't passed in fit_model, re‐fit here:
  if (is.null(fit_model)) {
    fit_model <- auto.arima(y_train, xreg = x_train, seasonal = TRUE,
                            stepwise = FALSE, approximation = FALSE)
  }
  
  # forecast over the test period
  fc <- forecast(fit_model, xreg = x_test, h = length(idx_test))
  preds <- as.numeric(fc$mean)
  
  # compute MAPE
  mape <- mean(abs((y_test - preds) / y_test)) * 100
  
  # return a tibble of test‐period results
  tibble(
    Date      = seq(from = as.Date(paste0(df$YEAR[ idx_test ][1],
                                          "-", df$month[ idx_test ][1],
                                          "-01")),
                    by = "month", length.out = length(idx_test)),
    Actual    = y_test,
    Forecast  = preds,
    District  = df$RDHS[1],
    MAPE      = mape
  )
}

# 2) Example: forecast for Colombo
df_col <- df_merged %>%
  filter(RDHS == "Colombo") %>%
  arrange(YEAR, month)

head(df_col)


train_rows <- 1:floor(0.8 * nrow(df_col))

summary(df_col[train_rows, c("Dengue_A_monthly", "PRECTOTCORR", "T2M")])



fit_col_uni <- auto.arima(
  ts(df_col$Dengue_A_monthly[train_rows], freq = 12,
     start = c(df_col$YEAR[1], df_col$month[1])),
  seasonal = TRUE
)
summary(fit_col_uni)




