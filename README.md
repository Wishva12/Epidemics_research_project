###Dengue Case Forecasting in Sri Lanka using XGBoost
##Project Overview
This project implements a machine learning approach to forecast monthly dengue fever cases in various districts (RDHS units) across Sri Lanka. The core of the solution is a two-part time series forecasting model built using the XGBoost algorithm.

The model is designed to handle the unique characteristics of dengue case data, such as zero-inflation (months with zero cases) and complex seasonal patterns. It integrates historical dengue case data with external weather variables to provide robust and data-driven forecasts.

##Core Methodology
The forecasting model employs a two-part strategy for each district:

Part 1: Classification Model: An XGBoost classification model predicts whether a given month will have any dengue cases at all. This addresses the challenge of zero-inflated data.

Part 2: Regression Model: For months predicted to have dengue cases, a separate XGBoost regression model forecasts the actual number of cases. The target variable is log-transformed to stabilize variance and improve model performance.

The models are trained on an extensive set of time series features, including:

Lagged Variables: Past values of dengue cases, rainfall, and temperature to capture their carryover effects.

Seasonal Features: Sine and cosine transformations (Fourier terms) of the month to model complex, non-linear seasonal patterns.

Time Trend: A continuous variable representing the passage of time to capture long-term trends.

Weather Features: Monthly total rainfall (PRECTOTCORR) and average temperature (T2M) from the NASA POWER API, along with their multi-month rolling averages to account for cumulative weather effects.
