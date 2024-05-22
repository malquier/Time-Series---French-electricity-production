#Erwann Hotellier - Mattéo Alquier
#Linear Time Series Projet
#Production Indice of electricity


# Loading libraries
library(readr)
library(forecast)
library(stargazer)
library(urca)
library(ggplot2)
library(zoo)
library(tseries)
library(ellipse)

# Defining the path to the CSV file
path_donnees <- "C:/Users/alqui/Downloads/serie_010768228_18052024/tmpZipSerieCsv2990033781039194211/valeurs_mensuelles.csv"

# Reading the CSV file
data1 <- read_delim(path_donnees, delim = ";", escape_double = FALSE, show_col_types = FALSE)

# Removing the 'Codes' column if necessary
data1$Codes <- NULL

# Selecting rows to use (adjust if necessary)
data1 <- data1[4:411, ]

# Renaming columns
colnames(data1) <- c("Date_string", "IPI_elec")

# Converting 'Date_string' column to Date format
data1$Date <- as.Date(paste(data1$Date_string, "01", sep = "-"), format = "%Y-%m-%d")
data1$Date_string <- NULL

# Extracting year and month for each date
data1$annee <- as.numeric(format(data1$Date, format = "%Y"))
data1$mois <- as.numeric(format(data1$Date, format = "%m"))

# Converting 'IPI_elec' to numeric
data1$IPI_elec <- as.numeric(gsub(",", ".", data1$IPI_elec))

# Sorting data by date
data1 <- data1[order(data1$Date), ]
rownames(data1) <- seq(length = nrow(data1))

# Displaying data structure
str(data1)

# Creating the time series object
Data1.ts <- ts(data1$IPI_elec, start = c(1991, 2), end = c(2023, 12), frequency = 12)

# Displaying the time series object
plot(Data1.ts)

# Testing for seasonality
seas_test <- nsdiffs(Data1.ts, test = "ocsb")
if (seas_test > 0) {
  cat("The series exhibits seasonality.\n")
} else {
  cat("The series does not exhibit seasonality.\n")
}

# Differencing the data
diff_ts_data1 <- diff(Data1.ts)

# Displaying the differenced data
plot(diff_ts_data1)

# Augmented Dickey-Fuller test on differenced data
adf_result_diff <- adf.test(diff_ts_data1)

# Displaying the test result
print(adf_result_diff)

# Function for significance tests of ARIMA coefficients
signif <- function(estim) {
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef / se
  pval <- (1 - pnorm(abs(t))) * 2
  return(rbind(coef, se, pval))
}

# Phillips-Perron test to determine the optimal number of differences
pp_test <- pp.test(diff_ts_data1)
if (pp_test$p.value < 0.05) {
  cat("The series requires additional differencing.\n")
} else {
  cat("The series does not require additional differencing.\n")
}

# Centering the differenced data
centered_data1 <- diff_ts_data1 - mean(diff_ts_data1)

# Displaying the centered data
plot(centered_data1)

# Autocorrelation and partial autocorrelation tests on centered data
acf(centered_data1)
pacf(centered_data1)

pmax <- 4
qmax <- 1

models <- list()

for (p in seq(0, pmax)) {
  for (q in seq(0, qmax)) {
    modele <- arima(centered_data1, c(p, 0, q))
    
    # Ljung-Box test
    ljung_box_test <- Box.test(modele$residuals, lag = 7, type = "Ljung-Box", fitdf = 5)
    
    # Significance tests
    significance <- signif(modele)
    
    # AIC and BIC
    aic_bic <- c("AIC" = AIC(modele), "BIC" = BIC(modele))
    
    # Storing the model
    model_name <- paste("ARIMA(", p, ", 0, ", q, ")", sep = "")
    models[[model_name]] <- list(model = modele, ljung_box_test = ljung_box_test, significance = significance, aic_bic = aic_bic)
  }
}

# Initializing AIC and BIC tables
tableau_aic <- matrix(NA, nrow = pmax + 1, ncol = qmax + 1)
tableau_bic <- matrix(NA, nrow = pmax + 1, ncol = qmax + 1)

for (p in seq(0, pmax)) {
  for (q in seq(0, qmax)) {
    modele <- arima(centered_data1, c(p, 0, q))
    
    # Retrieving AIC and BIC
    aic <- AIC(modele)
    bic <- BIC(modele)
    
    # Storing in tables
    tableau_aic[p + 1, q + 1] <- aic
    tableau_bic[p + 1, q + 1] <- bic
  }
}

# Displaying AIC and BIC tables
print("AIC Table:")
print(tableau_aic)

print("BIC Table:")
print(tableau_bic)

# Fitting the ARIMA model
modele <- arima(centered_data1, c(2, 0, 1))

# Summary of ARIMA model results
summary(modele)

# Extracting AIC and BIC values
AIC_value <- modele$aic
BIC_value <- modele$bic

# Displaying AIC and BIC values
print(paste("AIC:", AIC_value))
print(paste("BIC:", BIC_value))

# Plotting the autocorrelation function of residuals
par(mfrow = c(1, 1))
acf(modele$residuals, 50, main = "Autocorrelation Function of Residuals")

# Plotting the histogram of residuals
hist(modele$residuals, freq = FALSE, main = "Histogram of Residuals", xlab = "Residuals", ylab = "Density", breaks = 10)

# Plotting the normal density curve
curve(dnorm(x, mean = mean(modele$residuals), sd = sd(modele$residuals)), 
      from = min(modele$residuals), to = max(modele$residuals), 
      col = "blue", lwd = 2, add = TRUE)

# Legend
legend("topright", legend = c("Normal Density"), col = c("blue"), lwd = 2)

# Ljung-Box test for autocorrelation of residuals
portmanteau_test <- Box.test(modele$residuals, lag = 7, type = "Box-Pierce")

# Displaying the test result
print("Portmanteau test for residuals autocorrelation:")
print(portmanteau_test)

######################################################
#                Part 3: Forecasting                 #
######################################################
# Fitting the ARIMA model

data1_2015 <- subset(data1, Date >= as.Date("2015-01-01"))

# Creating the time series object from 2015
Data1_2015.ts <- ts(data1_2015$IPI_elec, start = c(2015, 1), end = c(2023, 12), frequency = 12)

modele <- arima(Data1_2015.ts, c(2, 1, 1))

# Forecasting for the next 4 months
Prevision <- predict(modele, n.ahead = 4)

# Calculating variances and covariance for forecasts
var_prev_1 <- modele$sigma2
var_prev_2 <- modele$sigma2 * (1 + modele$coef[1]^2)
Cov_prev_1_2 <- modele$sigma2 * modele$coef[1]

# Calculating upper and lower bounds of confidence intervals using theoretical formula
Bound_sup_formule <- rbind(Prevision$pred[1] + 1.96 * sqrt(var_prev_1), 
                           Prevision$pred[2] + 1.96 * sqrt(var_prev_2))
Bound_inf_formule <- rbind(Prevision$pred[1] - 1.96 * sqrt(var_prev_1), 
                           Prevision$pred[2] - 1.96 * sqrt(var_prev_2))

# Calculating upper and lower bounds of confidence intervals using package
Bound_sup_package <- rbind(Prevision$pred[1] + 1.96 * Prevision$se[1], 
                           Prevision$pred[2] + 1.96 * Prevision$se[2])
Bound_inf_package <- rbind(Prevision$pred[1] - 1.96 * Prevision$se[1], 
                           Prevision$pred[2] - 1.96 * Prevision$se[2])

# Displaying results of theoretical and package bounds
print("Bound_sup_formule:")
print(Bound_sup_formule)
print("Bound_inf_formule:")
print(Bound_inf_formule)
print("Bound_sup_package:")
print(Bound_sup_package)
print("Bound_inf_package:")
print(Bound_inf_package)

# Constructing the Sigma matrix with standard deviations and correlation coefficient
sigma_x <- sqrt(var_prev_1)
sigma_y <- sqrt(var_prev_2)
rho <- Cov_prev_1_2 / (sigma_x * sigma_y)
Sigma <- matrix(c(sigma_x^2, rho * sigma_x * sigma_y, rho * sigma_x * sigma_y, sigma_y^2), nrow = 2)

# Defining eigenvalues and eigenvectors of the Sigma matrix
eigenvalues <- eigen(Sigma)$values
eigenvectors <- eigen(Sigma)$vectors

# Plotting the bivariate confidence ellipse
ell <- ellipse(rho, scale = c(sigma_x, sigma_y), centre = Prevision$pred[1:2], level = 0.95, npoints = 10000)
pdf("Ellipse.pdf")
plot(ell, type = 'l', main = "95% Confidence Ellipse", xlab = "Forecast at T+1", ylab = "Forecast at T+2")
points(Prevision$pred[1], Prevision$pred[2], col = 'red', pch = 19)
dev.off()

# Graphical display of the confidence region
plot(forecast(modele, h = 4))


#Creation of the unit disc

arima201 <- arima(centered_data1, c(2,0,1))

adj_r2 <- function(model){
  p <- model$arma[1]
  q <- model$arma[2]
  n <- model$nobs-max(p,q)
  ss_res <- sum(model$residuals^2)
  ss_tot <- sum(centered_data1[-c(1:max(p,q))]^2)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}
adj_r2(arima201)