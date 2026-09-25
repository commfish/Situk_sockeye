# Situk Summary Figures
# authors: Sara Miller
# Last edited: September 2026

# load libraries
library(tidyverse)

Ahrn_data <- read_csv("data/Situk_sockeye_Ahrn_contrib.csv")

# fit regression through the origin
y <- Ahrn_data$Ahrnklin_prop
x <- Ahrn_data$age_0_diff
model <- lm(y ~ x - 1)  # '- 1' removes intercept

# view coefficients
coef(model)

# Calculate R^2 manually for regression through origin
y_pred <- fitted(model)
ss_res <- sum((y - y_pred)^2)       # Residual sum of squares
ss_tot <- sum(y^2)                  # Total sum of squares (no mean subtraction)
r2_origin <- 1 - ss_res / ss_tot

# Output results
cat("Slope:", coef(model), "\n")
cat("R^2 (through origin):", r2_origin, "\n")

