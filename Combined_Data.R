# Load required libraries
library(MASS)
library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(lmtest)  # For likelihood ratio test
library(cluster)     # For clustering algorithms
library(factoextra)  # For visualizing clusters
library(ggrepel)
library(mclust)

# Load dataset
data <- read.csv("~/Desktop/UC Berkeley/My Data/Inbred Lines - Life History:Resistance Assays/Signe's IL Data/Infection_Growth.csv")

# Consistently transform Dose (choose one: log10 or natural log)
data$Log_Dose <- log10(data$Dose)  # Base-10 log transformation (same in both models)

# Convert Population_ID to a factor
data$Population_ID <- as.factor(data$Population_ID)

# Define custom colors
custom_colors <- c(
  "Stock" = "black", "IL-1" = "#D73027", "IL-2" = "#FC8D59", "IL-4" = "#EFC000", 
  "IL-7" = "#91BFDB", "IL-9" = "#4575B4", "IL-10" = "#313695", "IL-14" = "#8E44AD", 
  "IL-17" = "#D73092", "IL-19" = "#666666", "IL-24" = "#1D91C0", "IL-29" = "#8C510A", 
  "IL-36" = "#A6DBA0"
)

# ---- FIT MODELS ----

# ✅ Logistic Regression (Binomial GLM)
glm_logit <- glm(cbind(Infection_Status, Total_Exposed - Infection_Status) ~ Log_Dose * Population_ID, 
                 family = binomial, 
                 data = data)

glm_logit
summary(glm_logit)
anova(glm_logit, test="Chisq")

# ✅ Negative Binomial GLM (accounts for overdispersion)
nb_glm <- glm.nb(Infection_Status ~ Log_Dose + Population_ID + offset(log(Total_Exposed)), 
                 data = data,
                 control = glm.control(maxit = 100))

nb_glm
summary(glm_logit)
anova(nb_glm, test="Chisq")

# ---- MODEL COMPARISON ----

# Log-Likelihood comparison
logLik_glm_logit <- logLik(glm_logit)
logLik_nb_glm <- logLik(nb_glm)

# AIC and BIC comparison
aic_glm_logit <- AIC(glm_logit)
bic_glm_logit <- BIC(glm_logit)

aic_nb_glm <- AIC(nb_glm)
bic_nb_glm <- BIC(nb_glm)

# ✅ Create a comparison table
model_comparison <- data.frame(
  Model = c("Logit GLM", "Negative Binomial GLM"),
  LogLik = c(logLik_glm_logit, logLik_nb_glm),
  AIC = c(aic_glm_logit, aic_nb_glm),
  BIC = c(bic_glm_logit, bic_nb_glm)
)

# Model    LogLik      AIC      BIC
# 1             Logit GLM -157.2767 366.5533 423.0874
# 2 Negative Binomial GLM -162.4594 354.9187 387.5345

# ---- LIKELIHOOD RATIO TEST ----

# Fit a reduced Negative Binomial model (without Population_ID)
nb_glm_reduced <- glm.nb(Infection_Status ~ Log_Dose + offset(log(Total_Exposed)), 
                         data = data, 
                         control = glm.control(maxit = 100))

# Compare full vs reduced Negative Binomial model
lr_nb <- lrtest(nb_glm, nb_glm_reduced)
print(lr_nb)

# Compare Logit vs Negative Binomial models
lr_glm_vs_nb <- lrtest(glm_logit, nb_glm)
print(lr_glm_vs_nb)

# ---- POST-HOC ANALYSIS ----

# Compare infection rates between populations using emmeans
population_comparisons <- emmeans(nb_glm, pairwise ~ Population_ID, type = "response")
print(population_comparisons)

# ---- VISUALIZATION ----

# Add predicted proportions and confidence intervals for Logit GLM
predictions <- predict(glm_logit, type = "link", se.fit = TRUE)

data <- data %>%
  mutate(
    Predicted_Proportion = predict(glm_logit, type = "response"),
    CI_Lower = plogis(predictions$fit - 1.96 * predictions$se.fit),  # Lower bound
    CI_Upper = plogis(predictions$fit + 1.96 * predictions$se.fit)   # Upper bound
  )

# ✅ Plot Dose-Response Curve (Logit GLM)
ggplot(data, aes(x = Dose, y = Predicted_Proportion, color = Population_ID, group = Population_ID)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Population_ID), alpha = 0.2, color = NA) +
  geom_point(aes(y = Infection_Status / Total_Exposed), size = 2, alpha = 0.7) +
  scale_x_log10() +  
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  scale_fill_manual(values = custom_colors) +   # Apply custom fills
  labs(
    title = "Dose-Response Curve (Logistic GLM)",
    x = "Dose",
    y = "Proportion Infected",
    color = "Population",
    fill = "Population"
  ) +
  theme_minimal()

# Calculate Proportion_Infected and its SEM
data <- data %>%
  mutate(Proportion_Infected = Infection_Status / Total_Exposed,
         SEM_Infected = sqrt((Proportion_Infected * (1 - Proportion_Infected)) / Total_Exposed))

# Generate predicted probabilities and confidence intervals
data <- data %>%
  mutate(
    Predicted_Prob = predict(nb_glm, type = "response"),
    Predicted_Prob_Capped = pmin(Predicted_Prob, Total_Exposed),
    CI_Lower = Predicted_Prob - 1.96 * predict(nb_glm, type = "response", se.fit = TRUE)$se.fit,
    CI_Upper = Predicted_Prob + 1.96 * predict(nb_glm, type = "response", se.fit = TRUE)$se.fit,
    CI_Upper_Capped = pmin(CI_Upper, Total_Exposed)
  )

# Reorder Population_ID for plotting
data$Population_ID <- factor(data$Population_ID, levels = c("IL-1", "IL-2", "IL-4", "IL-7", "IL-9", "IL-10", "IL-14", "IL-17", "IL-19", "IL-24", "IL-29", "IL-36", "Stock"))

# Combined observations and predictions plot
ggplot(data, aes(x = log10(Dose))) +
  geom_point(aes(y = Proportion_Infected, color = "Observed"), size = 2) +
  geom_line(aes(y = Predicted_Prob / Total_Exposed, color = "Predicted"), linetype = "dashed") +
  geom_ribbon(aes(ymin = CI_Lower / Total_Exposed, ymax = CI_Upper / Total_Exposed, fill = Population_ID), alpha = 0.2) +
  facet_wrap(~ Population_ID) +
  labs(title = "Observed vs Predicted Infection Rates by Population", x = "Log10(Dose)", y = "Proportion Infected") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

# ---- PLOT LOG-ODDS VS LOG10(DOSE) BY POPULATION ----

# Add predicted values to the dataset
data <- data %>%
  mutate(
    Predicted = predict(glm_logit, type = "response"),  # Predicted proportions# Observed proportions
    Residual = Observed - Predicted                    # Residuals
  )

# Summarize data by Population and Dose
results <- data %>%
  group_by(Population_ID, Dose) %>%
  summarize(
    Observed = mean(Observed, na.rm = TRUE),
    Predicted = mean(Predicted, na.rm = TRUE),
    Residual = mean(Residual, na.rm = TRUE),
    .groups = "drop"
  )

# Extract coefficients for slope and intercept
coefficients <- coef(glm_logit)
populations <- unique(data$Population_ID)

slope_intercept <- data.frame(
  Population = populations,
  Intercept = sapply(populations, function(pop) {
    intercept_term <- "(Intercept)"  # Global intercept
    pop_term <- paste0("Population_ID", pop)  # Specific population adjustment
    base_intercept <- ifelse(intercept_term %in% names(coefficients), coefficients[intercept_term], NA)
    pop_adjustment <- ifelse(pop_term %in% names(coefficients), coefficients[pop_term], 0)
    base_intercept + pop_adjustment  # Sum base intercept and population adjustment
  }),
  Slope = sapply(populations, function(pop) {
    dose_term <- "Log_Dose"  # Global dose coefficient
    interaction_term <- paste0("Log_Dose:Population_ID", pop)  # Interaction term for each population
    base_slope <- ifelse(dose_term %in% names(coefficients), coefficients[dose_term], NA)
    interaction_adjustment <- ifelse(interaction_term %in% names(coefficients), coefficients[interaction_term], 0)
    base_slope + interaction_adjustment  # Sum base slope and interaction adjustment
  })
)

# Generate dose sequence for plotting
dose_seq <- 10^seq(-6, -2, length.out = 100)

# Create data for plotting log-odds vs log10(dose)
plot_data <- slope_intercept %>%
  rowwise() %>%
  do({
    data.frame(
      Dose = dose_seq,
      LogOdds = .$Intercept + .$Slope * log10(dose_seq),
      Population = .$Population
    )
  }) %>%
  bind_rows()

# Plot log-odds vs log10(dose) for each population
ggplot(plot_data, aes(x = log10(Dose), y = LogOdds, color = Population, group = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Log-Odds vs Log10(Dose) by Population",
    x = "Log10(Dose)",
    y = "Log-Odds of Infection"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# ---- CLUSTER ANALYSIS WITH LOG-ODDS VS LOG-TRANSFORMED DOSES ----

# Summarize residual statistics for each population
# This creates a summary of model residuals for each population. Residuals represent the difference between 
# observed and predicted infection values. 

residual_summary <- data %>%
  group_by(Population_ID) %>%
  summarize(
    Mean_Residual = mean(Residual, na.rm = TRUE),    # Mean residual
    SD_Residual = sd(Residual, na.rm = TRUE),       # Standard deviation of residuals
    .groups = "drop"
  )

# Combine slope and intercept with residual statistics
# slope_intercept contains the slope and intercept from a model (likely a logistic or negative binomial GLM)
# for each population. This step merges slope and intercept estimates with the residual statistics, aligning 
# data for each Population_ID.
clustering_data <- slope_intercept %>%
  left_join(residual_summary, by = c("Population" = "Population_ID"))

# Standardize features for clustering
# scale() standardizes the selected columns (Slope, Intercept, Mean_Residual, and SD_Residual)
# Standardization ensures that all variables have mean = 0 and standard deviation = 1, preventing any variable 
# from dominating due to differences in scale.
scaled_features <- scale(clustering_data %>% select(Slope, Intercept, Mean_Residual, SD_Residual))

# Perform k-means clustering (e.g., 3 clusters based on initial insights)
# K-means clustering is applied to the standardized data. 3 clusters are specified and nstart = 25, which 
# ensures the algorithm runs 25 times with different initializations, selecting the best clustering based 
# on the lowest total within-cluster variance. set.seed(42) ensures results are reproducible.
set.seed(42)  # For reproducibility
kmeans_result <- kmeans(scaled_features, centers = 3, nstart = 25)

# Add cluster assignments to the data
# The cluster labels (1, 2, or 3) are assigned by kmeans() and stored in the clustering_data dataframe.
clustering_data$Cluster <- kmeans_result$cluster

# Clusters are merged back into the main dataset (data). Each row in data now has a Cluster column indicating
# which cluster the population belongs to based on infection response patterns.
data <- data %>%
  left_join(clustering_data %>% select(Population, Cluster), 
            by = c("Population_ID" = "Population"))

data <- data %>%
  mutate(Log_Odds = predict(glm_logit, type = "link"))  

#Plot log-odds vs. log-transformed dose for all populations (after clustering)
ggplot(data, aes(x = Log_Dose, y = Log_Odds, color = Population_ID, group = Population_ID)) +
  geom_line(size = 1) +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Log-Odds vs Log-Transformed Doses by Cluster",
    x = expression(Log[10](Dose)),
    y = "Log-Odds",
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "gray90", color = NA),  # Light gray background
    panel.grid.major = element_line(color = "white"),              # White major grid lines
    panel.grid.minor = element_line(color = "gray80")              # Subtle gray minor grid lines
  )

# Calculate proportion infected for correlation analysis
data <- data %>%
  mutate(Proportion_Infected = Infection_Status / Total_Exposed)

# Calculate correlation for each cluster
correlation_results <- data %>%
  group_by(Cluster) %>%
  summarize(Correlation = cor(Proportion_Infected, Growth_Rate, use = "complete.obs"))

print(correlation_results)

# Create a data frame for correlations
correlation_results <- data.frame(
  Cluster = c(2, 3, 1),  # Ensure this matches your cluster numbers
  Correlation = c(0.01054038, 0.06831828, 0.05020781)
)

# Ensure cluster column is consistent with the main data
correlation_results <- correlation_results %>%
  mutate(Cluster = as.integer(Cluster))

# Merge correlation values with the cluster data
data <- data %>%
  left_join(correlation_results, by = "Cluster")

# ---- CLUSTER ANALYSIS WITH PROPORTION INFECTED VS GROWTH RATE ----

# Plot proportion infected vs. growth rate by cluster
ggplot(data, aes(x = Growth_Rate, y = Proportion_Infected, color = Population_ID)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 1, color = "black") + # Set regression lines to black
  facet_wrap(~ Cluster, scales = "free") +
  geom_text(
    data = correlation_results,
    aes(
      x = 0.55,  # Adjusted to ensure labels are visible
      y = 0.9,   # Adjusted to position within each plot
      label = paste0("r = ", round(Correlation, 3))
    ),
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = custom_colors) +  # Apply IL-specific colors
  scale_x_continuous(limits = c(0.5, 0.7)) +    # Standardize Growth Rate axis
  labs(
    title = "Proportion Infected vs. Growth Rate by Cluster",
    x = "Growth Rate",
    y = "Proportion Infected",
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "gray90", color = NA),  # Light gray background
    panel.grid.major = element_line(color = "white"),              # White major grid lines
    panel.grid.minor = element_line(color = "gray80")              # Subtle gray minor grid lines
  )


# ---- CLUSTER ANALYSIS WITH LD50 ----
# Preprocess data for clustering
data <- data %>%
  mutate(
    Dose = ifelse(Dose == 0, 1e-6, Dose),  # Replace zero doses
    Log_Dose = log10(Dose + 1),            # Compute log-transformed dose
    Proportion_Infected = Infection_Status / Total_Exposed  # Calculate proportion infected
  )

# Perform correlation between LD50 and Growth Rate
correlation_all <- cor.test(data$LD50, data$Growth_Rate, method = "pearson")

# Print the results
print(correlation_all)

# Cluster populations based on LD50
ld50_clusters <- data %>%
  group_by(Population_ID) %>%
  summarize(Avg_LD50 = mean(LD50, na.rm = TRUE)) %>%
  ungroup()

# Perform k-means clustering (e.g., 3 clusters)
set.seed(42)
kmeans_result <- kmeans(ld50_clusters$Avg_LD50, centers = 3, nstart = 25)
ld50_clusters$Cluster <- as.factor(kmeans_result$cluster)

# Merge cluster information back to the main dataset
data <- data %>%
  left_join(ld50_clusters, by = "Population_ID")

data <- data %>%
  left_join(clustering_data %>% select(Population, Cluster), 
            by = c("Population_ID" = "Population"))

# Analyze correlation between Growth Rate and LD50 within each cluster
correlation_stats <- data %>%
  group_by(Cluster) %>%
  summarize(
    R = round(cor(Growth_Rate, LD50, method = "pearson"), 3),
    p_value = signif(cor.test(Growth_Rate, LD50)$p.value, 3)
  ) %>%
  mutate(
    label = paste0("bold('R = ", R, ", p = ", p_value, "')")  # Bold everything
  )

#----Graph LD50 vs. Growth by Cluster----

# Remove duplicate labels
unique_labels <- data %>%
  distinct(Population_ID, Growth_Rate, LD50, Cluster)

# Plot the data
ggplot(data, aes(x = Growth_Rate, y = LD50, color = as.factor(Cluster))) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 1) +  # Regression lines
  geom_text_repel(
    data = unique_labels,  # Use the unique labels dataset
    aes(label = Population_ID),
    size = 3, 
    max.overlaps = Inf,  # Show all labels dynamically
    box.padding = 0.3,   # Space around labels
    point.padding = 0.2, # Space between labels and points
    seed = 42            # Consistent positioning
  ) +
  annotate(
    "label",
    x = c(0.5, 0.5, 0.5),  # Approximate x-positions for each cluster
    y = c(0.05, 0.066, 0.082),  # Approximate y-positions for each cluster
    label = c(
      "R = -0.518, p = 0.00338",  # Correct for Cluster 1
      "R = 1, p = 9.43e-27",      # Correct for Cluster 2
      "R = -0.37, p = 0.108"      # Correct for Cluster 3
    ),
    fontface = "bold",  # Make text bold
    label.size = 0.5,   # Thickness of the box border
    fill = "white",     # Background color for the box
    color = c("#6A3D9A", "#1F78B4", "#33A02C"),  # Text color
    hjust = 0
  ) +
  scale_color_manual(values = c("#6A3D9A", "#1F78B4", "#33A02C"), name = "Cluster") +  # Custom colors
  labs(
    title = "LD50 vs. Growth Rate by Cluster",
    x = "Growth Rate",
    y = "LD50",
    color = "Cluster"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add borders
    legend.position = "right",  # Keep legend on the right
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12)
  )

# ---- PCA Data ----

## Add graph of PCA1 vs. LD50 here.

# Load the PCA dataset
pca_data <- read.csv("~/Desktop/UC Berkeley/My Data/Inbred Lines - Life History:Resistance Assays/Bioinformatics:Genomics/PCA_data.csv")

# Visualize PCA without clustering
ggplot(pca_data, aes(x = PC1, y = PC2, label = ID)) +
  geom_point(size = 3, alpha = 0.8, color = "black") +  # All points in blue
  geom_text_repel(size = 3, color = "black", box.padding = 0.2, seed = 42) +  # Add labels
  labs(
    title = "PCA of Populations",
    x = "Principal Component 1 (PC1) 21.72% of Variance",
    y = "Principal Component 2 (PC2) 13.22% of Variance"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",  # Remove legend since there's no clustering
    plot.title = element_text(hjust = 0.5)
  )
