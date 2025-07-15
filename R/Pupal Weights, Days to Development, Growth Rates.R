library(ggplot2)
library(dplyr)
library(readr)
library(multcompView)

# Load the data
data <- read.csv("~/Desktop/UC Berkeley/My Data/Inbred Lines - Life History:Resistance Assays/Signe's IL Data/Pupal Weights, Days to Development, Growth Rates.csv")

data$Line <- gsub("-", "", data$Line)

# Perform ANOVA and Tukey's HSD for Pupal Weight
anova_pupal <- aov(Weight ~ Line, data = data)
tukey_pupal <- TukeyHSD(anova_pupal)

# Extract Tukey's HSD results and create labels
tukey_labels_pupal <- multcompLetters4(anova_pupal, tukey_pupal)
labels_pupal <- as.data.frame.list(tukey_labels_pupal$`Line`)
labels_pupal$Line <- rownames(labels_pupal)

# Merge labels with the original data for plotting
max_pupal_weight <- max(data$Weight) + 1
labels_pupal <- labels_pupal %>%
  mutate(y = max_pupal_weight)

# Create the box plot
custom_colors <- c(
  "Stock" = "#FFFFFF", "IL1" = "#D73027", "IL2" = "#FC8D59", "IL4" = "#EFC000", 
  "IL7" = "#91BFDB", "IL9" = "#4575B4", "IL10" = "#313695", "IL14" = "#8E44AD", 
  "IL17" = "#D73092", "IL19" = "#666666", "IL24" = "#1D91C0", "IL29" = "#8C510A", 
  "IL36" = "#A6DBA0"
)

# Create the box plot for Pupal Weights with customized colors
plot1<-ggplot(data, aes(x = Line, y = Weight, fill = Line)) +
  geom_boxplot(color = "black", alpha = 1) +  # Set box outline to black and make box background opaque
  geom_jitter(color = "black", width = 0.2, size = 1, alpha = 0.7) +  # Jitter points are black
  geom_text(data = labels_pupal, aes(x = Line, y = y, label = Letters), vjust = -0.5) +
  scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 40, by = 5)) +  # Set y-axis from 0 to 40 with 5 mg increments
  labs(
    title = "Pupal Weight by Population", 
    y = "Pupal Weight (mg)", 
    x = "Population", 
    fill = "Population"  # Change legend title
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Show the legend on the right
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bolden the title
    axis.title.x = element_text(face = "bold", size = 12),  # Bolden the x-axis title
    axis.title.y = element_text(face = "bold", size = 12)   # Bolden the y-axis title
  )


# Perform ANOVA and Tukey's HSD for Days to Development
anova_days <- aov(Pupation ~ Line, data = data)
tukey_days <- TukeyHSD(anova_days)

# Extract Tukey's HSD results and create labels
tukey_labels_days <- multcompLetters4(anova_days, tukey_days)
labels_days <- as.data.frame.list(tukey_labels_days$`Line`)
labels_days$Line <- rownames(labels_days)

# Merge labels with the original data for plotting
min_days_to_development <- min(data$Pupation) - 1  # Adjust the y-position to be below the x-axis
labels_days <- labels_days %>%
  mutate(y = min_days_to_development)

# Define custom colors for each population, with Stock as white
custom_colors <- c(
  "Stock" = "#FFFFFF", "IL1" = "#D73027", "IL2" = "#FC8D59", "IL4" = "#EFC000", 
  "IL7" = "#91BFDB", "IL9" = "#4575B4", "IL10" = "#313695", "IL14" = "#8E44AD", 
  "IL17" = "#D73092", "IL19" = "#666666", "IL24" = "#1D91C0", "IL29" = "#8C510A", 
  "IL36" = "#A6DBA0"
)

# Create the box plot with customized colors, legend title, y-axis scale, and Tukey labels at the bottom
plot2<-ggplot(data, aes(x = Line, y = Pupation, fill = Line)) +
  geom_boxplot(color = "black", alpha = 1) +  # Set box outline to black and make box background opaque
  geom_jitter(color = "black", width = 0.2, size = 1, alpha = 0.7) +  # Jitter points are black
  geom_text(data = labels_days, aes(x = Line, y = y, label = Letters), vjust = 4) +  # Place letters below the x-axis
  scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
  scale_y_continuous(limits = c(10, 45), breaks = seq(10, 45, by = 5)) +  # Set y-axis from 0 to 40 with 5 mg increments
  labs(
    title = "Days to Development by Population", 
    y = "Days to Development (days)", 
    x = "Population", 
    fill = "Population"  # Change legend title
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # No legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bolden the title
    axis.title.x = element_text(face = "bold", size = 12),  # Bolden the x-axis title
    axis.title.y = element_text(face = "bold", size = 12)   # Bolden the y-axis title
  )

#----
# Perform ANOVA and Tukey's HSD for Growth Rate
anova_growth <- aov(Growth_Rate ~ Line, data = data)
tukey_growth <- TukeyHSD(anova_growth)

# Extract Tukey's HSD results and create labels
tukey_labels_growth <- multcompLetters4(anova_growth, tukey_growth)
labels_growth <- as.data.frame.list(tukey_labels_growth$`Line`)
labels_growth$Line <- rownames(labels_growth)

# Merge labels with the original data for plotting
min_growth <- min(data$Growth_Rate) - 1  # Adjust the y-position to be below the x-axis
labels_growth <- labels_growth %>%
  mutate(y=max(data$Growth_Rate) + 0.1)

# Define custom colors for each population, with Stock as white
custom_colors <- c(
  "Stock" = "#FFFFFF", "IL1" = "#D73027", "IL2" = "#FC8D59", "IL4" = "#EFC000", 
  "IL7" = "#91BFDB", "IL9" = "#4575B4", "IL10" = "#313695", "IL14" = "#8E44AD", 
  "IL17" = "#D73092", "IL19" = "#666666", "IL24" = "#1D91C0", "IL29" = "#8C510A", 
  "IL36" = "#A6DBA0"
)

# Create the box plot with customized colors, legend title, y-axis scale, and Tukey labels at the bottom
plot3<-ggplot(data, aes(x = Line, y = Growth_Rate, fill = Line)) +
  geom_boxplot(color = "black", alpha = 1) +  # Set box outline to black and make box background opaque
  geom_jitter(color = "black", width = 0.2, size = 1, alpha = 0.7) +  # Jitter points are black
  geom_text(data = labels_growth, aes(x = Line, y = y, label = Letters), vjust = 0) +  # Place letters below the x-axis
  scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1.3, by = 0.1)) +  # Set y-axis 
  labs(
    title = "Growth Rate by Population", 
    y = "Growth Rate", 
    x = "Population", 
    fill = "Population"  # Change legend title
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # No legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Center and bolden the title
    axis.title.x = element_text(face = "bold", size = 12),  # Bolden the x-axis title
    axis.title.y = element_text(face = "bold", size = 12),   # Bolden the y-axis title
  )

library(gridExtra)

# Combine the three plots into one
grid.arrange(plot1, plot2, plot3, ncol = 1)

# Combine relevant columns into a data frame (already present in 'data')
combined_data <- data.frame(
  Pupal_Weight = data$Weight,
  Days_to_Development = data$Pupation,
  Population_ID = data$Line
)

# Perform correlation test to get both r and p-value
cor_test <- cor.test(combined_data$Days_to_Development, combined_data$Pupal_Weight, use = "complete.obs")

# Extract r and p-value
r_value <- cor_test$estimate
p_value <- cor_test$p.value

# Print results
print(paste("Correlation (r):", round(r_value, 3)))
print(paste("p-value:", format.pval(p_value, digits = 3, scientific = TRUE)))

# Scatter plot with correlation, regression line, and CI shading
ggplot(combined_data, aes(x = Days_to_Development, y = Pupal_Weight, color = Population_ID)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot of the data points
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", fill = "gray70") +  # Add linear regression line with CI shading
  labs(
    title = "Correlation Between Pupal Weight and Days to Development",
    x = "Days to Development",
    y = "Pupal Weight",
    color = "Population"
  ) +
  annotate("text", x = max(combined_data$Days_to_Development) - 5, 
           y = max(combined_data$Pupal_Weight), 
           label = paste0("r = ", round(r_value, 3), "\n", "p = ", format.pval(p_value, digits = 3, scientific = TRUE)),
           hjust = 1, size = 5) +  # Display correlation coefficient
  scale_color_manual(values = c("Stock" = "#000000", "IL1" = "#D73027", "IL2" = "#FC8D59", "IL4" = "#EFC000", "IL7" = "#91BFDB", 
                                "IL9" = "#4575B4", "IL10" = "#313695", "IL14" = "#8E44AD", "IL17" = "#D73092", "IL19" = "#666666", 
                                "IL24" = "#1D91C0", "IL29" = "#8C510A", "IL36" = "#A6DBA0")) +  # Custom colors for populations
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold the title
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),  # Darker and larger x-axis label
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),  # Darker and larger y-axis label
    axis.text.x = element_text(color = "black", size = 10),  # Make x-axis tick labels darker
    axis.text.y = element_text(color = "black", size = 10)   # Make y-axis tick labels darker
  )

# Calculate means and standard errors for each population
summary_data <- data %>%
  group_by(Line) %>%
  summarise(
    Mean_Pupal_Weight = mean(Weight, na.rm = TRUE),
    SE_Pupal_Weight = sd(Weight, na.rm = TRUE) / sqrt(n()),
    Mean_Days_to_Development = mean(Pupation, na.rm = TRUE),
    SE_Days_to_Development = sd(Pupation, na.rm = TRUE) / sqrt(n())
  )

# Ensure Population_ID is a factor with the correct levels
summary_data$Line <- factor(summary_data$Line, levels = c("Stock", "IL1", "IL2", "IL4", "IL7", "IL9", "IL10", "IL14", "IL17", "IL19", "IL24", "IL29", "IL36"))

# Calculate correlation coefficient for mean values
mean_correlation <- cor(summary_data$Mean_Days_to_Development, summary_data$Mean_Pupal_Weight, use = "complete.obs")

# Scatter plot of mean values with error bars, regression line, CI shading, and correlation coefficient
ggplot(summary_data, aes(x = Mean_Days_to_Development, y = Mean_Pupal_Weight, fill = Line)) +
  geom_point(shape = 21, size = 4, color = "black") +  # Add points with custom fill and black outline for all
  geom_errorbar(aes(ymin = Mean_Pupal_Weight - SE_Pupal_Weight, ymax = Mean_Pupal_Weight + SE_Pupal_Weight), width = 0.1, color = "black") +  # Black error bars for Pupal Weight
  geom_errorbarh(aes(xmin = Mean_Days_to_Development - SE_Days_to_Development, xmax = Mean_Days_to_Development + SE_Days_to_Development), height = 0.1, color = "black") +  # Black error bars for Days to Development
  geom_smooth(aes(x = Mean_Days_to_Development, y = Mean_Pupal_Weight), method = "lm", se = TRUE, color = "black", linetype = "solid", fill = "gray70") +  # Regression line with CI shading
  labs(
    title = "Correlation Between Mean Pupal Weight and Mean Days to Development",
    x = "Mean Days to Development",
    y = "Mean Pupal Weight",
    fill = "Population"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(mean_correlation, 2)), hjust = 7.5, vjust = 5, size = 5) +  # Display correlation coefficient
  scale_fill_manual(values = c("Stock" = "#FFFFFF", "IL1" = "#D73027", "IL2" = "#FC8D59", "IL4" = "#EFC000", "IL7" = "#91BFDB", 
                               "IL9" = "#4575B4", "IL10" = "#313695", "IL14" = "#8E44AD", "IL17" = "#D73092", "IL19" = "#666666", 
                               "IL24" = "#1D91C0", "IL29" = "#8C510A", "IL36" = "#A6DBA0")) +  # Custom fill colors for populations
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold the title
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),  # Darker and larger x-axis label
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),  # Darker and larger y-axis label
    axis.text.x = element_text(color = "black", size = 10),  # Make x-axis tick labels darker
    axis.text.y = element_text(color = "black", size = 10)   # Make y-axis tick labels darker
  )
