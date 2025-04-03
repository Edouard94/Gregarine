# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)
library(tidyr)
library(car)
library(extrafont)
loadfonts()
setwd("C:/Users/eb826/OneDrive - University of Exeter/PhD_EdouardMicrosporidia/Thesis/Chapter4_Leidyana gryllorum/")

# Read the data from the Excel file
data <- read_excel("Rtable_Leidyana.xlsx", sheet = 1)

# Convert Species column to factor
data$Species <- as.factor(data$Species)

# Descriptive statistics
descriptive_stats <- data %>%
  group_by(Species) %>%
  summarise(across(starts_with(c("Total_length", "Protomerite_L", "Deutomerite_L", "Protomerite_WE", "Deutomerite_WE", "Protomerite_Ratio", "Deutomerite_Ratio")), list(mean = mean, se = ~sd(.x)/sqrt(length(.x))), .names = "{.col}_{.fn}"), na.rm = TRUE)

print(descriptive_stats)

# Reshape data to long format
data_long <- data %>%
  pivot_longer(cols = c(Total_length, Deutomerite_WE), names_to = "Trait", values_to = "Value")

# Plotting
Traits_boxplot <- ggplot(data_long, aes(x = Species, y = Value, color = Trait)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text = element_text(size = 12), 
    axis.text.x = element_text(face = "italic"), # Make x-axis labels italic
    axis.title = element_text(size = 12, vjust = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.direction = "vertical",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  labs(title = "Traits by Species", y = "Trait Values (µm)", color = "Trait") +
  scale_x_discrete(labels = c("Lerratica" = "L. erratica", "Lgryllorum" = "L. gryllorum")) +
  scale_color_discrete(labels = c("Total_length" = "Gamont total length", "Deutomerite_WE" = "Deutomerite width")) +
  theme(legend.position = "top", legend.direction = "horizontal")
Traits_boxplot
ggsave("Traits_boxplot.png", dpi = 600, width = 4, height = 4.5, units = "in")


# Perform Mann-Whitney U tests
mann_whitney_results <- data %>%
  gather(key = "Measurement", value = "Value", -Species) %>%
  group_by(Measurement) %>%
  do(tidy(wilcox.test(Value ~ Species, data = .))) %>%
  mutate(test = "Mann-Whitney")

# Print results
all_results <- mann_whitney_results
all_results %>%
  dplyr::ungroup() %>%
  dplyr::select(Measurement, statistic, p.value, test) %>%
  dplyr::arrange(Measurement, test) %>%
  print()

# Extract test statistic and p-value for annotations
annotation_data <- mann_whitney_results %>%
  dplyr::select(Measurement, statistic, p.value) %>%
  mutate(label = paste0("W = ", round(statistic, 2), ", p = ", signif(p.value, 3)))

# Ensure annotation_data includes the Trait column to match the plot
annotation_data <- annotation_data %>%
  mutate(Trait = case_when(
    Measurement == "Total_length" ~ "Total_length",
    Measurement == "Deutomerite_WE" ~ "Deutomerite_WE",
    TRUE ~ NA_character_ # Only include traits present in the plot
  )) %>%
  filter(!is.na(Trait)) # Remove rows for traits not in the plot

# Add annotations to the plot
Traits_boxplot <- ggplot(data_long, aes(x = Species, y = Value, color = Trait)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text = element_text(size = 12), 
    axis.text.x = element_text(face = "italic"), # Make x-axis labels italic
    axis.title = element_text(size = 12, vjust = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.direction = "vertical",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  labs(title = "Traits by Species", y = "Trait Values (µm)", color = "Trait") +
  scale_x_discrete(labels = c("Lerratica" = "L. erratica", "Lgryllorum" = "L. gryllorum")) +
  scale_color_discrete(labels = c("Total_length" = "Gamont total length", "Deutomerite_WE" = "Deutomerite width")) +
  theme(legend.position = "top", legend.direction = "horizontal") +
  geom_text(data = annotation_data, aes(x = 1.5, y = max(data_long$Value) + 5, label = label, color = Trait), inherit.aes = FALSE, size = 4, hjust = 0.5)

Traits_boxplot

# Save the updated plot
ggsave("Traits_boxplot_with_stats.svg", dpi = 600, width = 4, height = 4.5, units = "in")

# Assumption checks
# Normality
#normality_tests <- data %>%
#gather(key = "Measurement", value = "Value", -Species) %>%
#group_by(Measurement, Species) %>%
#do(tidy(shapiro.test(.$Value)))
#print(normality_tests)

# Equality of variances
#variance_tests <- data %>%
#gather(key = "Measurement", value = "Value", -Species) %>%
#group_by(Measurement) %>%
#do(tidy(car::leveneTest(Value ~ Species, data = .)))
#print(variance_tests)

# Calculate pooled standard deviations for each measurement
#std_devs <- data %>%
#gather(key = "Measurement", value = "Value", -Species) %>%
#group_by(Measurement) %>%
#summarise(sd = sqrt(((n() - 1) * var(Value)) / n()), .groups = "drop")

# Perform t-tests, join standard deviations, calculate Cohen's d, and round p-values
#t_test_results <- data %>%
#gather(key = "Measurement", value = "Value", -Species) %>%
#group_by(Measurement) %>%
#do(tidy(t.test(Value ~ Species, data = .))) %>%
#left_join(std_devs, by = "Measurement") %>%
#mutate(cohens_d = estimate / sd,
#p.value = round(p.value, 3),
#test = "t-test")

# Merge results
#all_results <- bind_rows(mann_whitney_results, t_test_results)

# Print results
#all_results %>%
#dplyr::ungroup() %>%
#dplyr::select(Measurement, estimate, statistic, p.value, conf.low, conf.high, cohens_d, test) %>%
#dplyr::arrange(Measurement, test) %>%
#print()

# Load the package
library(writexl)
# Write the data frame to an Excel file
write_xlsx(all_results, "all_results.xlsx")