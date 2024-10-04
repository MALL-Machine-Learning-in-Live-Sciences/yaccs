# Code to reproduce Figure 1b)
# Forest plot of Hazard ratio by feature: univariate model
# --------------------------
source("requirements.R")

# Inputpath
cox_inputpath <- "data/cox/tcga-train/cox_cvrts_yaccs.rds"
train_inputpath <- "data/cox/tcga-train/train_all.rds"

# Arguments
cpositions <- c(0.02, 0.22, 0.4)
fontsize <- 0.7

# Load data
train <- readRDS(train_inputpath)
names(train) <- make.names(names(train))

# Create an empty list to store Cox model results
cox_results_list <- list()

# List of predictor variables you want to analyze
predictor_vars <- names(train)
predictor_vars <- setdiff(predictor_vars, c("SurvTime", "vital_status"))

# Loop through the predictor variables and fit Cox models
for (var in predictor_vars) {
  # Create the formula for the Cox model
  formula <- as.formula(paste("Surv(SurvTime, vital_status) ~", var))
  
  # Fit the Cox model
  cox_model <- coxph(formula, data = train)
  
  # Extract HR and CI
  hr <- exp(coef(cox_model))
  ci <- exp(confint(cox_model))
  p_value <- summary(cox_model)$coef[5]
  
  res <- data.frame(
    var = var,
    hr = hr,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = p_value
  )
  
  # Store the Cox model results in the list
  cox_results_list[[var]] <- res
}

# Bind in a data.frame all the HR and CI of each feature
result <- rbindlist(cox_results_list)

# Format result
thresholds <- c(0.05, 0.01, 0.001)
result <- 
  result %>% 
  mutate(
    stars = dplyr::case_when(
      p_value < thresholds[3] ~ "***",
      p_value < thresholds[2] ~ "**",
      p_value < thresholds[1] ~ "*",
      TRUE ~ ""),
    p_value = round(p_value, 3),
    stars = paste(p_value, stars, sep = " "),
    stars = ifelse(stars == "0 ***", "<0.001 ***", stars),
    id = 1:nrow(result)
  ) %>% 
  dplyr::arrange(desc(id))
  

rangeb <- range(result$ci_lower, result$ci_upper)
breaks <- axisTicks(rangeb, log = F, nint = 8)
rangeplot <- rangeb
rangeplot[1] <- rangeplot[1]
rangeplot[2] <- rangeplot[2] + 0.5
main = "Univariate HR for each feature"
x_annotate <- seq_len(nrow(result))
width <- diff(rangeplot)
y_variable <- rangeplot[1] + cpositions[1] * width
annot_size_mm <- fontsize * as.numeric(convertX(unit(theme_get()$text$size,                                                     "pt"), "mm"))
y_nlevel <- rangeplot[1] + cpositions[2] * width
y_cistring <- rangeplot[1] + cpositions[3] * width
y_stars <- rangeb[2]

#setting up the basic plot
ggplot(
  result, 
  aes(seq_along(var), hr)
) + 
  geom_rect(
    aes(xmin = seq_along(var) - 0.5, 
        xmax = seq_along(var) + 0.5, 
        ymin = rangeplot[1], 
        ymax = rangeplot[2], 
        fill = ordered(seq_along(var)%%2 + 1))
  ) + 
  scale_fill_manual(
    values = c("#FFFFFF33", "#00000033"), 
    guide = "none"
  ) + 
  geom_point(
    pch = 15, 
    size = 4
  ) + 
  geom_errorbar(
    aes(ymin = ci_lower, 
        ymax = ci_upper), 
    width = .7
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_flip() +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    expand = c(0.4, 0.4),
    breaks = 1:7
  ) +
  theme_light() + 
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + 
  xlab("") + 
  annotate(
    geom = "text",
    x = x_annotate, 
    y = 0.35,
    label = result$var, 
    fontface = "bold",
    hjust = 1,
    size = annot_size_mm
  ) + 
  annotate(
    geom = "text",
    x = x_annotate, 
    y = 7.2,
    label = result$stars, 
    size = annot_size_mm, 
    hjust = 0, 
    fontface = "italic"
  )
