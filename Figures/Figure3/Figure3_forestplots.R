# Forest plot in external validation cohorts
# ---------------------
require(forestplot)
require(dplyr)
require(viridis)

# Inputpath
inputpath <- "data/cox/external/cindex_external_validations.rds"

# Load data
data <- readRDS(inputpath)

genes_cvrts_pos <- grep(" + cvrts", data$Study, fixed = TRUE)
cvrts_pos <- grep("^cvrts$", data$Study)

genes_cvrts <- data[genes_cvrts_pos, ]
cvrts <- 
    data[cvrts_pos, ] %>% 
    group_by(Subset) %>% 
    dplyr::slice(1) %>%
    mutate(
        signature = "cvrts"
    )
genes_cvrts <- rbind.data.frame(genes_cvrts, cvrts)

sig <- data[-c(genes_cvrts_pos, cvrts_pos), ]

# Plot genes + cvrts
View(genes_cvrts)
View(cvrts)

# Plot genes + cvrts
genes_cvrts |>
  as.data.frame() |>
  mutate(
    labeltext = Subset,
    mean = C,
    lower = ci_l,
    upper = ci_u
  ) |>
  group_by(signature) |>
  forestplot( 
    xlab = "C-Index",
    boxsize = 0.3,  # size of the boxes
    col = fpColors(box="darkgreen", line="darkblue", summary="royalblue")
    ) |>
  fp_add_lines("steelblue") |> 
  fp_add_header("Cohort") |> 
  fp_set_style(box = viridis(7))


# Plot genes
sig |>
  as.data.frame() |>
  mutate(
    labeltext = Subset,
    mean = C,
    lower = ci_l,
    upper = ci_u
  ) |>
  group_by(signature) |>
  forestplot( 
    xlab = "C-Index",
    boxsize = 0.3,  # size of the boxes
    col = fpColors(box="darkgreen", line="darkblue", summary="royalblue")
    ) |>
  fp_add_lines("steelblue") |> 
  fp_add_header("Cohort") |> 
  fp_set_style(box = viridis(6))
