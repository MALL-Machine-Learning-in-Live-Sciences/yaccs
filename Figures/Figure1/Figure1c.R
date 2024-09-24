# Scatter plot of C-Index in TCGA-COAD cohort
source("requirements.R")

# Inputpath
inputpath <- "data/cox/tcga/cindex_tcga.rds"

# Load data
data <- readRDS(inputpath)

# Plotting
xname <- expression(italic("C-Index"))
ggplot(data=data, aes(y=Study, x=C, xmin=ci_l, xmax=ci_u, fill = Subset, color = Subset)) + 
  
  geom_point(size = 5, position = position_dodge(width = .9)) +
  geom_errorbar(position = position_dodge2(width = .5, padding = 0.7)) +
  
  scale_color_manual(values = viridis(3)) +
  
  scale_x_continuous(limits=c(0,1), name=xname) +
  geom_vline(xintercept=0.5, color="black", linetype="dashed", alpha=.5) +
  
  theme_light() +
  theme(text=element_text(family="Times",size=15, color="black")) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(axis.title.y = element_blank())
