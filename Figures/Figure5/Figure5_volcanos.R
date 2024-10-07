# Plot volcano plots
source("requirements.R")

# Inputpath
fgsea_inputpath <- "data/repurposing/fgsea_by_signature.rds"

# Load data
fgseaRes <- readRDS(fgsea_inputpath)

pp <- list()
for (j in seq_along(fgseaRes)) {

    pp[[j]] <- ggscatter(
        fgseaRes[[j]],
        x = "NES",
        y = "log10adjpval",
        size = "size",
        label = "pathway",
        label.select = fgseaRes[[j]]$pathway[which(fgseaRes[[j]]$padj < 0.05)],
        repel = TRUE,
        title = names(fgseaRes)[j]
        ) +
    ylim(c(0, 6))
}

plot <- ggpubr::ggarrange(
    pp[[1]],
    pp[[2]],
    pp[[3]],
    pp[[4]],
    pp[[5]],
    pp[[6]],
    ncol = 2,
    nrow = 3,
    common.legend = TRUE
    )
annotate_figure(
    plot,
    top = text_grob("PRISM2",
    color = "black",
    size = 20))
