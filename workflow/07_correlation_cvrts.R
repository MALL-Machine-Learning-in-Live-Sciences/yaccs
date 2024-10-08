# Correlation with cvrts

data <- readRDS("extdata/tcga-coad/TCGA-COAD_clinvars_sigscores.rds")

sig <- as.character(unique(data$signature))
vars <- c('Molecular_Subtype', 'Stage', 'MSI.Status', 'clusters', 'kras', 'braf', 'tp53', 'apc', 'smad4')
res <- list()

for (i in seq_along(sig)) {
  x <- data[which(data$signature == sig[i]), ]
  res[[i]] <- data.frame(
    signatures = sig[i],
    Molecular_subtype = t.test(x[which(x[,vars[1]] == 'CIN'),]$value, x[which(x[,vars[1]] == 'noCIN'),]$value)$p.value,
    Stage = anova(lm(x$value ~ as.numeric(as.factor(x$Stage))))$`Pr(>F)`[1],
    MSI = t.test(x[which(x[,vars[3]] == 'MSI-H'),]$value, x[which(x[,vars[3]] == 'MSI-L/MSS'),]$value)$p.value,
    Consensus_clusters = kruskal.test(x$value, x$clusters)$p.value,
    kras = t.test(x[which(x[,vars[5]] == F),]$value, x[which(x[,vars[5]] == T),]$value)$p.value,
    braf = t.test(x[which(x[,vars[6]] == F),]$value, x[which(x[,vars[6]] == T),]$value)$p.value,
    tp53 = t.test(x[which(x[,vars[7]] == F),]$value, x[which(x[,vars[7]] == T),]$value)$p.value,
    apc = t.test(x[which(x[,vars[8]] == F),]$value, x[which(x[,vars[8]] == T),]$value)$p.value,
    smad4 = t.test(x[which(x[,vars[9]] == F),]$value, x[which(x[,vars[9]] == T),]$value)$p.value
  )
}

res <- as.data.frame(data.table::rbindlist(res))
rownames(res) <- res$signatures
res <- res[, -1]
res <- as.data.frame(t(res))
xtable::xtable(res, digits = -2, auto = T,
               caption = 'p-values')
