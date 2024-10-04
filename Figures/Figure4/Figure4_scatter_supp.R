require(ggpubr)
s1 = ggscatter(gsva, x = 'ours', y = 'G0.EarlyG1', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s2 = ggscatter(gsva, x = 'ours', y = 'G1Phase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s3 = ggscatter(gsva, x = 'ours', y = 'G1S.transition', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s4 = ggscatter(gsva, x = 'ours', y = 'S.phase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s5 = ggscatter(gsva, x = 'ours', y = 'G2.phase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s6 = ggscatter(gsva, x = 'ours', y = 'G2M.transition', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s7 = ggscatter(gsva, x = 'ours', y = 'M.prophase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s8 = ggscatter(gsva, x = 'ours', y = 'M.prometaphase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s9 = ggscatter(gsva, x = 'ours', y = 'M.metaphase.anaphase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s10 = ggscatter(gsva, x = 'ours', y = 'M.telophase.cytokinesis', add = 'reg.line', conf.int = T,
                add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')

ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, nrow = 5, ncol = 2)