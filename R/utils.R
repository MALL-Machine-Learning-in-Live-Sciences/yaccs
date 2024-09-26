# Fast Correlation Based Feature Selection
# -----------------------------------------
labelData = function(years, clinical){
  
  time = clinical$SurvTime
  event = clinical$vital_status
  
  limit = 365 * years
  
  for (i in 1:length(time)) {
    if (time[i] > limit) {
      time[i] = limit
      event[i] = 'Alive'
    }
  }
  
  clinical$SurvTime = time
  clinical$vital_status = event
  
  return(clinical)
}



fast.cor.FS = function(x, y, min_su){
  
  require(FCBF)
  
  dis = discretize_exprs(t(x))
  
  fcbf = fcbf(dis, y, verbose = T, minimum_su = min_su)
  res = x[,fcbf$index]
  
  return(res)
}

# Function to retrieve our signature from validation cohort
geneAnnot = function(clinical, platID, signature){
  
  # Retrieve platform data
  gpl = getGEO(platID)
  annot = gpl@dataTable@table
  annot = annot[, c('ID', 'ENTREZ_GENE_ID')]
  
  xx = annot[match(signature, annot$ENTREZ_GENE_ID), ]
  xx = xx[complete.cases(xx), ]
  
  notUnique = setdiff(signature, xx$ENTREZ_GENE_ID)
  l = list()
  for (i in 1:length(notUnique)) {
    l[[i]] = grep(notUnique[i], annot$ENTREZ_GENE_ID)
  }
  names(l) = notUnique
  
  ids = list()
  for (i in 1:length(l)) {
    if (length(l[[i]]) == 0) {
      print(paste0('This gene don´t have any sonda: ', names(l)[i]))
    } else{
      ids[[i]] = rep(names(l)[i], length(l[[i]]))
    }
  }
  
  yy = annot[unlist(l),]
  annotSig = rbind.data.frame(xx, yy)
  annotSig$genes = c(xx$ENTREZ_GENE_ID, unlist(ids))
  
  split = strsplit(annotSig$ENTREZ_GENE_ID, split = ' /// ')
  qos = list()
  for (i in 1:length(annotSig$genes)) {
    if(annotSig$genes[i] %in% split[[i]]){
      qos[[i]] = 'Yes'
    } else {
      qos[[i]] = 'No'
    }
  }
  
  annotSig = cbind.data.frame(annotSig ,qos = unlist(qos))
  annotSig = annotSig[which(annotSig$qos == 'Yes'), ]
  
  return(annotSig)
  
}


SelectSonda = function(platID, annotSig){
  
  if (platID == 'GPL570') {
    chip = 'hgu133plus2'
  } else if (platID == 'GPL96'){
    chip = 'hgu133a'
  }
  
  jetsetGenes = unique(annotSig$genes[which(duplicated(annotSig$genes) == T)])
  js = jscores(chip, eg = jetsetGenes)
  pre = annotSig[! annotSig$genes %in% jetsetGenes, ]
  
  if (length(js$nProbes) > 0) {
    pre = annotSig[! annotSig$genes %in% jetsetGenes, ]
    
    dd = unique(js$EntrezID[which(duplicated(js$EntrezID) == T)])
    
    if (length(dd) > 0) {
      uniqueSonda = js[match(setdiff(js$EntrezID, dd), js$EntrezID), ]
      uniqueSonda = annotSig[match(rownames(uniqueSonda), annotSig$ID), ]
      
      sondaJS = list()
      for (i in 1:length(dd)) {
        selSonda = js[which(js$EntrezID == dd[i]), ]
        selSonda = selSonda[which.max(selSonda$overall),]
        sondaJS[[i]] = rownames(selSonda)
      }
      ss = annotSig[match(unlist(sondaJS), annotSig$ID ), ]
      res = rbind.data.frame(pre, uniqueSonda, ss)
    }
  } else{
    res = pre
  }
  
  stopifnot(length(unique(rownames(res))) == dim(res)[1])
  
  return(res)
  
}

selectSonda = function(geoID, platID, chip, signature){
  
  # Packages
  require(biomaRt)
  require(jetset)
  require(GEOquery)
  
  # Get GEO
  gse = getGEO(geoID)[[1]]
  
  # ENSEMBL to affymetrix
  if (platID == 'GPL570'){
    martAffy = 'affy_hg_u133_plus_2'
  } else if (platID == 'GPL96'){
    martAffy = 'affy_hg_u133a'
  }
  
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  annot = getBM(filters= c('ensembl_gene_id'), 
                attributes= c('ensembl_gene_id', martAffy),
                values=signature,
                mart= mart)
  annot = annot[which(annot[, martAffy] != ''), ]
  
  # Retrieve counts
  counts = exprs(gse)
  
  # Checking!
  stopifnot(length(unique(annot$ensembl_gene_id)) == length(signature))
  
  # Selecting sonda
  toChoose = duplicated(annot$ensembl_gene_id)
  dups = unique(annot$ensembl_gene_id[which(toChoose == TRUE)])
  
  sonda = list()
  for (i in 1:length(dups)) {
    idx = grep(dups[i], annot$ensembl_gene_id)
    ps = annot[idx,][, martAffy]
    js = jscores(chip, probeset = ps)
    js = js[order(js$overall, decreasing = T),]
    if (rownames(js)[1] %in% rownames(counts) == TRUE) {
      sonda[[i]] = rownames(js)[1]
    } else {
      sonda[[i]] = rownames(js)[2]
    }
  }
  names(sonda) = dups
  
  noDup = setdiff(annot$ensembl_gene_id, names(sonda))
  annot2 = rbind.data.frame(annot[match(noDup, annot$ensembl_gene_id), ],
                            annot[match(unlist(sonda), annot[, martAffy]), ])
  annot2 = annot2[match(signature, annot2$ensembl_gene_id), ]
  
  # Checking!
  stopifnot(annot2$ensembl_gene_id == signature)
  ch1 = intersect(rownames(counts), annot2[, martAffy])
  stopifnot(length(ch1) == length(signature))
  
  # Subsetting in count data
  counts = as.data.frame(t(counts))
  counts = counts[, ch1]
  
  # Get clinical Data
  meta = phenoData(gse) 
  clinical = meta@data
  
  res = list(annotation = annot2, counts = counts, clinical = clinical)
  
  return(res)
}



predictTest = function(yearModel = c('5y', '3y', 'all'), algorithm = c('rf','glmnet', 'svm'), platform = c('gpl96', 'gpl570'), dirModels, outDir){
  
  require(mlr)
  
  dirData = 'd:/Users/jlinares/Documents/projects/SurvivalTCGA-COAD/data/'
  setwd(dirData)
  
  if(yearModel == '5y'){
    Modelfilename = 'FiveYearsModel.rds'
  } else if (yearModel == '3y'){
    Modelfilename = 'ThreeYearsModel.rds'
  } else if (yearModel == 'all'){
    Modelfilename = 'AllYearsModel.rds'
  }
  
  if(platform == 'gpl96'){
    plat = 1
  } else if (platform == 'gpl570'){
    plat = 2
  }
  
  cvrts = c('age_at_diagnosis', 'gender_male', 'Molecular_Subtype_noCIN', 'Stage_II', 'Stage_III', 'Stage_IV', 'MSI Status_MSI-L/MSS')
  
  # Read data according to model
  model_Y = readRDS(Modelfilename)
  
  # Create train task
  train = cbind.data.frame(model_Y[[plat]]$clinicalTrain[,cvrts],
                           model_Y[[plat]]$omicTrain_fcbf,
                           vital_status = model_Y[[plat]]$clinicalTrain[, 'vital_status'])
  names(train) = make.names(names(train))
  names(train)[ncol(train)] = 'target'
  train = makeClassifTask(data = train, target = 'target')
  train = normalizeFeatures(train) #return!
  
  # Create task for test data
  signature = names(model_Y[[plat]]$omicTrain_fcbf)
  test = cbind.data.frame(model_Y[[plat]]$clinicalTest[, cvrts],
                          model_Y[[plat]]$omicTest[, signature],
                          vital_status = model_Y[[plat]]$clinicalTest[, 'vital_status'])
  names(test) = make.names(names(test))
  names(test)[ncol(test)] = 'target'
  test = makeClassifTask(data = test, target = 'target')
  test = normalizeFeatures(test) #return!
  
  # Detect path to model
  setwd(dirModels)
  files = list.files(dirModels)
  files = files[grep(algorithm, files)]
  pathTOmodel = files[grep(yearModel, files)]
  pathTOmodel = pathTOmodel[grep(platform, pathTOmodel)]
  
  # Retrieve best model
  glmnet = readRDS(pathTOmodel)
  res = as.data.frame(glmnet)
  model = getBMRModels(glmnet)
  m = model[[1]][[1]]
  best = getLearnerModel(m[[which.max(res$auc)]]) #return!
  
  # Make prediction in train
  predTrain = predict(best, task = train) #return!
  probsTrain = getPredictionProbabilities(predTrain) #return!
  
  # Make prediction
  pred = predict(best, task = test) #return!
  confMatrix = calculateROCMeasures(pred) #return!
  perf = performance(pred, measures = list(acc, auc, mmce)) #return!
  probs = getPredictionProbabilities(pred) #return!
  
  result = list(data = list(train = train, test = test),
                bestModel = best,
                predictions = list(train = predTrain, test = pred),
                confMatrix = confMatrix,
                performances = perf,
                probabilities = list(ptrain = probsTrain, ptest = probs),
                info = list(modelReaded = pathTOmodel,
                            dataModelReaded = Modelfilename,
                            platform = platform))
  
  saveRDS(result, file = paste0(outDir, 'TestPredictions_', yearModel, platform, algorithm, '.rds'))
  return(result)
}


calculeIndx = function(test, train, ptest, ptrain){
  
  require(Hmisc)
  require(survminer)
  require(survAUC)
  
  c = rcorr.cens(ptest, Surv(test$SurvTime, test$vital_status))
  print(c)
  auc = AUC.uno(Surv(train$SurvTime, train$vital_status), Surv(test$SurvTime, test$vital_status), 1-ptest, seq(0, max(test$SurvTime), 2))
  brier = predErr(Surv(train$SurvTime, train$vital_status), Surv(test$SurvTime, test$vital_status), ptrain, ptest, seq(0, max(test$SurvTime), 2))
  
  # Plotting
  aucDF = data.frame(AUC = auc$auc, time = auc$times)
  brierDF = data.frame(Brier = brier$error, time = brier$times)
  
  gauc = ggplot(aucDF, aes(time, AUC)) + geom_line() +
    scale_color_manual(values = "#00AFBB") + theme(axis.title.x = element_blank())
  
  gbrier = ggplot(brierDF, aes(time, Brier)) + geom_line() +
    scale_color_manual(values = "#E7B800") + theme(axis.title.x = element_blank())
  
  print(ggarrange(gauc, gbrier, ncol = 2, nrow = 1))
  
  indx = list(c = c, auc = auc, brier = brier, plots = list(gauc = gauc, gbrier = gbrier))
  return(indx)
}


# == Data distribution
plotKM = function(filename = c('AllYearsModel.rds', 'ThreeYearsModel.rds', 'FiveYearsModel.rds'), platform = c('gpl96', 'gpl570')){
  
  require(survminer)
  require(survival)
  
  dirData = "d:/Users/jlinares/Documents/projects/SurvivalTCGA-COAD/data/"
  setwd(dirData)
  
  data = readRDS(filename)
  
  if (platform == 'gpl96') {
    plat = 1
  } else if (platform == 'gpl570'){
    plat = 2
  }
  
  x = data[[plat]]$clinicalTrain[, c('SurvTime', 'vital_status')]
  y = data[[plat]]$clinicalTest[, c('SurvTime', 'vital_status')]
  plotdf = rbind.data.frame(cbind(x, split = rep('train', nrow(x))),
                            cbind(y, split = rep('test', nrow(y))))
  plotdf$vital_status = ifelse(plotdf$vital_status == 'Dead', T, F)
  
  p = ggsurvplot(survfit(Surv(SurvTime, vital_status) ~ split, data = plotdf), data = plotdf, risk.table = T)
  p = p + ggtitle(paste0(gsub('.rds', '', filename), ' ', platform))
  
  print(p)
  
}


boxplotGenes = function(model = c('all', 'five', 'three'), platID = c('gpl96', 'gpl570')){
  
  require(data.table)
  require(ggplot2)
  require(ggpubr)
  require(viridis)
  require(biomaRt)
  
  dirData = 'd:/Users/jlinares/Documents/projects/SurvivalTCGA-COAD/data/'
  setwd(dirData)
  
  if (model == 'all') {
    modelPath = 'AllYearsModel.rds'
  } else if (model == 'three'){
    modelPath = 'ThreeYearsModel.rds'
  } else if (model == 'five'){
    modelPath = 'FiveYearsModel.rds'
  }
  
  if (platID == 'gpl96') {
    nplat = 1
  }else if (platID == 'gpl570'){
    nplat = 2
  }
  
  data = readRDS(file = modelPath)
  omic = data[[nplat]]$omicTrain_fcbf
  target = data[[nplat]]$clinicalTrain$vital_status
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  Sigconvert <- getBM(filters= c('ensembl_gene_id'),
                      attributes= c('ensembl_gene_id', 'hgnc_symbol'),
                      values=names(omic),
                      mart= mart)
  Sigconvert = Sigconvert[match(names(omic), Sigconvert$ensembl_gene_id), ]
  Sigconvert = Sigconvert$hgnc_symbol
  names(omic) = Sigconvert
  
  data = cbind.data.frame(omic, target)
  
  n = names(data)
  l = list()
  for (i in seq_along(n[1:length(n)-1])) {
    l[[i]] = cbind.data.frame(meta = rep(n[i], nrow(data)), level = data[, i], Vital_Status = data$target)
  }
  
  data = rbindlist(l)
  
  # Plot
  p = ggplot(data = data, aes(x = meta, y = level)) +
    geom_boxplot(aes(fill=Vital_Status))
  p = p + scale_fill_manual(values = viridis(2))
  p = p + facet_wrap( ~ meta, scales = 'free', strip.position = 'top', ncol = 5, nrow = 8) + stat_compare_means(aes(group = Vital_Status), method = 'wilcox.test',
                                                                                                                label = 'p.signif', label.y.npc = 0.80, hide.ns = T)
  p = p + theme_light()
  p = p + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  
  print(p)
}

Scaling = function(GeneSig, data){
  toScale = c('age_at_diagnosis', GeneSig)
  add = setdiff(names(data), toScale)
  
  sub = subset(data, select = toScale)
  res = as.data.frame(scale(sub, scale = F))
  
  res = cbind.data.frame(res, data[, add])
  
  return(res)
}

Scaling2 = function(GeneSig, data){
  toScale = GeneSig
  add = setdiff(names(data), toScale)
  
  sub = subset(data, select = toScale)
  res = as.data.frame(scale(sub, scale = F))
  
  res = cbind.data.frame(res, data[, add])
  
  return(res)
}

predIDX = function(cox, trData, valData){
  
  require(survival)
  require(survminer)
  require(survcomp)
  require(survAUC)
  
  ## C-Index
  pred = survConcordance(Surv(SurvTime, vital_status) ~ predict(cox, valData), valData)
  c = pred$concordance
  cUpper = c + 1.96 * pred$std.err
  cLower = c - 1.96 * pred$std.err
  
  ## Hazard Ratio
  lpnew = predict(cox, valData, type = 'lp')
  hr = hazard.ratio(x = lpnew, surv.time = valData$SurvTime, surv.event = valData$vital_status)
  hazarRatio = hr$hazard.ratio
  lower = hr$lower
  upper = hr$upper
  
  ## AUC
  lp = predict(cox, type = 'lp')
  Surv.rsp = Surv(trData$SurvTime, trData$vital_status)
  Surv.rsp.new = Surv(valData$SurvTime, valData$vital_status)
  times = seq(10, max(valData$SurvTime), 2)
  aucRes = AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
  auc = aucRes$iauc
  
  res = list(cIndex = list(C = c, CLower = cLower, CUpper = cUpper), HR = list(ratio = hazarRatio, lower = lower, upper = upper), AUC = auc)
  
  return(res)
}



grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}



