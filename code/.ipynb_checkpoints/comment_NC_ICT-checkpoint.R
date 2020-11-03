rm(list=ls())

library(rocc)
library(ggplot2)
library(Biobase)
library(edgeR)
library(limma)
library(biomaRt)
library(dplyr)
library(cancerclass)

######################## 
#calculate roc and draw a roc curve from Xiong et al. (https://github.com/donghaixiong/Immune_cells_analysis)
rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}


# single ROC plot
rocplot.single.V2 <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC = ",signif(auc, 2), " (95% CI: ", signif(ci.lower, 2), "-", signif(ci.upper, 2), ")", sep=""))
  }
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = "")) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_manual(labels = annotation, values = "#000000") +
    theme(
      plot.title = element_text(face="bold", size=14), 
      axis.text.x = element_text(face="bold", size=18),
      axis.text.y = element_text(face="bold", size=18),
      axis.title.x = element_text(face="bold", size=18),
      axis.title.y = element_text(face="bold", size=18, angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.justification=c(1,0), 
      legend.position=c(0.95,0.15),
      legend.text = element_text(size = 14),
      legend.title=element_blank(),
      legend.key = element_blank()
    )+
    labs(title=title)
  return(p)
}
               
               
############## 1. same experiment with a random gene set:
##### load ImmuneCell.sig and 4 dataset from Xiong et al. (https://github.com/donghaixiong/Immune_cells_analysis)
ImSig.genes <- readRDS('../data/ImSig.rds')

# GSE78220 28 samples
GSE78220_AltAnalyze <- readRDS('../data/GSE78220_expressionMatrix.rds')
GSE78220_PhenoInfo2 <- readRDS('../data/GSE78220_PhenoInfo2.rds')
GSE78220_AltAnalyze <- GSE78220_AltAnalyze[,c('Symbol',GSE78220_PhenoInfo2$sample)]
GSE78220_AltAnalyze <- GSE78220_AltAnalyze[rowSums(GSE78220_AltAnalyze[,-1]) > 1, ] # same operation with the original paper

# GSE91061 = BMS038 51 samples
BMS038.Pre.CountTable.normalized.log <- as.data.frame(readRDS('../data/BMS038.Pre.CountTable.normalized.log.rds'))
BMS038_phenoData <- readRDS('../data/BMS038_phenoData.rds')
BMS038.Pre.CountTable.normalized.log <- BMS038.Pre.CountTable.normalized.log[rowSums(BMS038.Pre.CountTable.normalized.log) > 10,] # same operation with the original paper, otherwise function fit() does not work 
BMS038.Pre.CountTable.normalized.log$Symbol <- rownames(BMS038.Pre.CountTable.normalized.log)
BMS038_PhenoInfo <- BMS038_phenoData@data

# CC_73samples from PRJEB23709 73 samples
CC_73samples_GE <- read.table('../data/DATASET-PRJEB23709_Pre_73samples.txt',sep="\t",header=T)
CC_73samples_pData <- readRDS('../data/PRJEB23709_Pre_73samples_phenoData.rds')
CC_73samples_GE <- CC_73samples_GE[CC_73samples_GE$Symbol != '',]
CC_73samples_GE_matrix <- CC_73samples_GE[,c('Symbol',CC_73samples_pData$sample)]


#### MGSP project: 103 samples
NatMed_103samples_GE_matrix <- readRDS('../data/NatMed_103samples_GE_matrix.rds')
NatMed_103samples_pData <- readRDS('../data/NatMed_103samples_pData.rds')


###################################################################################################################################################################
###################################################################################################################################################################

#### same function from Xiong et al. for each dataset
fitmodel <- function(expdata, features, phenoInfo, printngenes=FALSE){
  "
  expdata: col - sample, row - gene
  features: gene list (L1000 or nanostring)
  phenoInfo: dataframe: classs, sample
  "
  pData = data.frame(class=phenoInfo$class, sample=phenoInfo$sample,
                     row.names=phenoInfo$sample)
  phenoData <- new("AnnotatedDataFrame",data=pData)
  
  expdata_col_rearranged <- expdata[,c("Symbol",phenoInfo$sample)]
  expdata.sig <- expdata_col_rearranged[expdata_col_rearranged$Symbol %in% features,]
  expdata.sig <- expdata.sig[,-1]
  expdata.sig <- as.matrix(expdata.sig)
  #expdata.sig <- expdata.sig[rowSums(expdata.sig) > 1,]
  
  if (printngenes){
    print(dim(expdata.sig))
  }
  
  ExpSet_V5 <- ExpressionSet(assayData=as.matrix(expdata.sig),phenoData=phenoData)
  predictor_V5 <- fit(ExpSet_V5, method = "welch.test") # must remove 0 or low expression genes < 10, otherwise, fit function does not work!!!
  positive.class <- unique(pData(phenoData)[[1]])[2]
  negative.class <- unique(pData(phenoData)[[1]])[1]
  prediction_V5 <- predict(predictor_V5, ExpSet_V5, as.character(positive.class), ngenes=nrow(expdata.sig), dist = "cor")
  
  out_V5 <- as.factor(rep(c(1,2),c(table(pData(ExpSet_V5)[["class"]])[[negative.class]],table(pData(ExpSet_V5)[["class"]])[[positive.class]])))
  z_V5 <- as.numeric(prediction_V5@prediction[,'z'])
  Test_V5 <- cbind(out_V5,z_V5)
  colnames(Test_V5) <- c('grp','res')
  Test_V5 <- as.data.frame(Test_V5)
  
  #rocplot.single.V2(Test_GSE78220_V5$grp, Test_GSE78220_V5$res, title = "")
  
  p3_V5_data <- rocdata(Test_V5$grp, Test_V5$res)
  return(p3_V5_data$stats$auc)
}


#### fit multiple times
sampleSig <- function(expdata, genes, ngenes, phenoInfo, times){
  auclist <- c()
  genelist <- list()
  for (i in seq(1, times)){
    set.seed(i) # make result reproducible
    features <- sample(intersect(genes, expdata$Symbol), ngenes)
    a <- fitmodel(expdata,features,phenoInfo)
    print(i)
    auclist <- c(auclist, a)
    genelist[[i]] <- features
  }
  return(list(auclist=auclist, genelist = genelist))
}


################# 1.1 random select the same number genes with ImSig.genes from all genes in the expression dataset, use same pipeline and code with the original paper
############################################     data for Figure 1
iter <- 50
if (TRUE){
  # GSE78220：
  fitmodel(GSE78220_AltAnalyze, ImSig.genes, GSE78220_PhenoInfo2, printngenes = TRUE) # 103 signatures, AUC=0.98
  all.genes <- setdiff(GSE78220_AltAnalyze$Symbol, ImSig.genes)
  res1 <- sampleSig(expdata = GSE78220_AltAnalyze, genes = all.genes, ngenes = 103, phenoInfo = GSE78220_PhenoInfo2, iter)
  
  # GSE91061 = BMS038
  fitmodel(BMS038.Pre.CountTable.normalized.log, ImSig.genes, BMS038_PhenoInfo, printngenes = TRUE) # 103 signatures, AUC=0.96
  all.genes <- setdiff(BMS038.Pre.CountTable.normalized.log$Symbol, ImSig.genes)
  res2 <- sampleSig(expdata = BMS038.Pre.CountTable.normalized.log, genes = all.genes, ngenes = 103, phenoInfo = BMS038_PhenoInfo, iter)
  
  # CC_73samples from PRJEB23709
  fitmodel(CC_73samples_GE_matrix, ImSig.genes, CC_73samples_pData, printngenes = TRUE) # 101 signatures, AUC=0.86
  all.genes <- setdiff(CC_73samples_GE_matrix$Symbol, ImSig.genes)
  res3 <- sampleSig(expdata = CC_73samples_GE_matrix, genes = all.genes, ngenes = 101, phenoInfo = CC_73samples_pData, iter)
  
  #### MGSP project: 103 samples
  fitmodel(NatMed_103samples_GE_matrix, ImSig.genes, NatMed_103samples_pData, printngenes = TRUE) # 98 signatures, AUC=0.88
  all.genes <- setdiff(NatMed_103samples_GE_matrix$Symbol, ImSig.genes)
  res4 <- sampleSig(expdata = NatMed_103samples_GE_matrix, genes = all.genes, ngenes = 98, phenoInfo = NatMed_103samples_pData, iter)
  
  all.result <- data.frame(AUC=res1[['auclist']], Dataset=rep('GSE78220',length(res1[['auclist']])))
  all.result <- rbind(all.result, data.frame(AUC=res2[['auclist']], Dataset=rep('GSE91061',length(res2[['auclist']]))))
  all.result <- rbind(all.result, data.frame(AUC=res3[['auclist']], Dataset=rep('PRJEB23709',length(res3[['auclist']]))))
  all.result <- rbind(all.result, data.frame(AUC=res4[['auclist']], Dataset=rep('MGSP',length(res4[['auclist']]))))
}
ggboxplot(all.result, x="Dataset", y="AUC", color = "Dataset",
          palette = "jco", short.panel.labs = FALSE) + scale_x_discrete() + ylim(0,1.0) + xlab('') +
  theme(axis.title.y = element_text(size = 15, vjust = 0.5, hjust = 0.5)) + 
  geom_hline(aes(yintercept=0.7), colour="#990000", linetype="dashed")

group_by(all.result, Dataset) %>% summarise_each(funs(max))
group_by(all.result, Dataset) %>% summarise_each(funs(min))
group_by(all.result, Dataset) %>% summarise_each(funs(mean))
group_by(all.result, Dataset) %>% summarise_each(funs(which.max))
####################################################################################################################

################# 1.1 random select 100 genes in each dataset:
################# similar result with previous result
common_genes <- Reduce(intersect, list(GSE78220_AltAnalyze$Symbol, CC_73samples_GE_matrix$Symbol, BMS038.Pre.CountTable.normalized.log$Symbol, 
                                       NatMed_103samples_GE_matrix$Symbol))
length(common_genes)
#common_genes <- setdiff(common_genes, ImSig.genes) # no significant difference with or without ImSig in the total gene list
all.genes <- common_genes
iter <- 50
genenum <- 100
if (TRUE){
  # GSE78220：
  fitmodel(GSE78220_AltAnalyze, ImSig.genes, GSE78220_PhenoInfo2, printngenes = TRUE) # 103 signatures, AUC=0.98
  res1 <- sampleSig(expdata = GSE78220_AltAnalyze, genes = all.genes, ngenes = genenum, phenoInfo = GSE78220_PhenoInfo2, iter)
  
  # GSE91061 = BMS038
  fitmodel(BMS038.Pre.CountTable.normalized.log, ImSig.genes, BMS038_PhenoInfo, printngenes = TRUE) # 103 signatures, AUC=0.96
  res2 <- sampleSig(expdata = BMS038.Pre.CountTable.normalized.log, genes = all.genes, ngenes = genenum, phenoInfo = BMS038_PhenoInfo, iter)
  
  # CC_73samples from PRJEB23709
  fitmodel(CC_73samples_GE_matrix, ImSig.genes, CC_73samples_pData, printngenes = TRUE) # 101 signatures, AUC=0.86
  res3 <- sampleSig(expdata = CC_73samples_GE_matrix, genes = all.genes, ngenes = genenum, phenoInfo = CC_73samples_pData, iter)
  
  #### MGSP project: 103 samples
  fitmodel(NatMed_103samples_GE_matrix, ImSig.genes, NatMed_103samples_pData, printngenes = TRUE) # 98 signatures, AUC=0.88
  res4 <- sampleSig(expdata = NatMed_103samples_GE_matrix, genes = all.genes, ngenes = genenum, phenoInfo = NatMed_103samples_pData, iter)
  
  all.result <- data.frame(AUC=res1[['auclist']], Dataset=rep('GSE78220',length(res1[['auclist']])))
  all.result <- rbind(all.result, data.frame(AUC=res2[['auclist']], Dataset=rep('GSE91061',length(res2[['auclist']]))))
  all.result <- rbind(all.result, data.frame(AUC=res3[['auclist']], Dataset=rep('PRJEB23709',length(res3[['auclist']]))))
  all.result <- rbind(all.result, data.frame(AUC=res4[['auclist']], Dataset=rep('MGSP',length(res4[['auclist']]))))
}
ggboxplot(all.result, x="Dataset", y="AUC", color = "Dataset",
          palette = "jco", short.panel.labs = FALSE) + scale_x_discrete() + ylim(0,1.0) + xlab('') +
  theme(axis.title.y = element_text(size = 15, vjust = 0.5, hjust = 0.5)) + 
  geom_hline(aes(yintercept=0.7), colour="#990000", linetype="dashed")

group_by(all.result, Dataset) %>% summarise_each(funs(max))
min()group_by(all.result, Dataset) %>% summarise_each(funs(mean))
group_by(all.result, Dataset) %>% summarise_each(funs(which.max))


################# 2. use a resonable ML stratege: training dataset // test dataset (not used in the training dataset)
################# 2.1 use different dataset for train and test: 
table(GSE78220_phenoData@data$class)
table(BMS038_phenoData@data$class)
table(CC_73samples_phenoData@data$class)
table(NatMed_103samples_phenoData@data$class)

processData <- function(exp, phenoInfo, features, printngenes=FALSE){
  ### unify response name for cross dataset test
  phenoInfo$class <- as.character(phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'nonPD','NonResponder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'PD','Responder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'Nonresponder','NonResponder',phenoInfo$class)
  phenoInfo$class <- ifelse(phenoInfo$class == 'Progressor','NonResponder',phenoInfo$class)
  
  pData = data.frame(class=phenoInfo$class, sample=phenoInfo$sample,
                     row.names=phenoInfo$sample)
  phenoData <- new("AnnotatedDataFrame",data=pData)
  
  expdata_col_rearranged <- exp[,c("Symbol",phenoInfo$sample)]
  expdata.sig <- expdata_col_rearranged[expdata_col_rearranged$Symbol %in% features,]
  expdata.sig <- expdata.sig[,-1]
  expdata.sig <- as.matrix(expdata.sig)
  #expdata.sig <- expdata.sig[rowSums(expdata.sig) > 1,]
  rownames(expdata.sig) <- features
  
  if (printngenes){
    print(dim(expdata.sig))
  }
  ExpSet_V5 <- ExpressionSet(assayData=as.matrix(expdata.sig),phenoData=phenoData)
  return(ExpSet_V5)
  #return(list(exp=ExpSet_V5,genes=rownames(expdata.sig)))
}


trainANDtestModel <- function(traindata, testdata, features, phenoInfo_train, phenoInfo_test){
  "
  expdata: col - sample, row - gene
  features: gene list (L1000 or nanostring)
  phenoInfo: dataframe: classs, sample
  "
  common <- intersect(traindata$Symbol,testdata$Symbol)
  features <- intersect(common, features)
  
  ExpSet_train <- processData(exp = traindata, phenoInfo = phenoInfo_train, features = features)
  ExpSet_test <- processData(exp = testdata, phenoInfo = phenoInfo_test, features = features)
  
  predictor_V5 <- fit(ExpSet_train, method = "welch.test")
  positive.class <- unique(pData(ExpSet_test)$class)[2]
  negative.class <- unique(pData(ExpSet_test)$class)[1]
  #print(table(pData(ExpSet_test)$class))
  #print(length(features))
  prediction_V5 <- predict(predictor_V5, ExpSet_test, as.character(positive.class), ngenes=length(features), dist = "cor")
  
  out_V5 <- as.factor(rep(c(1,2),c(table(pData(ExpSet_test)[["class"]])[[negative.class]],table(pData(ExpSet_test)[["class"]])[[positive.class]])))
  z_V5 <- as.numeric(prediction_V5@prediction[,'z'])
  Test_V5 <- cbind(out_V5,z_V5)
  colnames(Test_V5) <- c('grp','res')
  Test_V5 <- as.data.frame(Test_V5)
  
  return(Test_V5)
}

#### change train dataset
testCrossDataset <- function(testdata, features, phenoInfo_test, testname='label'){
  test0 <- trainANDtestModel(traindata = GSE78220_AltAnalyze, testdata = testdata, 
                             features = features, phenoInfo_train = GSE78220_pData,phenoInfo_test = phenoInfo_test)
  colnames(test0) <- c(testname,'GSE78220')
  test1 <- trainANDtestModel(traindata = CC_73samples_GE_matrix, testdata = testdata, 
                             features = features, phenoInfo_train = CC_73samples_pData,phenoInfo_test = phenoInfo_test)
  colnames(test1) <- c(testname,'PRJEB23709')
  test2 <- trainANDtestModel(traindata = BMS038.Pre.CountTable.normalized.log, testdata = testdata, 
                             features = features, phenoInfo_train = BMS038_phenoData@data,phenoInfo_test = phenoInfo_test)
  colnames(test2) <- c(testname,'GSE91061')
  test3 <- trainANDtestModel(traindata = NatMed_103samples_GE_matrix, testdata = testdata, 
                             features = features, phenoInfo_train = NatMed_103samples_pData,phenoInfo_test = phenoInfo_test)
  colnames(test3) <- c(testname,'MGSP')
  test <- cbind(test0, test1['PRJEB23709'],test2['GSE91061'],test3['MGSP'])
  return(test)
}

#### fix traning dataset and change test dataset
testCrossDataset2 <- function(traindata, features, phenoInfo_train, testname='label'){
  test0 <- trainANDtestModel(traindata = traindata, testdata = GSE78220_AltAnalyze, 
                             features = features, phenoInfo_train = phenoInfo_train,phenoInfo_test = GSE78220_PhenoInfo2)
  colnames(test0) <- c('label','prediction')
  test0$Dataset <- 'GSE78220'
  
  test1 <- trainANDtestModel(traindata = traindata, testdata = CC_73samples_GE_matrix, 
                             features = features, phenoInfo_train = phenoInfo_train,phenoInfo_test = CC_73samples_pData)
  colnames(test1) <- c('label','prediction')
  test1$Dataset <- 'PRJEB23709'
  
  test2 <- trainANDtestModel(traindata = traindata, testdata = BMS038.Pre.CountTable.normalized.log, 
                             features = features, phenoInfo_train = phenoInfo_train,phenoInfo_test = BMS038_PhenoInfo)
  colnames(test2) <- c('label','prediction')
  test2$Dataset <- 'GSE91061'
  
  test3 <- trainANDtestModel(traindata = traindata, testdata = NatMed_103samples_GE_matrix, 
                             features = features, phenoInfo_train = phenoInfo_train, phenoInfo_test = NatMed_103samples_pData)
  colnames(test3) <- c('label','prediction')
  test3$Dataset <- 'MGSP'
  
  test <- rbind(test0, test1,test2,test3)
  return(test)
}


#### (gene siganature from Xiong et al.)
######################################### data for Figure 2
if (TRUE){
  test <- testCrossDataset2(traindata = GSE78220_AltAnalyze, features = ImSig.genes, phenoInfo_train = GSE78220_PhenoInfo2)
  write.csv(test, '../result/GSE78220_papersig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset2(traindata = CC_73samples_GE_matrix, features = ImSig.genes,phenoInfo_train = CC_73samples_pData)
  write.csv(test, '../result/PRJEB23709_papersig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset2(traindata = BMS038.Pre.CountTable.normalized.log, features = ImSig.genes,
                            phenoInfo_train = BMS038_phenoData@data)
  write.csv(test, '../result/GSE91061_papersig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset2(traindata = NatMed_103samples_GE_matrix, features = ImSig.genes,
                            phenoInfo_train = NatMed_103samples_pData)
  write.csv(test, '../result/MGSP_papersig_prediction.csv', quote = FALSE, row.names = FALSE)
}
####################################################################################################################

#### use a random sigature:
set.seed(17)
common_genes <- Reduce(intersect, list(GSE78220_AltAnalyze$Symbol, CC_73samples_GE_matrix$Symbol, BMS038.Pre.CountTable.normalized.log$Symbol, 
                                       NatMed_103samples_GE_matrix$Symbol))
length(common_genes)
common_genes <- setdiff(common_genes, ImSig.genes)
random.sig <- sample(common_genes,100)
if (TRUE){
  test <- testCrossDataset(testdata = GSE78220_AltAnalyze, features = random.sig,phenoInfo_test = GSE78220_PhenoInfo2)
  write.csv(test, '../result/GSE78220_randomsig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset(testdata = CC_73samples_GE_matrix, features = random.sig,phenoInfo_test = CC_73samples_pData)
  write.csv(test, '../result/PRJEB23709_randomsig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset(testdata = BMS038.Pre.CountTable.normalized.log, features = random.sig,
                           phenoInfo_test = BMS038_phenoData@data)
  write.csv(test, '../result/GSE91061_randomsig_prediction.csv', quote = FALSE, row.names = FALSE)
  
  test <- testCrossDataset(testdata = NatMed_103samples_GE_matrix, features = random.sig,
                           phenoInfo_test = NatMed_103samples_pData)
  write.csv(test, '../result/MGSP_randomsig_prediction.csv', quote = FALSE, row.names = FALSE)
}

################# 2.2 use the same dataset for train and test but split them into two datasets: 
table(GSE78220_PhenoInfo2$class) # 15:13 
table(CC_73samples_pData$class) # 27:46
table(BMS038_PhenoInfo$class) # 25:26 -> GSE91061
table(NatMed_103samples_pData$class) # 56:47 -> MGSP

library(caret)
testCV <- function(expdata, phenoInfo, features, num=5){
  set.seed(17)
  folds <- createFolds(y=phenoInfo$sample,k=num)
  auc_value<-as.numeric()
  test <- data.frame(label=c(), prediction=c(), nfold = c())
  for (i in 1:5){
    fold_test <- phenoInfo[folds[[i]],] #folds[[i]] for test
    fold_train <- phenoInfo[-folds[[i]],] # remaining data for train
    print(table(fold_train$class))
    print(table(fold_test$class))
    
    test0 <- trainANDtestModel(traindata = expdata, testdata = expdata, 
                               features = features, phenoInfo_train = fold_train, phenoInfo_test = fold_test)
    colnames(test0) <- c('label','prediction')
    auc_value <- append(auc_value, as.numeric(rocdata(test0$label, test0$prediction)$stats$auc))
    if (i==1){
      test <- test0
    }else{
      test <- rbind(test, test0)
    }
  }
  return(list(auc=auc_value, result=test))
}

test0 <- testCV(expdata = BMS038.Pre.CountTable.normalized.log, phenoInfo = BMS038_PhenoInfo, features = ImSig.genes)[['result']]
rocdata(test0$label, test0$prediction)
test1 <- testCV(expdata = BMS038.Pre.CountTable.normalized.log, phenoInfo = BMS038_PhenoInfo, features = random.sig)[['result']]
rocdata(test1$label, test1$prediction)
test <- cbind(test0, test1)