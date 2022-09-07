
## Run it with the agnostic then the mechanistic
# The desired outcome is to repeat the training process n times and calculate the AUC in the testing data for each n repeat
# Take the difference between the AUC in the test data using mechanistic and agnostic approach
# AUC_Test_Mech - AUC_Test_Agnostic
# Then calculate the Confidence interval of this difference


rm(list = ls())

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)
library(boot)
library(tidyverse)
library(patchwork)
library(qdapRegex)
library(data.table)

### Load expression and phenotype data
load("./Objs/progressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

#########
## Load the selected genes *(TF-MiR from Lotte)
TF_MiR <- load("/Users/mohamedomar/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/Genes/allTSPs.rda")

#TF_MiR <- load("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/Genes/TF_MiR.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#########
### TRAINING using restricted pairs
#######

### Set Feature number and max k
ktsp <- c(1:50) #12
featNo <- nrow(usedTrainMat)

# Combine the expression matrix with the phenotype in 1 data frame
Data <- cbind(t(usedTrainMat), usedTrainGroup)
Data <- as.data.frame(Data)
Data$usedTrainGroup <- as.factor(Data$usedTrainGroup)
levels(Data[, "usedTrainGroup"]) <- c("NoProgression", "Progression")


######################################################
## Mechanistic and agnostic using top 74 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=74)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  # get the pairs
  rownames(KTSP_Train_Agnostic$TSPs) <- gsub(',', '-', rownames(KTSP_Train_Agnostic$TSPs))
  rownames(KTSP_Train_Mech$TSPs) <- gsub(',', '-', rownames(KTSP_Train_Mech$TSPs))
  pairs_agnostic <- paste(rownames(KTSP_Train_Agnostic$TSPs), collapse = ',')
  pairs_Mech <- paste(rownames(KTSP_Train_Mech$TSPs), collapse = ',')
  # get the invidiual genes
  gene1_agnostic <- paste(KTSP_Train_Agnostic$TSPs[,1], collapse = ',')
  gene2_agnostic <- paste(KTSP_Train_Agnostic$TSPs[,2], collapse = ',')
  gene1_Mech <- paste(KTSP_Train_Mech$TSPs[,1], collapse = ',')
  gene2_Mech <- paste(KTSP_Train_Mech$TSPs[,2], collapse = ',')
  # predictions
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic, pairs_agnostic, pairs_Mech, gene1_agnostic, gene2_agnostic, gene1_Mech, gene2_Mech))
}

set.seed(333)
bootobject_74 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

######################################################
## Mechanistic and agnostic using top 100 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=100)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_100 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

######################################################
## Mechanistic and agnostic using top 200 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=200)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_200 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

######################################################
## Mechanistic and agnostic using top 500 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=500)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_500 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

############################################################
## Save all
save(bootobject_74, bootobject_100, bootobject_200, bootobject_500, file = "./Objs/KTSP/bootobjectKTSP.rda")

load("./Objs/KTSP/bootobjectKTSP.rda")


##############################################################
### Work with boot object 50  
All_74 <- bootobject_74$t
colnames(All_74) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic", "pairs_agnostic", "pairs_Mech", "gene1_agnostic", "gene2_agnostic", "gene1_Mech", "gene2_Mech")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_74 <- All_74[,"Diff_Agnostic"]
range(Diff_Agnostic_74)
quantile(Diff_Agnostic_74, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_74 <- All_74[,"Diff_Mechanistic"]
range(Diff_Mechanistic_74)
quantile(Diff_Mechanistic_74, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))

##########
# N of pairs
mean(All_74[,'N_Pairs_Agnostic'])
mean(All_74[,'N_Pairs_Mech'])

#########

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_74 <- data.frame(AUC = All_74[, "AUC_Train_Mech"])
AgnosticAUCTrain_74 <- data.frame(AUC = All_74[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_74$modelType <- "Mechanistic"
AgnosticAUCTrain_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_74 <- rbind(MechanisticAUCTrain_74, AgnosticAUCTrain_74)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_74 <- data.frame(AUC = All_74[, "AUC_Test_Mech"])
AgnosticAUCTest_74 <- data.frame(AUC = All_74[, "AUC_Test_Agnostic"])

MechanisticAUCTest_74$modelType <- "Mechanistic"
AgnosticAUCTest_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_74 <- rbind(MechanisticAUCTest_74, AgnosticAUCTest_74)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_74$data_type <- "Training"
ModelCompareAUCTest_74$data_type <- "Testing"

ModelCompareAUCTrain_74$NofFeatAgn <- "74_Genes"
ModelCompareAUCTest_74$NofFeatAgn <- "74_Genes"

###########################################################################
## Save for the main figure
ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_74, ModelCompareAUCTest_74)
ModelCompare_KTSP$algorithm <- "KTSP"
save(ModelCompare_KTSP, file = "./Objs/KTSP/ModelCompare_KTSP.rda")

################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_74 <- data.frame(NofPairs = All_74[, "N_Pairs_Mech"])
Agnostic_NofPairs_74 <- data.frame(NofPairs = All_74[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_74$modelType <- "Mechanistic"
Agnostic_NofPairs_74$modelType <- "Agnostic"

ModelCompare_NofPairs_74 <- rbind(Mechanistic_NofPairs_74, Agnostic_NofPairs_74)

ModelCompare_NofPairs_74$NofFeatAgn <- "74 Genes"


###############
# get the common pairs
pairs_agnostic <- strsplit(All_74[, 'pairs_agnostic'], ',')
pairs_mech <- strsplit(All_74[, 'pairs_Mech'], ',')

pairs_agnostic_df <- plyr::ldply(pairs_agnostic, rbind)
pairs_mech_df <- plyr::ldply(pairs_mech, rbind)

find_rep <- function(dat, feature){
  # NAs create problems in the function so we substitute that with "unknown"
  dat[is.na(dat)] <- "unknown"
  rep_rows <- sum(apply(dat, 1,  function(x) any(x == feature)))
  names(rep_rows) <- feature
  as.data.frame(rep_rows)
}

features_agnostic <- na.omit(unique(as.vector(as.matrix(pairs_agnostic_df))))
list_results_agnostic <- lapply(features_agnostic, find_rep, dat = pairs_agnostic_df)
sum_result_agnostic <- do.call(rbind, list_results_agnostic)
sum_result_agnostic$feature <- rownames(sum_result_agnostic)
sum_result_agnostic <- as.data.frame(sum_result_agnostic[order(sum_result_agnostic$rep_rows, decreasing = T), ])
#sum_result_agnostic_top <- sum_result_agnostic[c(1:20), ]

features_mech <- na.omit(unique(as.vector(as.matrix(pairs_mech_df))))
list_results_mech <- lapply(features_mech, find_rep, dat = pairs_mech_df)
sum_result_mech <- do.call(rbind, list_results_mech)
sum_result_mech$feature <- rownames(sum_result_mech)
sum_result_mech <- as.data.frame(sum_result_mech[order(sum_result_mech$rep_rows, decreasing = T), ])
#sum_result_mech_top <- sum_result_mech[c(1:20), ]


###########
## histogram
# # join the agnostic and mechanistic frequency
sum_result_agnostic_freq <- sum_result_agnostic
sum_result_agnostic_freq$type <- 'agnostic'
sum_result_agnostic_freq$frequency <- sum_result_agnostic_freq$rep_rows/1000
sum_result_agnostic_freq <- sum_result_agnostic_freq[c(1:93), ]
sum_result_agnostic_freq$TSP <- factor(1:nrow(sum_result_agnostic_freq))
levels(sum_result_agnostic_freq$TSP) <- paste0("TSP", levels(sum_result_agnostic_freq$TSP))

sum_result_mech_freq <- sum_result_mech
sum_result_mech_freq$type <- 'mechanistic'
sum_result_mech_freq$frequency <- sum_result_mech_freq$rep_rows/1000
sum_result_mech_freq <- sum_result_mech_freq[c(1:93), ]
sum_result_mech_freq$TSP <- factor(1:nrow(sum_result_mech_freq))
levels(sum_result_mech_freq$TSP) <- paste0("TSP", levels(sum_result_mech_freq$TSP))

sum_result_all <- rbind(sum_result_agnostic_freq, sum_result_mech_freq)
table(sum_result_all$type)
sum_result_all <- sum_result_all[order(sum_result_all$rep_rows, decreasing = T), ]

png('/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/bladder_pairs_frequency.png', width = 2000, height = 1500, res = 400)
ggplot(sum_result_all, aes(x= TSP, y = rep_rows, fill = type)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  ylab("frequency") + 
  theme(axis.text.x=element_text(angle=30,hjust=1, size = 0.3, colour = 'black'),
        axis.text.y=element_text(size = 8, colour = 'black'),
        axis.title.x=element_text(size = 10, colour = 'black'),
        axis.title.y=element_text(size = 10, colour = 'black'),
        axis.ticks.length=unit(0.05,"cm"), 
        panel.background = element_rect(fill = "transparent", color = NA_character_),
        axis.line = element_line(colour = 'black', size = 0.1), 
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10, colour = 'black'))+
        scale_y_continuous(expand = c(0, 0))
dev.off() 

#########
## get the individual genes for pairs of interest

# agnostic
agnostic_indvGns_good <- All_74 %>%
  as.data.frame() %>%
  dplyr::select(pairs_agnostic, gene1_agnostic, gene2_agnostic) %>%
  #mutate(tmp = strsplit(as.character(pairs_Mech),',')) %>%
  #unnest(tmp) %>%
  group_by(pairs_agnostic) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp) %>%
  #pivot_longer(cols = c(4:13), values_to = 'pair', values_drop_na = TRUE) %>%
  #select(-name) %>%
  #mutate(tmp_gene1 = strsplit(as.character(gene1_Mech),',')) %>%
  #group_by(pair) %>%
  #unnest(tmp_gene1, keep_empty = T) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp_gene1) %>%
  separate_rows(c(pairs_agnostic, gene1_agnostic, gene2_agnostic), sep = ',') %>%
  #select(-pairs_Mech) %>%
  dplyr::rename(gene1=gene1_agnostic, gene2=gene2_agnostic)

row.names(sum_result_agnostic) <- NULL

agnostic_genes_notGood <- sum_result_agnostic %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  column_to_rownames(var = 'feature') %>%
  select(-tmp) %>%
  relocate(rep_rows, .after = gene2)

agnostic_genes_freq <- data.frame(rep_rows = agnostic_genes_notGood$rep_rows, 
                                  pairs_agnostic = rownames(agnostic_genes_notGood),
                                  row.names = rownames(agnostic_genes_notGood))


# merge
agnostic_indvGns_good_unique <- agnostic_indvGns_good[!duplicated(agnostic_indvGns_good$pairs_agnostic), ]
agnostic_indvGns_good_clean <- merge(x = agnostic_indvGns_good_unique, y = agnostic_genes_freq, 
                                     by = 'pairs_agnostic', suffixes = colnames(agnostic_indvGns_good_unique))

agnostic_indvGns_good_clean <- agnostic_indvGns_good_clean[order(agnostic_indvGns_good_clean$rep_rows, decreasing = T), ] 
agnostic_indvGns_good_clean$pairs_agnostic <- NULL
#agnostic_indvGns_good_clean <- filter(agnostic_indvGns_good_clean, !grepl("///", gene1, ignore.case = TRUE))
#agnostic_indvGns_good_clean <- filter(agnostic_indvGns_good_clean, !grepl("///", gene2, ignore.case = TRUE))
agnostic_indvGns_good_clean <- aggregate(rep_rows~(gene1+ gene2),data=agnostic_indvGns_good_clean, FUN = sum)
agnostic_indvGns_good_clean <- agnostic_indvGns_good_clean[order(agnostic_indvGns_good_clean$rep_rows, decreasing = T), ] 


## matrix for vertices
agnostic_vertics <- agnostic_indvGns_good_clean %>%
  pivot_longer(cols = c('gene1', 'gene2'), values_to = 'gene') %>%
  rename(pair_frequency = rep_rows, order = name) %>%
  group_by(gene) %>%
  mutate(gene_frequency = n()) %>%
  relocate(gene, .before = pair_frequency) %>%
  relocate(order, .after = gene) %>%
  group_by(gene) %>%
  dplyr::filter(!duplicated(gene))

# get the top pairs only
agnostic_indvGns_good_clean <- agnostic_indvGns_good_clean[c(1:93), ]
agnostic_vertics <- agnostic_vertics[agnostic_vertics$gene %in% agnostic_indvGns_good_clean$gene1 | agnostic_vertics$gene %in% agnostic_indvGns_good_clean$gene2, ]

########
# mechanistic
mech_indvGns_good <- All_74 %>%
  as.data.frame() %>%
  dplyr::select(pairs_Mech, gene1_Mech, gene2_Mech) %>%
  #mutate(tmp = strsplit(as.character(pairs_Mech),',')) %>%
  #unnest(tmp) %>%
  group_by(pairs_Mech) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp) %>%
  #pivot_longer(cols = c(4:13), values_to = 'pair', values_drop_na = TRUE) %>%
  #select(-name) %>%
  #mutate(tmp_gene1 = strsplit(as.character(gene1_Mech),',')) %>%
  #group_by(pair) %>%
  #unnest(tmp_gene1, keep_empty = T) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp_gene1) %>%
  separate_rows(c(pairs_Mech, gene1_Mech, gene2_Mech), sep = ',') %>%
  #select(-pairs_Mech) %>%
  dplyr::rename(gene1=gene1_Mech, gene2=gene2_Mech)

row.names(sum_result_mech) <- NULL
mech_genes_notGood <- sum_result_mech %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  column_to_rownames(var = 'feature') %>%
  select(-tmp) %>%
  relocate(rep_rows, .after = gene2)

mech_genes_freq <- data.frame(rep_rows = mech_genes_notGood$rep_rows, 
                              pairs_Mech = rownames(mech_genes_notGood),
                              row.names = rownames(mech_genes_notGood))

# merge
mech_indvGns_good_unique <- mech_indvGns_good[!duplicated(mech_indvGns_good$pairs_Mech), ]
mech_indvGns_good_clean <- merge(x = mech_indvGns_good_unique, y = mech_genes_freq, 
                                 by = 'pairs_Mech', suffixes = colnames(mech_indvGns_good_unique))

mech_indvGns_good_clean <- mech_indvGns_good_clean[order(mech_indvGns_good_clean$rep_rows, decreasing = T), ] 
mech_indvGns_good_clean$pairs_Mech <- NULL

## matrix for vertices
mech_vertics <- mech_indvGns_good_clean %>%
  pivot_longer(cols = c('gene1', 'gene2'), values_to = 'gene') %>%
  rename(pair_frequency = rep_rows, order = name) %>%
  group_by(gene) %>%
  mutate(gene_frequency = n()) %>%
  relocate(gene, .before = pair_frequency) %>%
  relocate(order, .after = gene) %>%
  group_by(gene) %>%
  dplyr::filter(!duplicated(gene))

# get the top pairs only
#mech_indvGns_good_clean <- mech_indvGns_good_clean[c(1:100), ]
mech_vertics <- mech_vertics[mech_vertics$gene %in% mech_indvGns_good_clean$gene1 | mech_vertics$gene %in% mech_indvGns_good_clean$gene2, ]


#####################################
# bar plots

Mechfreq <- mech_indvGns_good_clean %>%
  mutate(feature = paste(gene1, gene2, sep = '-')) %>%
  select(-gene1, -gene2) %>%
  top_n(10, wt = rep_rows)

Agnosticfreq <- agnostic_indvGns_good_clean %>%
  mutate(feature = paste(gene1, gene2, sep = '-')) %>%
  select(-gene1, -gene2) %>%
  top_n(10, wt = rep_rows)

MechFreq_barplot <- ggplot(data=Mechfreq, aes(x=rep_rows, y=reorder(feature, rep_rows))) +
  geom_col(width=0.5, fill='#00BFC4') + 
  #coord_cartesian(xlim = c(1, 600))
  scale_x_continuous(limits = c(1,600), n.breaks =5, oob = scales::squish, expand = c(0, 0)) +
  labs(y = "Pair", x = "Frequency", title = "Mechanistic") +
  theme(plot.title = element_text(size=16, hjust=0.5, face = 'plain', color = 'black'), 
        axis.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(size=10, color = 'black'),
        axis.title.x = element_text(size=12, color = 'black'),
        axis.title.y = element_text(size=12, color = 'black')
  )

AgnosticFreq_barplot <- ggplot(data=Agnosticfreq, aes(x=rep_rows, y=reorder(feature, rep_rows))) +
  geom_col(width=0.5, fill='#F8766D') + 
  #coord_cartesian(xlim = c(1, 600))
  scale_x_continuous(limits = c(1,75), n.breaks =5, oob = scales::squish, expand = c(0, 0)) +
  labs(y = "Pair", x = "Frequency", title = "Agnostic") +
  theme(plot.title = element_text(size=16, hjust=0.5, face = 'plain', color = 'black'), 
        axis.text.y = element_text(size=10, color = 'black'),
        axis.text.x = element_text(size=10, color = 'black'),
        axis.title.x = element_text(size=12, color = 'black'),
        axis.title.y = element_text(size=12, color = 'black')
  )


tiff(filename = "/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/bladder_frequency.tiff", width = 3000, height = 2000, res = 300)
((AgnosticFreq_barplot + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) | 
    (MechFreq_barplot + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) 
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The 10 most frequent gene pairs returned by the k-TSPs model',
    tag_levels = c('A'),
    theme = theme(plot.title = element_text(size = 20, face = "plain", hjust=0.5))
  )
dev.off()

########################################
## network plots

library(igraph)
library(qgraph)

# make the graph
network_mech <- graph_from_data_frame(mech_indvGns_good_clean, 
                                      directed = T, 
                                      vertices = mech_vertics)

# edges parameters
E(network_mech)$width <- log2(E(network_mech)$rep_rows)/2
E(network_mech)$edge.color <- "black"
edge.start <- ends(network_mech, es=E(network_mech), names=F)[,1]

# vertices parameters
V(network_mech)$size <- V(network_mech)$gene_frequency*0.7
V(network_mech)$color <- ifelse(V(network_mech)$order == 'gene1', "tomato", "gold")

# filter edges
#hist(mech_indvGns_good_clean$rep_rows)
mean(mech_indvGns_good_clean$rep_rows)
sd(mech_indvGns_good_clean$rep_rows)
cut.off <- mean(mech_indvGns_good_clean$rep_rows) 
network_mech_FilEdges <- delete_edges(network_mech, E(network_mech)[rep_rows<cut.off])

l_mech <- qgraph.layout.fruchtermanreingold(get.edgelist(network_mech,names=FALSE),vcount=vcount(network_mech),
                                            area=8*(vcount(network_mech)^2),repulse.rad=(vcount(network_mech)^3.1))

###
## plot
set.seed(00000)
tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/gene_networks/bladder_mech.tiff', width =  6000, height = 6000, res = 600)
#l_mech <- layout_nicely(network_mech)
#l_mech <- norm_coords(l_mech, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5) #default -- scaled
plot(network_mech, 
     #vertex.size=degree(network_mech), 
     vertex.label.dist=0,
     #edge.arrow.width=0, 
     edge.arrow.size=0, 
     arrow.size = 0, 
     arrow.width = 0,
     vertex.label.cex=1,
     vertex.label.degree = -pi/2,
     vertex.label.color = "black",
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     rescale=T,
     layout = l_mech, 
     edge.curved=0.2,
     margin = -0.2,
     #main = 'Mechanistic'
     )
#legend(x=-1.5, y=-1.0, c("gene1","gene2"), pch=21,
#       col="#777777", pt.bg=c("tomato", "gold"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

#####
# cfg
cfg_mech <- cluster_fast_greedy(as.undirected(network_mech))
set.seed(00000)
tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/gene_networks/bladder_mech_cfg.tiff', width =  5000, height = 3000, res = 300)
plot(cfg_mech, network_mech, 
     #vertex.size=degree(network_mech), 
     vertex.label.dist=0,
     #edge.arrow.width=0, 
     edge.arrow.size=0, 
     arrow.size = 0, 
     arrow.width = 0,
     vertex.label.cex=0.7,
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     rescale=T,
     layout = l_mech, 
     edge.curved=0.2,
     #main = 'Mechanistic'
     )
#legend(x=-1.5, y=-1.0, c("gene1","gene2"), pch=21,
#       col="#777777", pt.bg=c("tomato", "gold"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

############################
## agnostic

# make the graph
network_agnostic <- graph_from_data_frame(agnostic_indvGns_good_clean, 
                                          directed = T, 
                                          vertices = agnostic_vertics)

# edges parameters
E(network_agnostic)$width <- log2(E(network_agnostic)$rep_rows)/2
E(network_agnostic)$edge.color <- "black"
edge.start <- ends(network_agnostic, es=E(network_agnostic), names=F)[,1]

# vertices parameters
V(network_agnostic)$size <- V(network_agnostic)$gene_frequency*0.1
V(network_agnostic)$color <- ifelse(V(network_agnostic)$order == 'gene1', "tomato", "gold")

# filter edges
#hist(agnostic_indvGns_good_clean$rep_rows)
mean(agnostic_indvGns_good_clean$rep_rows)
sd(agnostic_indvGns_good_clean$rep_rows)
cut.off <- mean(agnostic_indvGns_good_clean$rep_rows) 
network_agnostic_FilEdges <- delete_edges(network_agnostic, E(network_agnostic)[rep_rows<cut.off])

l_agnostic <- qgraph.layout.fruchtermanreingold(get.edgelist(network_agnostic,names=FALSE),vcount=vcount(network_agnostic),
                                                area=8*(vcount(network_agnostic)^2),repulse.rad=(vcount(network_agnostic)^3.1))

###
## plot
set.seed(00000)
tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/gene_networks/bladder_agnostic.tiff', width = 7000, height = 7000, res = 600)
#l_agnostic <- layout_nicely(network_agnostic)
#l_agnostic <- norm_coords(l_agnostic, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5) #default -- scaled
plot(network_agnostic, 
     #vertex.size=degree(network_agnostic), 
     vertex.label.dist=0,
     #vertex.label.degree = 0,
     vertex.label.cex=1,
     vertex.label.degree = -pi/2,
     vertex.label.color = "black",
     
     edge.arrow.width=0, 
     edge.arrow.size=0, 
     #edge.width = 0.5,
     arrow.size = 0, 
     arrow.width = 0,
     #label.cex=0.05,
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     layout = l_agnostic,
     edge_curved =0.2,
     margin = -0.1,
     #main = 'Agnostic' 
)
#legend(x=-1.5, y=-1.0, c("gene1","gene2"), pch=21,
#       col="#777777", pt.bg=c("tomato", "gold"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

#####
# cfg
cfg_agnostic <- cluster_fast_greedy(as.undirected(network_agnostic))
tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/gene_networks/bladder_agnostic_cfg.tiff', width =  5000, height = 3000, res = 300)
plot(cfg_agnostic, network_agnostic, 
     #vertex.size=degree(network_mech), 
     vertex.label.dist=0,
     #edge.arrow.width=0, 
     edge.arrow.size=0, 
     arrow.size = 0, 
     arrow.width = 0,
     #label.cex=0.05,
     
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     
     vertex.label.cex=0.7,
     layout = l_agnostic, 
     edge.curved=0.2,
     #main = 'Mechanistic gene pairs in predicting bladder cancer progression (clusters)'
     )
#legend(x=-1.5, y=-1.0, c("gene1","gene2"), pch=21,
#       col="#777777", pt.bg=c("tomato", "gold"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


##########################################
# Frequency per gene sets 
library(msigdbr)

genes_mech <- mech_vertics$gene
genes_agnostic <- agnostic_vertics$gene

# get the hallmarks gene sets
hallmarks_gs <- msigdbr(species = "human", category = "H")

hallmarks_gs_long <- hallmarks_gs %>%
  dplyr::select(gs_name, human_gene_symbol)

# if genes exist in mech genes, add 1 otherwise 0
hallmarks_gs_long_mech <- hallmarks_gs_long %>%
  dplyr::distinct() %>%
  dplyr::mutate(mech_existence = dplyr::if_else(
    human_gene_symbol %in% genes_mech, 
    1, 0
  ))

hallmarks_gs_long_agnostic <- hallmarks_gs_long %>%
  dplyr::distinct() %>%
  dplyr::mutate(agnostic_existence = dplyr::if_else(
    human_gene_symbol %in% genes_agnostic, 
    1, 0
  ))

# long format
hallmarks_gs_long_mech <- hallmarks_gs_long_mech %>%
  tidyr::pivot_wider(
    names_from = gs_name,
    values_from = mech_existence,
    values_fill = 0
  )

hallmarks_gs_long_agnostic <- hallmarks_gs_long_agnostic %>%
  tidyr::pivot_wider(
    names_from = gs_name,
    values_from = agnostic_existence,
    values_fill = 0
  )

# filter to mech genes
hallmarks_gs_long_mech <- hallmarks_gs_long_mech %>%
  dplyr::filter(
    human_gene_symbol %in% genes_mech
  )

# filter to agnostic genes
hallmarks_gs_long_agnostic <- hallmarks_gs_long_agnostic %>%
  dplyr::filter(
    human_gene_symbol %in% genes_agnostic
  ) 

# add missing genes
missing_mech <- genes_mech[! genes_mech %in% unique(hallmarks_gs_long_mech$human_gene_symbol)]
missing_agnostic <- genes_agnostic[! genes_agnostic %in% unique(hallmarks_gs_long_agnostic$human_gene_symbol)]

hallmarks_gs_long_mech <- hallmarks_gs_long_mech %>%
  dplyr::add_row(
    human_gene_symbol = missing_mech
  )

hallmarks_gs_long_agnostic <- hallmarks_gs_long_agnostic %>%
  dplyr::add_row(
    human_gene_symbol = missing_agnostic
  )

# replace na with 0
hallmarks_gs_long_mech[is.na(hallmarks_gs_long_mech)] <- 0
hallmarks_gs_long_agnostic[is.na(hallmarks_gs_long_agnostic)] <- 0

#####
# make genes as rownames
hallmarks_gs_long_mech <- as.data.frame(hallmarks_gs_long_mech)
rownames(hallmarks_gs_long_mech) <- hallmarks_gs_long_mech$human_gene_symbol
hallmarks_gs_long_mech$human_gene_symbol <- NULL

hallmarks_gs_long_agnostic <- as.data.frame(hallmarks_gs_long_agnostic)
rownames(hallmarks_gs_long_agnostic) <- hallmarks_gs_long_agnostic$human_gene_symbol
hallmarks_gs_long_agnostic$human_gene_symbol <- NULL

#####
# matrix then fix the names
hallmarks_gs_long_mech <- as.matrix(hallmarks_gs_long_mech)
colnames(hallmarks_gs_long_mech) <- gsub('HALLMARK_', '', colnames(hallmarks_gs_long_mech))
colnames(hallmarks_gs_long_mech) <- gsub('_', ' ', colnames(hallmarks_gs_long_mech))

hallmarks_gs_long_agnostic <- as.matrix(hallmarks_gs_long_agnostic)
colnames(hallmarks_gs_long_agnostic) <- gsub('HALLMARK_', '', colnames(hallmarks_gs_long_agnostic))
colnames(hallmarks_gs_long_agnostic) <- gsub('_', ' ', colnames(hallmarks_gs_long_agnostic))

####
## matrix for vertices
#mech_sets <- data.frame(names = colnames(hallmarks_gs_long_mech), type = 'geneset')
#mech_genes <- data.frame(names = rownames(hallmarks_gs_long_mech), type = 'gene')

#gs_mech_vertics <- rbind(mech_sets, mech_genes)

hallmarks_gs_long_mech_fil <- hallmarks_gs_long_mech[, colSums(hallmarks_gs_long_mech) !=0 ]
hallmarks_gs_long_mech_adj <- t(hallmarks_gs_long_mech_fil) %*% hallmarks_gs_long_mech_fil

hallmarks_gs_long_agnostic_fil <- hallmarks_gs_long_agnostic[, colSums(hallmarks_gs_long_agnostic) !=0 ]
hallmarks_gs_long_agnostic_adj <- t(hallmarks_gs_long_agnostic_fil) %*% hallmarks_gs_long_agnostic_fil

########################
# the network
network_gs_mech <- graph.adjacency(hallmarks_gs_long_mech_adj, diag = F)
network_gs_agnostic <- graph.adjacency(hallmarks_gs_long_agnostic_adj, diag = F)

network_gs_mech <- simplify(network_gs_mech)
network_gs_agnostic <- simplify(network_gs_agnostic)

l_mech_gs <- qgraph.layout.fruchtermanreingold(get.edgelist(network_gs_mech,names=FALSE),vcount=vcount(network_gs_mech),
                                       area=30*(vcount(network_gs_mech)^2),repulse.rad=(vcount(network_gs_mech)^3.1), niter=10000)

l_agnostic_gs <- qgraph.layout.fruchtermanreingold(get.edgelist(network_gs_agnostic,names=FALSE),vcount=vcount(network_gs_agnostic),
                                                   area=30*(vcount(network_gs_agnostic)^2), repulse.rad=(vcount(network_gs_agnostic)^3.1), niter=10000)

# CFG
gs_cfg_mech <- cluster_fast_greedy(as.undirected(network_gs_mech))
gs_cfg_agnostic <- cluster_fast_greedy(as.undirected(network_gs_agnostic))


# DEGREE: the number of genesets it overlaps with


# set.seed(333)
# tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/genesets_networks/bladder_mech_gs.tiff', width = 5000, height = 4500, res = 400)
# plot(network_gs_mech, 
#      vertex.size=log(degree(network_gs_mech, mode = 'out')+1)*4, 
#      vertex.label.dist=0,
#      edge.arrow.width=0, 
#      edge.arrow.size=0, 
#      vertex.label.cex=1,
#      vertex.label.degree = -pi/2,
#      vertex.label.color = "black",
#      
#      #edge.width = 0.5,
#      arrow.size = 0, 
#      arrow.width = 0,
#      #label.cex=0.05,
#      #vertex.frame.color = NA,
#      vertex.frame.width = 0.2,
#      layout = l_mech_gs,
#      edge_curved =0.2,
#      #main = 'Mechanistic',
#      margin = -0.1
# )
# 
# dev.off()

#l_agnostic_gs <- layout_nicely(network_gs_agnostic)
#l_agnostic_gs <- norm_coords(l_agnostic_gs, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5) #default -- scaled


# set.seed(333)
# tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/genesets_networks/bladder_agnostic_gs.tiff', width =10000, height =8000, res = 700)
# plot(network_gs_agnostic, 
#      vertex.size=log(degree(network_gs_agnostic, mode = 'out')+1)*5, 
#      vertex.label.dist=0,
#      edge.arrow.width=0.1, 
#      edge.arrow.size=0.1, 
#      #edge.width = 0.5,
#      arrow.size = 0, 
#      arrow.width = 0,
#      vertex.label.cex=0.8,
#      #vertex.frame.color = NA,
#      vertex.frame.width = 0.2,
#      layout = l_agnostic_gs,
#      edge_curved =0.2,
#      #main = 'Agnostic',
#      margin = -0.1
# )
# dev.off()

set.seed(333)
tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/genesets_networks/bladder_mech_cfg_gs.tiff', width =10000, height =8000, res = 700)
plot(gs_cfg_mech, network_gs_mech, 
     vertex.size=log(degree(network_gs_mech, mode = 'out')+1)*5, 
     vertex.label.dist=0,
     edge.arrow.width=0, 
     edge.arrow.size=0.0, 
     
     vertex.label.cex=1,
     vertex.label.degree = -pi/2,
     vertex.label.color = "black",
     
     edge.width = 0.2,
     arrow.size = 0, 
     arrow.width = 0,
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     layout = l_mech_gs,
     edge_curved =0.2,
     #main = 'Mechanistic',
     margin = -0.06
)
dev.off()

tiff(filename = '/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/genesets_networks/bladder_agnostic_cfg_gs.tiff', width =10000, height =10000, res = 800)
plot(gs_cfg_agnostic, network_gs_agnostic, 
     vertex.size=log(degree(network_gs_agnostic, mode = 'out')+1)*5, 
     vertex.label.dist=0,
     edge.arrow.width=0.0, 
     edge.arrow.size=0.0, 
     
     vertex.label.cex=1.2,
     vertex.label.degree = -pi/2,
     vertex.label.color = "black",
     
     edge.width = 0.2,
     arrow.size = 0, 
     arrow.width = 0,
     #vertex.frame.color = NA,
     vertex.frame.width = 0.2,
     layout = l_agnostic_gs,
     edge_curved =0.2,
     #main = 'Agnostic',
     margin = -0.04
)
dev.off()


####################################################################################
### Geneset enrichment on all genes

## organize the data:  mechanistic

# get the gene names from the k-tsp classifiers 
mech_forEnr <- All_74 %>%
  as.data.frame() %>%
  dplyr::select(pairs_Mech, gene1_Mech, gene2_Mech) %>%
  group_by(pairs_Mech) %>%
  separate_rows(c(pairs_Mech, gene1_Mech, gene2_Mech), sep = ',') %>%
  dplyr::rename(gene1=gene1_Mech, gene2=gene2_Mech)

# just for the frequency
mech_forEnr_freq <- sum_result_mech %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  dplyr::select(-tmp, -feature) %>%
  relocate(rep_rows, .after = gene2)

mech_forEnr_freq <- data.frame(rep_rows = mech_forEnr_freq$rep_rows, 
                                     pairs_Mech = rownames(mech_forEnr_freq),
                                     row.names = rownames(mech_forEnr_freq))

# merge this with that
mech_forEnr_unique <- mech_forEnr[!duplicated(mech_forEnr$pairs_Mech), ]
mech_forEnr_clean <- merge(x = mech_forEnr_unique, y = mech_forEnr_freq, 
                                        by = 'pairs_Mech', suffixes = colnames(mech_forEnr_unique))

mech_forEnr_clean <- mech_forEnr_clean[order(mech_forEnr_clean$rep_rows, decreasing = T), ] 
mech_forEnr_clean$pairs_Mech <- NULL

#############
## organize the data:  agnostic

# get the gene names from the k-tsp classifiers 
agnostic_forEnr <- All_74 %>%
  as.data.frame() %>%
  dplyr::select(pairs_agnostic, gene1_agnostic, gene2_agnostic) %>%
  group_by(pairs_agnostic) %>%
  separate_rows(c(pairs_agnostic, gene1_agnostic, gene2_agnostic), sep = ',') %>%
  dplyr::rename(gene1=gene1_agnostic, gene2=gene2_agnostic)

# just for the frequency
agnostic_forEnr_freq <- sum_result_agnostic %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  dplyr::select(-tmp, -feature) %>%
  relocate(rep_rows, .after = gene2)

agnostic_forEnr_freq <- data.frame(rep_rows = agnostic_forEnr_freq$rep_rows, 
                               pairs_agnostic = rownames(agnostic_forEnr_freq),
                               row.names = rownames(agnostic_forEnr_freq))

# merge this with that
agnostic_forEnr_unique <- agnostic_forEnr[!duplicated(agnostic_forEnr$pairs_agnostic), ]
agnostic_forEnr_clean <- merge(x = agnostic_forEnr_unique, y = agnostic_forEnr_freq, 
                           by = 'pairs_agnostic', suffixes = colnames(agnostic_forEnr_unique))

agnostic_forEnr_clean <- agnostic_forEnr_clean[order(agnostic_forEnr_clean$rep_rows, decreasing = T), ] 
agnostic_forEnr_clean$pairs_agnostic <- NULL


################################
## cleaning the mess
# if a gene is repeated in gene1 and gene2, keep the order in which it is more frequent 

# function to check how many times "a" (a length 1 atomic vector) occurs in "b":
f <- function(a, b) {
  a <- as.character(a)
  
  # make a lookup table a.k.a dictionary of values in b:
  b_freq <- table(b, useNA = "always")
  
  # if a is in b, return it's frequency:
  if (a %in% names(b_freq)) {
    return(b_freq[a])
  }
  
  # else (ie. a is not in b) return 0:
  return(0)
}

# vectorise that, enabling intake of any length of "a":
ff <- function(a, b) {
  purrr::map_dbl(.x = a, .f = f, b = b)
}

################
## clean mechanistic

# detect and count the flips
mech_forEnr_clean_legacy <- mech_forEnr_clean %>% 
  dplyr::select(-rep_rows) %>%
  group_by(gene1) %>%
  mutate(gene1_freq_as_gene1 = length(gene1 %in% gene1)) %>%
  ungroup() %>%
  mutate(
    gene1_freq_as_gene2 = ff(gene1, gene2)
  ) %>%
  group_by(gene2) %>%
  mutate(gene2_freq_as_gene2 = length(gene2 %in% gene2)) %>%
  ungroup() %>%
  mutate(
    gene2_freq_as_gene1 = ff(gene2, gene1)
  ) 

mech_forEnr_clean_fil <- mech_forEnr_clean %>% 
  dplyr::select(-rep_rows) %>%
  group_by(gene1) %>%
  mutate(gene1_freq_as_gene1 = length(gene1 %in% gene1)) %>%
  ungroup() %>%
  mutate(
    gene1_freq_as_gene2 = ff(gene1, gene2)
  ) %>%
  filter(!(gene1_freq_as_gene2 > gene1_freq_as_gene1)) %>%
  filter(!(gene1_freq_as_gene2 == gene1_freq_as_gene1)) %>%
  group_by(gene2) %>%
  mutate(gene2_freq_as_gene2 = length(gene2 %in% gene2)) %>%
  ungroup() %>%
  mutate(
    gene2_freq_as_gene1 = ff(gene2, gene1)
  ) %>%
  filter(!(gene2_freq_as_gene1 > gene2_freq_as_gene2)) %>%
  filter(!(gene2_freq_as_gene1 == gene2_freq_as_gene2))
  
  
  
## get the genes
mech_gene1 <- unique(mech_forEnr_clean_fil$gene1)
mech_gene2 <- unique(mech_forEnr_clean_fil$gene2)

# check
summary(mech_gene1 %in% mech_gene2)
summary(mech_gene2 %in% mech_gene1)

################
## clean agnostic

# detect and count the flips
agnostic_forEnr_clean_legacy <- agnostic_forEnr_clean %>% 
  dplyr::select(-rep_rows) %>%
  group_by(gene1) %>%
  mutate(gene1_freq_as_gene1 = length(gene1 %in% gene1)) %>%
  ungroup() %>%
  mutate(
    gene1_freq_as_gene2 = ff(gene1, gene2)
  ) %>%
  group_by(gene2) %>%
  mutate(gene2_freq_as_gene2 = length(gene2 %in% gene2)) %>%
  ungroup() %>%
  mutate(
    gene2_freq_as_gene1 = ff(gene2, gene1)
  ) 

agnostic_forEnr_clean_fil <- agnostic_forEnr_clean %>% 
  dplyr::select(-rep_rows) %>%
  group_by(gene1) %>%
  mutate(gene1_freq_as_gene1 = length(gene1 %in% gene1)) %>%
  ungroup() %>%
  mutate(
    gene1_freq_as_gene2 = ff(gene1, gene2)
  ) %>%
  filter(!(gene1_freq_as_gene2 > gene1_freq_as_gene1)) %>%
  filter(!(gene1_freq_as_gene2 == gene1_freq_as_gene1)) %>%
  group_by(gene2) %>%
  mutate(gene2_freq_as_gene2 = length(gene2 %in% gene2)) %>%
  ungroup() %>%
  mutate(
    gene2_freq_as_gene1 = ff(gene2, gene1)
  ) %>%
  filter(!(gene2_freq_as_gene1 > gene2_freq_as_gene2)) %>%
  filter(!(gene2_freq_as_gene1 == gene2_freq_as_gene2))


# get the top n pairs where n is the number of pairs in mechanistic (After cleaning)
agnostic_forEnr_clean_fil <- agnostic_forEnr_clean_fil[c(1:nrow(mech_forEnr_clean_fil)), ]

## get the genes
agnostic_gene1 <- unique(agnostic_forEnr_clean_fil$gene1)
agnostic_gene2 <- unique(agnostic_forEnr_clean_fil$gene2)

# check
summary(agnostic_gene1 %in% agnostic_gene2)
summary(agnostic_gene2 %in% agnostic_gene1)

################
## enrichment

library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2021", 'MSigDB_Oncogenic_Signatures')

###############
## mechanistic

# gobp
Enriched_mech_gene1 <- enrichr(genes = mech_gene1, databases = dbs)
Enriched_mech_gene2 <- enrichr(genes = mech_gene2, databases = dbs)

Enriched_mech_gene1_gobp <- Enriched_mech_gene1["GO_Biological_Process_2021"]$GO_Biological_Process_2021
Enriched_mech_gene1_gobp <- Enriched_mech_gene1_gobp[Enriched_mech_gene1_gobp$Adjusted.P.value <= 0.05, ]

Enriched_mech_gene2_gobp <- Enriched_mech_gene2["GO_Biological_Process_2021"]$GO_Biological_Process_2021
Enriched_mech_gene2_gobp <- Enriched_mech_gene2_gobp[Enriched_mech_gene2_gobp$Adjusted.P.value <= 0.05, ]

# fix the names and add variable for gene count and ratio
Enriched_mech_gene1_gobp$Term <- gsub(' \\(.+', '', Enriched_mech_gene1_gobp$Term)
Enriched_mech_gene2_gobp$Term <- gsub(' \\(.+', '', Enriched_mech_gene2_gobp$Term)

Enriched_mech_gene1_gobp <- Enriched_mech_gene1_gobp %>%
  mutate(Genes = strsplit(Genes, ';')) %>%
  mutate(gene_counts = as.numeric(lengths(Genes))) %>%
  mutate(gs_length = as.numeric(gsub('.+/', '',  Overlap))) %>%
  mutate(ratio = gene_counts/gs_length) %>%
  #arrange(gene_counts) %>%
  mutate(Term = factor(Term, levels = Term)) %>%
  mutate(gene_names = gsub('c', '', Genes)) %>%
  mutate(gene_names = gsub('[[:punct:]]', '', gene_names))
  
Enriched_mech_gene1_gobp$gene_names

Enriched_mech_gene2_gobp <- Enriched_mech_gene2_gobp %>%
  mutate(Genes = strsplit(Genes, ';')) %>%
  mutate(gene_counts = as.numeric(lengths(Genes))) %>%
  mutate(gs_length = as.numeric(gsub('.+/', '',  Overlap))) %>%
  mutate(ratio = gene_counts/gs_length) %>%
  #arrange(gene_counts) %>%
  mutate(Term = factor(Term, levels = Term)) %>%
  mutate(gene_names = gsub('c', '', Genes)) %>%
  mutate(gene_names = gsub('[[:punct:]]', '', gene_names))

## plot
#plot_mech_gene1 <- plotEnrich(Enriched_mech_gene1_gobp, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
#plot_mech_gene2 <- plotEnrich(Enriched_mech_gene2_gobp, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

# plot_mech_gene1 <- ggplot(Enriched_mech_gene1_gobp, aes(x = gene_counts, y = Term, fill = Adjusted.P.value)) +
#   geom_bar(stat = "identity") + 
#   theme_bw() +
#   scale_fill_continuous(low = "blue", high = "red") + 
#   guides(fill = guide_colorbar(title = "adjusted p-value", reverse = FALSE)) +
#   xlab("Gene count") + 
#   ylab("gene sets")
#   
# ggplot(Enriched_mech_gene1_gobp, aes(x = gene_counts, y = Term, fill = Adjusted.P.value, size = gene_counts, label = gene_names)) +
#   #geom_dotplot(binwidth = 1, binpositions = "all", position = "dodge", binaxis = "y") + 
#   geom_point(shape = 21, alpha=.5) +
#   theme_bw() +
#   scale_fill_continuous(low = "blue", high = "red") + 
#   guides(fill = guide_colorbar(title = "adjusted p-value", reverse = FALSE)) +
#   geom_text(show.legend = F, size = 2, nudge_x = 0, nudge_y = 0) +
#   xlab("Gene count") + 
#   ylab("gene sets")  
#############
## agnostic

# gobp
Enriched_agnostic_gene1 <- enrichr(genes = agnostic_gene1, databases = dbs)
Enriched_agnostic_gene2 <- enrichr(genes = agnostic_gene2, databases = dbs)

Enriched_agnostic_gene1_gobp <- Enriched_agnostic_gene1["GO_Biological_Process_2021"]$GO_Biological_Process_2021
Enriched_agnostic_gene1_gobp <- Enriched_agnostic_gene1_gobp[Enriched_agnostic_gene1_gobp$Adjusted.P.value <= 0.05, ]

Enriched_agnostic_gene2_gobp <- Enriched_agnostic_gene2["GO_Biological_Process_2021"]$GO_Biological_Process_2021
Enriched_agnostic_gene2_gobp <- Enriched_agnostic_gene2_gobp[Enriched_agnostic_gene2_gobp$Adjusted.P.value <= 0.05, ]

# fix the names
Enriched_agnostic_gene1_gobp$Term <- gsub(' \\(.+', '', Enriched_agnostic_gene1_gobp$Term)
Enriched_agnostic_gene2_gobp$Term <- gsub(' \\(.+', '', Enriched_agnostic_gene2_gobp$Term)


## plot
plot_agnostic_gene1 <- plotEnrich(Enriched_agnostic_gene1_gobp, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value", title = NULL)
plot_agnostic_gene2 <- plotEnrich(Enriched_agnostic_gene2_gobp, showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value", title = NULL)


#########
# all in one plot
# tiff(filename = "/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/bladder_enrichment.tiff", width = 3000, height = 2000, res = 300)
# ((plot_agnostic_gene1 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) | 
#     (plot_agnostic_gene1 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) / 
#     (plot_mech_gene1 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) |
#     (plot_mech_gene2 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10)))
# ) +
#   #plot_layout(widths = c(0.4, 1)) + 
#   plot_annotation(
#     #title = 'gene s',
#     tag_levels = c('A'),
#     theme = theme(plot.title = element_text(size = 20, face = "plain", hjust=0.5))
#   )
# dev.off()


############
## save

library(openxlsx)
list_of_tables <- list("agnostic_gene2" = Enriched_agnostic_gene2_gobp, "mechanistic_gene1" = Enriched_mech_gene1_gobp, "mechanistic_gene2" = Enriched_mech_gene2_gobp)
write.xlsx(list_of_tables, file = "./objs/bladder_enrichment.xlsx")

#####################################################################################
## venn diagram
library(VennDiagram)
# genes

mech_venn_lists <- apply(mech_indvGns_good_clean,1,as.list)
mech_venn_lists <- lapply(mech_venn_lists, function(x){
  x <- c(x[['gene1']], x[['gene2']])
  #name(x) <- paste(x[['gene1']], x[['gene2']], sep = '-')
  x
})

agnostic_venn_lists <- apply(agnostic_indvGns_good_clean[c(1:93), ],1,as.list)
agnostic_venn_lists <- lapply(agnostic_venn_lists, function(x){
  x <- c(x[['gene1']], x[['gene2']])
  #name(x) <- paste(x[['gene1']], x[['gene2']], sep = '-')
  x
})

mech_venn_lists <- unlist(mech_venn_lists, use.names = F)
agnostic_venn_lists <- unlist(agnostic_venn_lists, use.names = F)

list_venn <- list(agnostic_venn_lists, mech_venn_lists)
names(list_venn) <- c('Agnostic', 'Mechanistic')

###
## The venn diagram
# Prepare a palette of 5 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")[c(1,2)]

genes_venn_plot <- venn.diagram(
  x = list_venn,
  filename = NULL,
  #output=TRUE,
  
  # Output features
  imagetype="png" ,
  #height = 3000 , 
  #width = 2000 , 
  resolution = 600,
  disable.logging = T,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#F8766D", '#00BFC4'),
  # Numbers 
  cex=c(1,1,1),
  #fontface = "bold",
  fontfamily = "sans",

  ext.text = F,
  scaled = F, 
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  #cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #cat.default.pos = "text",
  #cat.pos = c(177, 177)
  margin = 0.05, 
  #cat.col = c("#66C2A5", '#FC8D62')
  #main = "Common genes"
)

# have a look at the default plot
grid.newpage()
grid.draw(genes_venn_plot)

# have a look at the names in the plot object v
lapply(genes_venn_plot,  names)
# We are interested in the labels
lapply(genes_venn_plot, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in agnostic only
#mech_venn_plot[[5]]$label  <- paste(setdiff(list_venn$Agnostic, list_venn$Mechanistic), collapse="\n")  
# in mech only
#mech_venn_plot[[6]]$label <- paste(setdiff(list_venn$Mechanistic, list_venn$Agnostic)  , collapse="\n")  
# intesection
genes_venn_plot[[7]]$label <- paste(intersect(list_venn$Agnostic, list_venn$Mechanistic), collapse="\n")  

# plot  
grid.newpage()
grid.draw(genes_venn_plot)

# save
ggsave(filename="/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/venn/Venn_genes_Bladder.tiff", plot=genes_venn_plot, device = 'tiff', width = 2500, height = 2500, units = 'px', dpi = 500)

####################
## genesets
gs_venn_lists <- list(Agnostic = colnames(hallmarks_gs_long_agnostic_fil), Mechanistic = colnames(hallmarks_gs_long_mech_fil))

## The venn diagram
# Prepare a palette of 5 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")[c(1,2)]

gs_venn_plot <- venn.diagram(
  x = gs_venn_lists,
  filename = NULL,
  #output=TRUE,
  
  # Output features
  imagetype="png" ,
  #height = 3000 , 
  #width = 2000 , 
  resolution = 600,
  disable.logging = T,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#F8766D", '#00BFC4'),
  # Numbers 
  cex=c(0.7,0.7,0.7),
  #fontface = "bold",
  fontfamily = "sans",
  
  ext.text = F,
  scaled = F, 
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  #cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #cat.pos = c(177, 177)
  margin = 0.05, 
  #cat.col = c("#66C2A5", '#FC8D62')
  #main = "Common genes"
)

# have a look at the default plot
grid.newpage()
grid.draw(gs_venn_plot)

# have a look at the names in the plot object v
lapply(gs_venn_plot,  names)
# We are interested in the labels
lapply(gs_venn_plot, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in agnostic only
gs_venn_plot[[5]]$label  <- paste(setdiff(gs_venn_lists$Agnostic, gs_venn_lists$Mechanistic), collapse="\n")  
# in mech only
gs_venn_plot[[6]]$label <- paste(setdiff(gs_venn_lists$Mechanistic, gs_venn_lists$Agnostic)  , collapse="\n")  
# intesection
gs_venn_plot[[7]]$label <- paste(intersect(gs_venn_lists$Agnostic, gs_venn_lists$Mechanistic), collapse="\n")  

# plot  
grid.newpage()
grid.draw(gs_venn_plot)

# save
ggsave(filename="/Users/mohamedomar/Library/CloudStorage/Box-Box/MechPaper/iScience/manuscript/main_figures/venn/Venn_gs_Bladder.tiff", plot=gs_venn_plot, device = 'tiff', width = 4000, height = 4000, units = 'px', dpi = 500)



###########################################
### Plot genes in the training set
## Which TSPs

mech_features <- gsub('-', ',', sum_result_mech_top)

i <- 1:nrow(sum_result_mech_top)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(str_split(x[i], '-', simplify = F))), 
              x=sum_result_mech_top$feature)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTrainMat, g=usedTrainGroup)


names(dfTspTrain) <- sum_result_mech_top$feature

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspTrain)) names(dfTspTrain[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=usedTrainMat, g=usedTrainGroup)

## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
##################################
# ## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=7)
# )
# 
# NofPairs_DistrHist_74 <- ggplot(ModelCompare_NofPairs_74, aes(NofPairs, fill = modelType)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(1, 35)) +
#   labs(title="Number of output pairs distribution of the mechanistic and agnostic (top 74 DEGs) K-TSPs models") + My_Theme


##############################################################
### Work with boot object 100  
All_100 <- bootobject_100$t
colnames(All_100) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_100 <- All_100[,"Diff_Agnostic"]
range(Diff_Agnostic_100)
quantile(Diff_Agnostic_100, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_100 <- All_100[,"Diff_Mechanistic"]
range(Diff_Mechanistic_100)
quantile(Diff_Mechanistic_100, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_100 <- data.frame(AUC = All_100[, "AUC_Train_Mech"])
AgnosticAUCTrain_100 <- data.frame(AUC = All_100[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_100$modelType <- "Mechanistic"
AgnosticAUCTrain_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain_100, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_100 <- data.frame(AUC = All_100[, "AUC_Test_Mech"])
AgnosticAUCTest_100 <- data.frame(AUC = All_100[, "AUC_Test_Agnostic"])

MechanisticAUCTest_100$modelType <- "Mechanistic"
AgnosticAUCTest_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest_100, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Genes"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Genes"


## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Mech"])
Agnostic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_100$modelType <- "Mechanistic"
Agnostic_NofPairs_100$modelType <- "Agnostic"

ModelCompare_NofPairs_100 <- rbind(Mechanistic_NofPairs_100, Agnostic_NofPairs_100)

ModelCompare_NofPairs_100$NofFeatAgn <- "100 Genes"

# #######
# ## Plots
# NofPairs_DistrHist_100 <- ggplot(ModelCompare_NofPairs_100, aes(NofPairs, fill = modelType)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(1, 50)) +
#   labs(title="Number of output pairs distribution of the mechanistic and agnostic (top 100 DEGs) K-TSPs models") + My_Theme

###############################
##############################################################
### Work with boot object 200  
All_200 <- bootobject_200$t
colnames(All_200) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_200 <- All_200[,"Diff_Agnostic"]
range(Diff_Agnostic_200)
quantile(Diff_Agnostic_200, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_200 <- All_200[,"Diff_Mechanistic"]
range(Diff_Mechanistic_200)
quantile(Diff_Mechanistic_200, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_200 <- data.frame(AUC = All_200[, "AUC_Train_Mech"])
AgnosticAUCTrain_200 <- data.frame(AUC = All_200[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_200$modelType <- "Mechanistic"
AgnosticAUCTrain_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_200 <- rbind(MechanisticAUCTrain_200, AgnosticAUCTrain_200)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_200 <- data.frame(AUC = All_200[, "AUC_Test_Mech"])
AgnosticAUCTest_200 <- data.frame(AUC = All_200[, "AUC_Test_Agnostic"])

MechanisticAUCTest_200$modelType <- "Mechanistic"
AgnosticAUCTest_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_200 <- rbind(MechanisticAUCTest_200, AgnosticAUCTest_200)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_200$data_type <- "Training"
ModelCompareAUCTest_200$data_type <- "Testing"

ModelCompareAUCTrain_200$NofFeatAgn <- "200_Genes"
ModelCompareAUCTest_200$NofFeatAgn <- "200_Genes"

## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Mech"])
Agnostic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_200$modelType <- "Mechanistic"
Agnostic_NofPairs_200$modelType <- "Agnostic"

ModelCompare_NofPairs_200 <- rbind(Mechanistic_NofPairs_200, Agnostic_NofPairs_200)

ModelCompare_NofPairs_200$NofFeatAgn <- "200 Genes"

# ##################################
# ## Plots

# NofPairs_DistrHist_200 <- ggplot(ModelCompare_NofPairs_200, aes(NofPairs, fill = modelType)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(1, 50)) +
#   labs(title="Number of output pairs distribution of the mechanistic and agnostic (top 200 DEGs) K-TSPs models") + My_Theme

###############################
##############################################################
### Work with boot object 500  
All_500 <- bootobject_500$t
colnames(All_500) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_500 <- All_500[,"Diff_Agnostic"]
range(Diff_Agnostic_500)
quantile(Diff_Agnostic_500, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_500 <- All_500[,"Diff_Mechanistic"]
range(Diff_Mechanistic_500)
quantile(Diff_Mechanistic_500, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_500 <- data.frame(AUC = All_500[, "AUC_Train_Mech"])
AgnosticAUCTrain_500 <- data.frame(AUC = All_500[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_500$modelType <- "Mechanistic"
AgnosticAUCTrain_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_500 <- rbind(MechanisticAUCTrain_500, AgnosticAUCTrain_500)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_500 <- data.frame(AUC = All_500[, "AUC_Test_Mech"])
AgnosticAUCTest_500 <- data.frame(AUC = All_500[, "AUC_Test_Agnostic"])

MechanisticAUCTest_500$modelType <- "Mechanistic"
AgnosticAUCTest_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_500 <- rbind(MechanisticAUCTest_500, AgnosticAUCTest_500)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_500$data_type <- "Training"
ModelCompareAUCTest_500$data_type <- "Testing"

ModelCompareAUCTrain_500$NofFeatAgn <- "500_Genes"
ModelCompareAUCTest_500$NofFeatAgn <- "500_Genes"

## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Mech"])
Agnostic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_500$modelType <- "Mechanistic"
Agnostic_NofPairs_500$modelType <- "Agnostic"

ModelCompare_NofPairs_500 <- rbind(Mechanistic_NofPairs_500, Agnostic_NofPairs_500)

ModelCompare_NofPairs_500$NofFeatAgn <- "500 Genes"

########
# plot
# NofPairs_DistrHist_500 <- ggplot(ModelCompare_NofPairs_500, aes(NofPairs, fill = modelType)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(1, 50)) +
#   labs(title="Number of output pairs distribution of the mechanistic and agnostic (top 500 DEGs) K-TSPs models") + My_Theme
# 

####################################################################################
# combine the number of pairs distribution

ModelCompare_KTSP_Npairs <- rbind(ModelCompare_NofPairs_74,
                                      ModelCompare_NofPairs_100,
                                      ModelCompare_NofPairs_200,
                                      ModelCompare_NofPairs_500
)

save(ModelCompare_KTSP_Npairs, file = "./Objs/KTSP/ModelCompare_KTSP_Npairs.rda")

####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_KTSP_DiffNoFeat <- rbind(ModelCompareAUCTrain_74,
                                    ModelCompareAUCTest_74,
                                    ModelCompareAUCTrain_100,
                                    ModelCompareAUCTest_100,
                                    ModelCompareAUCTrain_200,
                                    ModelCompareAUCTest_200,
                                    ModelCompareAUCTrain_500,
                                    ModelCompareAUCTest_500
)

save(ModelCompare_KTSP_DiffNoFeat, file = "./Objs/KTSP/ModelCompare_KTSP_DiffNoFeat.rda")

####################################################################################
####################################################################################
####################################################################################


## Plot the distributions of the N of Pairs from both methods
# Mechanistic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Mech"])
# Agnostic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Agnostic"])
# 
# Mechanistic_NofPairs_500$modelType <- "Mech"
# Agnostic_NofPairs_500$modelType <- "Agnostic"
# 
# ModelCompare_NofPairs_500 <- rbind(Mechanistic_NofPairs_500, Agnostic_NofPairs_500)
# 
# ##################################
# ## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=7)
# )
# 
# DiffHist_Agnostic_500 <- ggplot(as.data.frame(Diff_Agnostic_500), aes(Diff_Agnostic_500, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.4)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech_500 <- ggplot(as.data.frame(Diff_Mechanistic_500), aes(Diff_Mechanistic_500, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.4)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic (top 500 DEGs) and mechanistic KTSP models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic (top 500 DEGs) and mechanistic KTSP models in the testing data") + My_Theme
# 
# NofPairs_DistrHist_500 <- ggplot(ModelCompare_NofPairs_500, aes(NofPairs, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(1, 50)) +
#   labs(title="N of pairs distribution of the agnostic and mechanistic KTSP models") + My_Theme
# 
# png("./Figs/KTSP/KTSP_BS_AUC_500.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) | (DiffHist_Agnostic_500 / DiffHist_Mech_500 | NofPairs_DistrHist_500) + plot_layout(widths = c(1, 1.5)) 
# dev.off()
# 
# 
# ###########################################
# ###########################################
# My_Theme2 = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=9)
# )
# 
# #############################
# AUC_Train_DistrHist_50 <- ggplot(ModelCompareAUCTrain_50, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 50 DEGs) vs mechanistic KTSP models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_50 <- ggplot(ModelCompareAUCTest_50, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 50 DEGs) vs mechanistic KTSP models in the testing data") + My_Theme2
# #############################
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 100 DEGs) vs mechanistic KTSP models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 100 DEGs) vs mechanistic KTSP models in the testing data") + My_Theme2
# ##############################
# AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUCTrain_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 200 DEGs) vs mechanistic KTSP models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUCTest_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 200 DEGs) vs mechanistic KTSP models in the testing data") + My_Theme2
# ###############################
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 500 DEGs) vs mechanistic KTSP models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 500 DEGs) vs mechanistic KTSP models in the testing data") + My_Theme2
# 
# 
# # All in one figure
# png("./Figs/BootStrap_Diff_No_features/KTSP_BS_AUC_DiffFeatures.png", width = 3000, height = 1500, res = 150)
# (AUC_Train_DistrHist_50 / AUC_Test_DistrHist_50) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(1, 1)) 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
