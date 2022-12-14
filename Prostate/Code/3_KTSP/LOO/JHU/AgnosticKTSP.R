
rm(list = ls())

####################### 
##Load required packages
library(switchBox)
library(Biobase)
library(limma)
library(pROC)
library(caret)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(plotROC)
library(xtable)
library(mltools)

#####################################################################
### Load data
load("./Objs/LOO/MetastasisData_JHUOut.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes, ]
usedTestMat <- testMat[RGenes, ]

usedTestMat <- t(scale(t(usedTestMat), center = F, scale = TRUE))

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

###########################################################################
### TRAINING using all expressed genes

## Set Feature number and max K
featN <- nrow(usedTrainMat) 
ktsp <- c(3:25) 


###########
### Train a classifier using the default filter function
ktspPredictorUnRes <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = ktsp, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= featN)
ktspPredictorUnRes

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTrainUnRes$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, plot = FALSE, ci = TRUE, print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")
ROCTrain

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionTrain <- confusionMatrix(usedTrainPredictionUnRes, usedTrainGroup, positive = "Mets")
confusionTrain

MCC_Agnostic_Train <- mltools::mcc(pred = usedTrainPredictionUnRes, actuals = usedTrainGroup)
MCC_Agnostic_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Agnostic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

###############################################################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTestUnRes$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes$statistics, plot = F, ci = T, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")
ROCTest

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionTest <- confusionMatrix(usedTestPredictionUnRes, usedTestGroup, positive = "Mets", mode = "everything")
confusionTest

MCC_Agnostic_Test <- mltools::mcc(pred = usedTestPredictionUnRes, actuals = usedTestGroup)
MCC_Agnostic_Test

## Save the performance metrics
TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Agnostic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
JHU_Out_AgnosticPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(JHU_Out_AgnosticPerformance, file = "./Objs/KTSP/JHU_Out_AgnosticPerformance.rda")

