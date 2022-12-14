###################################################
## SVM

rm(list = ls())

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
# library(DMwR)
library(mltools)
library(boot)
library(patchwork)

##################################################
## SVM
# Mechanistic
## Load data
load("Objs/KTSP/TNBC_KTSP_STATs_Mechanistic_NotchAndMyc2.rda")
load("Objs/ChemoDataNew.rda")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Mechanistic)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Mechanistic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Mechanistic) <- make.names(names(Data_train_Mechanistic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 3, scale = 10, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Mech <- varImp(SVM, scale = TRUE)
  Importance_Mech <- Importance_Mech$importance
  Importance_Mech <- Importance_Mech[order(Importance_Mech$Resistant, decreasing = TRUE),]
  Importance_Mech <- Importance_Mech[Importance_Mech$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Mech)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainMechanistic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMechanistic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMechanistic$auc, ROCTestMechanistic$auc, N_ImportanVariables))
}

print("A")
set.seed(333)
bootobjectMech <- boot(data= Data_train_Mechanistic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

###################################################################################
### Agnostic

### 50 Random genes

## Load data
load("Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

Grid <- expand.grid(degree = 1, scale = 1, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  sel.genes = sample(setdiff(colnames(data), "usedTrainGroup"), 50)
  data = data[, sel.genes]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

print("B")
set.seed(333)
bootobjectAgnostic_50 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15)

########################################################################################
### Agnostic

### 100 Random genes

## Load data
load("Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

Grid <- expand.grid(degree = 1, scale = 1, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  sel.genes = sample(setdiff(colnames(data), "usedTrainGroup"), 100)
  data = data[, sel.genes]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

print("C")
set.seed(333)
bootobjectAgnostic_100 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### 200 Random genes

## Load data
load("Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

Grid <- expand.grid(degree = 1, scale = 1, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  sel.genes = sample(setdiff(colnames(data), "usedTrainGroup"), 200)
  data = data[, sel.genes]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

print("D")
set.seed(333)
bootobjectAgnostic_200 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### 500 Random genes

## Load data
load("Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

Grid <- expand.grid(degree = 1, scale = 1, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  sel.genes = sample(setdiff(colnames(data), "usedTrainGroup"), 500)
  data = data[, sel.genes]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

print("E")
set.seed(333)
bootobjectAgnostic_500 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
########################################################################################
## Save all Objects
save(bootobjectMech, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_200, bootobjectAgnostic_500, file= "Objs/SVM/SVMBootObjects_NotchAndMyc_RAND.rda")

## Load
# load("../../Objs/SVM/SVMBootObjects_NotchAndMyc_RAND.rda")
load("Objs/SVM/SVMBootObjects_NotchAndMyc_RAND.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 50 vs mechanistic

AUCs_SVM_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_SVM_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_50 <- AUCs_SVM_Agnostic_50[, "AUC_Train"] - AUCs_SVM_Agnostic_50[, "AUC_Test"]
quantile(DiffAgnostic_50, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_50$modelType <- "Random_genes"

ModelCompareAUCTrain_50 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_50$modelType <- "Random_genes"

ModelCompareAUCTest_50 <- rbind(MechanisticAUCTest, AgnosticAUCTest_50)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_50$data_type <- "Training"
ModelCompareAUCTest_50$data_type <- "Testing"

ModelCompareAUCTrain_50$NofFeatAgn <- "50_Genes"
ModelCompareAUCTest_50$NofFeatAgn <- "50_Genes"

save(ModelCompareAUCTrain_50, ModelCompareAUCTest_50, file = "Objs/SVM/ModelCompare_RAND_AUC_50.rda")

# ###########################################################################3
# ## Save for the main figure
ModelCompare_SVM <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTest_50)
ModelCompare_SVM$algorithm <- "SVM"
save(ModelCompare_SVM, file = "Objs/SVM/ModelCompare_SVM_RAND.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 100 vs mechanistic

AUCs_SVM_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_SVM_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_100 <- AUCs_SVM_Agnostic_100[, "AUC_Train"] - AUCs_SVM_Agnostic_100[, "AUC_Test"]
quantile(DiffAgnostic_100, c(0.025, 0.975))

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_100$modelType <- "Random_genes"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_100$modelType <- "Random_genes"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Genes"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Genes"

save(ModelCompareAUCTrain_100, ModelCompareAUCTest_100, file = "Objs/SVM/ModelCompare_RAND_AUC_100.rda")

###########################################################################3
## Save for the main figure
#ModelCompare_SVM <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTest_100)
#ModelCompare_SVM$algorithm <- "SVM"
#save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

AUCs_SVM_Agnostic_200 <- bootobjectAgnostic_200$t
colnames(AUCs_SVM_Agnostic_200) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_200 <- AUCs_SVM_Agnostic_200[, "AUC_Train"] - AUCs_SVM_Agnostic_200[, "AUC_Test"]
quantile(DiffAgnostic_200, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_200 <- data.frame(AUC = AUCs_SVM_Agnostic_200[, "AUC_Train"])

AgnosticAUCTrain_200$modelType <- "Random_genes"

ModelCompareAUCTrain_200 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_200)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_200 <- data.frame(AUC = AUCs_SVM_Agnostic_200[, "AUC_Test"])

AgnosticAUCTest_200$modelType <- "Random_genes"

ModelCompareAUCTest_200 <- rbind(MechanisticAUCTest, AgnosticAUCTest_200)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_200$data_type <- "Training"
ModelCompareAUCTest_200$data_type <- "Testing"

ModelCompareAUCTrain_200$NofFeatAgn <- "200_Genes"
ModelCompareAUCTest_200$NofFeatAgn <- "200_Genes"

save(ModelCompareAUCTrain_200, ModelCompareAUCTest_200, file = "Objs/SVM/ModelCompare_RAND_AUC_200.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 500 vs mechanistic

AUCs_SVM_Agnostic_500 <- bootobjectAgnostic_500$t
colnames(AUCs_SVM_Agnostic_500) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_500 <- AUCs_SVM_Agnostic_500[, "AUC_Train"] - AUCs_SVM_Agnostic_500[, "AUC_Test"]
quantile(DiffAgnostic_500, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_500 <- data.frame(AUC = AUCs_SVM_Agnostic_500[, "AUC_Train"])

AgnosticAUCTrain_500$modelType <- "Random_genes"

ModelCompareAUCTrain_500 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_500)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_500 <- data.frame(AUC = AUCs_SVM_Agnostic_500[, "AUC_Test"])

AgnosticAUCTest_500$modelType <- "Random_genes"

ModelCompareAUCTest_500 <- rbind(MechanisticAUCTest, AgnosticAUCTest_500)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_500$data_type <- "Training"
ModelCompareAUCTest_500$data_type <- "Testing"

ModelCompareAUCTrain_500$NofFeatAgn <- "500_Genes"
ModelCompareAUCTest_500$NofFeatAgn <- "500_Genes"

save(ModelCompareAUCTrain_500, ModelCompareAUCTest_500, file = "Objs/SVM/ModelCompare_RAND_AUC_500.rda")

###############################
## save all

ModelCompare_SVM_RAND_DiffNoFeat <- rbind(ModelCompareAUCTrain_50,
                                         ModelCompareAUCTest_50,
                                         ModelCompareAUCTrain_100,
                                         ModelCompareAUCTest_100,
                                         ModelCompareAUCTrain_200,
                                         ModelCompareAUCTest_200,
                                         ModelCompareAUCTrain_500,
                                         ModelCompareAUCTest_500
)

save(ModelCompare_SVM_RAND_DiffNoFeat, file = "./Objs/SVM/ModelCompare_SVM_RAND_DiffNoFeat.rda")

