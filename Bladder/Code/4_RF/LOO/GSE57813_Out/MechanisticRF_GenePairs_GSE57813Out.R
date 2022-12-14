########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: RF for predicting bladder cancer progression
# Mechanistic: gene pairs
# Corss validation: GSE57813

###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
#setwd("/Volumes/Macintosh/Research/Projects/Bladder")

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(varSelRF)


## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_GSE57813Out.rda")
load("./Objs/ProgressionDataGood_GSE57813Out.rda")


usedTrainGroup <- trainGroup
usedTestGroup <- testGroup


predictor_data <- t(KTSP_STATs_Train_Mechanistic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)



set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 100, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

set.seed(333)
RF_Mechanistic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 100, mtry = 6, proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
RF_Mechanistic_OnKTSP
#plot(RF_Mechanistic_OnKTSP)


#VAR <- varSelRF(xdata = predictor_data, Class = usedTrainGroup, ntree = 501, ntreeIterat = 20, whole.range = FALSE, keep.forest = TRUE)

#SelectedVars <- VAR$selected.vars

## RandomForest calculates an importance measures for each variable.
# rf_importances <- randomForest::importance(RF_Mechanistic_OnKTSP, scale=FALSE)
# rf_importances <- rf_importances[order(rf_importances[,4], decreasing = TRUE), ]
# ImportantVariables <- rownames(rf_importances)[1:50]
# 
# varImpPlot(RF_Mechanistic_OnKTSP, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
train_pred_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data, type="response")
confusionTrain <- confusionMatrix(train_pred_responses_Mechanistic_OnKTSP, usedTrainGroup, positive = "Progression")
confusionTrain

# Calculate Matthews correlation coefficient
MCC_Mechanistic_Train <- mcc(preds = train_pred_responses_Mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Mechanistic_Train

# ROC curve in training data
train_pred_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, newdata = predictor_data, type = "vote")
ROCTrain <- roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTrain

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Mechanistic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_Test_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusionTest <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "Progression")
confusionTest

# Calculate Matthews correlation coefficient
MCC_Mechanistic_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Mechanistic_Test

ROCTest <- roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Mechanistic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE57813_Out_RF_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(GSE57813_Out_RF_MechPerformance, file = "./Objs/RF/GSE57813_Out_RF_MechPerformance.rda")


############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_RF <- apply(rbind(
  ci(roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], levels = c("NoProgression", "Progression"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], levels = c("NoProgression", "Progression"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP[,2], levels = c("NoProgression", "Progression"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP[,2], levels = c("NoProgression", "Progression"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_RF <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup, levels=c("NoProgression", "Progression")),
  Agnostic.Training.Pairs = train_pred_votes_Agnostic_OnKTSP[,2],
  Mechanistic.Training.Pairs = train_pred_votes_Mechanistic_OnKTSP[,2]
))

### Change Colnames
colnames(datTrn_RF) <- c("Progression", "RF_type", "RF_sum")


### Testing
datTst_RF <- melt(data.frame(
  ## Testing group
  Testing=factor(usedTestGroup, levels=c("NoProgression", "Progression")),
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing.Pairs = RF_predictions_votes_Agnostic_OnKTSP[,2],
  Mechanistic.Testing.Pairs = RF_predictions_votes_Mechanistic_OnKTSP[,2]
))
### Change Colnames
colnames(datTst_RF) <- c("Progression", "RF_type", "RF_sum")

### Combine
dat_RF <- rbind(datTrn_RF, datTst_RF)
dat_RF$Progression <- as.numeric(dat_RF$Progression)-1

### Replace levels
levels(dat_RF$RF_type) <- gsub("\\.", "-", levels(dat_RF$RF_type))
levels(dat_RF$RF_type) <- paste(levels(dat_RF$RF_type), forLegend_RF[c(3,1,4,2)])

###################
### Plot Curve
png("./Figs/RF/CompareAUCggplot_GSE57813Out.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic RF (GSE57813 for validation))"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_RF_GSE57813Out <- ggplot(dat_RF, aes(d=Progression, m=RF_sum, color=RF_type,
                                               linetype = RF_type)) +
  geom_roc(cutoffs.at = seq(0.1,1,0.1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="bold", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_RF_GSE57813Out
### Close device
dev.off()


save(basicplot_RF_GSE57813Out, file = "./Objs/RF/BasicPlot_RF_GSE57813Out.rda")

