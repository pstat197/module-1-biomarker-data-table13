library(tidyverse)
library(rsample)
library(yardstick)
library(ranger)
library(glmnet)
library(pROC)

set.seed(42)

# ---- load data (choose one path) ----
load("data/biomarker-clean.RData")

split <- initial_split(biomarker_clean, prop = 0.7, strata = group)   
train <- training(split)
test  <- testing(split)

train$group <- as.factor(train$group)
test$group  <- as.factor(test$group)

feature_cols <- setdiff(names(biomarker_clean), "group")
X_train <- train[, feature_cols, drop = FALSE] |> mutate(across(everything(), as.numeric))
y_train <- train$group
X_test  <- test[,  feature_cols, drop = FALSE] |> mutate(across(everything(), as.numeric))
y_test  <- test$group

