library(tidyverse)
library(rsample)
library(yardstick)
library(glmnet)
library(pROC)

set.seed(42)

load("data/biomarker-clean.RData")
partitions <- initial_split(biomarker_clean, prop = 0.8, strata = group)
train <- training(partitions)
test  <- testing(partitions)

fit <- glm(class ~ ., 
           data = train, 
           family = binomial(link = "logit"))

