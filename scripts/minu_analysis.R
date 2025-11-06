library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

current_path <- getwd()
# ROOT <- dirname(current_path)
# load(file.path(ROOT, "data", "biomarker-clean.RData")) # loads data as biomarker_clean
load(file.path(current_path, "data", "biomarker-clean.RData"))

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE))) %>% 
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# predict on test data
predictions <- testing(biomarker_split) %>%
  mutate(pred = predict(fit, ., type = "response"),
         pred_class = factor(ifelse(pred > 0.5, TRUE, FALSE),
                             levels = c(FALSE, TRUE)))

# evaluate metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy)

bind_rows(
  class_metrics(predictions,
                truth = class,
                estimate = pred_class,
                event_level = "second"),
  roc_auc(predictions,
          truth = class,
          pred,
          event_level = "second")
)

# sensitivity: 0.812
# specificity: 0.867
# accuracy: 0.839
# roc_auc: 0.946

# Instead of selecting for just 10 proteins, we selected for 20. This improved the sensitivity by 6%, 
# specificity by 13%, accuracy by 9%, and AUC by 4%.

# QUESTION 3 PART 3 (fuzzy intersection)
proteins_s1_fuzzy <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

proteins_s2_fuzzy <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

proteins_sstar_fuzzy <- union(proteins_s1_fuzzy, proteins_s2_fuzzy)

biomarker_sstar_fuzzy <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_fuzzy)) %>%
  mutate(class = factor(group == 'ASD', levels = c(FALSE, TRUE))) %>% 
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_fuzzy <- biomarker_sstar_fuzzy %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit_fuzzy <- glm(class ~ ., 
           data = training(biomarker_split_fuzzy), 
           family = 'binomial')

# predict on test data
predictions_fuzzy <- testing(biomarker_split_fuzzy) %>%
  mutate(pred = predict(fit_fuzzy, ., type = "response"),
         pred_class = factor(ifelse(pred > 0.5, TRUE, FALSE),
                             levels = c(FALSE, TRUE)))

# evaluate metrics
class_metrics <- metric_set(sensitivity, specificity, accuracy)

bind_rows(
  class_metrics(predictions_fuzzy,
                truth = class,
                estimate = pred_class,
                event_level = "second"),
  roc_auc(predictions_fuzzy,
          truth = class,
          pred,
          event_level = "second")
)

# sensitivity: 0.562
# specificity: 0.733
# accuracy: 0.645
# AUC: 0.779

# Instead of doing a hard intersection, we looked at a fuzzy intersection of the proteins. This decreased the
#  sensitivity by 20%, specificity by 7%, accuracy by 10%, and AUC by 12%.