library(tidyverse)
library(rsample)
library(yardstick)
library(ranger)
library(glmnet)
library(pROC)
library(scales)

set.seed(42)

# ---- load data (choose one path) ----

load("data/biomarker-clean.RData")

colnames(biomarker_clean) <- make.names(colnames(biomarker_clean), unique = TRUE)


split <- initial_split(biomarker_clean, prop = 0.7, strata = group)   
train <- training(split)
test  <- testing(split)

topN <- 30

train$group <- factor(train$group, levels = c("TD", "ASD"))
test$group  <- factor(test$group,  levels = c("TD", "ASD"))

feature_cols <- setdiff(names(biomarker_clean), "group")

# numeric conversion
train[feature_cols] <- lapply(train[feature_cols], as.numeric)
test[feature_cols]  <- lapply(test[feature_cols],  as.numeric)

# drop features that are all-NA in TRAIN
all_na <- names(which(sapply(train[feature_cols], function(x) all(is.na(x)))))
if (length(all_na)) feature_cols <- setdiff(feature_cols, all_na)

# median impute on TRAIN; apply the same medians to TEST
train_medians <- sapply(train[feature_cols], function(x) median(x, na.rm = TRUE))

# if any column's median is NA (still all-NA), drop it
na_median <- names(train_medians)[is.na(train_medians)]
if (length(na_median)) {
  feature_cols <- setdiff(feature_cols, na_median)
  train_medians <- sapply(train[feature_cols], function(x) median(x, na.rm = TRUE))
}

for (f in feature_cols) {
  train[[f]][is.na(train[[f]])] <- train_medians[[f]]
  test[[f]][is.na(test[[f]])]   <- train_medians[[f]]
}

# drop zero-variance features in TRAIN
is_const <- sapply(train[feature_cols], function(x) {
  s <- sd(x)
  !is.finite(s) || s < 1e-8
})
if (any(is_const)) {
  feature_cols <- feature_cols[!is_const]
}

X_train <- train[, feature_cols, drop = FALSE] |> mutate(across(everything(), as.numeric))
y_train <- train$group
X_test  <- test[,  feature_cols, drop = FALSE] |> mutate(across(everything(), as.numeric))
y_test  <- test$group

# ---- method 1: Univariate AUC ----
uni_auc <- map_dfr(feature_cols, function(f) {
  v <- train[[f]]
  ok <- is.finite(v) & !is.na(v)
  if (sum(ok) < 10 || length(unique(v[ok])) < 3) return(tibble(feature = f, auc = NA_real_))
  out <- tryCatch({
    r <- roc(response = train$group[ok], predictor = v[ok], quiet = TRUE)
    tibble(feature = f, auc = as.numeric(auc(r)))
  }, error = function(e) tibble(feature = f, auc = NA_real_))
  out
})

top_auc <- uni_auc |>
  arrange(desc(auc)) |>
  filter(!is.na(auc)) |>
  slice_head(n = topN) |>
  pull(feature)

# ---- method 2: Random Forest importance ----
rf_fit <- ranger(
  group ~ .,
  data = train[, c("group", feature_cols)],
  probability = TRUE,
  num.trees = 1000,
  importance = "impurity",
  seed = 42
)

imp <- sort(rf_fit$variable.importance, decreasing = TRUE)
top_rf <- names(imp)[seq_len(min(topN, length(imp)))]

# ---- method 3: LASSO logistic ----
x_mat <- as.matrix(train[, feature_cols])
y_bin <- as.numeric(train$group) - 1  # TD=0, ASD=1

cvfit <- cv.glmnet(x_mat, y_bin, family = "binomial", alpha = 1, nfolds = 5)

coef_min <- coef(cvfit, s = "lambda.min")
sel <- rownames(coef_min)[which(coef_min[, 1] != 0)]
sel <- setdiff(sel, "(Intercept)")

if (length(sel) < topN) {
  beta <- as.vector(coef_min[-1, 1]); names(beta) <- colnames(x_mat)
  ranked <- names(sort(abs(beta), decreasing = TRUE))
  sel <- unique(c(sel, ranked))[1:topN]
} else {
  sel <- sel[1:topN]
}
top_lasso <- sel

# ---- fuzzy intersection ----
rankify <- function(all_feats, top_list) {
  r <- match(all_feats, top_list)
  score <- ifelse(is.na(r), 0,
                  rescale(-r, to = c(0,1), from = c(-length(top_list), -1)))
  tibble(feature = all_feats, score = score)
}

all_feats <- feature_cols

fuzzy <- rankify(all_feats, top_auc)   |>
  rename(score_auc = score)     |>
  left_join(rankify(all_feats, top_rf) |>
              rename(score_rf = score), by = "feature") |>
  left_join(rankify(all_feats, top_lasso) |>
              rename(score_lasso = score), by = "feature") |>
  mutate(total = score_auc + score_rf + score_lasso) |>
  arrange(desc(total))

panel_size <- 20
panel_fuzzy <- fuzzy |> slice_head(n = panel_size) |> pull(feature)

# Hard intersection for comparison
panel_hard <- Reduce(intersect, list(top_auc, top_rf, top_lasso))

# ---- Evaluate on test ----
eval_panel <- function(features, label) {
  if (length(features) < 1) return(tibble(panel = label, n = 0, auc = NA, accuracy = NA))
  df_tr <- train[, c("group", features)]
  df_te <- test[,  c("group", features)]
  
  x_tr <- as.matrix(df_tr[, features, drop = FALSE])
  y_tr <- as.numeric(df_tr$group) - 1
  x_te <- as.matrix(df_te[, features, drop = FALSE])
  
  cvfit <- cv.glmnet(x_tr, y_tr, family = "binomial", alpha = 0, nfolds = 5)
  probs <- as.numeric(predict(cvfit, newx = x_te, s = "lambda.min", type = "response"))
  
  preds <- factor(ifelse(probs >= 0.5, "ASD", "TD"), levels = c("TD","ASD"))
  auc_val <- tryCatch(as.numeric(auc(roc(test$group, probs, quiet = TRUE))), error = function(e) NA_real_)
  tibble(panel = label, n = length(features), auc = auc_val, accuracy = mean(preds == test$group))
}

results_tbl <- bind_rows(
  eval_panel(top_auc,   sprintf("Top-%d Univariate AUC", topN)),
  eval_panel(top_rf,    sprintf("Top-%d RandomForest",  topN)),
  eval_panel(top_lasso, sprintf("Top-%d LASSO",         topN)),
  eval_panel(panel_fuzzy, sprintf("Fuzzy union (top %d by score)", panel_size)),
  eval_panel(panel_hard,  sprintf("Hard intersection (|S|=%d)", length(panel_hard)))
)

print(results_tbl)