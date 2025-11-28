## =========================================================
##  二分类诊断模型 · 105 组合 · 多队列外部验证
##  训练: ZZ；验证: 其余 8 个
##  指标: AUC（列名沿用 Cindex）
## =========================================================
rm(list = ls()); options(stringsAsFactors = FALSE)

## ---- Packages ----
need <- c(
  "dplyr","tibble","tidyr","pROC","glmnet","ranger","gbm",
  "e1071","xgboost","mboost","ggplot2","ggsci","ggbreak"
)
for(pk in need){
  if(!requireNamespace(pk, quietly = TRUE)){
    message("请安装: install.packages('", pk, "')")
  }
}
library(dplyr); library(tibble); library(tidyr); library(pROC)
library(glmnet); library(ranger); library(gbm); library(e1071)
library(xgboost); library(mboost); library(ggplot2)

set.seed(91)

## ---- Data (与截图一致：第1列 Group，其余为基因特征) ----
## 你已有 ee；若尚未 load，请取消下一行注释并改路径
load('Output_data/14_gene_ee.rda')
rm(ID,rr_df,need,pk)

mm <- lapply(ee, function(x){
  stopifnot("Group" %in% colnames(x))
  x[,-match("Group", colnames(x))] <- scale(x[,-match("Group", colnames(x)), drop = FALSE])
  x
})

train_id <- "ZZ"
est_data <- mm[[train_id]]
val_list <- mm
pre_var <- setdiff(colnames(est_data), "Group")

## ---- Helpers ----
to01 <- function(v){
  vv <- suppressWarnings(as.numeric(v)); if(any(is.na(vv))) vv <- as.numeric(as.character(v)); vv
}
est_data$Group <- to01(est_data$Group)
val_list <- lapply(val_list, function(d){ d$Group <- to01(d$Group); d })

auc_safe <- function(y, p){
  if (length(unique(y[is.finite(y)])) < 2) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(y, p, levels = c(0,1), direction = "<")),
           error = function(e) NA_real_)
}
auc_from_rslist <- function(rs_list, label = "Group"){
  data.frame(Cindex = sapply(rs_list, function(df){ auc_safe(df[[label]], df$RS) })) %>%
    rownames_to_column("ID")
}
choose_k <- function(p){ max(5, min(30, round(p/3))) }

## -------- Feature selectors --------
sel_enet <- function(X, y, alpha){
  cv <- cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = 10)
  b  <- as.matrix(coef(cv, s = "lambda.min"))[-1, , drop = FALSE]
  nz <- rownames(b)[abs(b[,1]) > 0]
  if(length(nz) >= 2) return(nz)
  k <- choose_k(ncol(X)); rownames(b)[order(abs(b[,1]), decreasing = TRUE)][1:k]
}
sel_lasso <- function(X, y) sel_enet(X, y, 1)
sel_ridge <- function(X, y) sel_enet(X, y, 0)

sel_stepglm <- function(df, direction = c("both","backward","forward")){
  direction <- match.arg(direction)
  form <- reformulate(setdiff(colnames(df), "Group"), "Group")
  fit0 <- glm(form, data = df, family = binomial())
  fit  <- step(fit0, direction = direction, trace = 0)
  keep <- setdiff(names(coef(fit)), "(Intercept)")
  if(length(keep) >= 2) return(keep)
  sm <- summary(fit0)$coefficients[-1,,drop = FALSE]
  rownames(sm)[order(abs(sm[,"t value"]), decreasing = TRUE)][1:choose_k(nrow(sm))]
}
sel_rf <- function(df){
  fit <- ranger(reformulate(setdiff(colnames(df),"Group"), "Group"),
                data = df, probability = TRUE, num.trees = 1000,
                min.node.size = 4, importance = "permutation", seed = 91)
  vi  <- fit$variable.importance
  vi  <- (vi - min(vi)) / (max(vi) - min(vi) + 1e-12)
  rid <- names(vi)[vi > 0.2]
  if(length(rid) >= 2) return(rid)
  names(sort(vi, decreasing = TRUE))[1:choose_k(length(vi))]
}
sel_glmboost <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  df2 <- df
  df2$Group <- factor(df2$Group, levels = c(0,1))  # 关键：factor
  fit0 <- mboost::glmboost(reformulate(Xn, "Group"), data = df2,
                           family = mboost::Binomial(), center = TRUE)
  cv   <- mboost::cvrisk(fit0, folds = mboost::cv(model.weights(fit0), type = "kfold"))
  mboost::mstop(fit0) <- mboost::mstop(cv)
  vi <- mboost::varimp(fit0)
  if (length(vi) == 0) {
    co <- coef(fit0)
    vi <- abs(co[setdiff(names(co), "(Intercept)")])
  }
  names(sort(vi, decreasing = TRUE))[1:choose_k(length(vi))]
}
sel_gbm <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  fit <- gbm(reformulate(Xn,"Group"), data = df, distribution = "bernoulli",
             n.trees = 3000, interaction.depth = 3, n.minobsinnode = 10,
             shrinkage = 0.01, cv.folds = 5, verbose = FALSE)
  best <- gbm.perf(fit, method = "cv", plot.it = FALSE)
  si <- summary.gbm(fit, n.trees = best, plotit = FALSE)
  rid <- si$var[si$rel.inf >= 2]
  if(length(rid) >= 2) return(rid)
  si$var[order(si$rel.inf, decreasing = TRUE)][1:choose_k(nrow(si))]
}
sel_xgb <- function(X, y){
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = y)
  params <- list(objective = "binary:logistic", eval_metric = "auc",
                 max_depth = 3, eta = 0.1, subsample = 0.9, colsample_bytree = 0.9)
  cv <- xgb.cv(params = params, data = dtrain, nrounds = 500, nfold = 5,
               early_stopping_rounds = 20, verbose = 0)
  bst <- xgb.train(params = params, data = dtrain, nrounds = cv$best_iteration, verbose = 0)
  imp <- xgb.importance(model = bst)
  if(!is.null(imp) && nrow(imp) >= 2) return(head(imp$Feature, choose_k(ncol(X))))
  r <- apply(X, 2, function(z) suppressWarnings(cor(z, y)))
  names(sort(abs(r), decreasing = TRUE))[1:choose_k(ncol(X))]
}
sel_nb <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  sc <- sapply(Xn, function(v){
    xv <- df[[v]]
    if(all(is.finite(xv))){
      au <- tryCatch(as.numeric(pROC::auc(df$Group, xv, levels = c(0,1), direction = "<")),
                     error = function(e) NA_real_)
      if(is.na(au)) 0 else abs(au - 0.5)
    } else 0
  })
  names(sort(sc, decreasing = TRUE))[1:choose_k(length(sc))]
}
sel_svm <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  X  <- as.matrix(df[, Xn, drop = FALSE])
  y  <- factor(df$Group, levels = c(0,1))

  w <- NULL
  # 1) 首选 LiblineaR（若可用，最稳）
  if (requireNamespace("LiblineaR", quietly = TRUE)) {
    # type = 0: L2-regularized L2-loss SVM (primal)
    fit <- LiblineaR::LiblineaR(data = X, target = as.numeric(y) - 1,
                                type = 0, cost = 1, bias = TRUE, verbose = 0)
    w <- as.numeric(fit$W[seq_len(ncol(X))])
    names(w) <- Xn
  } else {
    # 2) 退而求其次：e1071 线性 SVM，w = t(coefs) %*% SV
    fit <- e1071::svm(x = X, y = y, kernel = "linear",
                      type = "C-classification", cost = 1, scale = FALSE)
    w_try <- tryCatch(as.numeric(t(fit$coefs) %*% as.matrix(fit$SV)),
                      error = function(e) NULL)
    if (!is.null(w_try) && length(w_try) == ncol(X)) {
      w <- w_try; names(w) <- Xn
    }
  }

  # 3) 兜底：用 |cor(feature, y)| 排序
  if (is.null(w)) {
    sc <- sapply(Xn, function(v) suppressWarnings(abs(cor(X[, v], as.numeric(y)))))
    w <- sc
  }

  names(sort(abs(w), decreasing = TRUE))[1:choose_k(length(Xn))]
}
## -------- Learners (返回预测概率函数) --------
learn_enet <- function(X, y, alpha){
  feat_names <- colnames(X)
  cv <- cv.glmnet(as.matrix(X), y, family = "binomial", alpha = alpha, nfolds = 10)
  function(newdf){
    newx <- as.matrix(newdf[, feat_names, drop = FALSE])
    as.numeric(predict(cv, newx = newx, s = "lambda.min", type = "response"))
  }
}
learn_lasso <- function(X, y) learn_enet(X, y, 1)
learn_ridge <- function(X, y) learn_enet(X, y, 0)

learn_stepglm <- function(df, direction = c("both","backward","forward")){
  direction <- match.arg(direction)
  Xn <- setdiff(colnames(df),"Group")
  fit <- step(glm(reformulate(Xn, "Group"), data = df, family = binomial()),
              direction = direction, trace = 0)
  function(newdf){
    newx <- newdf[, Xn, drop = FALSE]
    as.numeric(predict(fit, newdata = newx, type = "response"))
  }
}

learn_rf <- function(df){
  Xn <- setdiff(colnames(df),"Group")
  fit <- ranger(reformulate(Xn, "Group"),
                data = df, probability = TRUE,
                num.trees = 1000, min.node.size = 4, seed = 91)
  function(newdf){
    newx <- newdf[, Xn, drop = FALSE]
    as.numeric(predict(fit, data = newx)$predictions[,"1"])
  }
}

learn_glmboost <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  df2 <- df
  df2$Group <- factor(df2$Group, levels = c(0,1))  # 关键：factor
  fit0 <- mboost::glmboost(reformulate(Xn, "Group"), data = df2,
                           family = mboost::Binomial(), center = TRUE)
  cv   <- mboost::cvrisk(fit0, folds = mboost::cv(model.weights(fit0), type = "kfold"))
  mboost::mstop(fit0) <- mboost::mstop(cv)
  function(newdf){
    newx <- newdf[, Xn, drop = FALSE]
    as.numeric(predict(fit0, newdata = newx, type = "response"))  # 概率
  }
}

learn_gbm <- function(df){
  Xn <- setdiff(colnames(df),"Group")
  fit <- gbm(reformulate(Xn,"Group"), data = df, distribution = "bernoulli",
             n.trees = 8000, interaction.depth = 3, n.minobsinnode = 10,
             shrinkage = 0.005, cv.folds = 10, verbose = FALSE)
  best <- gbm.perf(fit, method = "cv", plot.it = FALSE)
  function(newdf){
    newx <- newdf[, Xn, drop = FALSE]
    as.numeric(predict(fit, newdata = newx, n.trees = best, type = "response"))
  }
}

learn_xgb <- function(X, y){
  feat_names <- colnames(X)
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = y)
  params <- list(objective = "binary:logistic", eval_metric = "auc",
                 max_depth = 3, eta = 0.1, subsample = 0.9, colsample_bytree = 0.9)
  cv <- xgb.cv(params = params, data = dtrain, nrounds = 1000, nfold = 5,
               early_stopping_rounds = 20, verbose = 0)
  bst <- xgb.train(params = params, data = dtrain, nrounds = cv$best_iteration, verbose = 0)
  function(newdf){
    newx <- xgb.DMatrix(as.matrix(newdf[, feat_names, drop = FALSE]))
    as.numeric(predict(bst, newdata = newx))
  }
}

learn_nb <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  x  <- df[, Xn, drop = FALSE]
  y  <- factor(df$Group, levels = c(0,1))

  ## 简单稳健化：缺失/非有限值用列中位数填，去零方差列
  for (v in Xn) {
    if (any(!is.finite(x[[v]]))) {
      col <- x[[v]]
      col[!is.finite(col)] <- median(col[is.finite(col)], na.rm = TRUE)
      x[[v]] <- col
    }
  }
  vars <- sapply(x, function(z) stats::sd(z, na.rm = TRUE))
  keep <- Xn[which(vars > 0)]
  if (length(keep) < 2) keep <- Xn[seq_len(min(2, length(Xn)))]

  fit <- e1071::naiveBayes(x = x[, keep, drop = FALSE], y = y)
  function(newdf){
    newx <- newdf[, keep, drop = FALSE]
    for (v in keep) {
      if (any(!is.finite(newx[[v]]))) {
        col <- newx[[v]]
        ## 用训练集该列的中位数做填充
        col[!is.finite(col)] <- median(x[, v], na.rm = TRUE)
        newx[[v]] <- col
      }
    }
    pr <- predict(fit, newdata = newx, type = "raw")
    cls1 <- if ("1" %in% colnames(pr)) "1" else colnames(pr)[2]
    as.numeric(pr[, cls1])
  }
}

learn_svm <- function(df){
  Xn <- setdiff(colnames(df),"Group")
  fit <- e1071::svm(x = df[,Xn,drop=FALSE], y = factor(df$Group, levels = c(0,1)),
                    type = "C-classification", kernel = "radial",
                    cost = 1, gamma = 1/length(Xn), probability = TRUE)
  function(newdf){
    newx <- newdf[, Xn, drop = FALSE]
    pr <- attr(predict(fit, newdata = newx, probability = TRUE), "probabilities")
    if(!is.null(pr) && ("1" %in% colnames(pr))) as.numeric(pr[,"1"]) else as.numeric(pr[,2])
  }
}

## -------- 解析/执行 --------
dispatch_selector <- function(name, param, train_df){
  X <- as.matrix(train_df[, setdiff(colnames(train_df),"Group"), drop = FALSE])
  y <- train_df$Group
  switch(name,
         "Enet"       = sel_enet(X, y, alpha = as.numeric(param$alpha)),
         "Lasso"      = sel_lasso(X, y),
         "Ridge"      = sel_ridge(X, y),
         "Stepglm"    = sel_stepglm(train_df, direction = match.arg(param$direction, c("both","backward","forward"))),
         "RF"         = sel_rf(train_df),
         "glmBoost"   = sel_glmboost(train_df),
         "GBM"        = sel_gbm(train_df),
         "XGBoost"    = sel_xgb(X, y),
         "NaiveBayes" = sel_nb(train_df),
         "SVM"        = sel_svm(train_df),     # ← 新增这一行
         stop("未知选择器: ", name)
  )
}
dispatch_learner <- function(name, param, train_df){
  X <- train_df[, setdiff(colnames(train_df),"Group"), drop = FALSE]
  y <- train_df$Group
  switch(name,
         "Enet"       = learn_enet(X, y, alpha = as.numeric(param$alpha)),
         "Lasso"      = learn_lasso(X, y),
         "Ridge"      = learn_ridge(X, y),
         "Stepglm"    = learn_stepglm(train_df, direction = match.arg(param$direction, c("both","backward","forward"))),
         "RF"         = learn_rf(train_df),
         "glmBoost"   = learn_glmboost(train_df),
         "NaiveBayes" = learn_nb(train_df),
         "GBM"        = learn_gbm(train_df),
         "XGBoost"    = learn_xgb(X, y),
         "SVM"        = learn_svm(train_df),
         stop("未知学习器: ", name)
  )
}
parse_one <- function(s){
  s0 <- gsub("\\s+", "", s)
  parts <- strsplit(s0, "\\+")[[1]]
  parse_token <- function(tok){
    m <- regexec("^([A-Za-z]+)(\\[(.*)\\])?$", tok)
    r <- regmatches(tok, m)[[1]]
    nm <- r[2]; pr <- r[4]; par <- list()
    if(!is.na(pr) && nzchar(pr)){
      if(grepl("=", pr)){
        kv <- strsplit(pr, ",")[[1]]
        for(k in kv){ kv2 <- strsplit(k, "=")[[1]]; par[[kv2[1]]] <- kv2[2] }
      }else{
        par$direction <- pr
      }
    }
    list(name = nm, param = par)
  }
  out <- list(learner = parse_token(parts[1]))
  if(length(parts) >= 2) out$selector <- parse_token(parts[2]) else out$selector <- NULL
  out$label <- s0; out
}
run_combos <- function(combos, result_list, pred_store){
  for (cmb in combos){
    pc <- parse_one(cmb)
    ## 选择特征
    feat <- pre_var
    if(!is.null(pc$selector)){
      rid <- dispatch_selector(pc$selector$name, pc$selector$param, est_data[, c("Group", pre_var)])
      feat <- intersect(pre_var, rid)
    }
    if(length(feat) < 2){
      cors <- sapply(pre_var, function(v) suppressWarnings(abs(cor(est_data[[v]], est_data$Group))))
      feat <- names(sort(cors, decreasing = TRUE))[1:2]
    }
    train_df <- est_data[, c("Group", feat), drop = FALSE]
    ## 学习器
    pred_fun <- dispatch_learner(pc$learner$name, pc$learner$param, train_df)
    ## 预测所有队列
    rs <- lapply(val_list, function(x){
      df <- x[, c("Group", feat), drop = FALSE]
      data.frame(Group = df$Group, RS = pred_fun(df), row.names = rownames(df))
    })
    ## AUC
    cc <- auc_from_rslist(rs); cc$Model <- pc$label
    result_list[[cmb]] <- cc
    pred_store[[cmb]] <- rs
    cat("Done ->", pc$label, "\n")
  }
  list(result_list = result_list, pred_store = pred_store)
}

## =========================================================
##  分块执行（和你风格一致，清晰标注）
## =========================================================
result <- list(); pred_store <- list()

#### 1. Lasso 系 ############################################
## 1-1 Lasso
## 1-2 Lasso+Stepglm[both]
## 1-3 Lasso+RF
## 1-4 Lasso+glmBoost
## 1-5 Lasso+NaiveBayes
## 1-6 Lasso+GBM
## 1-7 Lasso+XGBoost
## 1-8 Lasso+SVM
block1 <- c(
  "Lasso",
  "Lasso+Stepglm[both]",
  "Lasso+RF",
  "Lasso+glmBoost",
  "Lasso+NaiveBayes",
  "Lasso+GBM",
  "Lasso+XGBoost",
  "Lasso+SVM"
)
tmp <- run_combos(block1, result, pred_store)
result <- tmp$result_list
pred_store <- tmp$pred_store

#### 2. Enet (alpha=0.1~0.9) 系 #############################
## （每个 alpha 下按你列出的子组合执行）
block2 <- c(
  "Enet[alpha=0.1]",
  "Enet[alpha=0.1]+Stepglm[backward]",
  "Enet[alpha=0.1]+glmBoost",
  "Enet[alpha=0.1]+GBM",
  "Enet[alpha=0.1]+SVM",
  "Enet[alpha=0.1]+XGBoost",
  "Enet[alpha=0.2]",
  "Enet[alpha=0.2]+Stepglm[both]",
  "Enet[alpha=0.2]+SVM",
  "Enet[alpha=0.2]+XGBoost",
  "Enet[alpha=0.2]+glmBoost",
  "Enet[alpha=0.3]",
  "Enet[alpha=0.3]+Stepglm[backward]",
  "Enet[alpha=0.3]+RF",
  "Enet[alpha=0.3]+GBM",
  "Enet[alpha=0.4]",
  "Enet[alpha=0.4]+Stepglm[both]",
  "Enet[alpha=0.4]+RF",
  "Enet[alpha=0.5]",
  "Enet[alpha=0.5]+RF",
  "Enet[alpha=0.6]",
  "Enet[alpha=0.6]+Stepglm[backward]",
  "Enet[alpha=0.6]+Stepglm[forward]",
  "Enet[alpha=0.6]+RF",
  "Enet[alpha=0.7]",
  "Enet[alpha=0.7]+Stepglm[backward]",
  "Enet[alpha=0.7]+glmBoost",
  "Enet[alpha=0.8]",
  "Enet[alpha=0.8]+Stepglm[both]",
  "Enet[alpha=0.8]+glmBoost",
  "Enet[alpha=0.9]",
  "Enet[alpha=0.9]+RF"
)
tmp <- run_combos(block2, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 3. Ridge 与 Stepglm 系 ################################
block3 <- c(
  "Ridge",
  "Ridge+Stepglm[both]",
  "Stepglm[backward]",
  "Stepglm[backward]+RF",
  "Stepglm[backward]+glmBoost",
  "Stepglm[backward]+GBM",
  "Stepglm[backward]+Enet[alpha=0.1]",
  "Stepglm[backward]+Enet[alpha=0.2]",
  "Stepglm[backward]+Enet[alpha=0.3]",
  "Stepglm[backward]+Enet[alpha=0.5]",
  "Stepglm[backward]+Enet[alpha=0.6]",
  "Stepglm[backward]+Lasso",
  "Stepglm[both]+glmBoost",
  "Stepglm[both]+Enet[alpha=0.2]",
  "Stepglm[both]+Enet[alpha=0.3]",
  "Stepglm[both]+Enet[alpha=0.4]",
  "Stepglm[both]+Enet[alpha=0.6]",
  "Stepglm[both]+Ridge",
  "Stepglm[forward]+glmBoost",
  "Stepglm[forward]+Enet[alpha=0.1]",
  "Stepglm[forward]+Enet[alpha=0.3]",
  "Stepglm[forward]+Enet[alpha=0.5]",
  "Stepglm[forward]+Enet[alpha=0.6]",
  "Stepglm[forward]+Enet[alpha=0.8]",
  "Stepglm[forward]+Enet[alpha=0.9]"
)
tmp <- run_combos(block3, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 4. glmBoost 系 ########################################
block4 <- c(
  "glmBoost",
  "glmBoost+Enet[alpha=0.2]",
  "glmBoost+Lasso"
)
tmp <- run_combos(block4, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 5. SVM 系 #############################################
block5 <- c(
  "SVM",
  "SVM+Enet[alpha=0.1]",
  "SVM+Enet[alpha=0.3]",
  "SVM+Enet[alpha=0.4]",
  "SVM+Enet[alpha=0.5]",
  "SVM+Enet[alpha=0.6]",
  "SVM+Enet[alpha=0.8]",
  "SVM+Enet[alpha=0.9]",
  "SVM+NaiveBayes",
  "SVM+glmBoost",
  "SVM+GBM",
  "SVM+Lasso"
)
tmp <- run_combos(block5, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 6. GBM 系 #############################################
block6 <- c(
  "GBM+XGBoost",
  "GBM+glmBoost",
  "GBM+NaiveBayes",
  "GBM+Lasso",
  "GBM",
  "GBM+Enet[alpha=0.1]",
  "GBM+Enet[alpha=0.3]",
  "GBM+Enet[alpha=0.4]",
  "GBM+Enet[alpha=0.5]",
  "GBM+Enet[alpha=0.6]"
)
tmp <- run_combos(block6, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 7. Naive Bayes 系 #####################################
block7 <- c(
  "NaiveBayes",
  "NaiveBayes+Enet[alpha=0.1]",
  "NaiveBayes+Enet[alpha=0.2]",
  "NaiveBayes+Enet[alpha=0.5]",
  "NaiveBayes+Enet[alpha=0.6]",
  "NaiveBayes+Enet[alpha=0.7]",
  "NaiveBayes+Enet[alpha=0.8]",
  "NaiveBayes+Ridge"
)
tmp <- run_combos(block7, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

#### 8. XGBoost 系 #########################################
block8 <- c(
  "XGBoost+GBM",
  "XGBoost+SVM",
  "XGBoost+NaiveBayes",
  "XGBoost+Enet[alpha=0.2]",
  "XGBoost+Enet[alpha=0.3]",
  "XGBoost+Enet[alpha=0.5]",
  "XGBoost+Enet[alpha=0.9]"
)
tmp <- run_combos(block8, result, pred_store); result <- tmp$result_list; pred_store <- tmp$pred_store

## =========================================================
##  汇总与保存（与生存版风格一致）
## =========================================================
result2 <- do.call(rbind, result)



save(result2, pred_store, file = "Output_data/diagnostic_ZZ_105combo.rda")

## 排除训练集，计算外部验证平均 AUC
dd <- result2 %>%
  dplyr::filter(ID != train_id) %>%
  group_by(Model) %>%
  summarise(Cindex = mean(Cindex, na.rm = TRUE)) %>%
  arrange(desc(Cindex))

dd2 <- tidyr::pivot_wider(result2, names_from = "ID", values_from = "Cindex") %>% as.data.frame()
if(ncol(dd2) > 1) dd2[,-1] <- apply(dd2[,-1,drop=FALSE], 2, as.numeric)
aa <- dd2[dd2$Model=='Lasso+GBM',]
save(aa,dd, dd2, file = "Output_data/diagnostic_ZZ_105combo_summary.rda")

# -------------------------------------------------------------------------
load("Output_data/diagnostic_ZZ_105combo_summary.rda")

get_genes_for_model <- function(model_str){
  pc <- parse_one(model_str)
  # 有选择器：按选择器返回
  if (!is.null(pc$selector)){
    rid <- dispatch_selector(pc$selector$name, pc$selector$param,
                             est_data[, c("Group", pre_var), drop = FALSE])
    return(unique(intersect(pre_var, rid)))
  }
  # 无选择器：按学习器类型取“有效特征”
  nm <- pc$learner$name
  df <- est_data[, c("Group", pre_var), drop = FALSE]
  X  <- as.matrix(df[, pre_var, drop = FALSE]); y <- df$Group

  if (nm == "Enet"){
    a  <- as.numeric(pc$learner$param$alpha)
    cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = a, nfolds = 10)
    b  <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1, drop = FALSE]
    return(rownames(b)[abs(b[,1]) > 0])
  } else if (nm == "Lasso"){
    cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
    b  <- as.matrix(coef(cv, s = "lambda.min"))[-1, 1, drop = FALSE]
    return(rownames(b)[abs(b[,1]) > 0])
  } else if (nm == "Ridge"){
    return(pre_var)  # 岭回归不稀疏化
  } else if (nm == "Stepglm"){
    dir <- match.arg(pc$learner$param$direction, c("both","backward","forward"))
    fit <- step(glm(reformulate(pre_var, "Group"), data = df, family = binomial()),
                direction = dir, trace = 0)
    return(setdiff(names(coef(fit)), "(Intercept)"))
  } else if (nm == "glmBoost"){
    df2 <- df; df2$Group <- factor(df2$Group, levels = c(0,1))
    fit0 <- mboost::glmboost(reformulate(pre_var, "Group"), data = df2,
                             family = mboost::Binomial(), center = TRUE)
    cv   <- mboost::cvrisk(fit0, folds = mboost::cv(model.weights(fit0), type = "kfold"))
    mboost::mstop(fit0) <- mboost::mstop(cv)
    vi <- mboost::varimp(fit0)
    if (length(vi) == 0) {
      co <- coef(fit0); vi <- abs(co[setdiff(names(co), "(Intercept)")])
    }
    keep <- names(vi[vi != 0])
    if (length(keep) == 0) keep <- pre_var
    return(keep)
  } else {
    return(pre_var)  # RF/GBM/SVM/XGBoost/NaiveBayes：无显式稀疏
  }
}

## 计算每个模型的基因集合与数量
gene_sets <- lapply(dd$Model, get_genes_for_model)
gene_counts <- vapply(gene_sets, length, 1L)

dd_with_genes <- dd %>%
  mutate(Genes = gene_counts)

dd_with_genes %>% arrange(desc(Cindex)) %>% print(n = 50)


save(dd_with_genes,file = 'Output_data/dd_with_genes.rda')



load('Output_data/diagnostic_ZZ_105combo_summary.rda')
load('Output_data/dd_with_genes.rda')
aa <- dd2[dd2$Model=='Lasso+GBM',]




















