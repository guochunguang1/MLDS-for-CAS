## ==========================================
## 最终模型：Enet(alpha=0.3) + GBM（3基因）
## 训练：ZZ；验证：其余 8 个
## 输出：最终3基因+系数；各队列AUC
## ==========================================

rm(list = ls()); options(stringsAsFactors = FALSE)
set.seed(91)

## ---- Packages ----
need <- c("dplyr","tibble","glmnet","gbm","pROC")
for(pk in need){
  if(!requireNamespace(pk, quietly = TRUE)){
    stop("缺少包：", pk, "  请先 install.packages('", pk, "')")
  }
}
library(dplyr); library(tibble); library(glmnet); library(gbm); library(pROC)

## ---- Helpers ----
to01 <- function(v){
  vv <- suppressWarnings(as.numeric(v))
  if(any(is.na(vv))) vv <- as.numeric(as.character(v))
  vv
}
auc_safe <- function(y, p){
  if(length(unique(y[is.finite(y)])) < 2) return(NA_real_)
  tryCatch(as.numeric(pROC::auc(y, p, levels=c(0,1), direction="<")),
           error = function(e) NA_real_)
}
choose_k <- function(p) max(5, min(30, round(p/3)))

## GBM 选择器（与你之前一致的逻辑）
sel_gbm <- function(df){
  Xn <- setdiff(colnames(df), "Group")
  fit <- gbm(reformulate(Xn, "Group"), data = df, distribution = "bernoulli",
             n.trees = 3000, interaction.depth = 3, n.minobsinnode = 10,
             shrinkage = 0.01, cv.folds = 5, verbose = FALSE)
  best <- gbm.perf(fit, method = "cv", plot.it = FALSE)
  si <- summary.gbm(fit, n.trees = best, plotit = FALSE)
  rid <- si$var[si$rel.inf >= 2]                 # 阈值法
  if(length(rid) >= 2) return(rid)
  si$var[order(si$rel.inf, decreasing = TRUE)][1:choose_k(nrow(si))]  # 兜底Top-K
}

## ---- Data ----

load("Output_data/14_gene_ee.rda")

mm <- lapply(ee, function(x){
  stopifnot("Group" %in% colnames(x))
  x[,-match("Group", colnames(x))] <- scale(x[,-match("Group", colnames(x)), drop = FALSE])
  x$Group <- to01(x$Group)
  x
})
train_id <- "ZZ"
est_data <- mm[[train_id]]
val_list <- mm
pre_var  <- setdiff(colnames(est_data), "Group")

set.seed(91)

## ---- 1) ZZ上用GBM筛特征（与你框架一致）----
gbm_feats <- sel_gbm(est_data)
gbm_feats <- intersect(setdiff(colnames(est_data), "Group"), gbm_feats)
stopifnot(length(gbm_feats) >= 2)

## ---- 2) Enet(alpha=0.3) 用 cv.glmnet + lambda.min 建模（不强制=3）----
X <- as.matrix(est_data[, gbm_feats, drop = FALSE])
y <- est_data$Group

set.seed(91)  # 锁定CV折分
cvfit <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0.3, nfolds = 10)
lam_star <- "lambda.min"

## 系数（注意：若某个基因在 lambda.min 下被压成0，非零数可能<3；这与你前表一致）
b  <- as.matrix(coef(cvfit, s = lam_star))[-1, 1, drop = FALSE]
final_coef <- tibble::tibble(
  Gene = rownames(b),
  Beta = as.numeric(b[,1])
) |>
  dplyr::filter(Beta != 0) |>
  dplyr::arrange(dplyr::desc(abs(Beta))) |>
  dplyr::mutate(Direction = ifelse(Beta > 0, "Risky", "Protective"))

print(final_coef, n = nrow(final_coef))  # 这里就是最终“进入模型”的基因与系数

## ---- 3) 逐队列AUC（与你前表口径一致）----
predict_fun <- function(df){
  newx <- as.matrix(df[, gbm_feats, drop = FALSE])
  as.numeric(predict(cvfit, newx = newx, s = lam_star, type = "response"))
}
rs_list <- lapply(val_list, function(df){
  data.frame(Group = df$Group, RS = predict_fun(df), row.names = rownames(df))
})
auc_by_id <- tibble::tibble(
  ID  = names(rs_list),
  AUC = sapply(rs_list, function(dd) auc_safe(dd$Group, dd$RS))
)
print(auc_by_id, n = nrow(auc_by_id))

cat("\n=== 各队列AUC（含ZZ训练集行）===\n"); print(auc_by_id, n = nrow(auc_by_id))

## ---- 保存 ----
write.csv(final_coef, "Output_data/Enet_a0.3_GBM_final3genes_coef.csv", row.names = FALSE)
write.csv(auc_by_id, "Output_data/Enet_a0.3_GBM_per_cohort_auc.csv", row.names = FALSE)
save(rs_list,lam_star, gbm_feats, final_coef, auc_by_id,
     file = "Output_data/Enet_a0.3_GBM_final_model.rda")


# ROC曲线绘制 -----------------------------------------------------------------
rm(list = ls())
library(pROC)
library(ggplot2)
library(dplyr)
library(ggsci)
load('Output_data/Enet_a0.3_GBM_final_model.rda')
names(rs_list)
zz <- rs_list$GSE43292
map <- roc(zz$Group,as.numeric(zz$RS),#与单因素ROC比，不同的地方
           ylim=c(0,1), # Y轴的范围
           xlim=c(1,0), # X轴的范围
           # smooth=T, #绘制平滑曲线
           ci=TRUE, # 计算95%CI置信区间
           auc.polygon.col="#fff7f7",
           auc.polygon = TRUE,
           main='GSE43292', # 定义一个图的名字
           # print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
           col='darkorange',#线的颜色
           plot = T,
           lwd=2, #线的粗细
           legacy.axes=T,  ## 横坐标是“1-specificity”，从0到1
           newpage = F) #不重新开始一页（用于画图）
legend.name <- paste0("AUC = ",round(map$auc,3)) # 我们要输出的legend
legend( # 在图的右下角
  legend=legend.name, # 输出文字
  col = 'darkorange', # 颜色定义
  lwd = 2, #线的粗细
  bty="n",
  cex = 1,x = .55,y = 0.15) # 不绘制边框

library(export)
graph2pdf(file = 'Figure/ROC/GSE43292.pdf',width = 5,height = 5)


























