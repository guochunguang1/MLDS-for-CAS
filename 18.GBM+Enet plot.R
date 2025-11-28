library(gbm); library(dplyr); library(ggplot2)

load("Output_data/14_gene_ee.rda")
est_data <- ee$ZZ
# 重新拟合一次 GBM（与 sel_gbm 相同参数，保证可视化一致）
Xn <- setdiff(colnames(est_data), "Group")
set.seed(91)
sel_gbm <- function(df){
  Xn <- setdiff(names(df), "Group")
  fit <- gbm(Group ~ ., data=df, distribution="bernoulli",
             n.trees=3000, interaction.depth=3, n.minobsinnode=10,
             shrinkage=0.01, cv.folds=5, verbose=FALSE)
  best <- gbm.perf(fit, method="cv", plot.it=FALSE)
  si <- summary.gbm(fit, n.trees=best, plotit=FALSE)
  rid <- si$var[si$rel.inf >= 2]
  if (length(rid) >= 2) return(rid)
  choose_k <- function(p) max(5, min(30, round(p/3)))
  si$var[order(si$rel.inf, decreasing=TRUE)][1:choose_k(nrow(si))]
}

# 3) 生成 gbm_feats（与建模同一个随机种子）
gbm_feats <- sel_gbm(est_data)



fit_gbm <- gbm(reformulate(Xn, "Group"), data = est_data, distribution = "bernoulli",
               n.trees = 3000, interaction.depth = 3, n.minobsinnode = 10,
               shrinkage = 0.01, cv.folds = 5, verbose = FALSE)
best <- gbm.perf(fit_gbm, method = "cv", plot.it = FALSE)

# 1) VI 条形图（高亮被 sel_gbm 选中的基因）
vi <- summary.gbm(fit_gbm, n.trees = best, plotit = FALSE) %>%
  arrange(desc(rel.inf)) %>%
  mutate(selected = var %in% gbm_feats,
         var = factor(var, levels = rev(var)))

vi$var <- c('FDFT1','FYCO1','IRF1','MCH4','HMOX1','BCL2A1','FGR',
            'PRKD1','TCIRG1','AR','BID','CLTB','JAK3','SLC11A1')

vi <- vi[order(vi$rel.inf, decreasing = TRUE), ]

topN <- min(20, nrow(vi))   # 只画前20个，需全部就改成 nrow(vi)
p_vi <- ggplot(vi[1:topN,], aes(x = reorder(var, rel.inf), y = rel.inf, fill = selected)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE`="#D55E00", `FALSE`="grey80"), guide = "none") +
  labs(x = NULL, y = "Relative influence (%)",
       title = "GBM variable importance (CV-optimal trees)");p_vi
ggsave("Figure/GBM_VI_bar.png", p_vi, width = 7, height = 5)
ggsave("Figure/GBM_VI_bar.pdf", p_vi, width = 7, height = 5)






# 2) （可选）对进入最终模型的基因画 PDP
for(g in gbm_feats){
  pd <- plot(fit_gbm, i.var = g, n.trees = best, return.grid = TRUE)
  p_pd <- ggplot(pd, aes_string(x = g, y = "y")) +
    geom_line() +
    labs(x = g, y = "Partial dependence", title = paste("GBM PDP -", g)) +
    theme_classic(base_size = 12)
  ggsave(sprintf("Figure/GBM_PDP_%s.png", g), p_pd, width = 5, height = 4, dpi = 300)
}
# 2) Enet（alpha=0.3）CV 曲线 + 系数路径

library(glmnet); library(ggplot2); library(dplyr)
train_df <- est_data

## 2) 选择用于 Enet 的特征集合
all_feats <- setdiff(names(train_df), "Group")
feats <- if (exists("gbm_feats")) intersect(all_feats, gbm_feats) else all_feats
stopifnot(length(feats) >= 2)

## 3) 构造 X / y
X <- as.matrix(train_df[, feats, drop = FALSE])
y <- as.numeric(train_df$Group)
cvfit    <- cv.glmnet(X, y, family = "binomial", alpha = 0.3, nfolds = 10)
fit_path <- glmnet(X, y, family = "binomial", alpha = 0.3)

# 2) 系数正则化路径（横轴 log(lambda)），并标注 lambda.min / 1se
png("Figure/Enet_a0.3_coef_path.png", width = 1600, height = 1200, res = 150)
plot(fit_path, xvar = "lambda", label = FALSE)
abline(v = log(cvfit$lambda.min), lty = 2, col = "red")
abline(v = log(cvfit$lambda.1se), lty = 3, col = "red")
title("Elastic Net (alpha = 0.3) coefficient path")

# 一键把两条竖线对应的“最终基因+系数”列出来 --------------------------------------------------

# 2) 取在 lambda.min 下真正进入模型的系数（非零）
co <- as.matrix(coef(cvfit, s = "lambda.min"))[-1, 1]
co <- co[co != 0]                             # 只保留非零
genes <- names(co)
ys    <- as.numeric(co)

# 3) 在“最左侧 x”附近标注名字（向左偏一点），正负用不同颜色
par(xpd = NA)                                 # 允许在绘图区外标注
x_left <- min(log(fit_path$lambda))           # 左端的 log(lambda)
cols   <- ifelse(ys > 0, "firebrick", "steelblue")

text(x = x_left, y = ys, labels = genes,
     pos = 2,                               # 文本在点的左侧
     col = cols, cex = 0.9)
par(xpd = FALSE)

library(export)
graph2pdf(file = 'Figure/enet.pdf',width =7,height = 5)

dev.off()




