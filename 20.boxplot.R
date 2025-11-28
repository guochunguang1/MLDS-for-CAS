rm(list = ls())
library(data.table)
library(dplyr)
library(Hmisc)
library(DESeq2)
library(reshape2 )
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggplot2)
load('Output_data/14_gene_ee.rda')
ID <- c("ENSG00000125347","ENSG00000079459","ENSG00000163820")
Symbol <- c('IRF1','FDFT1','FYCO1')
for (i in names(ee)) {
  ee[[i]] <- ee[[i]][,c('Group',ID)]
  colnames(ee[[i]])[2:4] <- Symbol
  ee[[i]]$Group <- ifelse(ee[[i]]$Group==1,'CS','Con')
}


## 调色 & 输出目录
mycol <- c('#21b6af','#eeba4d')
out_dir <- "Figure/Boxplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#书写批量导出函数
plot_one <- function(dd, nm){
  stopifnot(is.data.frame(dd))
  if(!("Group" %in% names(dd))) {
    message("跳过 ", nm, "：无 Group 列")
    return(invisible(NULL))
  }
  dat <- melt(dd, id.vars = "Group",
              variable.name = "Gene", value.name = "Expression")
  # 分组顺序（没有的水平会被自动丢弃）
  dat$Group <- factor(dat$Group, levels = c("Con","CS"))

  p <- ggplot(dat, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.7, size = .3, outlier.size = 0.001) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          legend.position = "top") +
    xlab("") + ylab("Gene expression") +
    labs(fill = "Group") +
    scale_fill_manual(values = mycol) +
    ggtitle(nm)

  ## 仅在分组≥2时添加显著性
  if(length(na.omit(unique(dat$Group))) >= 2){
    p <- p + stat_compare_means(
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***","**","*"," ")),
      label = "p.signif"
    )
  }

  ## 导出（任选其一）
  ggsave(filename = file.path(out_dir, paste0(nm, ".pdf")),
         plot = p, width = 5, height = 5)

  # 若更偏好 export::graph2pdf：
  # print(p)
  # export::graph2pdf(file = file.path(out_dir, paste0(nm, ".pdf")),
  #                   width = 5, height = 5)
}

## 批量导出
invisible(lapply(names(ee), function(nm){
  try(plot_one(ee[[nm]], nm), silent = TRUE)
}))

message("完成：输出在 ", normalizePath(out_dir))


































