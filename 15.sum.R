rm(list = ls())
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
dd <- read_xlsx(path = '12种细胞死亡基因数据合并.xlsx')

symbols <- dd$`Gene List`
aa <- unique(symbols)


map <- bitr(symbols,
            fromType = "SYMBOL",
            toType   = "ENSEMBL",
            OrgDb    = org.Hs.eg.db)

CDD_gene <- map$ENSEMBL
save(CDD_gene,file = 'CDD_gene.rda')



# 合并数据 --------------------------------------------------------------------
rm(list = ls())
library(data.table)
library(dplyr)
library(Hmisc)
library(DESeq2)
load('all_data.rda')

#交集基因数
PID <- Reduce(intersect,list(rownames(ee_list$E_MTAB_2055),
                             rownames(ee_list$GSE100927),
                             rownames(ee_list$GSE111782),
                             rownames(ee_list$GSE118481),
                             rownames(ee_list$GSE163154),
                             rownames(ee_list$GSE21545),
                             rownames(ee_list$GSE23746),
                             rownames(ee_list$GSE24495),
                             rownames(ee_list$GSE28829),
                             rownames(ee_list$GSE43292),
                             rownames(ee_list$ZZ)
))
load('CDD_gene.rda')
ID <- intersect(PID,CDD_gene)

for (i in names(ee_list)) {
  ee_list[[i]] <- ee_list[[i]][ID, , drop = FALSE]
}


for (i in names(ee_list)) {
  ee_list[[i]] <- as.data.frame(t(ee_list[[i]]))
}
CS_data <- list()
for (i in names(ee_list)) {
  CS_data[[i]] <- merge(cc_list[[i]],ee_list[[i]],by.x=1,by.y=0)
  CS_data[[i]] <- tibble::column_to_rownames(CS_data[[i]],'Sample')
  CS_data[[i]]$Group <- ifelse(CS_data[[i]]$Group=='CS',1,0)
}

class(CS_data$E_MTAB_2055$Group)

save(CS_data,file = 'Output_data/CDD_cohort.rda')

# 单因素计算 -------------------------------------------------------------------
# 1) 小工具：给定一个数据集和一个基因，返回“带符号”的p值
signed_p <- function(df, gene){
  if(!"Group" %in% names(df)) return(NA_real_)
  dd <- df[, c("Group", gene)]
  dd <- dd[complete.cases(dd), ]
  if(length(unique(dd$Group)) < 2) return(NA_real_)  # 只有单一分组则跳过
  dd$Group <- as.numeric(as.character(dd$Group))  # 0/1

  fit <- try(glm(Group ~ x, data = setNames(dd, c("Group","x")), family = binomial()), silent = TRUE)
  if(inherits(fit, "try-error")) return(NA_real_)

  sm <- summary(fit)$coefficients
  beta <- sm["x","Estimate"]; p <- sm["x","Pr(>|z|)"]
  if(is.na(p)) p <- 1
  if(is.na(beta)) beta <- 0
  if(beta < 0) p <- -p  # 负号代表“保护”（表达越高，Group=1几率越低）
  return(as.numeric(p))
}

# 2) 针对所有基因、所有队列计算“带符号”的p值矩阵 rr
genes <- setdiff(colnames(CS_data[[1]]), "Group")  # 如果各队列基因相同；否则用 Reduce(intersect, lapply(..., colnames))
rr <- t(sapply(genes, function(g){
  sapply(CS_data, function(df) signed_p(df, g))
}))
rr <- as.data.frame(rr)
rr1 <- rr[,-c(1,7)]

# 3) 统计各基因在多少个队列里显著且方向一致（阈值仍用 0.1，可改）
RS <- rowSums(rr1 > 0 & abs(rr1) < 0.05, na.rm = TRUE)  # 风险（系数>0，p<0.1）
PS <- rowSums(rr1 < 0 & abs(rr1) < 0.05, na.rm = TRUE)  # 保护（系数<0，p<0.1）

rr_df <- as.data.frame(rr)
rr_df$RS <- RS
rr_df$PS <- PS



# 4) 按你原来的口径挑基因（例如至少 8 个队列方向一致且显著）
ID <- rownames(rr_df)[rr_df$PS >= 7 | rr_df$RS >= 7]
#根据后续的热图剔除些基因
ID <- ID[-c(3,6,10,13,14,16)]#将危险和保护因素共存的基因剔除
ID <- ID[c(1:4,6:10,5,11:14)]#调整顺序，热图好看，最后是14个基因


ee <- lapply(CS_data,function(x){
  x <- x[,c('Group',ID)]
  x[,-c(1)] <- scale(x[,-c(1)])
  return(x)
})
ee <- ee[!names(ee) %in% c("E_MTAB_2055",'GSE23746')]
ee <- c(list(ZZ = ee$ZZ), ee[setdiff(names(ee), "ZZ")])
names(ee[1]) <- 'ZZ_Cohort'


save(rr_df,ID,ee,file = 'Output_data/14_gene_ee.rda')
load('Output_data/14_gene_ee.rda')

rr2 <- rr1[ID,]

for (i in 1:9) {
  rr2[,i] <- ifelse(rr2[,i]>0&abs(rr2[,i])<0.05,'Risky',ifelse(rr2[,i]<0&abs(rr2[,i])<0.05,'Protective','Not Significant'))
}




library(ComplexHeatmap)
library(ggsci)
library(RColorBrewer)
library(scales)
colnames(rr2)[9] <- 'ZZ_Cohort'
rr2 <- rr2[,c(9,1:8)]
show_col(pal_npg()(10))
colnames(rr2)

top = HeatmapAnnotation(Cohort=colnames(rr2),show_legend = F,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 12),border = T,
                                                     title_gp = gpar(fontsize = 12,fontface = "bold"),
                                                     ncol=1),
                        border = T,which = 'column',
                        col=list(Cohort = c('ZZ_Cohort'=pal_npg(alpha = 0.8)(10)[3],'GSE43292'=pal_nejm(alpha = 0.8)(8)[3],
                                            'GSE100927'=pal_npg(alpha = 0.8)(10)[5],'GSE111782'=pal_nejm(alpha = 0.8)(8)[5],
                                            'GSE118481'=pal_npg(alpha = 0.8)(10)[7],'GSE163154'=pal_nejm(alpha = 0.8)(8)[7],
                                            'GSE21545'=pal_npg(alpha = 0.8)(10)[8],'GSE24495'=pal_nejm(alpha = 0.8)(8)[8],
                                            # 'GSE28829'=pal_d3(alpha = 0.8)(8)[7],'GSE24495'=pal_d3(alpha = 0.8)(8)[6],
                                            'GSE28829'=pal_d3(alpha = 0.8)(8)[1])),
                        show_annotation_name = F,
                        annotation_name_gp = gpar(fontsize = 12))


p <- Heatmap(rr2,col = c(pal_npg()(3)[14],pal_npg()(3)[2],pal_npg()(3)[1]),
        name = 'Type',border = T,
        rect_gp = gpar(col='black',lwd=2),
        row_names_side = 'left',
        top_annotation = top,
        column_names_side = 'top',
        column_split = 1:9,
        column_title = NULL,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 12), border = T,
                                  title_gp = gpar(fontsize = 12, fontface = "bold")));p


library(export)
library(ggplot2)
graph2pdf(file='Figure/9cohorts-Gene14-heatmap.pdf',width=7,height=7)



#将行名换位symbol
Symbol <- c('FGR','MCH4','BID' ,'SLC11A1','HMOX1','JAK3','TCIRG1','IRF1',
            'BCL2A1','FDFT1' ,'FYCO1' ,'AR' ,'CLTB' ,'PRKD1')
rownames(rr2) <- Symbol

top = HeatmapAnnotation(Cohort=colnames(rr2),show_legend = F,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 12),border = T,
                                                     title_gp = gpar(fontsize = 12,fontface = "bold"),
                                                     ncol=1),
                        border = T,which = 'column',
                        col=list(Cohort = c('ZZ_Cohort'=pal_npg(alpha = 0.8)(10)[3],'GSE43292'=pal_nejm(alpha = 0.8)(8)[3],
                                            'GSE100927'=pal_npg(alpha = 0.8)(10)[5],'GSE111782'=pal_nejm(alpha = 0.8)(8)[5],
                                            'GSE118481'=pal_npg(alpha = 0.8)(10)[7],'GSE163154'=pal_nejm(alpha = 0.8)(8)[7],
                                            'GSE21545'=pal_npg(alpha = 0.8)(10)[8],'GSE24495'=pal_nejm(alpha = 0.8)(8)[8],
                                            # 'GSE28829'=pal_d3(alpha = 0.8)(8)[7],'GSE24495'=pal_d3(alpha = 0.8)(8)[6],
                                            'GSE28829'=pal_d3(alpha = 0.8)(8)[1])),
                        show_annotation_name = F,
                        annotation_name_gp = gpar(fontsize = 12))


p <- Heatmap(rr2,col = c(pal_npg()(3)[14],pal_npg()(3)[2],pal_npg()(3)[1]),
             name = 'Type',border = T,
             rect_gp = gpar(col='black',lwd=2),
             row_names_side = 'left',
             top_annotation = top,
             column_names_side = 'top',
             column_split = 1:9,
             column_title = NULL,
             heatmap_legend_param=list(labels_gp = gpar(fontsize = 12), border = T,
                                       title_gp = gpar(fontsize = 12, fontface = "bold")));p


library(export)
library(ggplot2)
graph2pdf(file='Figure/9cohorts-Symbol-Gene14-heatmap.pdf',width=7,height=7)






