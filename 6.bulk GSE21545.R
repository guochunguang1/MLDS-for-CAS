rm(list = ls())
library(stringr)
library(dplyr)
library(limma)
library(data.table)
library(tibble)
library(BioEnricher)
library(stringr)
library(data.table)
library(tibble)
library(dplyr)
library(GEOquery)

eset <- getGEO(filename = 'GSE21545_series_matrix.txt.gz',getGPL = F)
expr <- as.data.frame(exprs(eset))
range(expr)
clin <- pData(eset)
clin <- clin[,c(1,2,38,41,42)]
group <- clin
rm(eset,clin)

colnames(group) <- c('Tissue','Sample','Age','Group','IS_time')
rownames(group) <- NULL
group <- tibble::column_to_rownames(.data = group,var = 'Sample')
table(group$Group)
group$Group <- ifelse(group$Group=='TRUE','IS','Con')


group$Tissue <- ifelse(grepl("pbmc", group$Tissue, ignore.case = TRUE),
                        "PBMC", "Plaque")
group$Age <- as.numeric(group$Age)
group$IS_time <- as.numeric(group$IS_time)

#数据转换
load('GPL570.rda')
ann <- ann[,c(1,3)]
nchar('ENSG00000223972')
ann$ENSG <- str_sub(ann$ENSG,1,15)
exp <- merge(ann,expr,by.x = 1,by.y = 0)[,-1]
#保留平均表达量最大的基因
exp1 <- exp %>%
  group_by(ENSG) %>%
  summarise(across(where(is.numeric), mean)) %>%  # 你也可以换成 max()
  ungroup()

rownames(exp1) <- NULL
exp1 <- tibble::column_to_rownames(exp1,'ENSG')
group1 <- tibble::rownames_to_column(group,'Sample')
exp2 <- exp1[,colnames(exp1)%in%group1$Sample]
exp3 <- apply(exp2,2,as.numeric)
exp3 <- as.data.frame(scale(exp3))
rownames(exp3) <- rownames(exp2)
range(exp3)

GSE21545_cc <- group1
GSE21545_ee <- exp3

save(GSE21545_cc,GSE21545_ee,file = 'CEA_data/GSE21545.rda')


# -------------------------------------------------------------------------
rm(list = ls())
load('CEA_data/GSE21545.rda')
