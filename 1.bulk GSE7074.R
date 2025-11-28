#颈动脉数据准备-----------
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

eset <- getGEO(filename = 'GSE7074_series_matrix.txt.gz',getGPL = F)
expr <- as.data.frame(exprs(eset))
range(expr)
clin <- pData(eset)
clin <- clin[,c(2,1)]
group <- clin
rm(eset,clin)
#仅保留颈动脉患者
group <- group[grepl("CA", group$title), ]
expr <- expr[,group$geo_accession]
colnames(group) <- c('Sample','Group')
rownames(group) <- NULL
group <- tibble::column_to_rownames(.data = group,var = 'Sample')
table(group$Group)
group$Group <- ifelse(group$Group=='Macroscopically intact tissue','Con','CS')

#数据转换
load('GPL4868.rda')
ann <- ann[,c(1,2)]
# nchar('ENSG00000223972')
# ann$ENSG <- str_sub(ann$ENSG,1,15)
exp <- merge(ann,expr,by.x = 1,by.y = 0)[,-1]
exp[,-1] <- apply(exp[,-1],2,as.numeric)

#保留平均表达量最大的基因
exp1 <- exp %>%
  group_by(SYM) %>%
  summarise(across(where(is.numeric), mean)) %>%  # 你也可以换成 max()
  ungroup()

rownames(exp1) <- NULL
exp1 <- tibble::column_to_rownames(exp1,'SYM')
group1 <- tibble::rownames_to_column(group,'Sample')
exp2 <- exp1[,colnames(exp1)%in%group1$Sample]
exp3 <- apply(exp2,2,as.numeric)
exp3 <- as.data.frame(scale(exp3))
rownames(exp3) <- rownames(exp2)
range(exp3,na.rm = T)

GSE7074_cc <- group1
GSE7074_ee <- exp3

save(GSE7074_cc,GSE7074_ee,file = 'CEA_data/GSE7074.rda')


# -------------------------------------------------------------------------
rm(list = ls())
load('CEA_data/GSE7074.rda')


