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

load('Raw data/zhengzhou/zhengzhou cohort.rda')
exp <- data3
rm(data3)
colnames(exp)[1] <- 'ENSG'
#保留平均表达量最大的基因
exp1 <- exp %>%
  group_by(ENSG) %>%
  summarise(across(where(is.numeric), mean)) %>%  # 你也可以换成 max()
  ungroup()
rownames(exp1) <- NULL
exp1 <- tibble::column_to_rownames(exp1,'ENSG')

group1 <- data.frame(Sample = colnames(exp1),
                     Group = rep('CS',46))

exp2 <- exp1[,colnames(exp1)%in%group1$Sample]
exp3 <- apply(exp2,2,as.numeric)
exp3 <- as.data.frame(scale(exp3))
rownames(exp3) <- rownames(exp2)
range(exp3)

ZZ_cc <- group1
ZZ_ee <- exp3

save(ZZ_cc,ZZ_ee,file = 'CEA_data/ZZ.rda')


# -------------------------------------------------------------------------
rm(list = ls())
load('CEA_data/ZZ.rda')












