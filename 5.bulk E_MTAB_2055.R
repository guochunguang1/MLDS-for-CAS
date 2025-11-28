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
eset <- read.delim('Data/Rawdata/Group_Probe_Profile_RAW.txt',
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

clin <- read.delim('Data/Rawdata/E-MTAB-2055.sdrf.txt',
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
# 提取探针ID
probe_id <- eset$PROBE_ID

# 保留包含 "AVG_Signal" 的列
expr_mat <- eset[, grepl("AVG_Signal", colnames(eset))]

# 设置行为 probe ID
rownames(expr_mat) <- probe_id
colnames(expr_mat) <- sub("\\..*$", "", colnames(expr_mat))

eset <- expr_mat

clin <- clin[,c(1,24)]
clin <- clin[-27,]
eset <- eset[,colnames(eset)%in%clin$`Source Name`]
expr <- eset

group <- clin
rm(eset,clin,expr_mat,probe_id)

colnames(group) <- c('Sample','Group')
rownames(group) <- NULL
group <- tibble::column_to_rownames(.data = group,var = 'Sample')
table(group$Group)
group$Group <- ifelse(group$Group=='ruptured','CS','Con')

#数据转换
load('GPL6947.rda')
ann <- ann[,c(1,4)]
nchar('ENSG00000223972')
ann$ENSG <- str_sub(ann$ENSG,1,15)
exp <- merge(ann,expr,by.x = 1,by.y = 0)[,-1]
#清洗数据格式，有的数据有两个点
# Step 1: 提取行名并保留表达部分
expr_char <- exp    # 例如你用 read.table 读入的原始数据
gene_ids <- expr_char[[1]]    # 第一列为行名（如ENSG）
expr_only <- expr_char[, -1]  # 去掉第一列，保留表达值部分

#  Step 2: 将字符格式转换为合法数值
expr_fixed <- as.data.frame(
  apply(expr_only, 2, function(x) {
    sapply(x, function(val) {
      # 如果是NA或空，直接返回NA
      if (is.na(val) || val == "") return(NA_real_)

      # 如果包含两个以上的点（即千位分隔形式），处理
      if (length(gregexpr("\\.", val)[[1]]) >= 2) {
        as.numeric(gsub("\\.", "", val)) / 1000
      } else {
        as.numeric(val)
      }
    })
  }),
  stringsAsFactors = FALSE
)

expr_fixed$ENSG <- gene_ids

exp <- expr_fixed[,c(48,1:47)]

class(exp$`100273-01`)

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

E_MTAB_2055_cc <- group1
E_MTAB_2055_ee <- exp3

save(E_MTAB_2055_cc,E_MTAB_2055_ee,file = 'CEA_data/E_MTAB_2055.rda')



# -------------------------------------------------------------------------
rm(list = ls())
load('CEA_data/E_MTAB_2055.rda')





