rm(list = ls())
library(ggplot2)
library(ggsci)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggbreak)
library(tibble)

load('Output_data/diagnostic_ZZ_105combo_summary.rda')
load('Output_data/dd_Cindex.rda')
dd <- dd_all_mean
#将筛选算法放到左边，建模算法放到右边
swap_plus <- function(x){
  i <- grepl("\\+", x)
  x[i] <- sub("^\\s*(.*?)\\s*\\+\\s*(.*?)\\s*$", "\\2+\\1", x[i])
  x
}

dd$Model <- swap_plus(dd$Model)
range(dd$Cindex)

dd$ll <- sprintf("%.3f",dd$Cindex)

dd%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width=0.7,stat = 'identity',fill='#1993CF',position = position_dodge(width = 0.9))+
  scale_x_continuous(breaks = c(0,0.65,1),labels = c(0,0.65,1),limits = c(0,1))+
  scale_x_break(c(0.05,0.4),scales = 20)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  ggtitle('Mean Cindex')+
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.6,vjust = -3),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())+
  geom_text(aes(label=ll),hjust=1.1,size=2.7,color='black')

ggsave(filename = 'Figure/all-model-bar.pdf',width = 4,height = 22)


# -------------------------------------------------------------------------
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)
dd3 <- dd2
dd3$Model <- swap_plus(dd3$Model)


dd3$Mean <- rowMeans(dd2[,-1])
dd3 <- dd3[order(dd3$Mean,decreasing = T),]
rownames(dd3) <- NULL
dd3 <- dd3%>%column_to_rownames('Model')


Cohort <- c(pal_npg(alpha = 0.8)(10)[3],
            pal_npg(alpha = 0.8)(10)[5],
            pal_nejm(alpha = 0.8)(8)[5],
            pal_npg(alpha = 0.8)(10)[7],
            pal_nejm(alpha = 0.8)(8)[7],
            pal_npg(alpha = 0.8)(10)[8],
            pal_nejm(alpha = 0.8)(8)[8],
            pal_d3(alpha = 0.8)(8)[1],
            pal_nejm(alpha = 0.8)(8)[3])


colnames(dd3)[1] <- 'ZZ_Cohort'
names(Cohort) <- colnames(dd3)[-10]
Top = HeatmapAnnotation(Cohort=colnames(dd3)[-10],
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 9,col='black'),
                                                     title_gp = gpar(fontsize = 9, fontface = "bold",col='black'),ncol=1),
                        gap=unit(1, "mm"),
                        border = T,
                        col=list(Cohort=Cohort),
                        show_annotation_name = TRUE,
                        annotation_name_side="right",
                        annotation_name_gp = gpar(fontsize = 10,col='black', fontface = "bold"))
dd3[,-10] <- round(dd3[-10],2)
range(dd3[,-10])
1.38/2

ht <- Heatmap(as.matrix(dd3[,-10]),col = colorRamp2(c(0.4,0.7,1),c('white','#21b6af','#eeba4d')),
        rect_gp = gpar(color='grey30',lwd=0.2),name = 'Cindex',
        cluster_rows = F,cluster_columns = F,column_names_gp = gpar(fontsize=13),
        row_names_side = 'left',column_names_side = 'top',column_names_rot = 0,column_names_centered = T,
        cell_fun = function(j,i,x,y,width,height,fill){grid.text(dd3[i,j],x,y,gp = gpar(fontsize = 9))},
        column_split = 1:9,top_annotation = Top,
        show_column_names = F,column_title = NULL,
        row_title = NULL
)

draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")


library(export)
graph2pdf(file='Figure/all_model_heatmap.pdf',width = 8,height = 22)
dev.off()
























