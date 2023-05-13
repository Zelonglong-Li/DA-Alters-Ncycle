
library(dplyr)
library(linkET)
library(ggplot2)

speciese <- read.csv("Gene_metabolite.csv",header = T,row.names = 1,encoding = "UTF-8")#读取物种数据
env <- read.csv("envfactor.csv",header = T,row.names = 1)#读取环境因子

mantel01 <- mantel_test(speciese, env,
                        spec_select = list(Genes = 1:13,Metabolite=14:16)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
mantel01

diag.node.size=100

set_corrplot_style(colours = c("OrangeRed", "DeepSkyBlue"))

qcorrplot(correlate(env), type = "upper", diag = FALSE) +
  
  geom_square() + 
  
  geom_mark( sig_thres = 0.05, only_mark = T, sig_level = 0.05, mark = "*", vjust = 0.65, size = 5, colour = "grey90")+ 
  
  geom_couple(aes(colour = pd, size = rd), data = mantel01, curvature = 0.1) +
  
  geom_diag_label(
    geom = "text",
    angle=0,
    nudge_x = -0.2,
    nudge_y = 0.3)+
  theme(
    axis.text.x = element_text(vjust = -0.2),
    axis.text.y = element_text(colour ="white")
    )+
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("OliveDrab", "burlywood", "azure3")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))


