##############
# 1. library #
##############

library(ggplot2)

############
# 2. input #
############

data_results <- read.delim('C:/Users/Camila/OneDrive/Área de Trabalho/Novos gráficos - All Samples/genot_data_.txt')
head(data_results)
table(data_results$computed_gender)

############
# Chr Plot #
############

ggplot(data_results, aes(x = data_results$cn.probe.chrXY.ratio_gender_meanX,
                         y = data_results$cn.probe.chrXY.ratio_gender_meanY,
                         colour = data_results$computed_gender)) +
  geom_point() +
  scale_colour_manual(values = c("darkgreen", "darkblue", "red")) +
  labs(x = "ChrX", y = "ChrY", 
       title = "Computed Gender", 
       colour = data_results$computed_gender )



