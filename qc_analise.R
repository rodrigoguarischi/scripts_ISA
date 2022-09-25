#####################################################
#################  CNV status #######################
#####################################################

setwd('C:/Users/Camila/OneDrive/Área de Trabalho/Análises/Results_novo/')
library(ggplot2)
library(dplyr)
library(rcompanion)
library (plyr)

qc_file <- read.delim("qc_table_cnv.txt")
head(qc_file)

boxplot(qc_file$MAPD,
        main = "Qualidade - MAPD",
        xlab = "Métrica de Qualidade",
        at = c(1),
        names = "MAPD",
        col = "darkred",
        border = "red",
        horizontal = TRUE,
        notch = TRUE)

boxplot(qc_file$WavinessSD,
        main = "Qualidade - Waviness",
        xlab = "Métrica de Qualidade",
        at = c(1),
        names = "Waviness",
        col = "darkblue",
        border = "blue",
        horizontal = TRUE,
        notch = TRUE)

plot(qc_file$Autosome...LOH,qc_file$MAPD,
     xlab = "%LOH Chr Aut", 
     ylab = "MAPD", 
     main= "% de LOH Chr por MAPD", 
     pch = 16, col="darkblue")


plotNormalHistogram(qc_file$MAPD,
                    main = "Histrograma MAPD", col = "deeppink4", 
                    linecol = "red", prob = TRUE, x = qc_file$MAPD,
                    xlab = "MAPD")

plotNormalHistogram(qc_file$WavinessSD,
                    main = "Histrograma Waviness", col = "darkgreen", 
                    linecol = "green", prob = TRUE, x = qc_file$WavinessSD,
                    xlab = "Waviness")


plotNormalHistogram(qc_file$Autosome...LOH,
                    main = "Histrograma LOH", col = "gold", 
                    linecol = "orange", prob = TRUE, x = qc_file$WavinessSD,
                    xlab = "LOH")


ggplot(qc_file, aes(x=reorder(CN.passes.MAPD, CN.passes.MAPD,  
                               function(x)-length(x)))) +
  geom_bar(fill='darkred', width=0.2, position = "dodge") +
  labs(x='MAPD Pass', y='Amostras (N)') +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(vjust=0.5), colour="black") +
  ggtitle("Status do MAPD")


ggplot(qc_file, aes(x=reorder(CN.passes.WavinessSD, CN.passes.WavinessSD,  
                              function(x)-length(x)))) +
  geom_bar(fill='orange', width=0.2, position = "dodge") +
  labs(x='WavinessSD Pass', y='Amostras (N)') +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(vjust=0.5), colour="black") +
  ggtitle("Status do WavinessSD")


#########################################################

library(ggplot2)
library("rcompanion")

dta_raw <- read.delim('cnv.bed')
summary(dta_raw)

samples <- read.delim('samples.txt')
count_calls <- read.delim('calls.txt')

calls_samples <- cbind.data.frame(samples,count_calls)
calls_samples_ord  <- data.frame(calls_samples[order(calls_samples$Calls),])
calls_samples_edit <- data.frame(calls_samples_ord$Sample, 
                                 calls_samples_ord$Calls)
head(calls_samples_edit)

barplot(calls_samples_edit$calls_samples_ord.Calls, 
        col="darkblue",xlab="Amostras", ylab="Número de CNVs", 
        main="Número de CNVs por amostra")

plotNormalHistogram(calls_samples_edit$calls_samples_ord.Calls,
                    main = "Número de CNVs por amostra", col = "darkred", 
                    linecol = "red", prob = FALSE )

summary(calls_samples_edit$calls_samples_ord.Calls)

boxplot(calls_samples_edit$calls_samples_ord.Calls, col = "darkgreen",
        main = "Número de CNVs por amostra", ylab = "Número de CNVs", 
        horizontal = TRUE)

write.csv(calls_samples_edit,
          file = "C:/Users/Camila/OneDrive/Área de Trabalho/Análises/Results_novo//calls_ord.txt",
          quote = FALSE)

dta_edit <- dplyr::filter(dta_raw,!grepl('X',dta_raw$Chromosome))
dta <- na.omit(dta_edit)

table(dta$Chromosome)

###########################################

gender <- read.delim('genero.txt')
summary(gender)

ggplot(data=gender, aes(x=Tipo, y=N)) +
  geom_bar(stat="identity", color=c("steelblue","deeppink4",
                                    "orange"), fill=c("steelblue","deeppink4",
                                                      "orange"),width = 0.3)+
  geom_text(aes(label=N), vjust=-1, color="black", size=3.5)+
  labs(title="Número de Individuos por Gênero", x="Tipo", 
       y="Número")+
  theme_minimal()

###############################################

segm <- read.delim('segmentos.txt')
summary(segm)

ggplot(data=segm, aes(x=Tipo, y=Calls)) +
  geom_bar(stat="identity", color=c("deeppink4","gold4","darkgreen","steelblue"), 
           fill=c("deeppink4","gold4","darkgreen","steelblue"),width = 0.3)+
  geom_text(aes(label=Calls), vjust=-1, color="black", size=3.5)+
  labs(title="Status de Número de Cópias por Segmento Cromossômico", 
       x="Tipos de Copy Number", 
       y="Número de Segmentos")+
  theme_minimal()

#################################################

chromo <- read.delim('chr_calls.txt')
summary(chromo)

color1 <- colorRampPalette(c("darkblue","darkred"))
barplot(chromo$Calls,
        col = color1(23), 
        main = "Quantidade de CNVs por cromossomo",
        ylab = "Numero de CNVs",
        xlab = "Cromossomos",
        ylim=c(0,4000),
        names.arg = c(chromo$Chr))

#################################################################

loh_pct <- data.frame(qc_file$Sample.Filename,qc_file$Autosome...LOH)
summary(loh_pct)

plotNormalHistogram(loh_pct$qc_file.Autosome...LOH,
                    main = "Porcentagem de LOH por amostra", 
                    col = rgb(0,0,1,0.6), 
                    linecol = "red", prob = TRUE, xlab="% Autosome LOH")

#################################################

tam <- read.delim('tamanho_cat.txt')
summary(tam)

ggplot(data=tam, aes(x=Cat, y=Call)) +
  geom_bar(stat="identity", color="deepskyblue4", 
           fill="deepskyblue4",width = 0.3)+
  geom_text(aes(label=Call), vjust=-1, color="black", size=3.5)+
  labs(title="Quantidade de CNV por Tamanho", x="Tamanho das CNVs", 
       y="Número de Chamadas")+
  theme_minimal()

###############################################

tipo_cnv <- read.delim('tipos_cnv.txt')
summary(tipo_cnv)

ggplot(data=tipo_cnv, aes(x=Tipo, y=Calls)) +
  geom_bar(stat="identity", color=c("steelblue","deepskyblue4"), 
           fill=c("steelblue","deepskyblue4"),width = 0.3)+
  geom_text(aes(label=Calls), vjust=-1, color="black", size=3.5)+
  labs(title="Número de chamadas de CNV por Tipo", x="Tipo de CNVs", 
       y="Número de Chamadas")+
  theme_minimal()


################################################

tam_loh <- read.delim('tamanho_cat_loh.txt')
summary(tam_loh)

ggplot(data=tam_loh, aes(x=Tipo, y=Calls)) +
  geom_bar(stat="identity", color="darkgreen", 
           fill="darkgreen",width = 0.3)+
  geom_text(aes(label=Calls), vjust=-1, color="black", size=3.5)+
  labs(title="Quantidade de LOH por Tamanho", x="Tamanho das LOHs", 
       y="Número de Chamadas")+
  theme_minimal()

################################################

all_size <- read.delim('cnv_filter_k.txt')
summary(all_size)

dta_chr <- data.frame(all_size$Chromosome,all_size$Size)
summary(dta_chr)

dta_chr_sorted <- dta_chr[order(dta_chr$all_size.Chromosome,
                                dta_chr$all_size.Size),]
summary(dta_chr_sorted)

ggplot(data = dta_chr_sorted, aes(x = all_size.Size)) +
  geom_histogram() + ggtitle('')+
  facet_wrap(~all_size.Chromosome)

stats <- tapply(dta_chr_sorted$all_size.Size, dta_chr$all_size.Chromosome,
       function(x) format(summary(x), scientific = FALSE))


df1 <- ldply (stats, data.frame)

write.csv(df1, "C:/Users/Camila/OneDrive/Área de Trabalho/chr_metrics.txt", 
          quote=FALSE)

#######################################################

all_size_loh <- read.delim('C:/Users/Camila/OneDrive/Área de Trabalho/Análises/loh_data_filter.txt')
summary(all_size_loh)

dta_chr_loh <- data.frame(all_size_loh$Chromosome,all_size_loh$Size)
summary(dta_chr_loh)

dta_chr_sorted_loh <- dta_chr_loh[order(dta_chr_loh$all_size_loh.Chromosome,
                                        dta_chr_loh$all_size_loh.Size),]
summary(dta_chr_sorted_loh)

ggplot(data = dta_chr_sorted_loh, aes(x = all_size_loh.Size)) +
  geom_histogram() + ggtitle('LOH')+
  facet_wrap(~all_size_loh.Chromosome)

stats2 <- tapply(dta_chr_sorted_loh$all_size_loh.Size, dta_chr_sorted_loh$all_size_loh.Chromosome,
                function(x) format(summary(x), scientific = FALSE))

df2 <- ldply (stats2, data.frame)

write.csv(df2, "C:/Users/Camila/OneDrive/Área de Trabalho/chr_metrics_loh.txt", 
          quote=FALSE)

##############################################################

summary(all_size)
dupdel_chr <- data.frame(all_size$Chromosome,all_size$State)
count_tipo <- tapply(dupdel_chr$all_size.State,dupdel_chr$all_size.Chromosome,
                 function(x) format(count(x), scientific = FALSE))

del_chr <- dplyr::filter(dupdel_chr, !grepl('3', dupdel_chr$all_size.State))
del_chr2 <- dplyr::filter(del_chr, !grepl('2', del_chr$all_size.State))
head(del_chr2)

p <- ggplot(data = del_chr2, aes(x=all_size.State)) +
  geom_histogram()+
  ggtitle('Deleções por Cromossomo')+
  facet_wrap(~all_size.Chromosome)

p + xlab("Categoria das deleções (0 e 1)") + ylab("Quantidade de chamadas") + 
  ylim(0,220)

###############################################################


