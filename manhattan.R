##############
# 1. library #
##############

library('qqman')
library("rio")

############
# 2. input #
############

data <- read.delim('C:/Users/Camila/OneDrive/�rea de Trabalho/stdout.txt')
head(data)

###############################################
# 3. Create a dataframe (load file in memory) #
###############################################

CHR <- c(data$Chromosome)
BP <- as.numeric(data$Physical.Position)
P <- c(data$H.W.p.Value)
SNP <- c(data$dbSNP.RS.ID)

CHR[CHR=="X"]=23
CHR[CHR=="Y"]=24
CHR[CHR=="MT"]=25

table(CHR)

dados <- data.frame(CHR,BP,P,SNP)
str(dados)

dados2 <- na.omit(dados)
dados2$CHR <- as.numeric(dados2$CHR)
str(dados2)
table(dados2$CHR)

dados3 <- data.frame(na.omit(dados2))

#############
# 4. Plots  #
#############

manhattan(dados3, col = c("gray10", "gray60"), main = "Manhattan Plot - ISA", 
          ylim = c(0, 11), annotateTop = TRUE, 
          suggestiveline = -log10(1e-06), genomewideline = FALSE)

qq(dados$P)

###################
# 5. Export Data  #
###################

export(dados2, file = "C:/Users/Camila/OneDrive/Area de Trabalho/dataR.txt")
