##############
# Libraries  #
##############

library("gdsfmt")
library("SNPRelate")
library("dplyr")
library("scatterplot3d")
library("ca")

###################################################

##################
# Open Databases #
##################

setwd('C:/Users/Camila/OneDrive/Documentos/Mestrado/mestrado-isa/')

# 707 samples and 1000g 

bed.fn <- "isamerge1000.bed"
bim.fn <- "isamerge1000.bim"
fam.fn <- "isamerge1000.fam"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn="isamerge1000.gds", verbose=TRUE)
genofile <- snpgdsOpen("isamerge1000.gds")

# ISA after LD filter

bed2.fn <- "isa_hg37.prunedld.bed"
bim2.fn <- "isa_hg37.prunedld.bim"
fam2.fn <- "isa_hg37.prunedld.fam"

snpgdsBED2GDS(bed2.fn, fam2.fn, bim2.fn, out.gdsfn="isa.filter.gds", verbose=TRUE)
genofile2 <- snpgdsOpen("isa.filter.gds")

# ISA before LD filter

bed3.fn <- "isa_hg37_filter2.bed"
bim3.fn <- "isa_hg37_filter2.bim"
fam3.fn <- "isa_hg37_filter2.fam"

snpgdsBED2GDS(bed3.fn, fam3.fn, bim3.fn, out.gdsfn="isafilter3.gds", verbose=TRUE)
genofile3 <- snpgdsOpen("isafilter3.gds")

# Related Samples 

bed4.fn <- "isacomm_excluidas_2.bed"
bim4.fn <- "isacomm_excluidas_2.bim"
fam4.fn <- "isacomm_excluidas_2.fam"

snpgdsBED2GDS(bed4.fn, fam4.fn, bim4.fn, out.gdsfn="isafilter4.gds", verbose=TRUE)
genofile4 <- snpgdsOpen("isafilter4.gds")

# Raw data - All SNPs

bed5.fn <- "isa_hg37.bed"
bim5.fn <- "isa_hg37.bim"
fam5.fn <- "isa_hg37.fam"

snpgdsBED2GDS(bed5.fn, fam5.fn, bim5.fn, out.gdsfn="isafilter5.gds", verbose=TRUE)
genofile5 <- snpgdsOpen("isafilter5.gds")

###################################################

##################
# SNP Frequency  #
##################

# After LD filter

snplist <- snpgdsSNPList(genofile2)
snplistchr <- data.frame(snplist$chromosome, snplist$snp.id)
freq <- as.data.frame(table(snplistchr$snplist.chromosome))
color1 <- colorRampPalette(c("darkblue","lightblue"))

sum(freq$Freq)
barplot(freq$Freq,
        col = color1(22), 
        main = "Quantidade de SNPs por cromossomo",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,60000),
        names.arg = c(freq$Var1))

# Before LD filter

snplist2 <- snpgdsSNPList(genofile3)
snplistchr2 <- data.frame(snplist2$chromosome, snplist2$snp.id)
freq2 <- as.data.frame(table(snplistchr2$snplist2.chromosome))
color2 <- colorRampPalette(c("darkgreen","lightgreen"))

sum(freq2$Freq)

barplot(freq2$Freq,
        col = color2(22), 
        main = "Quantidade de SNPs por cromossomo sem filtro de LD",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,60000),
        names.arg = c(freq2$Var1))

# All SNPs

snplist5 <- snpgdsSNPList(genofile5)
snplistchr5 <- data.frame(snplist5$chromosome, snplist5$snp.id)
freq5 <- as.data.frame(table(snplistchr5$snplist5.chromosome))
color5 <- colorRampPalette(c("darkgreen","lightgreen"))

sum(freq5$Freq)

barplot(freq5$Freq,
        col = "darkred", 
        main = "Quantidade de SNPs por cromossomo - Todos Marcadores",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,60000),
        names.arg = c(freq5$Var1))

###################################################

###################
# Calculating PCA #
###################

# 707 samples and 1000g

CP <- snpgdsPCA(genofile)
head(CP)
summary(CP)

# ISA after LD filter

CP2 <- snpgdsPCA(genofile2)

###################################################

#############################################
# Saving eigenvectors and eigenval on file  #
#############################################

write.table(CP$eigenvec, file="isamerge1000_eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$eigenvect, file="eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep="\t" , quote=FALSE)

write.table(CP$sample.id, file="samples.txt", col.names=TRUE,
            row.names=TRUE, sep="\t" , quote=FALSE)

write.table(CP$eigenval, file="isamerge1000_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$eigenval, file="isamerge1000_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$snp.id, file="snps.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP2$eigenvec, file="isa_eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep=",")

write.table(CP2$eigenval, file="isa_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep=",")

summary(CP$eigenvect)
summary(CP2$eigenvect)

###################################################

#########################################
#  Explanation percentage of components #
#########################################

# 707 samples and 1000g

pc.percent <- CP$eigenval[1:100] / sum(CP$eigenval[1:25])
write.table(pc.percent, file="percent_eigenval_isamerge1000.txt", sep=",", quote = FALSE)
cpplot <- c(pc.percent[1:6])
cp1cp2 <- cpplot[1] + cpplot[2]

barplot(cpplot,
        beside = TRUE, 
        width = 0.7, 
        ylim = c(0,0.6),
        col=1:32,
        legend = c("CP1","CP2","CP3","CP4","CP5","CP6"),
        main = "Explanation percentage of components - Isa and 1000g")

# ISA after LD filter

pc.percent2 <- CP2$eigenval[1:100] / sum(CP2$eigenval[1:25])
write.table(pc.percent2, file="percent_eigenval_isa.txt", sep=",")

cpplot2 <- c(pc.percent2[1:6])
cp1cp22 <- cpplot2[1] + cpplot2[2]

barplot(cpplot2,
        beside = TRUE, 
        width = 0.7, 
        ylim = c(0,0.6),
        col=1:8,
        legend = c("CP1","CP2","CP3","CP4","CP5","CP6"),
        main = "Explanation percentage of components - Isa")


#############################################################

#########################
# Building the PCA Plot #
#########################

par(mai=c(1.3, 1.1, .2, 0.2))
plot(CP$eigenvect[ , 2] ~ CP$eigenvect[ , 1],col=rgb(0,0,150,50,maxColorValue=255),las=1,
     pch=19,xlab="PC1",ylab="PC2", main = "Distribution - ISA and 1000g")

par(mai=c(1.3, 1.1, .2, 0.2))
plot(CP2$eigenvect[ , 2] ~ CP2$eigenvect[ , 1],col=rgb(0,0,155,155,maxColorValue=255),las=1,
     pch=19,xlab="PC1",ylab="PC2", main ="Distribution - ISA" )


##############################################################

#######################################################
# Plot of the relationships of the first 6 components #
#######################################################

names.PCA <- paste("PC", 1:6, "\n", format(pc.percent[1:6]), "%", sep="")
pairs(CP$eigenvec[ , 1:6], labels=names.PCA, 
      main = "Relationships of the first 6 components")

names.PCA2 <- paste("PC", 1:6, "\n", format(pc.percent2[1:6]), "%", sep="")
pairs(CP2$eigenvec[ , 1:6], labels=names.PCA2, 
      main = "Relationships of the first 6 components - ISA")

##############################################################

#############################################
# Correlation between eigenvectors and SNPs #
#     Obtaining the chromosome index        #
#############################################

# 707 and 1000g

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(CP, genofile, eig.which=1:6)

# CP1
plot(abs(CORR$snpcorr[1,]), ylim=c(0,1.0), xlab="SNP Index",
     ylab=paste("PC", 1), main="Correlation CP1 - chromosomes 1 a 22", 
     col=chr, pch="+")

# CP2
plot(abs(CORR$snpcorr[2,]), ylim=c(0,1.0), xlab="SNP Index",
     ylab=paste("PC", 2), main="Correlation CP2 - chromosomes 1 a 22", 
     col=chr, pch="+")

# ISA after LD filter

chr2 <- read.gdsn(index.gdsn(genofile2, "snp.chromosome"))
CORR2 <- snpgdsPCACorr(CP2, genofile2, eig.which=1:6)

# CP1
plot(abs(CORR2$snpcorr[1,]), ylim=c(0,1), xlab="SNP Index",
     ylab=paste("PC", 1), main="Correlation CP1 - chromosomes 1 a 22 - ISA", 
     col=chr2, pch="+")

# CP2
plot(abs(CORR2$snpcorr[2,]), ylim=c(0,1), xlab="SNP Index",
     ylab=paste("PC", 2), main="Correlation CP2 - chromosomes 1 a 22 - ISA", 
     col=chr2, pch="+")

##############################################################

#############################
# Plot Ancestry Proportions #
#############################

### 707 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group <- read.delim("isamerge1000g_2.psam")

tab <- data.frame(sample.id = samp.id, group = factor(group$SuperPop),
                  CP1 = CP$eigenvect[,1],    
                  CP2 = CP$eigenvect[,2],    
                  stringsAsFactors = FALSE)

plot(tab$CP1, tab$CP2, col=as.integer(tab$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations")
legend("topright", legend=levels(tab$group), col=1:6, cex = 0.5, pch = 16)

##############################################################

#############################################
# Plot Ancestry Proportions - Declared Race #
#############################################

### 707 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group2 <- read.delim("isamerge1000g-race_2.psam")

tab2 <- data.frame(sample.id = samp.id, group = factor(group2$SuperPop),
                   CP1 = CP$eigenvect[,1],    
                   CP2 = CP$eigenvect[,2],    
                   stringsAsFactors = FALSE)
summary(tab2)

# All colors

plot(tab2$CP1, tab2$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "darkviolet")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race")
legend("topright", legend=levels(tab2$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'darkviolet'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# ISA colored

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Without Brancos

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "gray","green2","green3","green4", "darkgreen",
                               "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "gray","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

##############################################################

##########################################
# Plot Ancestry Proportions - birthplace #
##########################################

### 707 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group3 <- read.delim("isamerge1000g-nascimento_2.psam")

tab3 <- data.frame(sample.id = samp.id, group = factor(group3$SuperPop),
                   CP1 = CP$eigenvect[,1],    
                   CP2 = CP$eigenvect[,2],    
                   stringsAsFactors = FALSE)

# All colors

plot(tab3$CP1, tab3$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", 
                               "darkviolet")[tab3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab3$group],
     main = "Ancestry Proportions - Local de Nascimento")
legend("topright", legend=levels(tab3$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'darkviolet'), 
       cex = 0.5, pch = c(16,16,16,16,8,17,16,15,4,16))

# ISA colored

plot(tab3$CP1, tab3$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", 
                               "gray")[tab3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab3$group],
     main = "Ancestry Proportions - Local de Nascimento")
legend("topright", legend=levels(tab3$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'gray'), 
       cex = 0.5, pch = c(16,16,16,16,8,17,16,15,4,16))

##############################################################

##########################################
# Plot Ancestry Proportions - Individual #
##########################################

### 707 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
groups <- list(AMR = samp.id[group$SuperPop == "AMR"],
               AFR = samp.id[group$SuperPop == "AFR"],
               EUR = samp.id[group$SuperPop == "EUR"])

prop <- snpgdsAdmixProp(CP, groups=groups,bound=TRUE)
prop_sort <- sort(prop)
propt <- t(prop)


### Save Admix Props on file

write.table(propt, file="admix-props_707samples.txt", col.names=TRUE,
            row.names=TRUE, sep="\t")

## All samples

BRpropt <- t(prop[tab$group=='BR - ISA',])
barplot(BRpropt, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",space = FALSE,
        ylim=c(0,1))
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Declared Race - Brancos

Branca <- (prop[tab2$group=='Isa - Branca',])
Brancat <- t(na.omit(Branca))

barplot(Brancat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Brancos-ISA", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

# Declared Race - Pardos

Parda <- (prop[tab2$group=='Isa - Parda',])
Pardat <- t(na.omit(Parda))

barplot(Pardat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Pardos-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

# Declared Race - Negros

Preta <- (prop[tab2$group=='Isa - Preta',])
Pretat <- t(na.omit(Preta))

barplot(Pretat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Pretos-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

###############################################################################

### Suggestion to organize the order of the figures

## All samples

# Test 1

BRpropt <- t(prop[tab$group=='BR - ISA',])
TBRpropt <- data.frame(t(BRpropt))
T1BRpropt_ord  <- TBRpropt[order(TBRpropt$AFR),]
T3BRpropt_ord  <- T1BRpropt_ord[order(T1BRpropt_ord$AMR),]
BRpropt_ord <- t(T3BRpropt_ord)

barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 2

BRpropt <- t(prop[tab$group=='BR - ISA',])
TBRpropt <- data.frame(t(BRpropt))
T1BRpropt_ord  <- TBRpropt[order(TBRpropt$AFR),]
T2BRpropt_ord  <- T1BRpropt_ord[order(T1BRpropt_ord$AMR),]
BRpropt_ord <- t(T2BRpropt_ord)

barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 3

BRpropt <- t(prop[tab$group=='BR - ISA',])
TBRpropt <- data.frame(t(BRpropt))
T2BRpropt_ord  <- TBRpropt_ord[order(TBRpropt_ord$EUR),]
T3BRpropt_ord  <- T2BRpropt_ord[order(T2BRpropt_ord$AMR),]
BRpropt_ord <- t(T3BRpropt_ord)

barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

## Declared Race - Brancos - Test 1

Branca <- (prop[tab2$group=='Isa - Branca',])
Brancat <- t(na.omit(Branca))
TBrancat <- data.frame(t(Brancat))
T1Brancat_ord  <- TBrancat[order(TBrancat$AFR),]
T2Brancat_ord  <- T1Brancat_ord[order(T1Brancat_ord$AMR),]
Brancapropt_ord <- t(T2Brancat_ord)

barplot(Brancapropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Brancos-ISA", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

## Declared Race - Pardos - Test 1

Parda <- (prop[tab2$group=='Isa - Parda',])
Pardat <- t(na.omit(Parda))
TPardat <- data.frame(t(Pardat))
T1Pardat_ord  <- TPardat[order(TPardat$AFR),]
T2Pardat_ord  <- T1Pardat_ord[order(T1Pardat_ord$AMR),]
Pardat_ord <- t(T2Pardat_ord)

barplot(Pardat_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Pardos-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

## Declared Race - Negros - Test 1

Preta <- (prop[tab2$group=='Isa - Preta',])
Pretat <- t(na.omit(Preta))
TPretat <- data.frame(t(Pretat))
T1Pretat_ord  <- TPretat[order(TPretat$AFR),]
T2Pretat_ord  <- T1Pretat_ord[order(T1Pretat_ord$AMR),]
Pretat_ord <- t(T2Pretat_ord)

barplot(Pretat_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Pretos-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

###############################################################################

################################
# Get the most correlated SNPs #
################################

## 707 ISA samples and 1000genomes

CORR_data_cp1 <- data.frame(CORR$snp.id,CORR$snpcorr[1,])
CORR_data_cp1_ord  <- CORR_data_cp1[order(-CORR_data_cp1$CORR.snpcorr.1...),]
head(CORR_data_cp1_ord)
write.table(CORR_data_cp1_ord, file="CORR_data_cp1_ord_707samplesAND1000g.csv", sep=","
            ,quote = FALSE)

CORR_data_cp2 <- data.frame(CORR$snp.id,CORR$snpcorr[2,])
CORR_data_cp2_ord  <- CORR_data_cp2[order(-CORR_data_cp2$CORR.snpcorr.2...),]
head(CORR_data_cp2_ord)
write.table(CORR_data_cp2_ord, file="CORR_data_cp2_ord_707samplesAND1000g.csv", sep=","
            ,quote = FALSE)

## Only 707 ISA samples

CORR2_data_cp1 <- data.frame(CORR2$snp.id,CORR2$snpcorr[1,])
CORR2_data_cp1_ord  <- CORR2_data_cp1[order(-CORR2_data_cp1$CORR2.snpcorr.1...),]
head(CORR2_data_cp1_ord)
write.table(CORR2_data_cp1_ord, file="CORR2_data_cp1_ord_707samples.csv", sep=","
            ,quote = FALSE)

CORR2_data_cp2 <- data.frame(CORR2$snp.id,CORR2$snpcorr[2,])
CORR2_data_cp2_ord  <- CORR2_data_cp2[order(-CORR2_data_cp2$CORR2.snpcorr.2...),]
head(CORR2_data_cp2_ord)
write.table(CORR2_data_cp2_ord, file="CORR2_data_cp2_ord_707samples.csv", sep=","
            ,quote = FALSE)

###############################################################################

###################################################
#              Joint Related Samples              #
# Calculate sample eigenvectors from SNP loadings #
###################################################

### Loading SNPs from 707 samples + 1000g  

SnpLoad <- snpgdsPCASNPLoading(CP, genofile)
names(SnpLoad)
dim(SnpLoad$snploading)

### PCA Loading - Related Samples

sample.id <- read.gdsn(index.gdsn(genofile4, "sample.id"))
samp_load <- snpgdsPCASampLoading(SnpLoad, genofile4, sample.id=sample.id)

# Check CP

diff <- CP$eigenvect[1:100,] - samp_load$eigenvect
summary(c(diff))

### Joint PCA Class

isamerge1000gfull <- list(
  sample.id = c(CP$sample.id, samp_load$sample.id),
  snp.id = CP$snp.id,
  eigenval = c(CP$eigenval, samp_load$eigenval),
  eigenvect = rbind(CP$eigenvect, samp_load$eigenvect),
  varprop = c(CP$varprop, samp_load$varprop),
  TraceXTX = CP$TraceXTX
)
class(isamerge1000gfull) <- "snpgdsPCAClass"

### New Database - All samples

isamerge1000gfull

###############################################################################

###########################################
# Plot Ancestry Proportions - All Samples #
###########################################

### Related Samples Colored

samp.id_all <- isamerge1000gfull$sample.id
group_all <- read.delim("isamerge1000g_all_relacionados.txt")

tab_all <- data.frame(sample.id = samp.id_all, group = factor(group_all$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)

plot(tab_all$CP1, tab_all$CP2, col=as.integer(tab_all$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations - All samples + 1000g")
legend("bottomright", legend=levels(tab_all$group), col=1:6, cex = 0.55, pch = 16)

### All Samples

group_all1 <- read.delim("isamerge1000g_all.txt")
tab_all1 <- data.frame(sample.id = samp.id_all, 
                       group = factor(group_all1$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)

plot(tab_all1$CP1, tab_all1$CP2, col=as.integer(tab_all1$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations -  All Samples")
legend("bottomright", legend=levels(tab_all1$group), col=1:6, cex = 0.55, pch = 16)

###############################################################################

###########################################################
# Plot Ancestry Proportions - Declared Race - All Samples #
###########################################################

group_all2 <- read.delim("isamerge1000g_all_race.txt")
tail(group_all2)

tab_all2 <- data.frame(sample.id = samp.id_all, group = factor(group_all2$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)
summary(tab_all2)


plot(tab_all2$CP1, tab_all2$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "darkviolet")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g")
legend("bottomright", legend=levels(tab_all2$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'darkviolet'), 
       cex = 0.42, pch = c(16,16,16,16,17,8,16,19,15,3,16))
    
# Only ISA colored

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.42, pch = c(16,16,16,16,17,8,16,19,15,3,16))


plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "gray","green2","green3","green4", "darkgreen",
                               "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "gray","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.42, pch = c(16,16,16,16,17,8,16,19,15,3,16))

###############################################################################

########################################################
# Plot Ancestry Proportions - birthplace - All Samples #
########################################################

group_all3 <- read.delim("isamerge1000g_all_nascimento.txt")

tab_all3 <- data.frame(sample.id = samp.id_all, group = factor(group_all3$SuperPop),
                   CP1 = isamerge1000gfull$eigenvect[,1],    
                   CP2 = isamerge1000gfull$eigenvect[,2],    
                   stringsAsFactors = FALSE)

plot(tab_all3$CP1, tab_all3$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", 
                               "darkviolet")[tab_all3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab_all3$group],
     main = "Ancestry Proportions - Local de Nascimento - All samples + 1000g")
legend("bottomright", legend=levels(tab_all3$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'darkviolet'), 
       cex = 0.40, pch = c(16,16,16,16,8,17,16,15,4,16))

# Only ISA colored

plot(tab_all3$CP1, tab_all3$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", 
                               "gray")[tab_all3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab_all3$group],
     main = "Ancestry Proportions - Local de Nascimento - All samples + 1000g")
legend("bottomright", legend=levels(tab_all3$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'gray'), 
       cex = 0.40, pch = c(16,16,16,16,8,17,16,15,4,16))

###############################################################################

########################################################
# Plot Ancestry Proportions - Individual - All Samples #
########################################################

groups4 <- list(AMR = samp.id_all[group$SuperPop == "AMR"],
               AFR = samp.id_all[group$SuperPop == "AFR"],
               EUR = samp.id_all[group$SuperPop == "EUR"])

prop4 <- snpgdsAdmixProp(isamerge1000gfull, groups=groups4,bound=TRUE)
prop_sort4 <- sort(prop4)
propt4 <- t(prop4)
head(prop4)

### Save Admix Prop - All samples

write.table(prop4, file="admix-props_all.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote = FALSE)

### All Samples

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
barplot(BRpropt4, col=c("red","green","blue"),xlab="Individual", 
        ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of BR-ISA - All samples + 1000g",space = FALSE,
        ylim=c(0,1))
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

### Declared Race - Brancos - All Samples

Branca4 <- (prop4[tab_all2$group=='Isa - Branca',])
Brancat4 <- t(na.omit(Branca4))

barplot(Brancat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Brancos-ISA - All samples + 1000g", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

### Declared Race - Pardos - All Samples

Parda4 <- (prop4[tab_all2$group=='Isa - Parda',])
Pardat4 <- t(na.omit(Parda4))

barplot(Pardat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Pardos-ISA - All samples + 1000g",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

### Declared Race - Negros - All Samples

Preta4 <- (prop4[tab_all2$group=='Isa - Preta',])
Pretat4 <- t(na.omit(Preta4))

barplot(Pretat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Pretos-ISA  - All samples + 1000g",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

###############################################################################

###############################
# Plot Simplex Representation #
###############################

### Using Prop4 Data frame - All Proportions 

prop4_frame <- data.frame(prop4)
prop4_frame2 <- tibble::rownames_to_column (prop4_frame, "Samples")
simplex_data <- dplyr::filter(prop4_frame2, 
                              grepl('a550778',prop4_frame2$Samples))

# Statistics 

dim(simplex_data)
colMeans(simplex_data[,2:4])
cov(simplex_data[,2:4])
sqrt(var(simplex_data[,2]))
sqrt(var(simplex_data[,3]))
sqrt(var(simplex_data[,4]))

# 5 trinominal on Simplex Representation

scatterplot3d(simplex_data[,2:4], 
              main ="Representacao das 5 trinomiais no simplex",
              angle = 55,color = "darkblue",pch = 16)

# Color by Group

group_all2_2 <- read.delim("isamerge1000g_all_race_2.txt")
simplex_data_2 <- cbind.data.frame(prop4_frame2,group_all2_2$SuperPop)
head(simplex_data_2)
table(simplex_data_2$`group_all2_2$SuperPop`)

colors <- c("black","deeppink3","blue","royalblue","seagreen","limegreen",
      "green2","green3","green4", "darkgreen","olivedrab4")
colors <- colors[as.factor(simplex_data_2$`group_all2_2$SuperPop`)]
s3d <- scatterplot3d(simplex_data_2[,2:4], color=colors,
                     main = "Representacao Simplex 
                     Proporcoes de Ancestralidade por Raca Declarada - All Samples",
                     type="h", pch=10,  angle = 120)
legend("topright",
       legend=levels(as.factor(simplex_data_2$`group_all2_2$SuperPop`)),
       col = c("black","deeppink3","blue","royalblue","seagreen","limegreen",
               "green2","green3","green4", "darkgreen","olivedrab4"), 
       pch = 16,cex = 0.55,inset = 0.06)

write.csv(simplex_data_2, file = "admix_prop_simplex.csv")

# Correspondence Analysis

fit.ca <- ca(simplex_data_2[,2:4])