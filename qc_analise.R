#####################################################
#################  CNV status #######################
#####################################################

setwd('C:/Users/Camila/OneDrive/Área de Trabalho/Results/')
library(ggplot2)
library(dplyr)
library("rcompanion") 

qc_file <- read.delim("all_samples_new_ref.txt")
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

source("http://bioconductor.org/biocLite.R")
library("CNTools")
library(stargazer)

samples = LETTERS[1:6]

a = c("A",5,10,25,25-10,5)
a = rbind(a,c("A",5,40,55,15, 2))
b = c("B",5,15,35,20,5)
c = c("C",5,5,45,40,5)
e = c("E",5,5,20,15,6)
f = c("F",5,45,55,10,1)


aux = as.data.frame(rbind(a,b,c,e,f))
aux
colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
aux$ID = as.character(aux$ID)
aux$chrom = as.numeric(as.character(aux$chrom))
aux$loc.start = as.numeric(as.character(aux$loc.start))
aux$loc.end = as.numeric(as.character(aux$loc.end))
aux$num.mark = as.numeric(as.character(aux$num.mark))
aux$seg.mean = as.numeric(as.character(aux$seg.mean))

aux

stargazer(aux, rownames = F, summary = F)

attach(aux)

n=length(samples)
aux <- data.frame(c(aux$ID,samples), c(aux$chrom,rep(unique(aux$chrom), n)), 
                  c(aux$loc.start,rep(1, n)), c(aux$loc.end,rep(max(aux$loc.end), n)), 
                  c(aux$num.mark,rep(1, n)), c(aux$seg.mean,rep(-1, n)))
colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
aux

aux$ID = as.character(aux$ID)
str(aux)
head(aux)
dim(aux)

seg <-  CNSeg(aux)
seg
rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
cnvA = rs(rsByregion)
cnvA

dim(cnvA)
cnvA = cbind(cnvA[,1:3],apply(cnvA[,-c(1:3)],2,arruma))
cnvA

stargazer(cnvA, rownames = F, summary = F)


#######################

########## Functions ##########

library("CNTools")
install.packages("kinship2")
library("kinship2")

arruma = function(x){
  x[which(x==1)] = 0
  x[which(x==2)] = 1
  x[which(x==-1)] = 2
  x[which(x==5)] = 3
  x[which(x==6)] = 4
  return(x)
}

cont_0 = function(x){
  length(which(x[-c(1,2,3)]==0))
}

cont_1 = function(x){
  length(which(x[-c(1,2,3)]==1))
}

cont_2 = function(x){
  length(which(x[-c(1,2,3)]==2))
}

cont_3 = function(x){
  length(which(x[-c(1,2,3)]==3))
}
cont_4 = function(x){
  length(which(x[-c(1,2,3)]==4))
}


grupos = function(a){
  t = sum(a)
  g = 5
  
  if(a[1]<t*maf)
    g = g-1
  if(a[2]<t*maf)
    g = g-1
  if(a[3]<t*maf)
    g = g-1
  if(a[4]<t*maf)
    g = g-1
  if(a[5]<t*maf)
    g = g-1
  
  return(g)
}



group = function(a){
  t = sum(a)
  g = 2
  
  if(a[1]<t*maf)
    g = g-1
  if(a[2]<t*maf)
    g = g-1
  
  return(g)
}


contaMut = function(x){
  return(length(which(x!=2))-3)
}


getwd()

########## APT/PennCNV outputs ##########

### Information about CNV regions
cnv = read.table("D:/ISA/", sep="\t", header = F)
colnames(cnv) = c("Chr","Start","End","Number","Length", "State", "CN", "Sample", "First Marker", "Last Marker")
head(cnv)
summary(cnv)
dim(cnv)
attach(cnv)

length(unique(Sample))



### Information about quality control PennCNV per sample
qc = read.table("/home/cicconella/DadosMestrado/tableQC", sep="\t", header = T)
qc = qc[order(qc$File),] 
head(qc)
summary(qc)
dim(qc)

########## Cleaning Bad Samples ##########
length(which(qc$LRR_SD > 0.35))

length(which(qc$BAF_mean > 0.6))+length(which(qc$BAF_mean < 0.4))

length(which(qc$BAF_drift > 0.01))

length(which(qc$WF > 0.04))+length(which(qc$WF < -0.04))


fail1 = intersect(which(qc$LRR_SD > 0.35), which(qc$BAF_drift > 0.01))
fail2 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF > 0.04))
fail3 = intersect(which(qc$LRR_SD > 0.35), which(qc$WF < -0.04))

head(cnv)
dim(cnv)

qc = qc[-c(fail2,fail3),]
dim(qc)
head(qc)
attach(qc)

cnv = merge(cnv, qc, by.x = "Sample", by.y = "File")
cnv = cnv[,-c(11:19)]
dim(cnv)
head(cnv)
summary(cnv)
attach(cnv)
ind = unique(Sample)
length(ind)

remove(qc)
remove(fail3)
remove(fail2)
remove(fail1)
#remove(ind)



########## Description ##########

dim(cnv)
head(cnv)
summary(cnv)

attach(cnv)


##### Samples #####
Sample
print(paste("Total de amostras:", length(unique(Sample))))

x = as.numeric(table(Sample))
summary(x)
sd(as.numeric(table(Sample)))^2

m=226

length(which(x<m))/length(x)

boxplot(x[which(x<m)], main = "Distribution of Number os CNVs for each sample", col = "lightblue")
hist(x[which(x<m)])
summary(x[which(x<m)])
sd(x[which(x<m)])

shapiro.test(x[which(x<m)])

shapiro.test(x)

png("../DadosMestrado/plots/numberCNVs.png")
hist(table(Sample), nc = 1000, main = "Number of CNVs per sample", xlab = "Number of CNVs")
dev.off()

table(Sample)[which(table(Sample)>500)]
paste(length(which(table(Sample)>500))/length(unique(Sample)), "são maiores que 500")
paste(length(which(table(Sample)>250))/length(unique(Sample)), "são maiores que 250")
paste(length(which(table(Sample)>100))/length(unique(Sample)), "são maiores que 100")

porcentagens = c()

for(i in seq(0,3000, by =25)){
  porcentagens = c(porcentagens, length(which(table(Sample)<i))/length(unique(Sample)))
}

porcentagens

# png("../DadosMestrado/plots/samplesSize.png")
# plot(seq(0,4000, by =25), porcentagens, pch="", ylab="% Samples", xlab = "Number of CNVs")
# lines(seq(0,4000, by =25), porcentagens, main = "Frequency of samples by number of CNVs")
# dev.off()
# 
bla = cbind(seq(0,3000, by =25), porcentagens)

stargazer(bla, digits = 2)

head(cnv)

barplot(table(CN)/sum(table(CN)), col = "blue", xlab = "Copy Number", 
        ylab = "Frequency (%)", main = "Frequency of Copy Number",
        ylim = c(0,0.7))

table(CN)/sum(table(CN))

which(table(Sample)<101)
names(which(table(Sample)<101))

bla = cnv[cnv$Sample %in% names(which(table(Sample)<101)),]

bla[1:100,]

barplot(table(bla$CN)/sum(table(bla$CN)), col = "blue", xlab = "Copy Number", 
        ylab = "Frequency (%)", main = "Frequency of Copy Number",
        ylim = c(0,0.7))

table(bla$CN)/sum(table(bla$CN))


##### Size #####
Size = Length
summary(Size)

hist(Size, nc=100)
barplot(table(Size))

barplot(c(1,5,6))

mean(Size)
summary(Size)

data = Size
hist(log10(data), nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", col="blue", main = "Histogram of CNV Length")
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))

head(cnv)

hist(log10(cnv[cnv$CN<2,5]), 
     nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", 
     col="blue", main = "Histogram of Deletions Length")
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))

hist(log10(cnv[cnv$CN>2,5]), 
     nc=100,xaxt="n", xlab = "Length of CNV", ylab = "Absolute Frequency", 
     col="blue", main = "Histogram of Deletions Length")
axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))


head(cnv)

for(i in 1:22){
  png(paste("../DadosMestrado/plots/length", i, ".png", sep=""))
  hist(log10(cnv$Length[cnv$Chr==3]), col = rgb(0,0,1), nc = 100, xaxt="n", 
       xlab = "Length of CNV", ylab = "Absolute Frequency", main = "Histogram of CNV Length")
  axis(1, at=1:7, labels = c("10bp", "100bp", "1kb", "10kb", "100kb", "1Mb", "10Mb"))
  dev.off()
}





########## BY chromosome ##########

# Find the minimal regions for each chromosome

chromosomes <- vector("list",22)

i=1

for(i in 1:22){
  
  cnvA = cnv[Chr==i,]
  dim(cnvA)
  head(cnvA)
  aux = cnvA[order(cnvA$Start),]
  dim(aux)
  head(aux)
  
  # Change matrix to data.frame and add an auxiliar region from the start of the 
  # chromossome
  
  aux <- data.frame(c(aux$Sample,ind), c(aux$Chr,rep(i, length(ind))), 
                    c(aux$Start,rep(1, length(ind))), c(aux$End,rep(max(aux$End), length(ind))), 
                    c(aux$Number,rep(1, length(ind))), c(aux$State,rep(-1, length(ind))))
  colnames(aux) = c("ID", "chrom", "loc.start", "loc.end", "num.mark",  "seg.mean")
  str(aux)
  head(aux)
  dim(aux)
  
  seg <-  CNSeg(aux)
  seg
  rsByregion <- getRS(seg, by = "region", imput = TRUE, XY = FALSE, what = "max")
  cnvA = rs(rsByregion)
  cnvA[1:10,1:10]
  
  dim(cnvA)
  cnvA = cbind(cnvA[,1:3],apply(cnvA[,-c(1:3)],2,arruma))
  
  chromosomes[[i]] = cnvA
  
  
}

remove(cnvA)
remove(aux)
remove(rsByregion)
remove(seg)

for(k in 300:600){
  print(table(as.numeric(chromosomes[[1]][k,-c(1:3)])))
}




##### Testes #####

chr =  6
chromosomes[[chr]][1,]

contaMut(chromosomes[[chr]][2,])

plotar = cbind(chromosomes[[chr]][,1:3], apply(chromosomes[[chr]], 1, contaMut))
plotar = as.data.frame(plotar)
colnames(plotar) = c("chrom", "start", "end", "score")
#plotar = plotar[-1,]
plotar[,1] = paste("chr",chr, sep = "")

plotar[1:10,]
tail(plotar)

chrom = plotar[1,1]
chromstart = chromosomes[[chr]][2,2]
chromend = chromosomes[[chr]][nrow(chromosomes[[chr]]),3]


# png("../DadosMestrado/Sushi/chr1.png")
plotBedgraph(plotar, chrom, chromstart,
             chromend, main = "CNVs per region (1019 samples)",
             colorbycol= SushiColors(5))
# 
labelgenome(chrom, chromstart,chromend,n=4,scale="Mb")
mtext("Number of CNVs",side=2,line=2.5,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
# dev.off()

plotar[which(plotar[,4]>400),]
dim(plotar[which(plotar[,4]>400),])
summary(plotar[which(plotar[,4]>400),])

# Teste para achar coisas

bla = merge(plotar[which(plotar[,4]>400),], chromosomes[[chr]], by="start")

bla[1:10,1:10]

dim(bla)

tabela = c()

for(i in 1:dim(bla)[1]){
  tabela = rbind(tabela,table(as.numeric(bla[i,-c(1:6)])))
}

tabela

mean(apply(tabela[1:29,1:2],1,sum))/mean(apply(tabela[,c(1:2,4:5)],1,sum))

mean(apply(tabela[,c(1:2,4:5)],1,sum))

cnv[Sample==1,]


##### CNVs, onde estao? #####

table(cnv$Chr)
barplot(table(cnv$Chr), xlab = "Chromosome", names.arg=c(1:22),col="blue", cex.names  = 0.75,
        main = "CNVs detected by Chromosome", ylab = "Absolute frequency")


tamanhos = read.table("../DadosMestrado/chromo", header = T, sep = "\t")
head(tamanhos)
tamanhos = tamanhos[,-c(3,4)]

prop = c()
prop2 = c()

for(i in 1:22){
  
  if(i!=2){
    plotar = cbind(chromosomes[[i]][,1:3], apply(chromosomes[[i]], 1, contaMut))  
  } else{
    plotar = cbind(chromosomes[[i]][,1:3], (apply(chromosomes[[i]], 1, contaMut)+1))    
  }
  
  
  plotar = as.data.frame(plotar)
  colnames(plotar) = c("chrom", "start", "end", "score")
  
  plotar[,1] = paste("chr",i, sep="")
  
  chrom = plotar[1,1]
  chromstart = chromosomes[[i]][1,2]
  chromend = chromosomes[[i]][nrow(chromosomes[[i]]),3]
  
  png(paste("../DadosMestrado/Sushi/chr",i,".png", sep=""))
  plotBedgraph(plotar, chrom, chromstart,
               chromend, main = "CNVs per region (1019 samples)",
               colorbycol= SushiColors(5))
  
  labelgenome(chrom, chromstart,chromend,n=4,scale="Mb")
  mtext("Number of CNVs",side=2,line=2.5,cex=1,font=2)
  axis(side=2,las=2,tcl=.2)
  dev.off()
  
  print(i)
  
  
  bla = plotar$end-plotar$start+1
  bla = sum(bla[which(plotar$score!=0)])
  
  prop = c(prop, bla/tamanhos$Total.length..bp.[i])
  
  bla = plotar$end-plotar$start+1
  bla = sum(bla[which(plotar$score>20)])
  
  prop2 = c(prop2, bla/tamanhos$Total.length..bp.[i])
  
}

png("../DadosMestrado/plots/prop.png")
barplot(prop*100, names.arg=c(1:22),col="blue", ylab = "Proportion of CNVs (%)", xlab = "Chromosome", 
        main = "Proportion of CNVs by Chromosome",cex.names  = 0.75, ylim = c(0,100))
dev.off()

png("../DadosMestrado/plots/prop2.png")
barplot(prop2*100, names.arg=c(1:22),col="blue", ylab = "Proportion of CNVs (%)", xlab = "Chromosome", 
        main = "Proportion of CNVs by Chromosome",cex.names  = 0.75, ylim=c(0,70))
dev.off()

png("../DadosMestrado/plots/props.png")
barplot(t(cbind(prop,prop2))*100, beside = T, col = c("blue", "green"), names.arg = c(1:22),
        ylab = "Proportion of CNVs (%)", xlab = "Chromosome", cex.names  = 0.75, main = "Proportion of CNVs by Chromosome")
legend("topleft",c(">1 Sample", ">20 Samples"),col = c("blue", "green"), pch = 15)
dev.off()


sizesT= c()
for(x in 1:22){
  y = dim(chromosomes[[x]])[1]
  sizesT = c(sizesT,y)
}

png("/Users/Ana/Google Drive/2016/before.png")
plot(1:22,sizesT, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions before Cleaning",pch=16)
dev.off()

remove(x)
remove(y)
remove(i)




########## Cleaning CNV Regions ##########

i=1

maf = 0.02

for(i in 1:22){
  mut_0 = apply(chromosomes[[i]],1,cont_0)
  mut_1 = apply(chromosomes[[i]],1,cont_1)
  mut_2 = apply(chromosomes[[i]],1,cont_2)
  mut_3 = apply(chromosomes[[i]],1,cont_3)
  mut_4 = apply(chromosomes[[i]],1,cont_4)
  
  mutations = cbind(mut_0,mut_1,mut_2,mut_3,mut_4)
  
  dim(chromosomes[[i]])
  dim(mutations)
  
  groupCNV = apply(mutations, 1, grupos)
  
  summary(as.factor(groupCNV))
  
  exclude = which(groupCNV==1)
  chromosomes[[i]] = chromosomes[[i]][-exclude,]
  dim(chromosomes[[i]])
  head(chromosomes[[i]][,1:15])
  print(i)
}

sizes = c()
for(i in 1:22){
  y = dim(chromosomes[[i]])[1]
  sizes = c(sizes,y)
}


########## Cleaning CNV Regions - Binary##########

maf = 0.05

chrBin = vector("list",22)

for(i in 1:22){
  chrBin[[i]] = chromosomes[[i]]
  chrBin[[i]][which(chrBin[[i]]!=2)] = 1
  chrBin[[i]][which(chrBin[[i]]==2)] = 0
  
  chrBin[[i]][,c(1,2,3)] = chromosomes[[i]][,c(1,2,3)]
}

head(chromosomes[[1]][,1:10])
head(chrBin[[1]][,1:10])


for(i in 1:22){
  mut_0 = apply(chrBin[[i]],1,cont_0)
  mut_1 = apply(chrBin[[i]],1,cont_1)
  
  mutations = cbind(mut_0,mut_1)
  
  dim(chrBin[[i]])
  dim(mutations)
  
  groupCNV = apply(mutations, 1, group)
  
  summary(as.factor(groupCNV))
  
  exclude = which(groupCNV==1)
  chrBin[[i]] = chrBin[[i]][-exclude,]
  dim(chrBin[[i]])
  head(chrBin[[i]][,1:15])
  print(i)
}

sizesB = c()
for(x in 1:22){
  y = dim(chrBin[[x]])[1]
  sizesB = c(sizesB,y)
}

png("/Users/Ana/Google Drive/2016/after-bin-maf-0.02.png")
plot(1:22,sizesB, xlab="Chromosomes", ylab="Number of CNV Regions", 
     main = "CNV Regions after Cleaning",pch=16)
dev.off()

png("/Users/Ana/Google Drive/2016/befofe-after-bin-0.02.png")
plot(sizesB,sizesT, xlab="Number of Filtered Regions", ylab = "Total of Regions",
     xlim=c(-20,2400),ylim=c(-20,11000),pch=16)
identify(sizesB,sizesT)
dev.off()

remove(mutations)
remove(exclude)
remove(groupCNV)
remove(i)
remove(maf)
remove(mut_0)
remove(mut_1)
remove(x)
remove(y)plotar[,2]>125000000

########## Individuals ##########

# Sample Data

info = read.table("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/Project/dados2", header = T, sep = ",")

head(info)  
summary(info)
dim(info)
str(info)
attach(info)

# Association between sample id and celfiles

ind = read.table("/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/Mestrado/Project/individuos")
colnames(ind) = c("cel", "IID")

class(ind)
head(ind)
dim(ind)
summary(ind)
str(ind)

# Associacao 

cel = IID 

for(i in 1:length(cel)){
  
  aux = which(ind$IID==cel[i])
  if(length(aux)!=0)
    cel[i] = max(ind$cel[aux])
  else
    cel[i] = NA
}

cel[1:10]

#  Info + cel

info = cbind(info,cel)

head(info)  
summary(info)
dim(info)
str(info)


# Correct CEL Files
info[which(info$IID==4919),9]=1444
info[which(info$IID==15908),9]=1455
info[which(info$IID==15911),9]=2113
info[which(info$IID==15921),9]=482
info[which(info$IID==32608),9]=1457
info[which(info$IID==96502),9]=2360

attach(info)

head(info)

positions = rep(0,nrow(info))
names = as.numeric(colnames(chromosomes[[1]])[-c(1,2,3)])

for (i in 1:length(IID)){
  k = which(info$cel[i]==names)
  if(length(k)==1)
    positions[i] = k
  else
    positions[i] = NA
}

remove(ind)
remove(aux)
remove(i)
remove(names)
remove(k)
remove(cel)

##### Genotipo e fenotipo #####

chromosomes[[1]][1:10,1:10]
table(as.numeric(chromosomes[[1]][4,-c(1:3)]))

head(info)

info[which(info$cel==1),]

chromosomes[[1]][1:10,1:10]

head(info)

chromosomes[[1]][1,1:50]

dim(chromosomes[[1]])

seq(1,1239, by=10)

## Seleciona o genotipo

for(c in seq(1,1239, by=10)){
  genotipo = as.numeric(chromosomes[[1]][c,-c(1:3)])
  
  genotipo = (cbind(colnames(chromosomes[[1]])[-c(1:3)], genotipo))
  
  colnames(genotipo) = c("celfiles", "cnv")
  genotipo = as.data.frame(genotipo)
  head(genotipo)
  
  head(info)
  
  final = merge(info, genotipo, by.x = "cel", by.y = "celfiles", all.x = T)
  final[1:100,]
  tail(final)
  
  head(final)
  
  comb = expand.grid(c(0:4),c(0:4),c(0:4))
  comb = data.frame(cbind(comb, rep(0,125)))
  comb[,4] = as.numeric(as.character(comb[,4]))
  
  colnames(comb) = c("P1", "P2", "OF", "CN")
  
  head(comb)
  head(final)
  
  for(i in 1:nrow(final)){
    final[i,]
    
    if(is.na(as.character(final[i,10]))){
      #print("next")
      next
    }
    if(final[i,4]!=0 & final[i,5]!=0){
      if(is.na(as.character(final[which(final[,3]==final[i,4]),10]))
         | is.na(as.character(final[which(final[,3]==final[i,5]),10]))){
        next
      }else{
        a = t(as.matrix(c(as.numeric(as.character(final[which(final[,3]==final[i,4]),10])),
                          as.numeric(as.character(final[which(final[,3]==final[i,5]),10])),
                          as.numeric(as.character(final[i,10])))))
        colnames(a) = colnames(comb)[-4]
        comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4] = comb[which(apply(comb, 1, function(x) identical(x[1:3], a[1,]))),4]+1
      }
    }else{
      #print("orfao")
    }
  }
  
  #png(paste(getwd(),"/trios.png", sep = ""), width = 1400, height = 460)
  #barplot(comb[,4],pch=16, names.arg = apply(comb[,1:3],1,paste,collapse = ""), las=2)
  #dev.off()
  
  for(i in 1:75){
    a = t(as.matrix(as.numeric(c(rev(comb[i, 1:2]),comb[i,3]))))
    colnames(a) = c("P1", "P2","OF")
    bla = which(apply(comb, 1, function(x) identical(x[1:3], a[1,])))
    if(bla!=i){
      comb[i,4] = comb[i,4]+comb[bla,4]
      comb = comb[-bla,]
    }
  }
  
  dim(comb)
  
  #png(paste(getwd(),"/trios2.png", sep = ""), width = 1400, height = 460)
  #barplot(comb[,4],pch=16, names.arg = apply(comb[,1:3],1,paste,collapse = ""), las=2)
  #dev.off()
  
  head(comb)
  
  if(c == 1){
    trios = comb
  }else{
    trios = cbind(trios,comb[,4])  
  }
  print(c)
}

head(trios)
trios[,1:20]

combinacoes = apply(trios[,-(1:3)], 1,sum)
names(combinacoes) = apply(trios[,1:3],1,paste,collapse="")
combinacoes_limpo = combinacoes[-which(names(combinacoes)=="222")] 

barplot(combinacoes_limpo, las=2, xlab = "Genotype", ylab = "Absolute frequency",
        col = "blue", main = "CNV Occurances in Trios")


trios[,1:20]
colnames(trios) = c(colnames(trios)[1:4], rep("CN", 123))

trios[,1:15]

x = cbind(apply(trios[,1:3],1,paste,collapse=""),trios[,4:15])
colnames(x) = c("CNs", colnames(x)[-1])

stargazer(x, summary = F)

?stargazer
##### Para ver o heredograma #####

kin = final[,c(2,3,4,5,6,10)]

head(kin)  

class(kin$cnv)

kin$cnv = as.numeric(as.character(kin$cnv))


kin$cnv[which(kin$cnv!=2)] = 1
kin$cnv[which(kin$cnv==2)] = 0
kin$cnv[is.na(kin$cnv)] = 0

head(kin)

kin <- with(kin, pedigree(kin$IID, kin$PAT, kin$MAT, kin$SEX, famid = kin$FID, affected = kin$cnv))
kin

round(8*kinship(kin))/4

unique(kin$famid)

kin1 = kin['110']
kin1
plot(kin1)

kin = info[1:48,c(1,2,3,4,5,8)]

head(kin)  

class(kin$altura)

kin$altura[which(kin$altura<161)] = 0
kin$altura[which(kin$altura>160)] = 1

class(kin$altura)
table(kin$altura)


kin <- with(kin, pedigree(kin$IID, kin$PAT, kin$MAT, kin$SEX, famid = kin$FID, affected = kin$altura))
kin

round(8*kinship(kin))/4

kin1 = kin['2']
kin1
plot(kin1)

##### Comparando os dados #####


head(info)

getwd()

dados = read.table()

########## Getting the files ped and phen ##########

# .ped

ped = cbind(info$IID,info$PAT,info$MAT,info$SEX,info$FID)
colnames(ped) = c("id","fa","mo","sex","fid")

head(ped)

write.table(ped, "/Users/Ana/Google Drive/2016/files/samples.ped", row.names = F, quote = F, sep = ",")

# .phen

phen = vector("list",22)

j=1

for(j in 1:22){
  fen = c()
  for(i in 1:nrow(chrBin[[j]])){  
    print(c("CHR ", j, "Reg ",i))
    cn = chrBin[[j]][i,-c(1,2,3)]
    cn = cn[positions]
    
    fen = cbind(fen,cn)
  }
  phen[[j]] = cbind(info,fen)
}

dim(phen[[1]])
dim(chrBin[[1]])

for(i in 1:22){
  
  phenotypes = phen[[i]][,-c(1,3,4,6,9)]
  
  name = paste("/Users/Ana/Google Drive/2016/files/phen",i,".phen",sep="")
  write.table(phenotypes,name,row.names = F, quote = F, sep = ",")
}

remove(fen)
remove(info)
remove(ped)
remove(phenotypes)
remove(chrBin) #ja esta organizado nos arquivos .phen
remove(cn)
remove(i)
remove(j)
remove(name)
remove(positions)
remove(sizesB)
remove(sizesT)
remove(phen)

########## Correctins files to Solar##########

for(i in 2:22){
  
  name = paste("/Users/Ana/Google Drive/2016/files/phen",i,".phen",sep="")
  chr = read.table(name, header = T, sep=",")
  
  chr[1:10,1:10]
  dim(chr)
  
  colnames(chr)[1:3] = c("id", "sexo", "idade")
  
  for(i in 1:nrow(chr)){
    chr[i,which(is.na(chr[i,]))] = ""
  }
  
  write.table(chr, name, row.names = F, quote = F, sep = ",")
}


rm(name)

# comb = expand.grid(c(0:4),c(0:4),c(0:4))
# comb = apply(comb, 1, paste, collapse="")
# comb = data.frame(cbind(comb, rep(0,125)))
# 
# comb[,2] = as.numeric(as.character(comb[,2]))
# 
# head(comb)
# 
# barplot(comb[1:50,2],pch=16)
# barplot(comb[51:100,2],pch=16)
# barplot(comb[101:125,2],pch=16)
# 
# head(final)
# 
# i=82
# 
# 
# for(i in 1:nrow(final)){
#   final[i,]
#   
#   if(is.na(as.character(final[i,10]))){
#     #print("next")
#     next
#   }
#   if(final[i,4]!=0 & final[i,5]!=0){
#     if(is.na(as.character(final[which(final[,3]==final[i,4]),10]))
#        | is.na(as.character(final[which(final[,3]==final[i,5]),10]))){
#       #print("nah")
#     }else{
#       a = paste(c(as.character(final[which(final[,3]==final[i,4]),10]),
#                     as.character(final[which(final[,3]==final[i,5]),10]),
#                     as.character(final[i,10])), collapse ="")
#       comb[which(comb[,1]==a),2] = comb[which(comb[,1]==a),2] + 1
#     }
#   }else{
#     #print("orfao")
#   }
# }
# comb
# 
# png(paste(getwd(),"/trios.png", sep = ""), width = 1400, height = 460)
# barplot(comb[,2],pch=16, names.arg = comb[,1], las=2)
# dev.off()
# 
# comb[1:20,]
