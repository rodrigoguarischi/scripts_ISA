
library(CNTools)
library(stargazer)
library(CNVRanger)

setwd('C:/Users/Camila/OneDrive/Área de Trabalho/Análises/Results_novo')

dta_raw <- read.delim('cnv_data.txt')
summary(dta_raw)

dta__ <- dplyr::filter(dta_raw,!grepl('X',dta_raw$Chromosome))

dta <- na.omit(dta__)
table(dta$Chromosome)

aux_f <- data.frame(dta$Sample.Filename,dta$Chromosome,
                  dta$Start.Position,dta$Stop.Position,dta$Marker.Count,
                  dta$MedianLog2Ratio)

write.table(aux_f, "C:/Users/Camila/OneDrive/Área de Trabalho/TESTE.txt", quote=FALSE,
            col.names = TRUE, row.names = FALSE)

aux <- read.delim("seg_chr1.txt")

samples <- read.delim("samples.txt")

aux$ID = as.character(aux$dta.Sample.Filename)
aux$chrom = as.numeric(as.character(aux$dta.Chromosome))
aux$loc.start = as.numeric(as.character(aux$dta.Start.Position))
aux$loc.end = as.numeric(as.character(aux$dta.Stop.Position))
aux$num.mark = as.numeric(as.character(aux$dta.Marker.Count))
aux$seg.mean = as.numeric(as.character(aux$dta.MedianLog2Ratio))

stargazer(aux, rownames = F, summary = F)

attach(aux)

n=length(samples)
aux <- data.frame(c(aux$ID,samples), 
                  c(aux$chrom,rep(unique(aux$chrom), n)), 
                  c(aux$loc.start,rep(1, n)), 
                  c(aux$loc.end,rep(max(aux$loc.end), n)), 
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

arruma = function(x){
  x[which(x==1)] = 0
  x[which(x==2)] = 1
  x[which(x==-1)] = 2
  x[which(x==5)] = 3
  x[which(x==6)] = 4
  return(x)
}

dim(cnvA)
cnvA = cbind(cnvA[,1:3],apply(cnvA[,-c(1:3)],2,arruma))
cnvA

stargazer(cnvA, rownames = F, summary = F)

#The following R code exemplifies the procedure to remove CNV with no or low mutations:
#Load the functions that count the number of samples of each group (CN = 0; :::; 4)
#and the function that counts the number of valid groups:

cont_0 = function ( x ) {
  length (which( x[c ( 1 , 2 , 3 ) ]==0) )
}

cont_1 = function ( x ) {
  length (which( x[c ( 1 , 2 , 3 ) ]==1) )
}

cont_2 = function ( x ) {
  length (which( x[c ( 1 , 2 , 3 ) ]==2) )
}

cont_3 = function ( x ) {
  length (which( x[c ( 1 , 2 , 3 ) ]==3) )
}

cont_4 = function ( x ) {
  length (which( x[c ( 1 , 2 , 3 ) ]==4) )
}

grupos = function ( a ) {
  t = sum( a )
  g = 5
  
  if ( a[1]<t*mcf )
    g = g-1
  
  if ( a [2]<t*mcf )
    g = g-1
  
  if ( a [3]<t*mcf )
    g = g-1
  
  if ( a [4]<t*mcf )
    g = g-1
  
  if ( a [5]<t*mcf )
    g = g-1
  
  return ( g )
}

#Load the dataset (we use the example from Appendix C) and define the mcf:
#> cnv
#chrom start end A B C D E F
#1 5 1 5 2 2 2 2 2 2
# 2 5 6 10 2 2 3 2 4 2
# 3 5 11 15 3 2 3 2 4 2
# 4 5 16 20 3 3 3 2 4 2
# 5 5 21 25 3 3 3 2 2 2
# 6 5 26 35 2 3 3 2 2 2
# 7 5 36 40 2 2 3 2 2 2
# 8 5 41 45 1 2 3 2 2 2
# 9 5 46 55 1 2 2 2 2 0
#
mcf = 0.02

#Calculate the number of groups of each CNV and remove the ones with only one group:

mut_0 = apply( cnvA , 1 , cont_0)
mut_1 = apply( cnvA , 1 , cont_1)
mut_2 = apply( cnvA , 1 , cont_2)
mut_3 = apply( cnvA , 1 , cont_3)
mut_4 = apply( cnvA , 1 , cont_4)

mutations = cbind(mut_0 ,mut_1 ,mut_2 ,mut_3 ,mut_4)

dim( cnvA )
dim(mutations )

groupCNV = apply(mutations , 1 , grupos )

summary( as.factor (groupCNV) )

exclude = which(groupCNV==1)
cnvA = cnvA[-exclude , ]
dim( cnvA )

# > cnv
# chrom start end A B C D E F
# 2 5 6 10 2 2 3 2 4 2
# 3 5 11 15 3 2 3 2 4 2
# 4 5 16 20 3 3 3 2 4 2
# 5 5 21 25 3 3 3 2 2 2
# 6 5 26 35 2 3 3 2 2 2
# 7 5 36 40 2 2 3 2 2 2
# 8 5 41 45 1 2 3 2 2 2
# 9 5 46 55 1 2 2 2 2 0