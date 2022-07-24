
#!/bin/bash

#-------------------------#
#     Databases - ISA     #
#-------------------------#

###########
# Sources #
###########

# apt = Axiom Affymetrix
# plink = 1.9
# plink2 = 2.0


########################
## APT convert to ped ##
########################


/mnt/c/Users/Camila/OneDrive/Área\ de\ Trabalho/mestrado-isa/apt_2.11.4_linux_64_bit_x86_binaries/bin/apt-format-result \
--calls-file AxiomGT1.calls.txt --annotation-file Axiom_PMRA.na35.annot.db --snp-list-file Recommended.ps \
--snp-identifier-column 'dbSNP_RS_ID' --sample-filter-file samplenames.txt --export-plink-file isahg37.ped

/mnt/c/Users/Camila/OneDrive/Área\ de\ Trabalho/mestrado-isa/apt_2.11.4_linux_64_bit_x86_binaries/bin/apt-format-result \
--calls-file AxiomGT1.calls.txt --annotation-file Axiom_PMRA.na36.r1.a1.annot.db --snp-list-file Recommended.ps \
--snp-identifier-column 'dbSNP_RS_ID' --sample-filter-file samplenames.txt --export-plink-file isahg38.ped


#######################
## 1- Convert to bed ##
#######################

./plink --file isahg37 --no-fid --no-parents --no-sex --no-pheno --make-bed
./plink --file isahg38 --no-fid --no-parents --no-sex --no-pheno --make-bed


######################################################################################
## 2 - Filter by ChrAut, Without palindrômicas/HLA regions, MAF 0,0001, HWE 0,00001 ##
######################################################################################

# Palindromics/HLA regions  

./plink2 --bfile isahg37 --exclude palind-hla-onlyrs.txt --make-bed --out isahg37_tmp
./plink2 --bfile isahg38 --exclude palind-hla-onlyrs.txt --make-bed --out isahg38_tmp

# Filter maf, hwe and chrsex

./plink2 --bfile isahg37_tmp --autosome --make-bed --maf 0.0001 --hwe 0.000001 --out isahg37_B
./plink2 --bfile isahg38_tmp --autosome --make-bed --maf 0.0001 --hwe 0.000001 --out isahg38_B

# LD com r2 de 0.5

./plink2 --bfile isahg37_B --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.5 --out isahg37_B_tmp
./plink2 --bfile isahg37_B --extract isahg37_B_tmp.prune.in --make-bed --out isahg37_C

./plink2 --bfile isahg38_B --exclude range high-LD-Regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.5 --out isahg38_B_tmp
./plink2 --bfile isahg38_B --extract isahg38_B_tmp.prune.in --make-bed --out isahg38_C

# LD com r2 de 0.11

./plink2 --bfile isahg37_B --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.11 --out isahg37_B_tmp2
./plink2 --bfile isahg37_B --extract isahg37_B_tmp2.prune.in --make-bed --out isahg37_D

./plink2 --bfile isahg38_B --exclude range high-LD-Regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.11 --out isahg38_B_tmp2
./plink2 --bfile isahg38_B --extract isahg38_B_tmp2.prune.in --make-bed --out isahg38_D

# MAF 0,01, Eq-HW 0,000001

./plink2 --bfile isahg37_tmp --autosome --make-bed --maf 0.01 --hwe 0.000001 --out isahg37_F
./plink2 --bfile isahg38_tmp --autosome --make-bed --maf 0.01 --hwe 0.000001 --out isahg38_F


# MAF 0,01 - LD com r2 de 0.5

./plink2 --bfile isahg37_F --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.5 --out isahg37_F_tmp
./plink2 --bfile isahg37_F --extract isahg37_F_tmp.prune.in --make-bed --out isahg37_G

./plink2 --bfile isahg38_F --exclude range high-LD-Regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.5 --out isahg38_F_tmp
./plink2 --bfile isahg38_F --extract isahg38_F_tmp.prune.in --make-bed --out isahg38_G

# MAF 0,01 - LD com r2 de 0.11

./plink2 --bfile isahg37_F --exclude range high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.11 --out isahg37_F_tmp2
./plink2 --bfile isahg37_F --extract isahg37_F_tmp2.prune.in --make-bed --out isahg37_H

./plink2 --bfile isahg38_F --exclude range high-LD-Regions-hg38-GRCh38.txt --indep-pairwise 50 5 0.11 --out isahg38_F_tmp2
./plink2 --bfile isahg38_F --extract isahg38_F_tmp2.prune.in --make-bed --out isahg38_H


###############################
## 3 - Convert to Raw Format ##
###############################


for bed in $(ls *.bed | cut -f1 -d".");
do
	echo ${bed};
	./plink --bfile ${bed} --recodeAD --out ${bed}_raw --no-fid --no-parents --no-sex --no-pheno
done