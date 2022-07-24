#!/bin/bash

#----------------------------#
#     Database selected      #
#    isa_rec_rs_aut_hg37     #
#----------------------------#

###########
# Sources #
###########

# plink = 1.9
# plink2 = 2.0


#########
# Steps #
#########

# Convert Bed

./plink --file old_analysis/isa_rec_rs_aut_hg37 --no-fid --no-parents --no-sex --no-pheno --make-bed --out isa_hg37



# Filtered by Chr, MAF and HWE

./plink2 --bfile isa_hg37 --autosome --make-bed --maf 0.0001 --hwe 0.000001 --out isa_hg37_filter



# Filtered by Palind,HLA

./plink2 --bfile isa_hg37_filter --exclude old_analysis/palind-hla-onlyrs.txt --make-bed --out isa_hg37_filter2



# Filtered by LD

./plink2 --bfile isa_hg37_filter2 --exclude range old_analysis/high-LD-regions-hg19-GRCh37.txt --indep-pairwise 50 5 0.5 --out isa_hg37.prunedLD


# Extract by LD - ISA

./plink2 --bfile isa_hg37_filter2 --extract isa_hg37.prunedLD.prune.in --make-bed --out isa_hg37.prunedld


# Extract by LD - 1000g

./plink2 --bfile old_analysis/1000gp --extract isa_hg37.prunedLD.prune.in --make-bed --out 1000gp.prunedld


# Extract 1000g

./plink2 --bfile 1000gp.prunedld --keep allref.txt --make-bed --out 1000gp_filter.prunedld

# Filter SNPs

cut -f2 isa_hg37.prunedld.bim > isa.snp.txt
cut -f2 1000gp_filter.prunedld.bim > 1000gp.snps.txt
grep -f isa.snp.txt 1000gp.snps.txt > comuns-isa-1000g.txt


# Extract SNPs

./plink2 --bfile isa_hg37.prunedld --extract comuns-isa-1000g.txt --max-alleles 2 --make-bed --out isacomm

./plink2 --bfile 1000gp_filter.prunedld --extract comuns-isa-1000g.txt --max-alleles 2 --make-bed --out 1000gpcomm


# Extract Samples ISA - Remove 

./plink2 --bfile isacomm --remove incluidas.txt --make-bed --out isacomm_excluidas

# Extract Samples ISA - Keep

./plink2 --bfile isacomm --keep incluidas.txt --make-bed --out isacomm_2

# Separando os multialélicos

./plink --bfile isacomm_excluidas --bmerge 1000gpcomm.bed 1000gpcomm.bim 1000gpcomm.fam --make-bed --out isacomm_excluidas_1

# Extract multialelicos - ISA - 134

./plink --bfile isacomm_excluidas_1 --exclude  isacomm_excluidas_1-merge.missnp --make-bed --out isacomm_excluidas_2

# Extract multialelicos - 1000g - 134

./plink --bfile 1000gpcomm --exclude isacomm_excluidas_1-merge.missnp --make-bed --out 1000gpcomm2

# Merge ISA 134 and 1000g

./plink --bfile isamerge1000ANDExcluidas --bmerge 1000gpcomm2 --make-bed --out isamerge1000ANDExcluidas

# Separando os multialélicos

./plink --bfile isacomm_2 --bmerge 1000gpcomm.bed 1000gpcomm.bim 1000gpcomm.fam --make-bed --out isamerge1000_1

# Extract multialelicos - 1000g - 707

./plink --bfile 1000gpcomm --exclude isamerge1000_1-merge.missnp --make-bed --out 1000gpcomm2_2

# Extract multialelicos - ISA - 707

./plink --bfile isacomm_2 --exclude isamerge1000_1-merge.missnp --make-bed --out isamerge1000g

# Merge ISA 707 and 1000g

./plink --bfile isamerge1000g --bmerge 1000gpcomm2_2 --make-bed --out isamerge1000