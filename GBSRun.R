#!/usr/bin/env Rscript

gform <- "ANGSDcounts"   # uneak (default), Tassel or chip

genofile <- "../angsd/raso26.NRscafs.min100reads.counts.counts.gz"

source("~/programs/KGD/GBS-Chip-Gmatrix.R")

# Gfull      <- calcG()
# GHWdgm.05  <- calcG(which(HWdis > -0.05), "HWdgm.05")  # recalculate using Hardy-Weinberg disequilibrium cut-off at -0.05
# G_HW05_MAF4  <- calcG(which(HWdis > -0.05 & maf < 0.4), "HWdgm.05_MAF.4")  # recalculate using Hardy-Weinberg disequilibrium cut-off at -0.05 and a MAF frequency of 0.4
# GHW0  <- calcG(which(HWdis > 0), "HW0")  # recalculate using Hardy-Weinberg disequilibrium cut-off at 0
# G_HW0_MAF3  <- calcG(which(HWdis > 0 & maf<0.3), "HW0_MAF.3")  # recalculate using Hardy-Weinberg disequilibrium cut-off at 0
# GHWp2  <- calcG(which(l10pstar < 2), "HWp2")  # recalculate using Hardy-Weinberg log p value below 2
# G_HW_0_snpdepth_6  <- calcG(which(HWdis > 0 & snpdepth <= 6), "HW_0_snpdepth_5_20")  # recalculate using Hardy-Weinberg log p value below 2
G_HW0_1_raso <- calcG(which(HWdis > 0 & HWdis < 0.1), "HW0_1")  # recalculate using Hardy-Weinberg log p value below 2

genofile <- "../angsd/skyN10.NRscafs.min50reads.counts.counts.gz"

source("~/programs/KGD/GBS-Chip-Gmatrix.R")

G_HW0_1_skyN <- calcG(which(HWdis > 0 & HWdis < 0.1), "HW0_1")  # recalculate using Hardy-Weinberg log p value below 2




###export matrices
#fetch sample names
bam_names_raso <- read.table("../angsd/raso26.bam.NRscafs.list", as.is=T)[,1]
names_raso <- gsub("_1alignmentonly.NRscafs.bam", "", gsub("../bam/","",bam_names_raso))
rownames(G_HW0_1_raso$G5) <- names_raso
colnames(G_HW0_1_raso$G5) <- names_raso

bam_names_skyN <- read.table("../angsd/skylarkNL10.bam.NRscafs.list", as.is=T)[,1]
names_skyN <- gsub("_1alignmentonly.NRscafs.bam", "", gsub("../bam/","",bam_names_skyN))
rownames(G_HW0_1_skyN$G5) <- names_skyN
colnames(G_HW0_1_skyN$G5) <- names_skyN

#export
write.table(G_HW0_1_raso$G5, file="raso26_KGD_G5_HW0_1.tsv", sep="\t",quote=F)

write.table(G_HW0_1_skyN$G5, file="skyN10_KGD_G5_HW0_1.tsv", sep="\t",quote=F)

