## Maddalena Ardissino
# June 2023
# Transcriptome wide MR for pre-eclampsia

# bonf for preec 0.05/116943 = 4.28e-07
# bonf for gesthtn 0.05/ 132841 = 3.76e-07

#### Load packages ####
library(data.table)
library(dplyr)
library(tidyr)
library(foreign)
library(tibble)
library(TwoSampleMR)
library(MendelianRandomization)
library(meta)
library(metafor)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)
library(R.utils)
library(stringr)
library(survival)
library(survminer)
library(forestplot)
library(RadialMR)
library(MRPRESSO)
library(TwoSampleMR)
library(openxlsx)
library(MVMR)
library(coloc)
library(forestplot)
library(PredictABEL)
library(strex)
library(readr)
library(beepr)
library(PheWAS)
library(RColorBrewer)
library(vroom)
library(bread)
library(ieugwasr)
library(ggrepel)
library(SeqArray)
library(RColorBrewer)
library(ggmanh)
setwd('~/desktop/preec')

#### Working notes ####

#### LIST of tissue types ####
# # Adipose_Subcutaneous 
# # Adipose_Visceral_Omentum 
# # Adrenal_Gland 
# Artery_Aorta 
# Artery_Coronary 
# Artery_Tibial 
# Brain_Amygdala 
# Brain_Anterior_cingulate_cortex_BA24 --
# Brain_Caudate_basal_ganglia --
# Brain_Cerebellar_Hemisphere --
# Brain_Cerebellum --
# Brain_Cortex 
# Brain_Frontal_Cortex_BA9 --
# Brain_Hippocampus --
# Brain_Hypothalamus 
# Brain_Nucleus_accumbens_basal_ganglia --
# Brain_Putamen_basal_ganglia --
# Brain_Spinal_cord_cervical_c-1 --
# Brain_Substantia_nigra --
# Breast_Mammary_Tissue 
# Cells_Cultured_fibroblasts
# Cells_EBV-transformed_lymphocytes
# Colon_Sigmoid 
# Colon_Transverse
# Esophagus_Gastroesophageal_Junction
# Esophagus_Mucosa
# Esophagus_Muscularis
# Heart_Atrial_Appendage
# Heart_Left_Ventricle
# Kidney_Cortex
# Liver
# Lung
# Minor_Salivary_Gland
# Muscle_Skeletal
# Nerve_Tibial
# Ovary
# Pancreas --
# Pituitary --
# Prostate --
# Skin_Not_Sun_Exposed_Suprapubic --
# Skin_Sun_Exposed_Lower_leg --
# Small_Intestine_Terminal_Ileum --
# Spleen --
# Stomach --
# Thyroid --
# Uterus --
# Vagina --

#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------ GTEX ---------------------------------------------- ####
#### ------------------------------------------ preec --------------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####

#### Format GTEX annotation document - ONCE ONLY ####

# # # Merge chr pos rsid from gtex annotation file - ONLY RUN THIS ONCE!
# setwd('/volumes/maddy2/gtex')
# # length(count.fields('annotate.txt')) # 46569705 rows!
# # fread('/volumes/maddy2/gtex/annotate.txt', sep = "\t", nrow = 10)
# meta<-bmeta(file='/volumes/maddy2/gtex/annotate.txt')
# bfile_split( file='/volumes/maddy2/gtex/annotate.txt', by_columns = 'chr')
# 
# rm(list=ls())

# #### Format instruments document - ONCE ONLY ####
# # Load merged data 
# setwd('~/desktop/preec')
# eqtls <- as.data.frame(fread('gtex/gtex-eqtls.csv')) # the chr and pos here are in build 37!!
# head(eqtls)
# eqtls$phenotype <- str_c(eqtls$gene_name, ':', eqtls$tissue)
# dput(colnames(eqtls))
# eqtls <- eqtls[,c("phenotype", "Chromosome", "Position", 
#                   "effect_allele", "other_allele", "beta.exposure", "se.exposure", "pval.exposure")] 
# colnames(eqtls) <- c('phenotype', 'chr.exposure', 'pos.exposureb37', 'effect_allele.exposure', 'other_allele.exposure', 
#                      'beta.exposure', 'se.exposure', 'pval.exposureb37') # ? need to pull eaf from mapper file? not done for now
# table(eqtls$chr.exposure) # up to 22, no X
# eqtls$chrpos37 <-  str_c(eqtls$chr.exposure, ':', eqtls$pos.exposure)
# 
# # Chr 1
# eqtls_chr1 <- filter(eqtls, eqtls$chr.exposure == 1)
# annot_chr1 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr1.csv'))
# annot_chr1$chrpos38 <- str_c('1:', annot_chr1$variant_pos)
# annot_chr1 <- annot_chr1[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr1$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr1$variant_id_b37))
# eqtls_chr1_annot <- merge(eqtls_chr1, annot_chr1, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr1_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr1_annot.csv')
# rm(eqtls_chr1, annot_chr1, eqtls_chr1_annot)
# 
# # Chr 2
# eqtls_chr2 <- filter(eqtls, eqtls$chr.exposure == 2)
# annot_chr2 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr2.csv'))
# annot_chr2$chrpos38 <- str_c('2:', annot_chr2$variant_pos)
# annot_chr2 <- annot_chr2[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr2$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr2$variant_id_b37))
# eqtls_chr2_annot <- merge(eqtls_chr2, annot_chr2, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr2_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr2_annot.csv')
# rm(eqtls_chr2, annot_chr2, eqtls_chr2_annot)
# 
# 
# # Chr 3
# eqtls_chr3 <- filter(eqtls, eqtls$chr.exposure == 3)
# annot_chr3 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr3.csv'))
# annot_chr3$chrpos38 <- str_c('3:', annot_chr3$variant_pos)
# annot_chr3 <- annot_chr3[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr3$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr3$variant_id_b37))
# eqtls_chr3_annot <- merge(eqtls_chr3, annot_chr3, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr3_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr3_annot.csv')
# rm(eqtls_chr3, annot_chr3, eqtls_chr3_annot)
# 
# # Chr 4
# eqtls_chr4 <- filter(eqtls, eqtls$chr.exposure == 4)
# annot_chr4 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr4.csv'))
# annot_chr4$chrpos38 <- str_c('4:', annot_chr4$variant_pos)
# annot_chr4 <- annot_chr4[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr4$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr4$variant_id_b37))
# eqtls_chr4_annot <- merge(eqtls_chr4, annot_chr4, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr4_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr4_annot.csv')
# rm(eqtls_chr4, annot_chr4, eqtls_chr4_annot)
# 
# # Chr 5
# eqtls_chr5 <- filter(eqtls, eqtls$chr.exposure == 5)
# annot_chr5 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr5.csv'))
# annot_chr5$chrpos38 <- str_c('5:', annot_chr5$variant_pos)
# annot_chr5 <- annot_chr5[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr5$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr5$variant_id_b37))
# eqtls_chr5_annot <- merge(eqtls_chr5, annot_chr5, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr5_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr5_annot.csv')
# rm(eqtls_chr5, annot_chr5, eqtls_chr5_annot)
# 
# # Chr 6
# eqtls_chr6 <- filter(eqtls, eqtls$chr.exposure == 6)
# annot_chr6 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr6.csv'))
# annot_chr6$chrpos38 <- str_c('6:', annot_chr6$variant_pos)
# annot_chr6 <- annot_chr6[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr6$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr6$variant_id_b37))
# eqtls_chr6_annot <- merge(eqtls_chr6, annot_chr6, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr6_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr6_annot.csv')
# rm(eqtls_chr6, annot_chr6, eqtls_chr6_annot)
# 
# # Chr 7
# eqtls_chr7 <- filter(eqtls, eqtls$chr.exposure == 7)
# annot_chr7 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr7.csv'))
# annot_chr7$chrpos38 <- str_c('7:', annot_chr7$variant_pos)
# annot_chr7 <- annot_chr7[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr7$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr7$variant_id_b37))
# eqtls_chr7_annot <- merge(eqtls_chr7, annot_chr7, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr7_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr7_annot.csv')
# rm(eqtls_chr7, annot_chr7, eqtls_chr7_annot)
# 
# # Chr 8
# eqtls_chr8 <- filter(eqtls, eqtls$chr.exposure == 8)
# annot_chr8 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr8.csv'))
# annot_chr8$chrpos38 <- str_c('8:', annot_chr8$variant_pos)
# annot_chr8 <- annot_chr8[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr8$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr8$variant_id_b37))
# eqtls_chr8_annot <- merge(eqtls_chr8, annot_chr8, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr8_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr8_annot.csv')
# rm(eqtls_chr8, annot_chr8, eqtls_chr8_annot)
# 
# # Chr 9
# eqtls_chr9 <- filter(eqtls, eqtls$chr.exposure == 9)
# annot_chr9 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr9.csv'))
# annot_chr9$chrpos38 <- str_c('9:', annot_chr9$variant_pos)
# annot_chr9 <- annot_chr9[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr9$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr9$variant_id_b37))
# eqtls_chr9_annot <- merge(eqtls_chr9, annot_chr9, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr9_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr9_annot.csv')
# rm(eqtls_chr9, annot_chr9, eqtls_chr9_annot)
# 
# # Chr 10
# eqtls_chr10 <- filter(eqtls, eqtls$chr.exposure == 10)
# annot_chr10 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr10.csv'))
# annot_chr10$chrpos38 <- str_c('10:', annot_chr10$variant_pos)
# annot_chr10 <- annot_chr10[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr10$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr10$variant_id_b37))
# eqtls_chr10_annot <- merge(eqtls_chr10, annot_chr10, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr10_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr10_annot.csv')
# rm(eqtls_chr10, annot_chr10, eqtls_chr10_annot)
# 
# # Chr 11
# eqtls_chr11 <- filter(eqtls, eqtls$chr.exposure == 11)
# annot_chr11 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr11.csv'))
# annot_chr11$chrpos38 <- str_c('11:', annot_chr11$variant_pos)
# annot_chr11 <- annot_chr11[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr11$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr11$variant_id_b37))
# eqtls_chr11_annot <- merge(eqtls_chr11, annot_chr11, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr11_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr11_annot.csv')
# rm(eqtls_chr11, annot_chr11, eqtls_chr11_annot)
# 
# # Chr 12
# eqtls_chr12 <- filter(eqtls, eqtls$chr.exposure == 12)
# annot_chr12 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr12.csv'))
# annot_chr12$chrpos38 <- str_c('12:', annot_chr12$variant_pos)
# annot_chr12 <- annot_chr12[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr12$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr12$variant_id_b37))
# eqtls_chr12_annot <- merge(eqtls_chr12, annot_chr12, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr12_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr12_annot.csv')
# rm(eqtls_chr12, annot_chr12, eqtls_chr12_annot)
# 
# # Chr 13
# eqtls_chr13 <- filter(eqtls, eqtls$chr.exposure == 13)
# annot_chr13 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr13.csv'))
# annot_chr13$chrpos38 <- str_c('13:', annot_chr13$variant_pos)
# annot_chr13 <- annot_chr13[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr13$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr13$variant_id_b37))
# eqtls_chr13_annot <- merge(eqtls_chr13, annot_chr13, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr13_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr13_annot.csv')
# rm(eqtls_chr13, annot_chr13, eqtls_chr13_annot)
# 
# # Chr 14
# eqtls_chr14 <- filter(eqtls, eqtls$chr.exposure == 14)
# annot_chr14 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr14.csv'))
# annot_chr14$chrpos38 <- str_c('14:', annot_chr14$variant_pos)
# annot_chr14 <- annot_chr14[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr14$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr14$variant_id_b37))
# eqtls_chr14_annot <- merge(eqtls_chr14, annot_chr14, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr14_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr14_annot.csv')
# rm(eqtls_chr14, annot_chr14, eqtls_chr14_annot)
# 
# # Chr 15
# eqtls_chr15 <- filter(eqtls, eqtls$chr.exposure == 15)
# annot_chr15 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr15.csv'))
# annot_chr15$chrpos38 <- str_c('15:', annot_chr15$variant_pos)
# annot_chr15 <- annot_chr15[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr15$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr15$variant_id_b37))
# eqtls_chr15_annot <- merge(eqtls_chr15, annot_chr15, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr15_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr15_annot.csv')
# rm(eqtls_chr15, annot_chr15, eqtls_chr15_annot)
# 
# # Chr 16
# eqtls_chr16 <- filter(eqtls, eqtls$chr.exposure == 16)
# annot_chr16 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr16.csv'))
# annot_chr16$chrpos38 <- str_c('16:', annot_chr16$variant_pos)
# annot_chr16 <- annot_chr16[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr16$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr16$variant_id_b37))
# eqtls_chr16_annot <- merge(eqtls_chr16, annot_chr16, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr16_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr16_annot.csv')
# rm(eqtls_chr16, annot_chr16, eqtls_chr16_annot)
# 
# # Chr 17
# eqtls_chr17 <- filter(eqtls, eqtls$chr.exposure == 17)
# annot_chr17 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr17.csv'))
# annot_chr17$chrpos38 <- str_c('17:', annot_chr17$variant_pos)
# annot_chr17 <- annot_chr17[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr17$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr17$variant_id_b37))
# eqtls_chr17_annot <- merge(eqtls_chr17, annot_chr17, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr17_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr17_annot.csv')
# rm(eqtls_chr17, annot_chr17, eqtls_chr17_annot)
# 
# # Chr 18
# eqtls_chr18 <- filter(eqtls, eqtls$chr.exposure == 18)
# annot_chr18 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr18.csv'))
# annot_chr18$chrpos38 <- str_c('18:', annot_chr18$variant_pos)
# annot_chr18 <- annot_chr18[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr18$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr18$variant_id_b37))
# eqtls_chr18_annot <- merge(eqtls_chr18, annot_chr18, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr18_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr18_annot.csv')
# rm(eqtls_chr18, annot_chr18, eqtls_chr18_annot)
# 
# # Chr 19
# eqtls_chr19 <- filter(eqtls, eqtls$chr.exposure == 19)
# annot_chr19 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr19.csv'))
# annot_chr19$chrpos38 <- str_c('19:', annot_chr19$variant_pos)
# annot_chr19 <- annot_chr19[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr19$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr19$variant_id_b37))
# eqtls_chr19_annot <- merge(eqtls_chr19, annot_chr19, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr19_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr19_annot.csv')
# rm(eqtls_chr19, annot_chr19, eqtls_chr19_annot)
# 
# # Chr 20
# eqtls_chr20 <- filter(eqtls, eqtls$chr.exposure == 20)
# annot_chr20 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr20.csv'))
# annot_chr20$chrpos38 <- str_c('20:', annot_chr20$variant_pos)
# annot_chr20 <- annot_chr20[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr20$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr20$variant_id_b37))
# eqtls_chr20_annot <- merge(eqtls_chr20, annot_chr20, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr20_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr20_annot.csv')
# rm(eqtls_chr20, annot_chr20, eqtls_chr20_annot)
# 
# # Chr 21
# eqtls_chr21 <- filter(eqtls, eqtls$chr.exposure == 21)
# annot_chr21 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr21.csv'))
# annot_chr21$chrpos38 <- str_c('21:', annot_chr21$variant_pos)
# annot_chr21 <- annot_chr21[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr21$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr21$variant_id_b37))
# eqtls_chr21_annot <- merge(eqtls_chr21, annot_chr21, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr21_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr21_annot.csv')
# rm(eqtls_chr21, annot_chr21, eqtls_chr21_annot)
# 
# # Chr 22
# eqtls_chr22 <- filter(eqtls, eqtls$chr.exposure == 22)
# annot_chr22 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr22.csv'))
# annot_chr22$chrpos38 <- str_c('22:', annot_chr22$variant_pos)
# annot_chr22 <- annot_chr22[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr22$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr22$variant_id_b37))
# eqtls_chr22_annot <- merge(eqtls_chr22, annot_chr22, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr22_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr22_annot.csv')
# rm(eqtls_chr22, annot_chr22, eqtls_chr22_annot)
# 
# setwd('~/desktop/preec')
# library(readr)
# files <- list.files('/volumes/maddy2/gtex/temp_annot', 
#                     pattern = ".csv$", recursive = TRUE, full.names = TRUE)
# eqtls_annotated <- as.data.frame(read_csv(files, col_names = c("...1", "chrpos37", "phenotype", "chr.exposure", "pos.exposureb37", 
#                                                                "effect_allele.exposure", "other_allele.exposure", "beta.exposure", 
#                                                                "se.exposure", "pval.exposureb37", "chrpos38", "rs_id_dbSNP151_GRCh38p7", 
#                                                                "variant_id_b37"), 
#                                           col_types = c("-", "c", "c", "n", "n", 
#                                                         "c", "c", "n", 
#                                                         "n", "n", "c", "c", 
#                                                         "c")) %>% bind_rows())
# eqtls_annotated <- eqtls_annotated[-1,]
# head(eqtls_annotated)
# write.csv(eqtls_annotated, file = 'gtex/eqtls_instruments.csv')
# 
# rm(list=ls())
# 
#### Format outcome data ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('gtex/eqtls_instruments.csv', header = TRUE))[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7')]
head(pqtls)

preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt',  # in HG38!
                             drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos38 <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos38', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos38')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos38'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos38'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, 'outgtex/preec_out.csv')
rm(list=ls())

#### Separate unadjusted instruments, select only available in preeclampsia and save ####

setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('gtex/eqtls_instruments.csv', header = TRUE))[,-1]
preec_out <- as.data.frame(fread('outgtex/preec_out.csv', header = TRUE))[,c('SNP', 'chrpos38')]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(preec_out$chrpos38)),]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(preec_out$chrpos38)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec/*")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())


#### ---------------------------  MR --------------------------- ####
#### MR - Adipose_Subcutaneous ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Adipose_Subcutaneous.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adipose_Subcutaneous.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Subcutaneous")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Subcutaneous")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Subcutaneous/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adipose_Subcutaneous.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Subcutaneous", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres



mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Adipose_Subcutaneous.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]x 
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adipose_Subcutaneous.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Adipose_Subcutaneous.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adipose_Subcutaneous.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adipose_Subcutaneous', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adipose_Subcutaneous"),
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp) +  ylim(0, 16)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Adipose_Subcutaneous.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())


#### MR - Adipose_Visceral_Omentum ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Adipose_Visceral_Omentum.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adipose_Visceral_Omentum.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Visceral_Omentum")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Visceral_Omentum")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Visceral_Omentum/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adipose_Visceral_Omentum.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adipose_Visceral_Omentum", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Adipose_Visceral_Omentum.csv", row.names = FALSE)

mergedres <- rbind(filter(mergedres, mergedres$method == 'Inverse variance weighted'), filter(mergedres, mergedres$method == 'Wald ratio'))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adipose_Visceral_Omentum.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Adipose_Visceral_Omentum.csv", row.names = FALSE)
rm(list=ls())


# Plot results

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adipose_Visceral_Omentum.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adipose_Visceral_Omentum', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adipose_Visceral_Omentum"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Adipose_Visceral_Omentum.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Adrenal_Gland ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Adrenal_Gland.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adrenal_Gland.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adrenal_Gland")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adrenal_Gland")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adrenal_Gland/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adrenal_Gland.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Adrenal_Gland", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Adrenal_Gland.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adrenal_Gland.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Adrenal_Gland.csv", row.names = FALSE)
rm(list=ls())


# Plot results# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Adrenal_Gland.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adrenal_Gland', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adrenal_Gland"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Adrenal_Gland.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Artery_Aorta ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Artery_Aorta.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Aorta.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Aorta")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Aorta")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Aorta/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Aorta.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Aorta", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Artery_Aorta.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Aorta.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Artery_Aorta.csv", row.names = FALSE)
rm(list=ls())


# Plot results

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Aorta.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Aorta', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Aorta"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Artery_Aorta.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Artery_Coronary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Artery_Coronary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Coronary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Coronary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Coronary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Coronary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Coronary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Coronary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Artery_Coronary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Coronary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Artery_Coronary.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Coronary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Coronary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Coronary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Artery_Coronary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())








#### MR - Artery_Tibial ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Artery_Tibial.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Tibial.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Tibial")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Tibial")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Tibial/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Tibial.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Artery_Tibial", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Artery_Tibial.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Tibial.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Artery_Tibial.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Artery_Tibial.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Tibial', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Tibial"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Artery_Tibial.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())








#### MR - Brain_Amygdala ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Brain_Amygdala.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Amygdala.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Amygdala")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Amygdala")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Amygdala/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Amygdala.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Amygdala", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Brain_Amygdala.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Amygdala.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Brain_Amygdala.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Amygdala.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Amygdala', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Amygdala"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Brain_Amygdala.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Brain_Cortex ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Brain_Cortex.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Cortex.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Cortex")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Cortex")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Cortex/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Cortex.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Cortex", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Brain_Cortex.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Cortex.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Brain_Cortex.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Cortex.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Cortex', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Cortex"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Brain_Cortex.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Brain_Hypothalamus ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Brain_Hypothalamus.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Hypothalamus.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Hypothalamus")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Hypothalamus")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Hypothalamus/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Hypothalamus.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Brain_Hypothalamus", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Brain_Hypothalamus.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Hypothalamus.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Brain_Hypothalamus.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Brain_Hypothalamus.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Hypothalamus', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Hypothalamus"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Brain_Hypothalamus.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Breast_Mammary_Tissue ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Breast_Mammary_Tissue.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Breast_Mammary_Tissue.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Breast_Mammary_Tissue")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Breast_Mammary_Tissue")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Breast_Mammary_Tissue/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Breast_Mammary_Tissue.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Breast_Mammary_Tissue", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Breast_Mammary_Tissue.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Breast_Mammary_Tissue.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Breast_Mammary_Tissue.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Breast_Mammary_Tissue.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Breast_Mammary_Tissue', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Breast_Mammary_Tissue"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Breast_Mammary_Tissue.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())







# #### MR - Cells_Cultured_fibroblasts ####
# # Import all IVs
# setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
# files = list.files(pattern="*Cells_Cultured_fibroblasts.csv")   #make list of all csv names
# data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
# names(data_list) <- gsub(".csv","",
#                          list.files(pattern="*Cells_Cultured_fibroblasts.csv", full.names = FALSE),
#                          fixed = TRUE)   #make names for GE - withoout the .csv
# invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
# genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
# length(genelist)
# 
# rm(files, data_list)
# setwd("~/Desktop/preec")
# 
# # Format IVs as exposure data
# formfunc<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
#                     beta_col = "beta.exposure",
#                     se_col = "se.exposure",
#                     pval_col="pval.exposureb37",
#                     # eaf_col = "eaf.exposure",
#                     effect_allele_col = "effect_allele.exposure",
#                     other_allele_col = "other_allele.exposure", 
#                     phenotype_col =  "phenotype")
# }
# 
# join_list<-lapply(genelist, formfunc)
# names(join_list) <- genelist
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# rm(formfunc, join_list)
# 
# 
# # Import outcome association estimates
# outex<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <- read_outcome_data(snps = dat$SNP, 
#                            filename = "~/Desktop/preec/outgtex/preec_out.csv", 
#                            sep = ",", 
#                            snp_col = "SNP",
#                            beta_col = "beta.outcome",
#                            se_col = "se.outcome",
#                            effect_allele_col = "effect_allele.outcome",
#                            other_allele_col = "other_allele.outcome",
#                            eaf_col = "eaf.outcome",
#                            pval_col = "pval.outcome")
# }
# join_list<-lapply(genelist, outex)
# names(join_list) <- str_c("out_",genelist) 
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# 
# outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets
# 
# length(genelist)
# length(outlist)
# 
# rm(outex, join_list)
# 
# 
# # Harmonise data
# hrm<-function(exp, out){
#   exp <- get(exp, envir = .GlobalEnv)
#   out <- get(out, envir = .GlobalEnv)
#   dat <- harmonise_data(exp, out, action=2)
# }
# join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
# names(join_list) <- str_c("har_",names(join_list)) 
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm, hrm, join_list)
# genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets
# 
# rm(list=ls()[!(ls() %in% genelist)]) 
# genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets
# 
# # Perform MR
# mrfunc2<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   res<- mr(dat)
# }
# mr_table2 <- lapply(genelist, mrfunc2)
# names(mr_table2) <- gsub("har_","",genelist) 
# names(mr_table2) <- str_c(names(mr_table2), '_res')
# invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))
# 
# reslist <- names(mr_table2)
# rm(list=ls()[!(ls() %in% reslist)]) 
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# 
# reslist<-ls(pattern = "_res", mget(ls()))
# length(reslist)
# 
# dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_Cultured_fibroblasts")
# setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_Cultured_fibroblasts")
# unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_Cultured_fibroblasts/*")
# 
# files <- mget(reslist)
# for (i in 1:length(files)){
#   write.csv(files[[i]], paste(names(files[i]), "_Cells_Cultured_fibroblasts.csv", sep = ""))
# }
# rm(list=ls())
# 
# # Merge results
# setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
# files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_Cultured_fibroblasts", 
#                     pattern = ".csv", recursive = TRUE, full.names = TRUE)
# files <- files[which(file.info(files)$size>3)]
# 
# alldat <- lapply(setNames(nm = files), read.csv)
# mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
# rownames(mergedres) <- NULL
# mergedres$X <- NULL
# mergedres$filename <- NULL
# mergedres
# 
# mergedres$or <- exp(mergedres$b)
# mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
# mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
# mergedres[,c(1:3)] <- NULL
# mergedres$outcome <- 'Pre-eclampsia'
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# 
# mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
# mergedres <- mergedres[!duplicated(mergedres$exposure),]
# mergedres <- mergedres[order(mergedres$exposure),]
# mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
# mergedres <- mergedres[order(mergedres$pval),]
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# 
# mergedres <- filter(mergedres, mergedres$padj < 0.05)
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# rm(list=ls())
# 
# 
# # Plot results
# # Full manhattan plot
# mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Cells_Cultured_fibroblasts.csv'))
# mergedres <- mergedres[order(mergedres$exposure),]
# mergedres$exposure <- gsub(':Cells_Cultured_fibroblasts', '', mergedres$exposure)
# mergedres$Gene <- mergedres$exposure
# sigp <- 0.05/nrow(mergedres)
# originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
# originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
# originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
# originaldata <- originaldata[!is.na(originaldata$Gene),]
# originaldata <- originaldata[,c('Gene', "chr", "pos")]
# originaldata <- originaldata[!duplicated(originaldata$Gene),]
# mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)
# 
# mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
# mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
# mergedres$logp <- -log10(mergedres$pval)
# labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
# mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
# mergedres$pos <- as.numeric(mergedres$pos)
# g <- manhattan_plot(x = mergedres, 
#                     pval.colname = "pval", 
#                     chr.colname = "chromosome", 
#                     pos.colname = "pos", 
#                     plot.title = gsub("_", ' ', "Cells_Cultured_fibroblasts"), 
#                     label.colname = "siglabel", label.font.size = 2, 
#                     chr.order = c(1:22), signif = sigp)
# pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Cells_Cultured_fibroblasts.pdf', width = 9, height = 5)
# g
# dev.off()
# rm(list=ls())
# 
# 
# 
# 
# 
#### MR - Cells_EBV-transformed_lymphocytes ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Cells_EBV-transformed_lymphocytes.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Cells_EBV-transformed_lymphocytes.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_EBV-transformed_lymphocytes")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_EBV-transformed_lymphocytes")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_EBV-transformed_lymphocytes/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Cells_EBV-transformed_lymphocytes.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Cells_EBV-transformed_lymphocytes", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)
rm(list=ls())


# Plot results

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Cells_EBV-transformed_lymphocytes.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Cells_EBV-transformed_lymphocytes', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Cells_EBV-transformed_lymphocytes"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Cells_EBV-transformed_lymphocytes.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Colon_Sigmoid  ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Colon_Sigmoid.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Colon_Sigmoid.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Sigmoid")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Sigmoid")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Sigmoid/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Colon_Sigmoid.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Sigmoid", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Colon_Sigmoid.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Colon_Sigmoid.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Colon_Sigmoid.csv", row.names = FALSE)
rm(list=ls())


# Plot results
# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Colon_Sigmoid.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Colon_Sigmoid', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Colon_Sigmoid"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Colon_Sigmoid.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Colon_Transverse####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Colon_Transverse.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Colon_Transverse.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Transverse")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Transverse")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Transverse/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Colon_Transverse.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Colon_Transverse", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Colon_Transverse.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Colon_Transverse.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Colon_Transverse.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Colon_Transverse.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Colon_Transverse', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Colon_Transverse"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Colon_Transverse.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Esophagus_Gastroesophageal_Junction ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Esophagus_Gastroesophageal_Junction.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Esophagus_Gastroesophageal_Junction.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Gastroesophageal_Junction")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Gastroesophageal_Junction")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Gastroesophageal_Junction/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Esophagus_Gastroesophageal_Junction.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Gastroesophageal_Junction", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Esophagus_Gastroesophageal_Junction.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Esophagus_Gastroesophageal_Junction', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Esophagus_Gastroesophageal_Junction"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Esophagus_Gastroesophageal_Junction.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Esophagus_Mucosa ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Esophagus_Mucosa.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Esophagus_Mucosa.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Mucosa")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Mucosa")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Mucosa/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Esophagus_Mucosa.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Esophagus_Mucosa", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Esophagus_Mucosa.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Esophagus_Mucosa.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Esophagus_Mucosa.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Esophagus_Mucosa.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Esophagus_Mucosa', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Esophagus_Mucosa"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Esophagus_Mucosa.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Heart_Atrial_Appendage ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Heart_Atrial_Appendage.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Heart_Atrial_Appendage.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Atrial_Appendage")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Atrial_Appendage")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Atrial_Appendage/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Heart_Atrial_Appendage.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Atrial_Appendage", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Heart_Atrial_Appendage.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Heart_Atrial_Appendage.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Heart_Atrial_Appendage.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Heart_Atrial_Appendage.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Heart_Atrial_Appendage', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Heart_Atrial_Appendage"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Heart_Atrial_Appendage.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Heart_Left_Ventricle ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Heart_Left_Ventricle.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Heart_Left_Ventricle.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Left_Ventricle")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Left_Ventricle")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Left_Ventricle/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Heart_Left_Ventricle.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Heart_Left_Ventricle", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Heart_Left_Ventricle.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Heart_Left_Ventricle.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Heart_Left_Ventricle.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Heart_Left_Ventricle.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Heart_Left_Ventricle', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Heart_Left_Ventricle"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Heart_Left_Ventricle.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Kidney_Cortex ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Kidney_Cortex.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Kidney_Cortex.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Kidney_Cortex")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Kidney_Cortex")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Kidney_Cortex/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Kidney_Cortex.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Kidney_Cortex", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Kidney_Cortex.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Kidney_Cortex.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Kidney_Cortex.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Kidney_Cortex.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Kidney_Cortex', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Kidney_Cortex"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Kidney_Cortex.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())







#### MR - Liver ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Liver.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Liver.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Liver")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Liver")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Liver/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Liver.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Liver", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Liver.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Liver.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Liver.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Liver.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Liver', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Liver"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Liver.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Lung ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Lung.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Lung.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Lung")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Lung")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Lung/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Lung.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Lung", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Lung.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Lung.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Lung.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Lung.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Lung', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Lung"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Lung.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Minor_Salivary_Gland ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Minor_Salivary_Gland.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Minor_Salivary_Gland.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Minor_Salivary_Gland")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Minor_Salivary_Gland")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Minor_Salivary_Gland/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Minor_Salivary_Gland.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Minor_Salivary_Gland", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Minor_Salivary_Gland.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Minor_Salivary_Gland.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Minor_Salivary_Gland.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Minor_Salivary_Gland.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Minor_Salivary_Gland', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Minor_Salivary_Gland"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Minor_Salivary_Gland.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Muscle_Skeletal ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Muscle_Skeletal.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Muscle_Skeletal.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Muscle_Skeletal")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Muscle_Skeletal")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Muscle_Skeletal/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Muscle_Skeletal.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Muscle_Skeletal", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Muscle_Skeletal.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Muscle_Skeletal.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Muscle_Skeletal.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Muscle_Skeletal.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Muscle_Skeletal', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Muscle_Skeletal"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Muscle_Skeletal.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Nerve_Tibial ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Nerve_Tibial.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Nerve_Tibial.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Nerve_Tibial")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Nerve_Tibial")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Nerve_Tibial/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Nerve_Tibial.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Nerve_Tibial", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Nerve_Tibial.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Nerve_Tibial.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Nerve_Tibial.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Nerve_Tibial.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Nerve_Tibial', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Nerve_Tibial"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Nerve_Tibial.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Ovary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Ovary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Ovary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Ovary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Ovary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Ovary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Ovary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Ovary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Ovary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Ovary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Ovary.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Ovary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Ovary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Ovary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Ovary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())







#### MR - Pancreas ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Pancreas.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Pancreas.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pancreas")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pancreas")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pancreas/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Pancreas.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pancreas", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Pancreas.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Pancreas.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Pancreas.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Pancreas.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Pancreas', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Pancreas"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Pancreas.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Pituitary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Pituitary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Pituitary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pituitary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pituitary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pituitary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Pituitary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Pituitary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Pituitary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Pituitary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Pituitary.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Pituitary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Pituitary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Pituitary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Pituitary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Skin_Not_Sun_Exposed_Suprapubic ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Skin_Not_Sun_Exposed_Suprapubic.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Skin_Not_Sun_Exposed_Suprapubic.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Not_Sun_Exposed_Suprapubic")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Not_Sun_Exposed_Suprapubic")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Not_Sun_Exposed_Suprapubic/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Skin_Not_Sun_Exposed_Suprapubic.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Not_Sun_Exposed_Suprapubic", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Skin_Not_Sun_Exposed_Suprapubic.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Skin_Not_Sun_Exposed_Suprapubic', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Skin_Not_Sun_Exposed_Suprapubic"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Skin_Not_Sun_Exposed_Suprapubic.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Skin_Sun_Exposed_Lower_leg ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Skin_Sun_Exposed_Lower_leg.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Skin_Sun_Exposed_Lower_leg.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Sun_Exposed_Lower_leg")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Sun_Exposed_Lower_leg")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Sun_Exposed_Lower_leg/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Skin_Sun_Exposed_Lower_leg.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Skin_Sun_Exposed_Lower_leg", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Skin_Sun_Exposed_Lower_leg.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Skin_Sun_Exposed_Lower_leg', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Skin_Sun_Exposed_Lower_leg"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Skin_Sun_Exposed_Lower_leg.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Small_Intestine_Terminal_Ileum ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Small_Intestine_Terminal_Ileum.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Small_Intestine_Terminal_Ileum.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Small_Intestine_Terminal_Ileum")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Small_Intestine_Terminal_Ileum")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Small_Intestine_Terminal_Ileum/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Small_Intestine_Terminal_Ileum.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Small_Intestine_Terminal_Ileum", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Small_Intestine_Terminal_Ileum.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Small_Intestine_Terminal_Ileum', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Small_Intestine_Terminal_Ileum"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Small_Intestine_Terminal_Ileum.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Spleen ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Spleen.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Spleen.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list) 
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Spleen")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Spleen")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Spleen/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Spleen.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Spleen", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Spleen.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Spleen.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Spleen.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Spleen.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Spleen', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Spleen"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Spleen.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Stomach ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Stomach.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Stomach.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Stomach")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Stomach")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Stomach/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Stomach.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Stomach", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Stomach.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Stomach.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Stomach.csv", row.names = FALSE)
rm(list=ls())

rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Stomach.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Stomach', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Stomach"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Stomach.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Thyroid ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Thyroid.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Thyroid.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Thyroid")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Thyroid")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Thyroid/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Thyroid.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Thyroid", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Thyroid.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Thyroid.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Thyroid.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Thyroid.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Thyroid', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Thyroid"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Thyroid.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Uterus ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Uterus.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Uterus.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Uterus")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Uterus")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Uterus/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Uterus.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Uterus", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Uterus.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Uterus.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Uterus.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Uterus.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Uterus', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Uterus"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Uterus.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Vagina  ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Vagina.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Vagina.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Vagina")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Vagina")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Vagina/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Vagina.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Vagina", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Vagina.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Vagina.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Vagina.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Vagina.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Vagina', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Vagina"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Vagina.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Whole_Blood ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_preec")
files = list.files(pattern="*Whole_Blood.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Whole_Blood.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outgtex/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Whole_Blood")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Whole_Blood")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Whole_Blood/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Whole_Blood.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Whole_Blood", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Whole_Blood.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Whole_Blood.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Whole_Blood.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Whole_Blood.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Whole_Blood', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Whole_Blood"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Whole_Blood.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### Plot all genes that have at least one significant results in heatmap ####
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/", 
                    pattern = "^fdrsig_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>80)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)

# res <- mergedres[,c('Gene', 'Tissue', 'b')]
# res <- as.data.frame(pivot_wider(res, names_from = Tissue, values_from = b))
# rownames(res) <- res$Gene

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=b)) +
  geom_tile() +  scale_fill_gradient(low="yellow", high="blue") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# geom_text(aes(Tissue, Gene, label = b), color = "black", size = 4)
pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/gtex_sigonly_heatmap.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()

mergedres$posdirect <- ifelse(mergedres$b > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=posdirect)) +
  geom_tile() + scale_fill_gradient(low="gold", high="blue") +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/gtex_heatmap_sigonly_negpos.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()

## Heatmaps with all bets (not significant only)

genelist <- mergedres[!duplicated(mergedres$Gene),]
genelist <- genelist[,c('Gene', 'Tissue')]

setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>1)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)

mergedres <- mergedres[which(mergedres$Gene %in% genelist$Gene),]
mergedres$signif <- ifelse(mergedres$padj <0.05, '*', '')

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=b)) +
  geom_tile() +  scale_fill_gradient(low="gold", high="blue") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "none", panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/gtex_all_heatmap.pdf", width = 9, height = 16)
cor 
dev.off()

mergedres$posdirect <- ifelse(mergedres$b > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=posdirect)) +
  geom_tile() + scale_fill_gradient(low="gold", high="blue") +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/gtex_heatmap_all_negpos.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()
rm(list=ls())

#### Specific gene plots for scRNAseq project ####

## Heatmaps with specific genes

setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>1)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)


mergedres <- filter(mergedres, mergedres$Gene == 'ADAMTS10' | mergedres$Gene == 'ATAD3B' | mergedres$Gene == 'ATXN7L3' 
                    | mergedres$Gene == 'GABPB2' |mergedres$Gene == 'LATS2' | mergedres$Gene == 'METAP1' 
                    | mergedres$Gene == 'PARD3B' |mergedres$Gene == 'POLM' | mergedres$Gene == 'RSBN1' 
                    | mergedres$Gene == 'S100A12' |mergedres$Gene == 'S100A8' | mergedres$Gene == 'S100A9' 
                    | mergedres$Gene == 'SNX18' |mergedres$Gene == 'UBALD2' | mergedres$Gene == 'SUMF2' 
                    | mergedres$Gene == 'XPNPEP1' |mergedres$Gene == 'POLM')
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', n=length(mergedres$pval))

mergedres$signif <- ifelse(mergedres$pval <1,'ns',  '')
mergedres$signif<- ifelse(mergedres$pval <0.05,'*', mergedres$signif)
mergedres$signif<- ifelse(mergedres$pval <0.001,'***', mergedres$signif)

mergedres <- mergedres %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres$signif!='ns'  ~ "Direct",
  b<0 & mergedres$signif!='ns' ~ "Inverse",
  b>0 & mergedres$signif=='ns'  ~ "No association",
  b<0 & mergedres$signif=='ns' ~ "No association",
  is.na(b) ~ as.character(NA)))
mergedres$signif <- gsub('ns', '', mergedres$signif)

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Direct association", "Inverse association", "No significant association", 'Missing')) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/specific_genes_mch/heatmap_alltissues.pdf", width = 10, height = 4)
cor
dev.off()
write.csv(mergedres[,1:12], "~/Desktop/preec/specific_genes_mch/gtex_mr_alltissues.csv")

# only artery tissues
mergedres <- filter(mergedres, mergedres$Tissue == 'Artery Aorta' | mergedres$Tissue == 'Artery Tibial' )
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', n=length(mergedres$pval))

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = 'Pastel1', labels = c("Direct association",  "No significant association", 'Missing')) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/specific_genes_mch/heatmap_arteries.pdf", width =  4, height = 3)
cor
dev.off()
write.csv(mergedres[,1:12], "~/Desktop/preec/specific_genes_mch/gtex_mr_arteries.csv")

rm(list=ls())

#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------ GTEX ---------------------------------------------- ####
#### ----------------------------------------- geshtn --------------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####
# #### Format GTEX annotation document - ONCE ONLY ####
# 
# # # # Merge chr pos rsid from gtex annotation file - ONLY RUN THIS ONCE!
# # setwd('/volumes/maddy2/gtex')
# # # length(count.fields('annotate.txt')) # 46569705 rows!
# # # fread('/volumes/maddy2/gtex/annotate.txt', sep = "\t", nrow = 10)
# # meta<-bmeta(file='/volumes/maddy2/gtex/annotate.txt')
# # bfile_split( file='/volumes/maddy2/gtex/annotate.txt', by_columns = 'chr')
# # 
# # rm(list=ls())
# 
# #### Format instruments document ####
# # Load merged data (SuppTab2 from CHD paper)
# setwd('~/desktop/preec')
# eqtls <- as.data.frame(fread('gtex/gtex-eqtls.csv')) # the chr and pos here are in build 37!!
# head(eqtls)
# eqtls$phenotype <- str_c(eqtls$gene_name, ':', eqtls$tissue)
# dput(colnames(eqtls))
# eqtls <- eqtls[,c("phenotype", "Chromosome", "Position", 
#                   "effect_allele", "other_allele", "beta.exposure", "se.exposure", "pval.exposure")] 
# colnames(eqtls) <- c('phenotype', 'chr.exposure', 'pos.exposureb37', 'effect_allele.exposure', 'other_allele.exposure', 
#                      'beta.exposure', 'se.exposure', 'pval.exposureb37') # ? need to pull eaf from mapper file? not done for now
# table(eqtls$chr.exposure) # up to 22, no X
# eqtls$chrpos37 <-  str_c(eqtls$chr.exposure, ':', eqtls$pos.exposure)
# 
# # Chr 1
# eqtls_chr1 <- filter(eqtls, eqtls$chr.exposure == 1)
# annot_chr1 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr1.csv'))
# annot_chr1$chrpos38 <- str_c('1:', annot_chr1$variant_pos)
# annot_chr1 <- annot_chr1[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr1$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr1$variant_id_b37))
# eqtls_chr1_annot <- merge(eqtls_chr1, annot_chr1, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr1_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr1_annot.csv')
# rm(eqtls_chr1, annot_chr1, eqtls_chr1_annot)
# 
# # Chr 2
# eqtls_chr2 <- filter(eqtls, eqtls$chr.exposure == 2)
# annot_chr2 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr2.csv'))
# annot_chr2$chrpos38 <- str_c('2:', annot_chr2$variant_pos)
# annot_chr2 <- annot_chr2[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr2$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr2$variant_id_b37))
# eqtls_chr2_annot <- merge(eqtls_chr2, annot_chr2, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr2_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr2_annot.csv')
# rm(eqtls_chr2, annot_chr2, eqtls_chr2_annot)
# 
# 
# # Chr 3
# eqtls_chr3 <- filter(eqtls, eqtls$chr.exposure == 3)
# annot_chr3 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr3.csv'))
# annot_chr3$chrpos38 <- str_c('3:', annot_chr3$variant_pos)
# annot_chr3 <- annot_chr3[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr3$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr3$variant_id_b37))
# eqtls_chr3_annot <- merge(eqtls_chr3, annot_chr3, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr3_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr3_annot.csv')
# rm(eqtls_chr3, annot_chr3, eqtls_chr3_annot)
# 
# # Chr 4
# eqtls_chr4 <- filter(eqtls, eqtls$chr.exposure == 4)
# annot_chr4 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr4.csv'))
# annot_chr4$chrpos38 <- str_c('4:', annot_chr4$variant_pos)
# annot_chr4 <- annot_chr4[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr4$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr4$variant_id_b37))
# eqtls_chr4_annot <- merge(eqtls_chr4, annot_chr4, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr4_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr4_annot.csv')
# rm(eqtls_chr4, annot_chr4, eqtls_chr4_annot)
# 
# # Chr 5
# eqtls_chr5 <- filter(eqtls, eqtls$chr.exposure == 5)
# annot_chr5 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr5.csv'))
# annot_chr5$chrpos38 <- str_c('5:', annot_chr5$variant_pos)
# annot_chr5 <- annot_chr5[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr5$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr5$variant_id_b37))
# eqtls_chr5_annot <- merge(eqtls_chr5, annot_chr5, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr5_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr5_annot.csv')
# rm(eqtls_chr5, annot_chr5, eqtls_chr5_annot)
# 
# # Chr 6
# eqtls_chr6 <- filter(eqtls, eqtls$chr.exposure == 6)
# annot_chr6 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr6.csv'))
# annot_chr6$chrpos38 <- str_c('6:', annot_chr6$variant_pos)
# annot_chr6 <- annot_chr6[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr6$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr6$variant_id_b37))
# eqtls_chr6_annot <- merge(eqtls_chr6, annot_chr6, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr6_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr6_annot.csv')
# rm(eqtls_chr6, annot_chr6, eqtls_chr6_annot)
# 
# # Chr 7
# eqtls_chr7 <- filter(eqtls, eqtls$chr.exposure == 7)
# annot_chr7 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr7.csv'))
# annot_chr7$chrpos38 <- str_c('7:', annot_chr7$variant_pos)
# annot_chr7 <- annot_chr7[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr7$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr7$variant_id_b37))
# eqtls_chr7_annot <- merge(eqtls_chr7, annot_chr7, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr7_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr7_annot.csv')
# rm(eqtls_chr7, annot_chr7, eqtls_chr7_annot)
# 
# # Chr 8
# eqtls_chr8 <- filter(eqtls, eqtls$chr.exposure == 8)
# annot_chr8 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr8.csv'))
# annot_chr8$chrpos38 <- str_c('8:', annot_chr8$variant_pos)
# annot_chr8 <- annot_chr8[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr8$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr8$variant_id_b37))
# eqtls_chr8_annot <- merge(eqtls_chr8, annot_chr8, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr8_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr8_annot.csv')
# rm(eqtls_chr8, annot_chr8, eqtls_chr8_annot)
# 
# # Chr 9
# eqtls_chr9 <- filter(eqtls, eqtls$chr.exposure == 9)
# annot_chr9 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr9.csv'))
# annot_chr9$chrpos38 <- str_c('9:', annot_chr9$variant_pos)
# annot_chr9 <- annot_chr9[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr9$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr9$variant_id_b37))
# eqtls_chr9_annot <- merge(eqtls_chr9, annot_chr9, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr9_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr9_annot.csv')
# rm(eqtls_chr9, annot_chr9, eqtls_chr9_annot)
# 
# # Chr 10
# eqtls_chr10 <- filter(eqtls, eqtls$chr.exposure == 10)
# annot_chr10 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr10.csv'))
# annot_chr10$chrpos38 <- str_c('10:', annot_chr10$variant_pos)
# annot_chr10 <- annot_chr10[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr10$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr10$variant_id_b37))
# eqtls_chr10_annot <- merge(eqtls_chr10, annot_chr10, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr10_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr10_annot.csv')
# rm(eqtls_chr10, annot_chr10, eqtls_chr10_annot)
# 
# # Chr 11
# eqtls_chr11 <- filter(eqtls, eqtls$chr.exposure == 11)
# annot_chr11 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr11.csv'))
# annot_chr11$chrpos38 <- str_c('11:', annot_chr11$variant_pos)
# annot_chr11 <- annot_chr11[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr11$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr11$variant_id_b37))
# eqtls_chr11_annot <- merge(eqtls_chr11, annot_chr11, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr11_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr11_annot.csv')
# rm(eqtls_chr11, annot_chr11, eqtls_chr11_annot)
# 
# # Chr 12
# eqtls_chr12 <- filter(eqtls, eqtls$chr.exposure == 12)
# annot_chr12 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr12.csv'))
# annot_chr12$chrpos38 <- str_c('12:', annot_chr12$variant_pos)
# annot_chr12 <- annot_chr12[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr12$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr12$variant_id_b37))
# eqtls_chr12_annot <- merge(eqtls_chr12, annot_chr12, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr12_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr12_annot.csv')
# rm(eqtls_chr12, annot_chr12, eqtls_chr12_annot)
# 
# # Chr 13
# eqtls_chr13 <- filter(eqtls, eqtls$chr.exposure == 13)
# annot_chr13 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr13.csv'))
# annot_chr13$chrpos38 <- str_c('13:', annot_chr13$variant_pos)
# annot_chr13 <- annot_chr13[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr13$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr13$variant_id_b37))
# eqtls_chr13_annot <- merge(eqtls_chr13, annot_chr13, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr13_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr13_annot.csv')
# rm(eqtls_chr13, annot_chr13, eqtls_chr13_annot)
# 
# # Chr 14
# eqtls_chr14 <- filter(eqtls, eqtls$chr.exposure == 14)
# annot_chr14 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr14.csv'))
# annot_chr14$chrpos38 <- str_c('14:', annot_chr14$variant_pos)
# annot_chr14 <- annot_chr14[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr14$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr14$variant_id_b37))
# eqtls_chr14_annot <- merge(eqtls_chr14, annot_chr14, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr14_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr14_annot.csv')
# rm(eqtls_chr14, annot_chr14, eqtls_chr14_annot)
# 
# # Chr 15
# eqtls_chr15 <- filter(eqtls, eqtls$chr.exposure == 15)
# annot_chr15 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr15.csv'))
# annot_chr15$chrpos38 <- str_c('15:', annot_chr15$variant_pos)
# annot_chr15 <- annot_chr15[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr15$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr15$variant_id_b37))
# eqtls_chr15_annot <- merge(eqtls_chr15, annot_chr15, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr15_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr15_annot.csv')
# rm(eqtls_chr15, annot_chr15, eqtls_chr15_annot)
# 
# # Chr 16
# eqtls_chr16 <- filter(eqtls, eqtls$chr.exposure == 16)
# annot_chr16 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr16.csv'))
# annot_chr16$chrpos38 <- str_c('16:', annot_chr16$variant_pos)
# annot_chr16 <- annot_chr16[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr16$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr16$variant_id_b37))
# eqtls_chr16_annot <- merge(eqtls_chr16, annot_chr16, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr16_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr16_annot.csv')
# rm(eqtls_chr16, annot_chr16, eqtls_chr16_annot)
# 
# # Chr 17
# eqtls_chr17 <- filter(eqtls, eqtls$chr.exposure == 17)
# annot_chr17 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr17.csv'))
# annot_chr17$chrpos38 <- str_c('17:', annot_chr17$variant_pos)
# annot_chr17 <- annot_chr17[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr17$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr17$variant_id_b37))
# eqtls_chr17_annot <- merge(eqtls_chr17, annot_chr17, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr17_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr17_annot.csv')
# rm(eqtls_chr17, annot_chr17, eqtls_chr17_annot)
# 
# # Chr 18
# eqtls_chr18 <- filter(eqtls, eqtls$chr.exposure == 18)
# annot_chr18 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr18.csv'))
# annot_chr18$chrpos38 <- str_c('18:', annot_chr18$variant_pos)
# annot_chr18 <- annot_chr18[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr18$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr18$variant_id_b37))
# eqtls_chr18_annot <- merge(eqtls_chr18, annot_chr18, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr18_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr18_annot.csv')
# rm(eqtls_chr18, annot_chr18, eqtls_chr18_annot)
# 
# # Chr 19
# eqtls_chr19 <- filter(eqtls, eqtls$chr.exposure == 19)
# annot_chr19 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr19.csv'))
# annot_chr19$chrpos38 <- str_c('19:', annot_chr19$variant_pos)
# annot_chr19 <- annot_chr19[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr19$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr19$variant_id_b37))
# eqtls_chr19_annot <- merge(eqtls_chr19, annot_chr19, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr19_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr19_annot.csv')
# rm(eqtls_chr19, annot_chr19, eqtls_chr19_annot)
# 
# # Chr 20
# eqtls_chr20 <- filter(eqtls, eqtls$chr.exposure == 20)
# annot_chr20 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr20.csv'))
# annot_chr20$chrpos38 <- str_c('20:', annot_chr20$variant_pos)
# annot_chr20 <- annot_chr20[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr20$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr20$variant_id_b37))
# eqtls_chr20_annot <- merge(eqtls_chr20, annot_chr20, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr20_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr20_annot.csv')
# rm(eqtls_chr20, annot_chr20, eqtls_chr20_annot)
# 
# # Chr 21
# eqtls_chr21 <- filter(eqtls, eqtls$chr.exposure == 21)
# annot_chr21 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr21.csv'))
# annot_chr21$chrpos38 <- str_c('21:', annot_chr21$variant_pos)
# annot_chr21 <- annot_chr21[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr21$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr21$variant_id_b37))
# eqtls_chr21_annot <- merge(eqtls_chr21, annot_chr21, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr21_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr21_annot.csv')
# rm(eqtls_chr21, annot_chr21, eqtls_chr21_annot)
# 
# # Chr 22
# eqtls_chr22 <- filter(eqtls, eqtls$chr.exposure == 22)
# annot_chr22 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr22.csv'))
# annot_chr22$chrpos38 <- str_c('22:', annot_chr22$variant_pos)
# annot_chr22 <- annot_chr22[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'variant_id_b37')]
# annot_chr22$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr22$variant_id_b37))
# eqtls_chr22_annot <- merge(eqtls_chr22, annot_chr22, by='chrpos37', all.x=FALSE, all.y=FALSE)
# write.csv(eqtls_chr22_annot, '/volumes/maddy2/gtex/temp_annot/eqtls_chr22_annot.csv')
# rm(eqtls_chr22, annot_chr22, eqtls_chr22_annot)
# 
# setwd('~/desktop/preec')
# library(readr)
# files <- list.files('/volumes/maddy2/gtex/temp_annot', 
#                     pattern = ".csv$", recursive = TRUE, full.names = TRUE)
# eqtls_annotated <- as.data.frame(read_csv(files, col_names = c("...1", "chrpos37", "phenotype", "chr.exposure", "pos.exposureb37", 
#                                                                "effect_allele.exposure", "other_allele.exposure", "beta.exposure", 
#                                                                "se.exposure", "pval.exposureb37", "chrpos38", "rs_id_dbSNP151_GRCh38p7", 
#                                                                "variant_id_b37"), 
#                                           col_types = c("-", "c", "c", "n", "n", 
#                                                         "c", "c", "n", 
#                                                         "n", "n", "c", "c", 
#                                                         "c")) %>% bind_rows())
# eqtls_annotated <- eqtls_annotated[-1,]
# head(eqtls_annotated)
# write.csv(eqtls_annotated, file = 'gtex/eqtls_instruments.csv')
# 
# rm(list=ls())
# 
#### Format outcome data ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('gtex/eqtls_instruments.csv', header = TRUE))[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7')]
head(pqtls)

geshtn <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_geshtn_European_allBiobanks_omitNone_1.txt',  # in HG38!
                              drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
geshtn$chrpos38 <- str_c(geshtn$Chromosome, ':', geshtn$Position)
geshtn_out <- merge(pqtls, geshtn, by='chrpos38', all.x=FALSE, all.y = FALSE)
geshtn_out <- geshtn_out[,c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos38')]
geshtn_out$Allele1 <- toupper(geshtn_out$Allele1)
geshtn_out$Allele2 <- toupper(geshtn_out$Allele2)
setnames(geshtn_out, old=c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos38'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos38'))
geshtn_out$phenotype <- 'geshtn'
geshtn_out <- geshtn_out[!duplicated(geshtn_out$SNP),]
write.csv(geshtn_out, 'outgtex/geshtn_out.csv')
rm(list=ls())

#### Separate unadjusted instruments, select only available in geshtn and save ####

setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('gtex/eqtls_instruments.csv', header = TRUE))[,-1]
geshtn_out <- as.data.frame(fread('outgtex/geshtn_out.csv', header = TRUE))[,c('SNP', 'chrpos38')]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(geshtn_out$chrpos38)),]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(geshtn_out$chrpos38)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn/*")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())


#### --------------------  MR ------------------------- ####

#### MR - Adipose_Subcutaneous ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Adipose_Subcutaneous.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adipose_Subcutaneous.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Subcutaneous")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Subcutaneous")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Subcutaneous/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adipose_Subcutaneous.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Subcutaneous", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Adipose_Subcutaneous.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adipose_Subcutaneous.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Adipose_Subcutaneous.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adipose_Subcutaneous.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adipose_Subcutaneous', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adipose_Subcutaneous"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Adipose_Subcutaneous.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Adipose_Visceral_Omentum ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Adipose_Visceral_Omentum.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adipose_Visceral_Omentum.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Visceral_Omentum")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Visceral_Omentum")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Visceral_Omentum/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adipose_Visceral_Omentum.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adipose_Visceral_Omentum", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Adipose_Visceral_Omentum.csv", row.names = FALSE)

mergedres <- rbind(filter(mergedres, mergedres$method == 'Inverse variance weighted'), filter(mergedres, mergedres$method == 'Wald ratio'))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adipose_Visceral_Omentum.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Adipose_Visceral_Omentum.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adipose_Visceral_Omentum.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adipose_Visceral_Omentum', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adipose_Visceral_Omentum"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Adipose_Visceral_Omentum.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())







#### MR - Adrenal_Gland ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Adrenal_Gland.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Adrenal_Gland.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adrenal_Gland")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adrenal_Gland")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adrenal_Gland/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Adrenal_Gland.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Adrenal_Gland", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Adrenal_Gland.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adrenal_Gland.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Adrenal_Gland.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Adrenal_Gland.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Adrenal_Gland', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Adrenal_Gland"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Adrenal_Gland.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())








#### MR - Artery_Aorta ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Artery_Aorta.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Aorta.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Aorta")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Aorta")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Aorta/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Aorta.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Aorta", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Artery_Aorta.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Aorta.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Artery_Aorta.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Aorta.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Aorta', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Aorta"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Artery_Aorta.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Artery_Coronary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Artery_Coronary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Coronary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Coronary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Coronary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Coronary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Coronary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Coronary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Artery_Coronary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Coronary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Artery_Coronary.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Coronary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Coronary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Coronary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Artery_Coronary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Artery_Tibial ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Artery_Tibial.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Artery_Tibial.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Tibial")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Tibial")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Tibial/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Artery_Tibial.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Artery_Tibial", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Artery_Tibial.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Tibial.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Artery_Tibial.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Artery_Tibial.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Artery_Tibial', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Artery_Tibial"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Artery_Tibial.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Brain_Amygdala ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Brain_Amygdala.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Amygdala.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Amygdala")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Amygdala")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Amygdala/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Amygdala.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Amygdala", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Brain_Amygdala.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Amygdala.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Brain_Amygdala.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Amygdala.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Amygdala', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Amygdala"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Brain_Amygdala.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Brain_Cortex ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Brain_Cortex.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Cortex.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Cortex")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Cortex")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Cortex/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Cortex.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Cortex", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Brain_Cortex.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Cortex.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Brain_Cortex.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Cortex.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Cortex', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Cortex"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Brain_Cortex.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Brain_Hypothalamus ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Brain_Hypothalamus.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Brain_Hypothalamus.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Hypothalamus")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Hypothalamus")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Hypothalamus/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Brain_Hypothalamus.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Brain_Hypothalamus", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Brain_Hypothalamus.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Hypothalamus.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Brain_Hypothalamus.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Brain_Hypothalamus.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Brain_Hypothalamus', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Brain_Hypothalamus"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Brain_Hypothalamus.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Breast_Mammary_Tissue ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Breast_Mammary_Tissue.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Breast_Mammary_Tissue.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Breast_Mammary_Tissue")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Breast_Mammary_Tissue")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Breast_Mammary_Tissue/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Breast_Mammary_Tissue.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Breast_Mammary_Tissue", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Breast_Mammary_Tissue.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Breast_Mammary_Tissue.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Breast_Mammary_Tissue.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Breast_Mammary_Tissue.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Breast_Mammary_Tissue', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Breast_Mammary_Tissue"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Breast_Mammary_Tissue.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




# #### MR - Cells_Cultured_fibroblasts ####
# # Import all IVs
# setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
# files = list.files(pattern="*Cells_Cultured_fibroblasts.csv")   #make list of all csv names
# data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
# names(data_list) <- gsub(".csv","",
#                          list.files(pattern="*Cells_Cultured_fibroblasts.csv", full.names = FALSE),
#                          fixed = TRUE)   #make names for GE - withoout the .csv
# invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
# genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
# length(genelist)
# 
# rm(files, data_list)
# setwd("~/desktop/preec")
# 
# # Format IVs as exposure data
# formfunc<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
#                     beta_col = "beta.exposure",
#                     se_col = "se.exposure",
#                     pval_col="pval.exposureb37",
#                     # eaf_col = "eaf.exposure",
#                     effect_allele_col = "effect_allele.exposure",
#                     other_allele_col = "other_allele.exposure", 
#                     phenotype_col =  "phenotype")
# }
# 
# join_list<-lapply(genelist, formfunc)
# names(join_list) <- genelist
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# rm(formfunc, join_list)
# 
# 
# # Import outcome association estimates
# outex<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <- read_outcome_data(snps = dat$SNP, 
#                            filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
#                            sep = ",", 
#                            snp_col = "SNP",
#                            beta_col = "beta.outcome",
#                            se_col = "se.outcome",
#                            effect_allele_col = "effect_allele.outcome",
#                            other_allele_col = "other_allele.outcome",
#                            eaf_col = "eaf.outcome",
#                            pval_col = "pval.outcome")
# }
# join_list<-lapply(genelist, outex)
# names(join_list) <- str_c("out_",genelist) 
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# 
# outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets
# 
# length(genelist)
# length(outlist)
# 
# rm(outex, join_list)
# 
# 
# # Harmonise data
# hrm<-function(exp, out){
#   exp <- get(exp, envir = .GlobalEnv)
#   out <- get(out, envir = .GlobalEnv)
#   dat <- harmonise_data(exp, out, action=2)
# }
# join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
# names(join_list) <- str_c("har_",names(join_list)) 
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm, hrm, join_list)
# genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets
# 
# rm(list=ls()[!(ls() %in% genelist)]) 
# genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets
# 
# # Perform MR
# mrfunc2<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   # dat$mr_keep <- TRUE
#   res<- mr(dat)
# }
# mr_table2 <- lapply(genelist, mrfunc2)
# names(mr_table2) <- gsub("har_","",genelist) 
# names(mr_table2) <- str_c(names(mr_table2), '_res')
# invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))
# 
# reslist <- names(mr_table2)
# rm(list=ls()[!(ls() %in% reslist)]) 
# 
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# 
# reslist<-ls(pattern = "_res", mget(ls()))
# length(reslist)
# 
# dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_Cultured_fibroblasts")
# setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_Cultured_fibroblasts")
# unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_Cultured_fibroblasts/*")
# 
# files <- mget(reslist)
# for (i in 1:length(files)){
#   write.csv(files[[i]], paste(names(files[i]), "_Cells_Cultured_fibroblasts.csv", sep = ""))
# }
# rm(list=ls())
# 
# # Merge results
# setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
# files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_Cultured_fibroblasts", 
#                     pattern = ".csv", recursive = TRUE, full.names = TRUE)
# 
# files <- files[which(file.info(files)$size>3)]
# alldat <- lapply(setNames(nm = files), read.csv)
# mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
# rownames(mergedres) <- NULL
# mergedres$X <- NULL
# mergedres$filename <- NULL
# mergedres
# 
# mergedres$or <- exp(mergedres$b)
# mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
# mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
# mergedres[,c(1:3)] <- NULL
# mergedres$outcome <- 'Gestational hypertension'
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# 
# mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
# mergedres <- mergedres[!duplicated(mergedres$exposure),]
# mergedres <- mergedres[order(mergedres$exposure),]
# mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
# mergedres <- mergedres[order(mergedres$pval),]
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# 
# mergedres <- filter(mergedres, mergedres$padj < 0.05)
# write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Cells_Cultured_fibroblasts.csv", row.names = FALSE)
# rm(list=ls())
# 
# 
# # Full manhattan plot
# mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Cells_Cultured_fibroblasts.csv'))
# mergedres <- mergedres[order(mergedres$exposure),]
# mergedres$exposure <- gsub(':Cells_Cultured_fibroblasts', '', mergedres$exposure)
# mergedres$Gene <- mergedres$exposure
# sigp <- 0.05/nrow(mergedres)
# originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
# originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
# originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
# originaldata <- originaldata[!is.na(originaldata$Gene),]
# originaldata <- originaldata[,c('Gene', "chr", "pos")]
# originaldata <- originaldata[!duplicated(originaldata$Gene),]
# mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)
# 
# mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
# mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
# mergedres$logp <- -log10(mergedres$pval)
# labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
# mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
# mergedres$pos <- as.numeric(mergedres$pos)
# g <- manhattan_plot(x = mergedres, 
#                     pval.colname = "pval", 
#                     chr.colname = "chromosome", 
#                     pos.colname = "pos", 
#                     plot.title = gsub("_", ' ', "Cells_Cultured_fibroblasts"), 
#                     label.colname = "siglabel", label.font.size = 2, 
#                     chr.order = c(1:22), signif = sigp)
# pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Cells_Cultured_fibroblasts.pdf', width = 9, height = 5)
# g
# dev.off()
# rm(list=ls())
# 
# 
# 
# 
#### MR - Cells_EBV-transformed_lymphocytes ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Cells_EBV-transformed_lymphocytes.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Cells_EBV-transformed_lymphocytes.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_EBV-transformed_lymphocytes")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_EBV-transformed_lymphocytes")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_EBV-transformed_lymphocytes/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Cells_EBV-transformed_lymphocytes.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Cells_EBV-transformed_lymphocytes", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Cells_EBV-transformed_lymphocytes.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Cells_EBV-transformed_lymphocytes.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Cells_EBV-transformed_lymphocytes', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Cells_EBV-transformed_lymphocytes"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Cells_EBV-transformed_lymphocytes.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Colon_Sigmoid  ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Colon_Sigmoid.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Colon_Sigmoid.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Sigmoid")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Sigmoid")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Sigmoid/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Colon_Sigmoid.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Sigmoid", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Colon_Sigmoid.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Colon_Sigmoid.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Colon_Sigmoid.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Colon_Sigmoid.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Colon_Sigmoid', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Colon_Sigmoid"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Colon_Sigmoid.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())






#### MR - Colon_Transverse####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Colon_Transverse.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Colon_Transverse.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Transverse")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Transverse")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Transverse/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Colon_Transverse.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Colon_Transverse", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Colon_Transverse.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Colon_Transverse.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Colon_Transverse.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Colon_Transverse.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Colon_Transverse', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Colon_Transverse"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Colon_Transverse.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Esophagus_Gastroesophageal_Junction ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Esophagus_Gastroesophageal_Junction.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Esophagus_Gastroesophageal_Junction.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Gastroesophageal_Junction")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Gastroesophageal_Junction")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Gastroesophageal_Junction/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Esophagus_Gastroesophageal_Junction.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Gastroesophageal_Junction", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Esophagus_Gastroesophageal_Junction.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Esophagus_Gastroesophageal_Junction.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Esophagus_Gastroesophageal_Junction', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Esophagus_Gastroesophageal_Junction"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Esophagus_Gastroesophageal_Junction.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Esophagus_Mucosa ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Esophagus_Mucosa.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Esophagus_Mucosa.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Mucosa")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Mucosa")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Mucosa/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Esophagus_Mucosa.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Esophagus_Mucosa", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Esophagus_Mucosa.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Esophagus_Mucosa.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Esophagus_Mucosa.csv", row.names = FALSE)
rm(list=ls())



# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Esophagus_Mucosa.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Esophagus_Mucosa', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Esophagus_Mucosa"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Esophagus_Mucosa.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Heart_Atrial_Appendage ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Heart_Atrial_Appendage.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Heart_Atrial_Appendage.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Atrial_Appendage")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Atrial_Appendage")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Atrial_Appendage/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Heart_Atrial_Appendage.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Atrial_Appendage", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Heart_Atrial_Appendage.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Heart_Atrial_Appendage.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Heart_Atrial_Appendage.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Heart_Atrial_Appendage.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Heart_Atrial_Appendage', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Heart_Atrial_Appendage"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Heart_Atrial_Appendage.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Heart_Left_Ventricle ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Heart_Left_Ventricle.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Heart_Left_Ventricle.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Left_Ventricle")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Left_Ventricle")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Left_Ventricle/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Heart_Left_Ventricle.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Heart_Left_Ventricle", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Heart_Left_Ventricle.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Heart_Left_Ventricle.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Heart_Left_Ventricle.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Heart_Left_Ventricle.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Heart_Left_Ventricle', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Heart_Left_Ventricle"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Heart_Left_Ventricle.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Kidney_Cortex ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Kidney_Cortex.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Kidney_Cortex.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Kidney_Cortex")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Kidney_Cortex")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Kidney_Cortex/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Kidney_Cortex.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Kidney_Cortex", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Kidney_Cortex.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Kidney_Cortex.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Kidney_Cortex.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Kidney_Cortex.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Kidney_Cortex', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Kidney_Cortex"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Kidney_Cortex.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())





#### MR - Liver ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Liver.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Liver.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Liver")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Liver")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Liver/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Liver.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Liver", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Liver.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Liver.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Liver.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Liver.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Liver', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Liver"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Liver.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Lung ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Lung.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Lung.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Lung")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Lung")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Lung/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Lung.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Lung", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Lung.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Lung.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Lung.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Lung.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Lung', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Lung"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Lung.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Minor_Salivary_Gland ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Minor_Salivary_Gland.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Minor_Salivary_Gland.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Minor_Salivary_Gland")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Minor_Salivary_Gland")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Minor_Salivary_Gland/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Minor_Salivary_Gland.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Minor_Salivary_Gland", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Minor_Salivary_Gland.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Minor_Salivary_Gland.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Minor_Salivary_Gland.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Minor_Salivary_Gland.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Minor_Salivary_Gland', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Minor_Salivary_Gland"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Minor_Salivary_Gland.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Muscle_Skeletal ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Muscle_Skeletal.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Muscle_Skeletal.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Muscle_Skeletal")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Muscle_Skeletal")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Muscle_Skeletal/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Muscle_Skeletal.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Muscle_Skeletal", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Muscle_Skeletal.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Muscle_Skeletal.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Muscle_Skeletal.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Muscle_Skeletal.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Muscle_Skeletal', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Muscle_Skeletal"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Muscle_Skeletal.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Nerve_Tibial ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Nerve_Tibial.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Nerve_Tibial.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Nerve_Tibial")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Nerve_Tibial")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Nerve_Tibial/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Nerve_Tibial.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Nerve_Tibial", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Nerve_Tibial.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Nerve_Tibial.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Nerve_Tibial.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Nerve_Tibial.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Nerve_Tibial', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Nerve_Tibial"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Nerve_Tibial.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Ovary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Ovary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Ovary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Ovary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Ovary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Ovary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Ovary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Ovary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Ovary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Ovary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Ovary.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Ovary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Ovary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Ovary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Ovary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Pancreas ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Pancreas.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Pancreas.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pancreas")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pancreas")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pancreas/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Pancreas.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pancreas", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Pancreas.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Pancreas.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Pancreas.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Pancreas.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Pancreas', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Pancreas"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Pancreas.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Pituitary ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Pituitary.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Pituitary.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pituitary")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pituitary")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pituitary/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Pituitary.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Pituitary", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Pituitary.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Pituitary.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Pituitary.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Pituitary.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Pituitary', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Pituitary"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Pituitary.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Skin_Not_Sun_Exposed_Suprapubic ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Skin_Not_Sun_Exposed_Suprapubic.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Skin_Not_Sun_Exposed_Suprapubic.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Not_Sun_Exposed_Suprapubic")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Not_Sun_Exposed_Suprapubic")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Not_Sun_Exposed_Suprapubic/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Skin_Not_Sun_Exposed_Suprapubic.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Not_Sun_Exposed_Suprapubic", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Skin_Not_Sun_Exposed_Suprapubic.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Skin_Not_Sun_Exposed_Suprapubic.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Skin_Not_Sun_Exposed_Suprapubic', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Skin_Not_Sun_Exposed_Suprapubic"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Skin_Not_Sun_Exposed_Suprapubic.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Skin_Sun_Exposed_Lower_leg ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Skin_Sun_Exposed_Lower_leg.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Skin_Sun_Exposed_Lower_leg.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Sun_Exposed_Lower_leg")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Sun_Exposed_Lower_leg")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Sun_Exposed_Lower_leg/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Skin_Sun_Exposed_Lower_leg.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Skin_Sun_Exposed_Lower_leg", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Skin_Sun_Exposed_Lower_leg.csv", row.names = FALSE)
rm(list=ls())


# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Skin_Sun_Exposed_Lower_leg.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Skin_Sun_Exposed_Lower_leg', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Skin_Sun_Exposed_Lower_leg"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Skin_Sun_Exposed_Lower_leg.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Small_Intestine_Terminal_Ileum ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Small_Intestine_Terminal_Ileum.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Small_Intestine_Terminal_Ileum.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Small_Intestine_Terminal_Ileum")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Small_Intestine_Terminal_Ileum")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Small_Intestine_Terminal_Ileum/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Small_Intestine_Terminal_Ileum.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Small_Intestine_Terminal_Ileum", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Small_Intestine_Terminal_Ileum.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Small_Intestine_Terminal_Ileum.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Small_Intestine_Terminal_Ileum', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Small_Intestine_Terminal_Ileum"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Small_Intestine_Terminal_Ileum.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Spleen ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Spleen.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Spleen.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list) 
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Spleen")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Spleen")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Spleen/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Spleen.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Spleen", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Spleen.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Spleen.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Spleen.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Spleen.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Spleen', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Spleen"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Spleen.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Stomach ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Stomach.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Stomach.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Stomach")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Stomach")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Stomach/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Stomach.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Stomach", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Stomach.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Stomach.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Stomach.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Stomach.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Stomach', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Stomach"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Stomach.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Thyroid ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Thyroid.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Thyroid.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Thyroid")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Thyroid")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Thyroid/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Thyroid.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Thyroid", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Thyroid.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Thyroid.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Thyroid.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Thyroid.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Thyroid', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Thyroid"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Thyroid.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Uterus ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Uterus.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Uterus.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Uterus")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Uterus")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Uterus/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Uterus.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Uterus", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Uterus.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Uterus.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Uterus.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Uterus.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Uterus', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Uterus"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Uterus.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### MR - Vagina  ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Vagina.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Vagina.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Vagina")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Vagina")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Vagina/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Vagina.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Vagina", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Vagina.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Vagina.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Vagina.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Vagina.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Vagina', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Vagina"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Vagina.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())



#### MR - Whole_Blood ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsgtex_geshtn")
files = list.files(pattern="*Whole_Blood.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Whole_Blood.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)

rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposureb37",
                    # eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outgtex/geshtn_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets

length(genelist)
length(outlist)

rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets

rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Perform MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(genelist, mrfunc2)
names(mr_table2) <- gsub("har_","",genelist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Whole_Blood")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Whole_Blood")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Whole_Blood/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Whole_Blood.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/Whole_Blood", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Gestational hypertension'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/all_Whole_Blood.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Whole_Blood.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/fdrsig_Whole_Blood.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/main_merged_Whole_Blood.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Whole_Blood', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)
originaldata <- as.data.frame(fread('~/desktop/preec/gtex/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[!is.na(originaldata$Gene),]
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Whole_Blood"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_geshtn/full_manh_geshtn_Whole_Blood.pdf', width = 9, height = 5)
g
dev.off()
rm(list=ls())




#### Plot all genes that have at least one significant results in heatmap ####
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/", 
                    pattern = "^fdrsig_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)

# res <- mergedres[,c('Gene', 'Tissue', 'b')]
# res <- as.data.frame(pivot_wider(res, names_from = Tissue, values_from = b))
# rownames(res) <- res$Gene

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=b)) +
  geom_tile() +  scale_fill_gradient(low="yellow", high="blue") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# geom_text(aes(Tissue, Gene, label = b), color = "black", size = 4)
pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/gtex_sigonly_heatmap.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()

mergedres$posdirect <- ifelse(mergedres$b > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=posdirect)) +
  geom_tile() + scale_fill_gradient(low="gold", high="blue") +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/gtex_heatmap_sigonly_negpos.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()

## Heatmaps with all bets (not significant only)

genelist <- mergedres[!duplicated(mergedres$Gene),]
genelist <- genelist[,c('Gene', 'Tissue')]

setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)

mergedres <- mergedres[which(mergedres$Gene %in% genelist$Gene),]
mergedres$signif <- ifelse(mergedres$padj <0.05, '*', '')

cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=b)) +
  geom_tile() +  scale_fill_gradient(low="gold", high="blue") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "none", panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/gtex_all_heatmap.pdf", width = 9, height = 16)
cor 
dev.off()

mergedres$posdirect <- ifelse(mergedres$b > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=posdirect)) +
  geom_tile() + scale_fill_gradient(low="gold", high="blue") +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/gtex_heatmap_all_negpos.pdf", width = 9, height = 15)
cor + theme(legend.position = "none")
dev.off()

rm(list=ls())


#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------  PLACENTAL EQTL ---------------------------------- ####

# Placenta eQTLs summary statistics.
# The statistics include Beta, SE, Pvalue, gene_id, CHROM:POS:NON_EFF:EFF, and EAF. 
# SE, standard error of beta; POS, the SNP co-ordinates based on hg19; NON_EFF, the non-effective allele;
# EFF, the effective allele; EAF, the effective allele frequency.
# The cis-window is defined as [TSS - 500kb, TES + 500kb], where TSS denotes transcription start site and TES denotes transcription end site.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6329610/


#### --------------------------------------------------------------------------------------------- ####
#### Format eQTL data ####
setwd('~/desktop/preec')

#chr1
eqtls_chr1 <- as.data.frame(fread('placentaeqtl/S1_Data_chr1_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr1$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr1[,5])
annot_chr1 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr1.csv'))
annot_chr1$chrpos38 <- str_c('1:', annot_chr1$variant_pos)
annot_chr1$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr1$variant_id_b37))
annot_chr1 <- annot_chr1[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr1_annot <- merge(eqtls_chr1, annot_chr1, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr1_annot <- merge(eqtls_chr1_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr1_annot$phenotype <- str_c(eqtls_chr1_annot$gene_symbol, ':Placenta')
eqtls_chr1_annot$effect_allele.exposure <- str_match(eqtls_chr1_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr1_annot$other_allele.exposure <- str_match(eqtls_chr1_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr1_annot <- eqtls_chr1_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr1_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr1_annot, 'placentaeqtl/temp_annot/eqtls_chr1_annot.csv')
rm(list=ls())

#chr2
eqtls_chr2 <- as.data.frame(fread('placentaeqtl/S1_Data_chr2_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr2$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr2[,5])
annot_chr2 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr2.csv'))
annot_chr2$chrpos38 <- str_c('2:', annot_chr2$variant_pos)
annot_chr2$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr2$variant_id_b37))
annot_chr2 <- annot_chr2[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr2_annot <- merge(eqtls_chr2, annot_chr2, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr2_annot <- merge(eqtls_chr2_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr2_annot$phenotype <- str_c(eqtls_chr2_annot$gene_symbol, ':Placenta')
eqtls_chr2_annot$effect_allele.exposure <- str_match(eqtls_chr2_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr2_annot$other_allele.exposure <- str_match(eqtls_chr2_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr2_annot <- eqtls_chr2_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr2_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr2_annot, 'placentaeqtl/temp_annot/eqtls_chr2_annot.csv')
rm(list=ls())

#chr3
eqtls_chr3 <- as.data.frame(fread('placentaeqtl/S1_Data_chr3_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr3$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr3[,5])
annot_chr3 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr3.csv'))
annot_chr3$chrpos38 <- str_c('3:', annot_chr3$variant_pos)
annot_chr3$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr3$variant_id_b37))
annot_chr3 <- annot_chr3[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr3_annot <- merge(eqtls_chr3, annot_chr3, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr3_annot <- merge(eqtls_chr3_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr3_annot$phenotype <- str_c(eqtls_chr3_annot$gene_symbol, ':Placenta')
eqtls_chr3_annot$effect_allele.exposure <- str_match(eqtls_chr3_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr3_annot$other_allele.exposure <- str_match(eqtls_chr3_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr3_annot <- eqtls_chr3_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr3_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr3_annot, 'placentaeqtl/temp_annot/eqtls_chr3_annot.csv')
rm(list=ls())


#chr4
eqtls_chr4 <- as.data.frame(fread('placentaeqtl/S1_Data_chr4_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr4$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr4[,5])
annot_chr4 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr4.csv'))
annot_chr4$chrpos38 <- str_c('4:', annot_chr4$variant_pos)
annot_chr4$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr4$variant_id_b37))
annot_chr4 <- annot_chr4[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr4_annot <- merge(eqtls_chr4, annot_chr4, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr4_annot <- merge(eqtls_chr4_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr4_annot$phenotype <- str_c(eqtls_chr4_annot$gene_symbol, ':Placenta')
eqtls_chr4_annot$effect_allele.exposure <- str_match(eqtls_chr4_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr4_annot$other_allele.exposure <- str_match(eqtls_chr4_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr4_annot <- eqtls_chr4_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr4_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr4_annot, 'placentaeqtl/temp_annot/eqtls_chr4_annot.csv')
rm(list=ls())

#chr5
eqtls_chr5 <- as.data.frame(fread('placentaeqtl/S1_Data_chr5_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr5$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr5[,5])
annot_chr5 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr5.csv'))
annot_chr5$chrpos38 <- str_c('5:', annot_chr5$variant_pos)
annot_chr5$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr5$variant_id_b37))
annot_chr5 <- annot_chr5[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr5_annot <- merge(eqtls_chr5, annot_chr5, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr5_annot <- merge(eqtls_chr5_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr5_annot$phenotype <- str_c(eqtls_chr5_annot$gene_symbol, ':Placenta')
eqtls_chr5_annot$effect_allele.exposure <- str_match(eqtls_chr5_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr5_annot$other_allele.exposure <- str_match(eqtls_chr5_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr5_annot <- eqtls_chr5_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr5_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr5_annot, 'placentaeqtl/temp_annot/eqtls_chr5_annot.csv')
rm(list=ls())

#chr6
eqtls_chr6 <- as.data.frame(fread('placentaeqtl/S1_Data_chr6_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr6$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr6[,5])
annot_chr6 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr6.csv'))
annot_chr6$chrpos38 <- str_c('6:', annot_chr6$variant_pos)
annot_chr6$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr6$variant_id_b37))
annot_chr6 <- annot_chr6[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr6_annot <- merge(eqtls_chr6, annot_chr6, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr6_annot <- merge(eqtls_chr6_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr6_annot$phenotype <- str_c(eqtls_chr6_annot$gene_symbol, ':Placenta')
eqtls_chr6_annot$effect_allele.exposure <- str_match(eqtls_chr6_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr6_annot$other_allele.exposure <- str_match(eqtls_chr6_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr6_annot <- eqtls_chr6_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr6_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr6_annot, 'placentaeqtl/temp_annot/eqtls_chr6_annot.csv')
rm(list=ls())

#chr7
eqtls_chr7 <- as.data.frame(fread('placentaeqtl/S1_Data_chr7_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr7$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr7[,5])
annot_chr7 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr7.csv'))
annot_chr7$chrpos38 <- str_c('7:', annot_chr7$variant_pos)
annot_chr7$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr7$variant_id_b37))
annot_chr7 <- annot_chr7[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr7_annot <- merge(eqtls_chr7, annot_chr7, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr7_annot <- merge(eqtls_chr7_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr7_annot$phenotype <- str_c(eqtls_chr7_annot$gene_symbol, ':Placenta')
eqtls_chr7_annot$effect_allele.exposure <- str_match(eqtls_chr7_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr7_annot$other_allele.exposure <- str_match(eqtls_chr7_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr7_annot <- eqtls_chr7_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr7_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr7_annot, 'placentaeqtl/temp_annot/eqtls_chr7_annot.csv')
rm(list=ls())


#chr8
eqtls_chr8 <- as.data.frame(fread('placentaeqtl/S1_Data_chr8_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr8$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr8[,5])
annot_chr8 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr8.csv'))
annot_chr8$chrpos38 <- str_c('8:', annot_chr8$variant_pos)
annot_chr8$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr8$variant_id_b37))
annot_chr8 <- annot_chr8[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr8_annot <- merge(eqtls_chr8, annot_chr8, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr8_annot <- merge(eqtls_chr8_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr8_annot$phenotype <- str_c(eqtls_chr8_annot$gene_symbol, ':Placenta')
eqtls_chr8_annot$effect_allele.exposure <- str_match(eqtls_chr8_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr8_annot$other_allele.exposure <- str_match(eqtls_chr8_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr8_annot <- eqtls_chr8_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr8_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr8_annot, 'placentaeqtl/temp_annot/eqtls_chr8_annot.csv')
rm(list=ls())


#chr9
eqtls_chr9 <- as.data.frame(fread('placentaeqtl/S1_Data_chr9_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr9$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr9[,5])
annot_chr9 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr9.csv'))
annot_chr9$chrpos38 <- str_c('9:', annot_chr9$variant_pos)
annot_chr9$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr9$variant_id_b37))
annot_chr9 <- annot_chr9[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr9_annot <- merge(eqtls_chr9, annot_chr9, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr9_annot <- merge(eqtls_chr9_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr9_annot$phenotype <- str_c(eqtls_chr9_annot$gene_symbol, ':Placenta')
eqtls_chr9_annot$effect_allele.exposure <- str_match(eqtls_chr9_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr9_annot$other_allele.exposure <- str_match(eqtls_chr9_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr9_annot <- eqtls_chr9_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                        'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr9_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr9_annot, 'placentaeqtl/temp_annot/eqtls_chr9_annot.csv')
rm(list=ls())

#chr10
eqtls_chr10 <- as.data.frame(fread('placentaeqtl/S1_Data_chr10_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr10$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr10[,5])
annot_chr10 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr10.csv'))
annot_chr10$chrpos38 <- str_c('10:', annot_chr10$variant_pos)
annot_chr10$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr10$variant_id_b37))
annot_chr10 <- annot_chr10[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr10_annot <- merge(eqtls_chr10, annot_chr10, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr10_annot <- merge(eqtls_chr10_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr10_annot$phenotype <- str_c(eqtls_chr10_annot$gene_symbol, ':Placenta')
eqtls_chr10_annot$effect_allele.exposure <- str_match(eqtls_chr10_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr10_annot$other_allele.exposure <- str_match(eqtls_chr10_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr10_annot <- eqtls_chr10_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr10_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr10_annot, 'placentaeqtl/temp_annot/eqtls_chr10_annot.csv')
rm(list=ls())

#chr11
eqtls_chr11 <- as.data.frame(fread('placentaeqtl/S1_Data_chr11_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr11$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr11[,5])
annot_chr11 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr11.csv'))
annot_chr11$chrpos38 <- str_c('11:', annot_chr11$variant_pos)
annot_chr11$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr11$variant_id_b37))
annot_chr11 <- annot_chr11[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr11_annot <- merge(eqtls_chr11, annot_chr11, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr11_annot <- merge(eqtls_chr11_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr11_annot$phenotype <- str_c(eqtls_chr11_annot$gene_symbol, ':Placenta')
eqtls_chr11_annot$effect_allele.exposure <- str_match(eqtls_chr11_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr11_annot$other_allele.exposure <- str_match(eqtls_chr11_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr11_annot <- eqtls_chr11_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr11_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr11_annot, 'placentaeqtl/temp_annot/eqtls_chr11_annot.csv')
rm(list=ls())


#chr12
eqtls_chr12 <- as.data.frame(fread('placentaeqtl/S1_Data_chr12_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr12$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr12[,5])
annot_chr12 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr12.csv'))
annot_chr12$chrpos38 <- str_c('12:', annot_chr12$variant_pos)
annot_chr12$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr12$variant_id_b37))
annot_chr12 <- annot_chr12[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr12_annot <- merge(eqtls_chr12, annot_chr12, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr12_annot <- merge(eqtls_chr12_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr12_annot$phenotype <- str_c(eqtls_chr12_annot$gene_symbol, ':Placenta')
eqtls_chr12_annot$effect_allele.exposure <- str_match(eqtls_chr12_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr12_annot$other_allele.exposure <- str_match(eqtls_chr12_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr12_annot <- eqtls_chr12_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr12_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr12_annot, 'placentaeqtl/temp_annot/eqtls_chr12_annot.csv')
rm(list=ls())

#chr13
eqtls_chr13 <- as.data.frame(fread('placentaeqtl/S1_Data_chr13_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr13$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr13[,5])
annot_chr13 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr13.csv'))
annot_chr13$chrpos38 <- str_c('13:', annot_chr13$variant_pos)
annot_chr13$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr13$variant_id_b37))
annot_chr13 <- annot_chr13[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr13_annot <- merge(eqtls_chr13, annot_chr13, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr13_annot <- merge(eqtls_chr13_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr13_annot$phenotype <- str_c(eqtls_chr13_annot$gene_symbol, ':Placenta')
eqtls_chr13_annot$effect_allele.exposure <- str_match(eqtls_chr13_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr13_annot$other_allele.exposure <- str_match(eqtls_chr13_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr13_annot <- eqtls_chr13_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr13_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr13_annot, 'placentaeqtl/temp_annot/eqtls_chr13_annot.csv')
rm(list=ls())

#chr14
eqtls_chr14 <- as.data.frame(fread('placentaeqtl/S1_Data_chr14_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr14$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr14[,5])
annot_chr14 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr14.csv'))
annot_chr14$chrpos38 <- str_c('14:', annot_chr14$variant_pos)
annot_chr14$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr14$variant_id_b37))
annot_chr14 <- annot_chr14[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr14_annot <- merge(eqtls_chr14, annot_chr14, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr14_annot <- merge(eqtls_chr14_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr14_annot$phenotype <- str_c(eqtls_chr14_annot$gene_symbol, ':Placenta')
eqtls_chr14_annot$effect_allele.exposure <- str_match(eqtls_chr14_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr14_annot$other_allele.exposure <- str_match(eqtls_chr14_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr14_annot <- eqtls_chr14_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr14_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr14_annot, 'placentaeqtl/temp_annot/eqtls_chr14_annot.csv')
rm(list=ls())

#chr15
eqtls_chr15 <- as.data.frame(fread('placentaeqtl/S1_Data_chr15_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr15$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr15[,5])
annot_chr15 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr15.csv'))
annot_chr15$chrpos38 <- str_c('15:', annot_chr15$variant_pos)
annot_chr15$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr15$variant_id_b37))
annot_chr15 <- annot_chr15[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr15_annot <- merge(eqtls_chr15, annot_chr15, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr15_annot <- merge(eqtls_chr15_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr15_annot$phenotype <- str_c(eqtls_chr15_annot$gene_symbol, ':Placenta')
eqtls_chr15_annot$effect_allele.exposure <- str_match(eqtls_chr15_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr15_annot$other_allele.exposure <- str_match(eqtls_chr15_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr15_annot <- eqtls_chr15_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr15_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr15_annot, 'placentaeqtl/temp_annot/eqtls_chr15_annot.csv')
rm(list=ls())


#chr16
eqtls_chr16 <- as.data.frame(fread('placentaeqtl/S1_Data_chr16_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr16$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr16[,5])
annot_chr16 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr16.csv'))
annot_chr16$chrpos38 <- str_c('16:', annot_chr16$variant_pos)
annot_chr16$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr16$variant_id_b37))
annot_chr16 <- annot_chr16[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr16_annot <- merge(eqtls_chr16, annot_chr16, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr16_annot <- merge(eqtls_chr16_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr16_annot$phenotype <- str_c(eqtls_chr16_annot$gene_symbol, ':Placenta')
eqtls_chr16_annot$effect_allele.exposure <- str_match(eqtls_chr16_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr16_annot$other_allele.exposure <- str_match(eqtls_chr16_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr16_annot <- eqtls_chr16_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr16_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr16_annot, 'placentaeqtl/temp_annot/eqtls_chr16_annot.csv')
rm(list=ls())


#chr17
eqtls_chr17 <- as.data.frame(fread('placentaeqtl/S1_Data_chr17_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr17$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr17[,5])
annot_chr17 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr17.csv'))
annot_chr17$chrpos38 <- str_c('17:', annot_chr17$variant_pos)
annot_chr17$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr17$variant_id_b37))
annot_chr17 <- annot_chr17[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr17_annot <- merge(eqtls_chr17, annot_chr17, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr17_annot <- merge(eqtls_chr17_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr17_annot$phenotype <- str_c(eqtls_chr17_annot$gene_symbol, ':Placenta')
eqtls_chr17_annot$effect_allele.exposure <- str_match(eqtls_chr17_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr17_annot$other_allele.exposure <- str_match(eqtls_chr17_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr17_annot <- eqtls_chr17_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr17_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr17_annot, 'placentaeqtl/temp_annot/eqtls_chr17_annot.csv')
rm(list=ls())

#chr18
eqtls_chr18 <- as.data.frame(fread('placentaeqtl/S1_Data_chr18_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr18$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr18[,5])
annot_chr18 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr18.csv'))
annot_chr18$chrpos38 <- str_c('18:', annot_chr18$variant_pos)
annot_chr18$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr18$variant_id_b37))
annot_chr18 <- annot_chr18[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr18_annot <- merge(eqtls_chr18, annot_chr18, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr18_annot <- merge(eqtls_chr18_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr18_annot$phenotype <- str_c(eqtls_chr18_annot$gene_symbol, ':Placenta')
eqtls_chr18_annot$effect_allele.exposure <- str_match(eqtls_chr18_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr18_annot$other_allele.exposure <- str_match(eqtls_chr18_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr18_annot <- eqtls_chr18_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr18_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr18_annot, 'placentaeqtl/temp_annot/eqtls_chr18_annot.csv')
rm(list=ls())


#chr19
eqtls_chr19 <- as.data.frame(fread('placentaeqtl/S1_Data_chr19_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr19$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr19[,5])
annot_chr19 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr19.csv'))
annot_chr19$chrpos38 <- str_c('19:', annot_chr19$variant_pos)
annot_chr19$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr19$variant_id_b37))
annot_chr19 <- annot_chr19[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr19_annot <- merge(eqtls_chr19, annot_chr19, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr19_annot <- merge(eqtls_chr19_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr19_annot$phenotype <- str_c(eqtls_chr19_annot$gene_symbol, ':Placenta')
eqtls_chr19_annot$effect_allele.exposure <- str_match(eqtls_chr19_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr19_annot$other_allele.exposure <- str_match(eqtls_chr19_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr19_annot <- eqtls_chr19_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr19_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr19_annot, 'placentaeqtl/temp_annot/eqtls_chr19_annot.csv')
rm(list=ls())

#chr20
eqtls_chr20 <- as.data.frame(fread('placentaeqtl/S1_Data_chr20_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr20$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr20[,5])
annot_chr20 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr20.csv'))
annot_chr20$chrpos38 <- str_c('20:', annot_chr20$variant_pos)
annot_chr20$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr20$variant_id_b37))
annot_chr20 <- annot_chr20[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr20_annot <- merge(eqtls_chr20, annot_chr20, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr20_annot <- merge(eqtls_chr20_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr20_annot$phenotype <- str_c(eqtls_chr20_annot$gene_symbol, ':Placenta')
eqtls_chr20_annot$effect_allele.exposure <- str_match(eqtls_chr20_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr20_annot$other_allele.exposure <- str_match(eqtls_chr20_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr20_annot <- eqtls_chr20_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr20_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr20_annot, 'placentaeqtl/temp_annot/eqtls_chr20_annot.csv')
rm(list=ls())

#chr21
eqtls_chr21 <- as.data.frame(fread('placentaeqtl/S1_Data_chr21_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr21$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr21[,5])
annot_chr21 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr21.csv'))
annot_chr21$chrpos38 <- str_c('21:', annot_chr21$variant_pos)
annot_chr21$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr21$variant_id_b37))
annot_chr21 <- annot_chr21[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr21_annot <- merge(eqtls_chr21, annot_chr21, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr21_annot <- merge(eqtls_chr21_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr21_annot$phenotype <- str_c(eqtls_chr21_annot$gene_symbol, ':Placenta')
eqtls_chr21_annot$effect_allele.exposure <- str_match(eqtls_chr21_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr21_annot$other_allele.exposure <- str_match(eqtls_chr21_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr21_annot <- eqtls_chr21_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr21_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr21_annot, 'placentaeqtl/temp_annot/eqtls_chr21_annot.csv')
rm(list=ls())

#chr22
eqtls_chr22 <- as.data.frame(fread('placentaeqtl/S1_Data_chr22_eQTL.txt.gz')) # the chr and pos here are in build 37!!
eqtls_chr22$chrpos37 <-  sub("^(([^:]*:){1}[^:]*).*", "\\1", eqtls_chr22[,5])
annot_chr22 <- as.data.frame(fread('/volumes/maddy2/gtex/annotate_chr22.csv'))
annot_chr22$chrpos38 <- str_c('22:', annot_chr22$variant_pos)
annot_chr22$chrpos37 <-  gsub('_', ':', sub("^(([^_]*_){1}[^_]*).*", "\\1", annot_chr22$variant_id_b37))
annot_chr22 <- annot_chr22[,c('chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'chrpos37')]
eqtls_chr22_annot <- merge(eqtls_chr22, annot_chr22, by='chrpos37', all.x=FALSE, all.y=FALSE)
annotgene <- as.data.frame(fread('placentaeqtl/annotate.gz'))
eqtls_chr22_annot <- merge(eqtls_chr22_annot, annotgene, by='gene_id', all.x=FALSE, all.y=FALSE)
eqtls_chr22_annot$phenotype <- str_c(eqtls_chr22_annot$gene_symbol, ':Placenta')
eqtls_chr22_annot$effect_allele.exposure <- str_match(eqtls_chr22_annot[,6], '([^:]+)(?::[^:]+){0}$')[,2]
eqtls_chr22_annot$other_allele.exposure <- str_match(eqtls_chr22_annot[,6], '([^:]+)(?::[^:]+){1}$')[,2]
eqtls_chr22_annot <- eqtls_chr22_annot[,c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                          'effect_allele.exposure', 'other_allele.exposure', 'EAF', 'Beta', 'SE', 'Pvalue')]
colnames(eqtls_chr22_annot) <- c('chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure')
write.csv(eqtls_chr22_annot, 'placentaeqtl/temp_annot/eqtls_chr22_annot.csv')
rm(list=ls())



# Merge
setwd('~/desktop/preec')
library(readr)
files <- list.files('~/desktop/preec/placentaeqtl/temp_annot',
                    pattern = ".csv$", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]

eqtls_annotated <- as.data.frame(read_csv(files, col_names = c('...1', 'chrpos37', 'chrpos38', 'rs_id_dbSNP151_GRCh38p7', 'phenotype', 
                                                               'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure',
                                                               'beta.exposure', 'se.exposure', 'pval.exposure'),
                                          col_types = c("-", "c", "c", "c", "c", "c",
                                                        "c", "n", "n",
                                                        "n", "n")) %>% bind_rows())
eqtls_annotated <- eqtls_annotated[-1,]
head(eqtls_annotated)
eqtls_annotated$pval.exposure <- as.numeric(eqtls_annotated$pval.exposure)
eqtls_annotated <- filter(eqtls_annotated, eqtls_annotated$pval.exposure < 5e-8)
write.csv(eqtls_annotated, file = 'placentaeqtl/eqtls_instruments.csv')

rm(list=ls())
# 
#### ------------------- ** preec ** --------------- ####
#### Format outcome data preec ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('placentaeqtl/eqtls_instruments.csv', header = TRUE))[,c('chrpos38', 'chrpos37', 'rs_id_dbSNP151_GRCh38p7')]
head(pqtls)

preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt',  # in HG38!
                             drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos38 <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos38', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos38')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('rs_id_dbSNP151_GRCh38p7', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos38'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos38'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, 'outplacentaeqtl/preec_out.csv')
rm(list=ls())

#### Separate unadjusted instruments, preec ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('placentaeqtl/eqtls_instruments.csv', header = TRUE))[,-1]
preec_out <- as.data.frame(fread('outplacentaeqtl/preec_out.csv', header = TRUE))[,c('SNP', 'chrpos38')]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(preec_out$chrpos38)),]
pqtls <- pqtls[which(paste(pqtls$chrpos38) %in% paste(preec_out$chrpos38)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsplacentaeqtl_preec")
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivsplacentaeqtl_preec/*")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())
#### MR - Placenta preec  ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsplacentaeqtl_preec")
files = list.files(pattern="*Placenta.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*Placenta.csv", full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
length(genelist)
rm(files, data_list)
setwd("~/desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rs_id_dbSNP151_GRCh38p7",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposure",
                    eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "phenotype")
}

join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/desktop/preec/outplacentaeqtl/preec_out.csv", 
                           sep = ",", 
                           snp_col = "SNP",
                           beta_col = "beta.outcome",
                           se_col = "se.outcome",
                           effect_allele_col = "effect_allele.outcome",
                           other_allele_col = "other_allele.outcome",
                           eaf_col = "eaf.outcome",
                           pval_col = "pval.outcome")
}
join_list<-lapply(genelist, outex)
names(join_list) <- str_c("out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
outlist<-ls(pattern = "^out_", mget(ls())) #make list with names of all outcome sets
length(genelist)
length(outlist)
rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list)) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove all iv sets with 0 snps only, if any
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm, hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) #update list with names of all IV sets
rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) #re make list with names of all IV sets

# Clump locally v3
genelist<-ls(pattern = "har_", mget(ls()))
rsid<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat$rsid <- dat$SNP
  dat$pval <- dat$pval.exposure
  dat$id <- dat$exposure
  return(dat) } # format for local clumping
iv_list<- sapply(genelist, rsid, simplify = FALSE)
names(iv_list) <- genelist
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
rm(iv_list)
try_clp <- function(dat) {
  out <- tryCatch(
    {
      dat <- get(dat, envir = .GlobalEnv)
      dat <- ld_clump(dat,
                      plink_bin = genetics.binaRies::get_plink_binary(),
                      bfile = "/volumes/maddy2/mrdata/1kg.v3/EUR")
    },
    error=function(cond) {
      message(paste("Unable to clump:", dat$id))
      message("Here's the original error message:")
      message(cond)
      return(dat[which.min(dat$pval.exposure),])
    },
    warning=function(cond) {
      message(paste("Clumping caused a warning:", dat$id))
      message("Here's the original warning message:")
      message(cond)
      return(dat[which.min(dat$pval.exposure),])
    },
    finally={
      message(paste("Clumped:", dat$id))
    }
  )
} 
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1))# remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # 

unlink("/Volumes/MADDY2/datasets/preeclampsia/placentaeqtl_harm_clump_unadj_cis_preec/*")
setwd("/Volumes/MADDY2/datasets/preeclampsia/placentaeqtl_harm_clump_unadj_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/placentaeqtl_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
rm(files, data_list)
clplist<-ls(pattern = "clp_", mget(ls()))

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("^clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))
reslist <- names(mr_table2)
rm(list=ls()[!(ls() %in% reslist)]) 

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empty results
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)

reslist<-ls(pattern = "_res", mget(ls()))
length(reslist)

dir.create("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Placenta")
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Placenta")
unlink("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Placenta/*")

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), "_Placenta.csv", sep = ""))
}
rm(list=ls())

# Merge results
setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/Placenta", 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)

files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/all_Placenta.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[!duplicated(mergedres$exposure),]
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
mergedres <- mergedres[order(mergedres$pval),]
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Placenta.csv", row.names = FALSE)

mergedres <- filter(mergedres, mergedres$padj < 0.05)
write.csv(mergedres, "/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/fdrsig_Placenta.csv", row.names = FALSE)
rm(list=ls())

# Full manhattan plot
mergedres <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/main_merged_Placenta.csv'))
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$exposure <- gsub(':Placenta', '', mergedres$exposure)
mergedres$Gene <- mergedres$exposure
sigp <- 0.05/nrow(mergedres)

originaldata <- as.data.frame(fread('~/desktop/preec/placentaeqtl/eqtls_instruments.csv', header = TRUE))[,-1]
originaldata <- separate(originaldata, phenotype, c("Gene", "Tissue"), sep = ':')
originaldata <- separate(originaldata, chrpos38, c("chr", "pos"), sep = ':')
originaldata <- originaldata[,c('Gene', "chr", "pos")]
originaldata <- originaldata[!duplicated(originaldata$Gene),]
mergedres <- merge(mergedres, originaldata, by='Gene', all.x=FALSE, all.y=FALSE)

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
mergedres$siglabel <- ifelse(mergedres$padj < 0.05, mergedres$Gene, '')
mergedres$logp <- -log10(mergedres$pval)
labels <- filter(mergedres, mergedres$gws == 1)[,c('Gene')]
mergedres$chromosome <- factor(mergedres$chr,  c(1:22))
mergedres$pos <- as.numeric(mergedres$pos)
g <- manhattan_plot(x = mergedres, 
                    pval.colname = "pval", 
                    chr.colname = "chromosome", 
                    pos.colname = "pos", 
                    plot.title = gsub("_", ' ', "Placenta"), 
                    label.colname = "siglabel", label.font.size = 2, 
                    chr.order = c(1:22), signif = sigp)
pdf('~/desktop/preec/resgtex_preec/full_manh_preec_Placenta.pdf', width = 9, height = 5)
g
dev.off()


rm(list=ls())







#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------  PLOTTING ---------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####
#### GTEx + Placenta heatmaps - preec 116943 tests fdr ####

setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_preec/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
mergedres <- filter(mergedres, mergedres$Tissue != 'Placenta')

mergedres$padj <- p.adjust(mergedres$pval, 'fdr')
rm(files, alldat)
head(mergedres)

genelist <- filter(mergedres, mergedres$padj < 0.05)
genelist <- genelist[!duplicated(genelist$Gene),]
genelist <- genelist[,c('Gene', 'Tissue')]


write.csv(mergedres, '~/desktop/preec/restranscriptomic/transcriptomic_mergedres_main_preec.csv')

mergedres <- mergedres[which(mergedres$Gene %in% genelist$Gene),]
mergedres$signif <- ifelse(mergedres$pval <1,'ns',  '')
mergedres$signif<- ifelse(mergedres$pval <0.05,'*', mergedres$signif)
mergedres$signif<- ifelse(mergedres$pval <0.001,'**', mergedres$signif)
mergedres$signif<- ifelse(mergedres$padj <0.05, '***', mergedres$signif)

mergedres <- mergedres %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres$signif!='ns'  ~ "Direct",
  b<0 & mergedres$signif!='ns' ~ "Inverse",
  b>0 & mergedres$signif=='ns'  ~ "No association",
  b<0 & mergedres$signif=='ns' ~ "No association",
  is.na(b) ~ as.character(NA)))
mergedres$signif <- gsub('ns', '', mergedres$signif)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Direct", "Inverse", "No association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("~/desktop/preec/restranscriptomic/gtex_heatmap_all_preec.pdf", width = 9, height = 15)
cor 
dev.off()


# Split heatmap
genelist1 <- genelist[1:(nrow(genelist)/2),]
genelist2 <- genelist[((nrow(genelist)/2)+1):nrow(genelist),]

mergedres1 <- mergedres[which(mergedres$Gene %in% genelist1$Gene),]
mergedres2 <- mergedres[which(mergedres$Gene %in% genelist2$Gene),]

cor <- ggplot(data = mergedres1, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Direct", "Inverse", "No association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("~/desktop/preec/restranscriptomic/gtex_heatmap_all_preec_1.pdf", width = 10, height = 9)
cor 
dev.off()

cor <- ggplot(data = mergedres2, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Direct", "Inverse", "No association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("~/desktop/preec/restranscriptomic/gtex_heatmap_all_preec_2.pdf", width = 10, height = 9)
cor 
dev.off()




rm(list=ls())


#### GTEx + Placenta heatmaps - geshtn 132841 tests fdr ####

setwd("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/resgtex_geshtn/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
dim(mergedres)
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
mergedres$padj <- p.adjust(mergedres$pval, 'fdr')
rm(files, alldat)
head(mergedres)

genelist <- filter(mergedres, mergedres$padj < 0.05)
genelist <- genelist[!duplicated(genelist$Gene),]
genelist <- genelist[,c('Gene', 'Tissue')]

write.csv(mergedres, '~/desktop/preec/restranscriptomic/transcriptomic_mergedres_main_geshtn.csv')


mergedres <- mergedres[which(mergedres$Gene %in% genelist$Gene),]
mergedres$signif <- ifelse(mergedres$padj <0.05, '*', '')


mergedres$posdirect <- ifelse(mergedres$b > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=b)) +
  geom_tile() + scale_fill_gradient2(low="cyan3", high="coral", mid = 'white', midpoint = 0) +  
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("~/desktop/preec/restranscriptomic/gtex_heatmap_all_geshtn.pdf", width = 9, height = 7)
cor + theme(legend.position = "none")
dev.off()

rm(list=ls())


#### Immune cells heatmap - preec ####

setwd("/Volumes/MADDY2/datasets/preeclampsia/ressoskic_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/ressoskic_preec/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure
mergedres <- separate(mergedres, exposure, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
mergedres$padj <- p.adjust(mergedres$pval, 'fdr')
rm(files, alldat)
head(mergedres)
write.csv(mergedres, '~/desktop/preec/restranscriptomic/transcriptomic_mergedres_immunecells_preec.csv')

mergedres <- filter(mergedres, mergedres$padj < 0.05)
genelist <- mergedres[!duplicated(mergedres$Gene),]
genelist <- genelist[,c('Gene', 'Tissue')]

setwd("/Volumes/MADDY2/datasets/preeclampsia/ressoskic_preec")
files <- list.files("/Volumes/MADDY2/datasets/preeclampsia/ressoskic_preec/", 
                    pattern = "^main_", recursive = TRUE, full.names = TRUE)
files <- files[which(file.info(files)$size>3)]
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL
mergedres$X <- NULL
mergedres$filename <- NULL
mergedres$genetissue <- mergedres$exposure_res
mergedres <- separate(mergedres, exposure_res, c("Gene", "Tissue"), sep = ':')
mergedres$Tissue <- gsub('_', ' ', mergedres$Tissue)
rm(files, alldat)
head(mergedres)

mergedres <- mergedres[which(mergedres$Gene %in% genelist$Gene),]
mergedres$signif <- ifelse(mergedres$padj <0.05, '*', '')

mergedres$posdirect <- ifelse(mergedres$b_res > 0, 1, 0)
cor <- ggplot(data = mergedres, aes(x=Tissue, y=Gene, fill=posdirect)) +
  geom_tile() + scale_fill_gradient(low="gold", high="blue") +  
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("/Volumes/MADDY2/datasets/preeclampsia/restranscriptomic/soskic_heatmap_immunecells_preec.pdf", width = 9, height = 35)
cor + theme(legend.position = "none")
dev.off()

rm(list=ls())

