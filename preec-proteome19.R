## Maddalena Ardissino
# April 2023
# Proteome wide MR for pre-eclampsia

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
library(strex)
library(readr)
library(RColorBrewer)
library(vroom)
library(bread)
library(ieugwasr)
library(ggrepel)
library(ggmanh)
library(SeqArray)
library(RColorBrewer)
library(PredictABEL)
library(beepr)
library(PheWAS)
library(polars)
library(biomaRt)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(MRutils)
library(BSgenome)
library(mr.raps)
library(phenoscanner)

setwd('~/desktop/preec')

#### Working notes ####
# Notes: 
# File can be downloaded using token and then link e.g., 
# # # as.data.frame(fread("https://download.decode.is/s3/download?token=bce23d3b-6819-474a-8aaa-310773a21730&file=10089_7_ASMTL_ASML.txt.gz")
# However one single download = 28 minutes, therefore would take ~70d to do 4906 ptoeins 

# Proofs: 
# # Proof for function'inv_neglog10p' function:
# inv_neglog10p <- function(neglog10p) {
# logp <- - neglog10p / log(exp(1), base = 10)
# t_sq <- qchisq(logp, df = 1, lower.tail=FALSE, log.p = TRUE)
# sqrt(t_sq)
# }
# t <- 1.5
# p <- pchisq(t^2, df=1, lower.tail=FALSE)
# p <- 2*(1-pnorm(t))
# neglog10p <- - log(p, base=10)
# t_2 <- inv_neglog10p(neglog10p) # = 1.5, function works

## NB PREECLAMPSIA DATA IS IN HG38

#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------ DECODE ------------------------------------------ ####
#### --------------------------------------------------------------------------------------------- ####
#### ------------------- **** preec **** --------------- ####
#### Extract protein list ####

# Make merged file with protein list
proteins <- as.data.frame(fread('~/Desktop/preec/decode/decodelinks.txt'))
proteins <- grep(pattern = '.gz', x = proteins[,1], value = TRUE)
names <- sub(".txt.gz.*", "", proteins) 
names <- sub(".*&file=", "", names) 
names <- str_c('pr', names)
proteins <- as.data.frame(cbind(proteins, names))
colnames(proteins) <- c('link', 'name')
write.csv(proteins, 'decode/proteinlist.csv') # 4906 proteins
rm(list=ls())

# make long vs short name file
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('decode/decode-pqtls.csv')) 
dput(colnames(pqtls))
pqtls <- pqtls[,c("variant", "shortname\n(prot.)", "longname\n(prot.)","chr\n(var.)", "pos\n(var.)","beta\n(adj.)", 
                  "beta\n(unadj.)", "-Log10(P)\n(adj.)", "-Log10(P)\n(unadj.)", "Amin", "Amaj")] # only min/maj/maf - need to merge EA etc from reference file 
pqtls <- pqtls[,c("shortname\n(prot.)", "longname\n(prot.)")]
colnames(pqtls) <- c('shortname', 'longname')
pqtls <- pqtls[!duplicated(pqtls$longname),]
pqtls$phenotype <- pqtls$longname
pqtls$phenotype <- str_c('p_', pqtls$phenotype)
pqtls$phenotype <- gsub('/', '_', pqtls$phenotype)
pqtls$phenotype <- gsub(',', '_', pqtls$phenotype)
pqtls$phenotype <- gsub("\\[|\\]", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\\(|\\)", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\'", "", pqtls$phenotype)
pqtls$phenotype <- gsub("\\.", "", pqtls$phenotype)
pqtls$phenotype <- gsub('__', '_', pqtls$phenotype)
write.csv(pqtls, 'decode/namemap_decode.csv')
rm(list=ls())

#### Format instruments document ####

# Load merged data (Supp Tab 2)
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('decode/decode-pqtls.csv')) 
dput(colnames(pqtls))
pqtls <- pqtls[,c("variant", "shortname\n(prot.)", "longname\n(prot.)","chr\n(var.)", "pos\n(var.)","cis/\ntrans", "beta\n(adj.)", 
                  "beta\n(unadj.)", "-Log10(P)\n(adj.)", "-Log10(P)\n(unadj.)", "Amin", "Amaj", "MAF\n(%)", "variance\nexpl\n(var.)", "var.\nexpl\n(cumul.)")] # only min/maj/maf - need to merge EA etc from reference file 

# # Check MinA / EA
# reference <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/assocvariants.annotated.txt'))
# reference <- reference[,c('Chrom', 'Pos', 'effectAllele', 'otherAllele', 'effectAlleleFreq')]
# reference$chrpos <- str_c(reference$Chrom, ':', reference$Pos)
# pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)
# pqtls <- merge(pqtls, reference, by = 'chrpos', all.x=TRUE, all.y = FALSE) # 28191 variants, but returns 42546 due to duplicates
# pqtls_1 <- filter(pqtls, pqtls$effectAllele == pqtls$min) # keep only SNPs where EA is either the minor or major -minor here
# pqtls_2 <- filter(pqtls, pqtls$effectAllele == pqtls$maj) # keep only SNPs where EA is either the minor or major -majoe here
# pqtls <- rbind(pqtls_1, pqtls_2) # back to exactly 28191 SNPs

inv_neglog10p <- function(neglog10p) {
  logp <- - neglog10p / log(exp(1), base = 10)
  t_sq <- qchisq(logp, df = 1, lower.tail=FALSE, log.p = TRUE)
  sqrt(t_sq)
}

pqtls$tstat_adjusted <- inv_neglog10p(pqtls$logp_adjusted)
pqtls$se_adjusted <-  abs(pqtls$beta_adjusted) / pqtls$tstat_adjusted
pqtls$pval_adjusted <- 10^(- pqtls$logp_adjusted)

pqtls$tstat_unadjusted <- inv_neglog10p(pqtls$logp_unadjusted)
pqtls$se_unadjusted <-  abs(pqtls$beta_unadjusted) / pqtls$tstat_unadjusted
pqtls$pval_unadjusted <- 10^(- pqtls$logp_unadjusted)

pqtls$phenotype <- gsub(' ', '_', pqtls$longname)
pqtls$phenotype <- gsub('-', '_', pqtls$phenotype)

pqtls_adjusted <-copy(pqtls)
pqtls_unadjusted <- copy(pqtls)

setnames(pqtls_adjusted, old = c('beta_adjusted', 'se_adjusted', 'pval_adjusted',  'effectAllele', 'otherAllele', 'effectAlleleFreq', 'phenotype'),
         new = c('beta.exposure', 'se.exposure', 'pval.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'phenotype'))
pqtls_adjusted <- pqtls_adjusted[,c('SNP', 'shortname', 'phenotype', 'chr.exposure', 'pos.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure','beta.exposure', 'se.exposure', 'pval.exposure')]
write.csv(pqtls_adjusted, 'decode/pqtls_adjusted_instruments.csv')

setnames(pqtls_unadjusted, old = c('beta_unadjusted', 'se_unadjusted', 'pval_unadjusted',  'effectAllele', 'otherAllele', 'effectAlleleFreq', 'phenotype'),
         new = c('beta.exposure', 'se.exposure', 'pval.exposure', 'effect_allele.exposure','other_allele.exposure', 'eaf.exposure', 'phenotype'))
pqtls_unadjusted_trans <- 
  pqtls_unadjusted[,c('SNP', 'shortname', 'phenotype', 'chr.exposure', 'pos.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure','beta.exposure', 'se.exposure', 'pval.exposure')]
write.csv(pqtls_unadjusted_trans, 'decode/pqtls_unadjusted_instruments.csv')

pqtls_unadjusted_cis <- filter(pqtls_unadjusted, pqtls_unadjusted$cistrans == 'cis')
pqtls_unadjusted_cis <- 
  pqtls_unadjusted_cis[,c('SNP', 'shortname', 'phenotype', 'chr.exposure', 'pos.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure','beta.exposure', 'se.exposure', 'pval.exposure')]
write.csv(pqtls_unadjusted_cis, 'decode/pqtls_unadjusted_instruments_cis.csv')

rm(list=ls())

#### Format outcome data ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('decode/pqtls_unadjusted_instruments.csv'))[,c('SNP', 'chr.exposure', 'pos.exposure')]
pqtls$chr.exposure <- substring(pqtls$chr.exposure, 4)
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)

preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt', 
                             drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, 'outdecode/preec_out_unadj.csv')
rm(list=ls())

#### Separate, select p<1.8e-9 and save ####

# Trans
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('decode/pqtls_unadjusted_instruments.csv'))
pqtls$phenotype <- str_c('p_', pqtls$phenotype)
pqtls$phenotype <- gsub('/', '_', pqtls$phenotype)
pqtls$phenotype <- gsub(',', '_', pqtls$phenotype)
pqtls$phenotype <- gsub("\\[|\\]", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\\(|\\)", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\'", "", pqtls$phenotype)
pqtls$phenotype <- gsub("\\.", "", pqtls$phenotype)
pqtls$phenotype <- gsub('__', '_', pqtls$phenotype)
pqtls <- filter(pqtls, pqtls$pval.exposure < 1.8e-9)
pqtls <- filter(pqtls, pqtls$effect_allele.exposure != '!')
pqtls <- filter(pqtls, pqtls$other_allele.exposure != '!')
pqtls$eaf.exposure <- pqtls$eaf.exposure/100
preec_out <- as.data.frame(fread('outdecode/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls$chr.exposure <- substring(pqtls$chr.exposure, 4)
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)
pqtls <- pqtls[which(paste(pqtls$chrpos) %in% paste(preec_out$chrpos)),]
pqtls <- pqtls[which(paste(pqtls$SNP) %in% paste(preec_out$SNP)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_trans_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_trans_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_trans_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())

# cis
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('decode/pqtls_unadjusted_instruments_cis.csv'))[,-1]
pqtls$phenotype <- str_c('p_', pqtls$phenotype)
pqtls$phenotype <- gsub('/', '_', pqtls$phenotype)
pqtls$phenotype <- gsub(',', '_', pqtls$phenotype)
pqtls$phenotype <- gsub("\\[|\\]", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\\(|\\)", "_", pqtls$phenotype)
pqtls$phenotype <- gsub("\'", "", pqtls$phenotype)
pqtls$phenotype <- gsub("\\.", "", pqtls$phenotype)
pqtls$phenotype <- gsub('__', '_', pqtls$phenotype)
pqtls <- filter(pqtls, pqtls$pval.exposure < 1.8e-9)
pqtls <- filter(pqtls, pqtls$effect_allele.exposure != '!')
pqtls <- filter(pqtls, pqtls$other_allele.exposure != '!')
pqtls$eaf.exposure <- pqtls$eaf.exposure/100
preec_out <- as.data.frame(fread('outdecode/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls$chr.exposure <- substring(pqtls$chr.exposure, 4)
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)
pqtls <- pqtls[which(paste(pqtls$chrpos) %in% paste(preec_out$chrpos)),]
pqtls <- pqtls[which(paste(pqtls$SNP) %in% paste(preec_out$SNP)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_cis_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())

#### MR - cis ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    pval_col="pval.exposure",
                    eaf_col = "eaf.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure", 
                    phenotype_col =  "shortname")
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outdecode/preec_out_unadj.csv", 
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
names(join_list) <- gsub("p_","out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
outlist <- gsub("p_","out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- gsub("p_","har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

unlink("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_unadj_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
rm(files, data_list)
clplist<-ls(pattern = "clp_", mget(ls()))


# rsq and save
nsnp <- data.frame(do.call("rbind", lapply(ls(), function(x) {
  obj = get(x)
  if (class(obj) == "data.frame")
    c(protein = x, nsnp = nrow(obj))
})))
rsq <- function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  sum(2*dat$eaf.exposure*(1-dat$eaf.exposure)*(dat$beta.exposure^2))
}
rsqlist<-as.data.frame(sapply(clplist, rsq, simplify = FALSE))
rsqlist<-as.data.frame(t(rsqlist))
rsqlist$protein <- rownames(rsqlist)
rsqlist <- merge(rsqlist, nsnp, by = 'protein')
colnames(rsqlist) <- c('protein', 'rsq', 'nsnp')
rsqlist$nsnp <- as.numeric(rsqlist$nsnp)
rsqlist$fstat <- ((35559-rsqlist$nsnp-1)/rsqlist$nsnp) * (rsqlist$rsq/(1-rsqlist$rsq))
rsqlist <- rsqlist[order(rsqlist$fstat),]
write.csv(rsqlist, "~/desktop/preec/resdecode_preec/cis_rsqfstat.csv")
rm(rsqlist, rsq, nsnp)

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}
rm(list=ls())

#### Merge results - cis ####
setwd("~/desktop/preec/resdecode_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))


mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "~/Desktop/preec/resdecode_preec/all_merged_cis_preec.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[order(mergedres$pval),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
write.csv(mergedres, "main_cis_preec.csv", row.names = FALSE)


mergedres <- filter(mergedres, mergedres$padj< 0.05)
write.csv(mergedres, "fdrsig_cis_preec.csv", row.names = FALSE)

rm(list=ls())




#### Plot results - cis ####

# Manhattan plot NB GWS = 1.5e-5
mergedres <- as.data.frame(fread('main_cis_preec.csv'))
sigp <- 0.05/nrow(mergedres)
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$gws <- ifelse(mergedres$pval < sigp, 1, 0)
mergedres$logp <- -log10(mergedres$pval)
cols <- c("0" = "grey32","1"= "coral")
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]

p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) +
  geom_hline(yintercept = -log10(sigp), color = "gray32", linetype = "dashed", linewidth = 1, alpha = 0.5)
pdf('decode_cis_preec.pdf', width = 12, height = 6)
p1
dev.off()

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]
p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) 
pdf('proteome_cis_preec_fdrsig.pdf', width = 16, height = 6)
p1
dev.off()

rm(list=ls())


#### --------------------------------------------------------------------------------------------- ####
#### -------------------------------------------- Zheng ------------------------------------------ ####
#### ---------------------------------------------Tier 1/2  -------------------------------------- ####
#### ------------------- **** preec **** --------------- ####
#### Format outcome data ####

setwd('~/desktop/preec')
# pqtls <- as.data.frame(fread('zheng/pqtls_unadjusted_instruments_cis.csv')) # only run once 
# annot <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/assocvariants.annotated.txt'))
# setnames(annot, old = 'rsids', new = 'SNP')
# annot <- annot[,c('Chrom', 'Pos', 'Name', 'SNP')]
# pqtls <- merge(pqtls, annot, by='SNP', all.x=FALSE, all.y=FALSE)
# pqtls$chrpos <- str_c(gsub('chr', '', pqtls$Chrom), ':', pqtls$Pos)
# write.csv(pqtls, 'zheng/pqtls_unadjusted_instruments_cis_chrpos.csv')
# rm(list=ls())

pqtls <- as.data.frame(fread('zheng/pqtls_unadjusted_instruments_cis_chrpos.csv'))
preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt', drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, 'outzheng/preec_out_unadj.csv')
beep(2)
rm(list=ls())

#### Separate, select and save ####
# cis
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('zheng/instruments_zheng_2.csv'))
pqtls <- pqtls[,c("phenotype", "Chrom",  "Pos", "Name",  "SNP", "Type", 'Cis/Trans', 
                  "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
pqtls <- filter(pqtls, pqtls$Type != 'Tier3')
pqtls <- filter(pqtls, pqtls$`Cis/Trans` != 'trans')
pqtls$phenotype <- str_c('p_', pqtls$phenotype)
pqtls$pval.exposure <- NULL # errors w pvals, will infer later
pqtls$chrpos <- gsub('chr', '', pqtls$Name)
preec_out <- as.data.frame(fread('outzheng/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls <- pqtls[which(paste(pqtls$SNP) %in% paste(preec_out$SNP)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_cis_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())

#trans
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('zheng/instruments_zheng_2.csv'))
pqtls <- pqtls[,c("phenotype", "Chrom",  "Pos", "Name",  "SNP", "Type", 'Cis/Trans', 
                  "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
pqtls <- filter(pqtls, pqtls$Type != 'Tier3')
pqtls$phenotype <- str_c('p_', pqtls$phenotype)
pqtls$pval.exposure <- NULL # errors w pvals, will infer later
pqtls$chrpos <- gsub('chr', '', pqtls$Name)
preec_out <- as.data.frame(fread('outzheng/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls <- pqtls[which(paste(pqtls$SNP) %in% paste(preec_out$SNP)),] # same
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_trans_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_trans_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_trans_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())
#### MR - cis ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivszheng_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    # pval_col="pval.exposure",
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
                           filename = "~/Desktop/preec/outzheng/preec_out_unadj.csv", 
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
names(join_list) <- gsub("p_","out_",genelist) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
outlist <- gsub("p_","out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- gsub("p_","har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # 
unlink("/Volumes/MADDY2/datasets/preeclampsia/zheng_harm_clump_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/zheng_harm_clump_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/zheng_harm_clump_unadj_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/zheng_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
rm(files, data_list)
clplist<-ls(pattern = "clp_", mget(ls()))


# rsq and save
nsnp <- data.frame(do.call("rbind", lapply(ls(), function(x) {
  obj = get(x)
  if (class(obj) == "data.frame")
    c(protein = x, nsnp = nrow(obj))
})))
rsq <- function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  sum(2*dat$eaf.exposure*(1-dat$eaf.exposure)*(dat$beta.exposure^2))
}
rsqlist<-as.data.frame(sapply(clplist, rsq, simplify = FALSE))
rsqlist<-as.data.frame(t(rsqlist))
rsqlist$protein <- rownames(rsqlist)
rsqlist <- merge(rsqlist, nsnp, by = 'protein')
colnames(rsqlist) <- c('protein', 'rsq', 'nsnp')
rsqlist$nsnp <- as.numeric(rsqlist$nsnp)
rsqlist$fstat <- ((20330-rsqlist$nsnp-1)/rsqlist$nsnp) * (rsqlist$rsq/(1-rsqlist$rsq))
rsqlist <- rsqlist[order(rsqlist$fstat),]
write.csv(rsqlist, "~/desktop/preec/reszheng_preec/cis_rsqfstat.csv")
rm(rsqlist, rsq, nsnp)

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/reszheng_preec/individual_results_cis_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/reszheng_preec/individual_results_cis_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/reszheng_preec/individual_results_cis_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}
rm(list=ls())

#### Merge results - cis ####
setwd("~/Desktop/preec/reszheng_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/reszheng_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
mergedres$exposure <- gsub('p_', '', mergedres$exposure )

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'
write.csv(mergedres, "~/Desktop/preec/reszheng_preec/all_merged_cis_preec.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[order(mergedres$pval),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
write.csv(mergedres, "main_cis_preec.csv", row.names = FALSE)


mergedres <- filter(mergedres, mergedres$padj< 0.05)
write.csv(mergedres, "fdrsig_cis_preec.csv", row.names = FALSE)

rm(list=ls())




#### Plot results - cis ####

# Manhattan plot NB GWS = 1.5e-5
mergedres <- as.data.frame(fread('main_cis_preec.csv'))
sigp <- 0.05/nrow(mergedres)
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$gws <- ifelse(mergedres$pval < sigp, 1, 0)
mergedres$logp <- -log10(mergedres$pval)
cols <- c("0" = "grey32","1"= "coral")
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]

p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) +
  geom_hline(yintercept = -log10(sigp), color = "gray32", linetype = "dashed", linewidth = 1, alpha = 0.5)
pdf('zheng_cis_preec.pdf', width = 12, height = 6)
p1
dev.off()

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]
p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) 
pdf('proteome_cis_preec_fdrsig.pdf', width = 16, height = 6)
p1
dev.off()

rm(list=ls())


#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------ UKB  --------------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####
#### ------------------- **** preec **** --------------- ####
#### Extract protein list ####

# make long vs short name file
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('sunukb/pqtls_all.csv')) 
dput(colnames(pqtls))

pqtls$pval <- pqtls$log10p

pqtls$chrpos_hg37_a0_a1 <- gsub(':imp:v1', '', pqtls$chrpos_hg37_a0_a1)
names <- t(as.data.frame(strsplit(pqtls$chrpos_hg37_a0_a1, ":")))
pqtls$effectAllele <- names[,4]
pqtls$otherAllele <- names[,3]
pqtls$rsid <- ifelse(pqtls$rsid == '-', str_c(pqtls$chr, ':', pqtls$pos_hg38), pqtls$rsid)


setnames(pqtls, old = c('rsid', 'BETA', 'SE', 'pval',  'effectAllele', 'otherAllele', 'A1FREQ', 'protein', 'chr', 'pos_hg38'),
         new = c('SNP', 'beta.exposure', 'se.exposure', 'pval.exposure', 'effect_allele.exposure','other_allele.exposure', 'eaf.exposure', 'phenotype', 'chr.exposure', 'pos.exposure'))
pqtls_unadjusted_trans <- 
  pqtls[,c('SNP', 'phenotype', 'chr.exposure', 'pos.exposure', 'effect_allele.exposure', 'other_allele.exposure', 
           'eaf.exposure','beta.exposure', 'se.exposure', 'pval.exposure')]
write.csv(pqtls_unadjusted_trans, 'sunukb/pqtls_unadjusted_instruments.csv')

pqtls_unadjusted_cis <- filter(pqtls, pqtls$`cis/trans` == 'cis')
pqtls_unadjusted_cis <- 
  pqtls_unadjusted_cis[,c('SNP', 'phenotype', 'chr.exposure', 'pos.exposure', 'effect_allele.exposure', 'other_allele.exposure', 
                          'eaf.exposure','beta.exposure', 'se.exposure', 'pval.exposure')]
write.csv(pqtls_unadjusted_cis, 'sunukb/pqtls_unadjusted_instruments_cis.csv')

rm(list=ls())

#### Format outcome data ####
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('sunukb/pqtls_unadjusted_instruments.csv'))[,c('SNP', 'chr.exposure', 'pos.exposure')]
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)

preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt', 
                             drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))

preec$chrpos <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos'))
preec_out$phenotype <- 'preec'
preec_out$SNP <- ifelse(preec_out$SNP == '-',preec_out$chrpos, preec_out$SNP)
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, 'outsunukb/preec_out_unadj.csv')
rm(list=ls())

#### Separate and save - NB significance P<1.7×10−11 ####

# Trans
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('sunukb/pqtls_unadjusted_instruments.csv'))[,-1]
preec_out <- as.data.frame(fread('outsunukb/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)
pqtls <- pqtls[which(paste(pqtls$chrpos) %in% paste(preec_out$chrpos)),]
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_trans_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_trans_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_trans_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())

# cis
setwd('~/desktop/preec')
pqtls <- as.data.frame(fread('sunukb/pqtls_unadjusted_instruments_cis.csv'))[,-1]
preec_out <- as.data.frame(fread('outsunukb/preec_out_unadj.csv'))[,c('SNP', 'chrpos')]
pqtls$chrpos <- str_c(pqtls$chr.exposure, ':', pqtls$pos.exposure)
pqtls <- pqtls[which(paste(pqtls$chrpos) %in% paste(preec_out$chrpos)),]
pqlts_list <- split(pqtls, f = pqtls$phenotype)
unlink("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_cis_preec")
sapply(names(pqlts_list), 
       function (x) write.csv(pqlts_list[[x]], file=paste(x, "csv", sep=".") ))
rm(list=ls())

#### MR - cis ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
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
                           filename = "~/Desktop/preec/outsunukb/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

unlink("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_unadj_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_unadj_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_unadj_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
clplist<-ls(pattern = "clp_", mget(ls()))

# rsq and save
nsnp <- data.frame(do.call("rbind", lapply(ls(), function(x) {
  obj = get(x)
  if (class(obj) == "data.frame")
    c(protein = x, nsnp = nrow(obj))
})))
rsq <- function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  sum(2*dat$eaf.exposure*(1-dat$eaf.exposure)*(dat$beta.exposure^2))
}
rsqlist<-as.data.frame(sapply(clplist, rsq, simplify = FALSE))
rsqlist<-as.data.frame(t(rsqlist))
rsqlist$protein <- rownames(rsqlist)
rsqlist <- merge(rsqlist, nsnp, by = 'protein')
colnames(rsqlist) <- c('protein', 'rsq', 'nsnp')
rsqlist$nsnp <- as.numeric(rsqlist$nsnp)
rsqlist$fstat <- ((34557-rsqlist$nsnp-1)/rsqlist$nsnp) * (rsqlist$rsq/(1-rsqlist$rsq))
rsqlist <- rsqlist[order(rsqlist$fstat),]
write.csv(rsqlist, "~/desktop/preec/ressunukb_preec/cis_rsqfstat.csv")
rm(rsqlist, rsq, nsnp)

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}
beep(2)
rm(list=ls())

#### Merge results - cis ####
setwd("~/desktop/preec/ressunukb_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))


mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres[,c(1:3)] <- NULL
mergedres$outcome <- 'Pre-eclampsia'

mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec/', '', mergedres$filename))
mergedres$filename <- NULL

write.csv(mergedres, "~/Desktop/preec/ressunukb_preec/all_merged_cis_preec.csv", row.names = FALSE)

mergedres <- rbind((filter(mergedres, mergedres$method == 'Inverse variance weighted')), (filter(mergedres, mergedres$method == 'Wald ratio')))
mergedres <- mergedres[order(mergedres$pval),]
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr')
write.csv(mergedres, "main_cis_preec.csv", row.names = FALSE)


mergedres <- filter(mergedres, mergedres$padj< 0.05)
mergedres
write.csv(mergedres, "fdrsig_cis_preec.csv", row.names = FALSE)

rm(list=ls())




#### Plot results - cis ####
setwd("~/desktop/preec/ressunukb_preec")
# Manhattan plot NB GWS = 1.5e-5
mergedres <- as.data.frame(fread('main_cis_preec.csv'))
sigp <- 0.05/nrow(mergedres)
mergedres <- mergedres[order(mergedres$exposure),]
mergedres$gws <- ifelse(mergedres$pval < sigp, 1, 0)
mergedres$logp <- -log10(mergedres$pval)
cols <- c("0" = "grey32","1"= "coral")
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]

p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) +
  geom_hline(yintercept = -log10(sigp), color = "gray32", linetype = "dashed", linewidth = 1, alpha = 0.5)
pdf('sunukb_cis_preec.pdf', width = 12, height = 6)
p1
dev.off()

mergedres$gws <- ifelse(mergedres$padj < 0.05, 1, 0)
labels <- filter(mergedres, mergedres$gws == 1)[,c('exposure')]
p1 <- ggplot(mergedres, aes(x = exposure, y = logp, colour = gws)) +
  geom_point() + 
  theme_classic() + 
  scale_discrete_manual(cols) +
  scale_x_discrete(breaks = labels) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 10),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(color = "exposure", size = "Effect size", x = "", y = "-log10(p-value)") +
  ylim(c(0, 12)) 
pdf('proteome_cis_preec_fdrsig.pdf', width = 16, height = 6)
p1
dev.off()

rm(list=ls())


#### --------------------------------------------------------------------------------------------- ####
#### ------------------------------------------ PLOTTING ----------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####
#### HEATMAP - Preec cis&trans ####

### HEATMAP - Main & sensitivity ###
setwd("~/Desktop/preec/resdecode_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
decodecis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
decodecis_preec$analysis <- 'deCODE\n(SomaScan)'
decodecis_preec$padj <- p.adjust(decodecis_preec$pval, method = 'fdr')

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec', pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
sunukbcis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rm(files, alldat)
sunukbcis_preec$analysis <- 'UK Biobank\n(Olink)'
sunukbcis_preec$padj <- p.adjust(sunukbcis_preec$pval, method = 'fdr')
sunukbcis_preec$exposure <- gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec/', '', sunukbcis_preec$filename)
sunukbcis_preec$exposure <- gsub('_res.csv', '', sunukbcis_preec$exposure)


decodecis_preec_sig <- filter(decodecis_preec, decodecis_preec$padj <0.05 & decodecis_preec$method == 'Inverse variance weighted' |
                                decodecis_preec$padj <0.05 & decodecis_preec$method == 'Wald ratio')
sunukbcis_preec_sig <- filter(sunukbcis_preec, sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Inverse variance weighted' |
                                sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Wald ratio')

genelist <- rbind(decodecis_preec_sig, sunukbcis_preec_sig)
genelist <- genelist[!duplicated(genelist$exposure),]

mergedres <- rbind(decodecis_preec,  sunukbcis_preec)
rownames(mergedres) <- NULL
mergedres$filename <- NULL
mergedres$X <- NULL
mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted' | mergedres$method == 'Wald ratio')

mergedres$exposure <- gsub('p_', '', mergedres$exposure)
mergedres$exposure <- gsub('Pregnancy zone protein', 'PZP', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('MIC-1', 'GDF15', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('LPH', 'LCT', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('EBI3_IL27', 'IL27', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('Protease nexin I', 'SERPINE2', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('ATS13', 'ADAMTS13', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('NPPB', 'BNP', mergedres$exposure, fixed = TRUE)
genelist$exposure <- gsub('p_', '', genelist$exposure)
genelist$exposure <- gsub('Pregnancy zone protein', 'PZP', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('MIC-1', 'GDF15', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('LPH', 'LCT', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('EBI3_IL27', 'IL27', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('Protease nexin I', 'SERPINE2', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('ATS13', 'ADAMTS13', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('NPPB', 'BNP', genelist$exposure, fixed = TRUE)

write.csv(mergedres, '~/Desktop/preec/resproteomic/proteomic_mergedres_main_preec.csv')

mergedres <- mergedres[which(mergedres$exposure %in% genelist$exposure),]

mergedres$OR <- round(exp(mergedres$b), 2)
mergedres$LCI <- round(exp(mergedres$b - 1.96*mergedres$se), 2)
mergedres$UCI <- round(exp(mergedres$b + 1.96*mergedres$se), 2)
mergedres$AdjPround <- round(mergedres$padj, 3)
write.csv(mergedres, '~/Desktop/preec/resproteomic/proteomic_mergedres_mainandreplication_preec.csv')

# mergedres <- as.data.frame(fread('~/Desktop/preec/resproteomic/proteomic_mergedres_mainandreplication_preec.csv'))

mergedres$Protein <- mergedres$exposure
mergedres$Analysis <- mergedres$analysis

mergedres$signif <- ifelse(mergedres$pval <1,'ns',  '')
mergedres$signif<- ifelse(mergedres$pval <0.05,'*', mergedres$signif)
mergedres$signif<- ifelse(mergedres$pval <0.001,'**', mergedres$signif)
mergedres$signif<- ifelse(mergedres$padj <0.05, '***', mergedres$signif)

mergedres <- mergedres %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres$signif!='ns'  ~ "1 Positive",
  b<0 & mergedres$signif!='ns' ~ "2 Negative",
  b>0 & mergedres$signif=='ns'  ~ "3 No significant association",
  b<0 & mergedres$signif=='ns' ~ "3 No significant association",
  is.na(b) ~ as.character(NA)))
mergedres$signif <- gsub('ns', '', mergedres$signif)

cor <- ggplot(data = mergedres, aes(x=Analysis, y=Protein, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Positive", "Negative", "No significant association", 'Missing')) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/resproteomic/heatmap_preec_main_sens.pdf", width = 5, height = 5.2)
cor
dev.off()

#### HEATMAP wih sensitivity ####

### HEATMAP - Main & sensitivity ###
setwd("~/Desktop/preec/resdecode_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
decodecis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
decodecis_preec$analysis <- 'deCODE\n(SomaScan)'
decodecis_preec$padj <- p.adjust(decodecis_preec$pval, method = 'fdr')

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/reszheng_preec/individual_results_cis_preec', pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
zhengcis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rm(files, alldat)
zhengcis_preec$analysis <- 'Zheng et al\n(SomaScan)'
zhengcis_preec$padj <- p.adjust(zhengcis_preec$pval, method = 'fdr')
zhengcis_preec$exposure <- gsub('p_', '', zhengcis_preec$exposure)

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec', pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
sunukbcis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rm(files, alldat)
sunukbcis_preec$analysis <- 'UK Biobank\n(Olink)'
sunukbcis_preec$padj <- p.adjust(sunukbcis_preec$pval, method = 'fdr')
sunukbcis_preec$exposure <- gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec/', '', sunukbcis_preec$filename)
sunukbcis_preec$exposure <- gsub('_res.csv', '', sunukbcis_preec$exposure)

decodecis_preec_sig <- filter(decodecis_preec, decodecis_preec$padj <0.05 & decodecis_preec$method == 'Inverse variance weighted' |
                                decodecis_preec$padj <0.05 & decodecis_preec$method == 'Wald ratio')
zhengcis_preec_sig <- filter(zhengcis_preec, zhengcis_preec$padj <0.05 & zhengcis_preec$method == 'Inverse variance weighted' |
                               zhengcis_preec$padj <0.05 & zhengcis_preec$method == 'Wald ratio')
sunukbcis_preec_sig <- filter(sunukbcis_preec, sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Inverse variance weighted' |
                                sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Wald ratio')

genelist <- rbind(decodecis_preec_sig, sunukbcis_preec_sig)
genelist <- genelist[!duplicated(genelist$exposure),]

mergedres <- rbind(decodecis_preec, zhengcis_preec,  sunukbcis_preec)
rownames(mergedres) <- NULL
mergedres$filename <- NULL
mergedres$X <- NULL
mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted' | mergedres$method == 'Wald ratio')

mergedres$exposure <- gsub('p_', '', mergedres$exposure)
mergedres$exposure <- gsub('Pregnancy zone protein', 'PZP', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('MIC-1', 'GDF15', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('LPH', 'LCT', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('EBI3_IL27', 'IL27', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('Protease nexin I', 'SERPINE2', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('ATS13', 'ADAMTS13', mergedres$exposure, fixed = TRUE)
mergedres$exposure <- gsub('NPPB', 'BNP', mergedres$exposure, fixed = TRUE)
genelist$exposure <- gsub('p_', '', genelist$exposure)
genelist$exposure <- gsub('Pregnancy zone protein', 'PZP', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('MIC-1', 'GDF15', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('LPH', 'LCT', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('EBI3_IL27', 'IL27', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('Protease nexin I', 'SERPINE2', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('ATS13', 'ADAMTS13', genelist$exposure, fixed = TRUE)
genelist$exposure <- gsub('NPPB', 'BNP', genelist$exposure, fixed = TRUE)

mergedres <- mergedres[which(mergedres$exposure %in% genelist$exposure),]

mergedres$OR <- round(exp(mergedres$b), 2)
mergedres$LCI <- round(exp(mergedres$b - 1.96*mergedres$se), 2)
mergedres$UCI <- round(exp(mergedres$b + 1.96*mergedres$se), 2)
mergedres$AdjPround <- round(mergedres$padj, 3)

mergedres$Protein <- mergedres$exposure
mergedres$Analysis <- mergedres$analysis

mergedres$signif <- ifelse(mergedres$pval <1,'ns',  '')
mergedres$signif<- ifelse(mergedres$pval <0.05,'*', mergedres$signif)
mergedres$signif<- ifelse(mergedres$pval <0.001,'**', mergedres$signif)
mergedres$signif<- ifelse(mergedres$padj <0.05, '***', mergedres$signif)

mergedres <- mergedres %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres$signif!='ns'  ~ "1 Positive",
  b<0 & mergedres$signif!='ns'~ "2 Negative",
  b>0 & mergedres$signif=='ns' ~ "3 No significant association",
  b<0 & mergedres$signif=='ns'~ "3 No significant association",
  is.na(b) ~ as.character(NA)))
mergedres$signif <- gsub('ns', '', mergedres$signif)

cor <- ggplot(data = mergedres, aes(x=Analysis, y=Protein, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Positive", "Negative", "No significant association", 'Missing')) +
  geom_text(aes(label = signif), color = "black", size = 4) + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/resproteomic/sens_withzheng.pdf", width = 6.5, height = 5.2)
cor
dev.off()

#### HEATMAP - Main & Coloc - decode ####
head(mergedres)
coloc <- as.data.frame(fread('~/Desktop/preec/resdecode_preec/decode_coloc_sighits.csv'))
coloc <- coloc[,c("nsnps", "PP.H1.abf",  "PP.H2.abf",  "PP.H3.abf",  "PP.H4.abf", "exposure")]
coloc$exposure
decodecis_preec_sig$exposure
coloc$exposure <- gsub('ALDHE2', 'ALDH-E2', coloc$exposure, fixed = TRUE)
coloc$exposure <- gsub('JUND', 'jun-D', coloc$exposure, fixed = TRUE)
decodecis_preec_sig$exposure <- gsub('Pregnancy zone protein', 'PZP', decodecis_preec_sig$exposure, fixed = TRUE)
decodecis_preec_sig$exposure <- gsub('MIC-1', 'GDF15', decodecis_preec_sig$exposure, fixed = TRUE)
decodecis_preec_sig$exposure <- gsub('Protease nexin I', 'SERPINE2', decodecis_preec_sig$exposure, fixed = TRUE)
decodecis_preec_sig$exposure <- gsub('ATS13', 'ADAMTS13', decodecis_preec_sig$exposure, fixed = TRUE)

colocmerge <- merge(decodecis_preec_sig[,c('exposure', 'b', 'se', 'pval', 'padj')], coloc, by='exposure', all.x=TRUE, all.y=FALSE)
colocmerge$OR <- round(exp(colocmerge$b), 2)
colocmerge$LCI <- round(exp(colocmerge$b - 1.96*colocmerge$se), 2)
colocmerge$UCI <- round(exp(colocmerge$b + 1.96*colocmerge$se), 2)
colocmerge$AdjPround <- round(colocmerge$padj, 3)
colocmerge$`PP.H4.abf/(PP.H3+PP.H4)` <- colocmerge$`PP.H4.abf`/(colocmerge$`PP.H3.abf`+ colocmerge$`PP.H4.abf`)
write.csv(colocmerge, '~/Desktop/preec/resproteomic/mainres_significant_coloc_decode.csv')

# Plot colocalization results
colocmerge <- colocmerge[,c("PP.H1.abf",  "PP.H2.abf",  "PP.H3.abf",  "PP.H4.abf", "PP.H4.abf/(PP.H3+PP.H4)", "exposure")]
colnames(colocmerge) <- c("Causal for Trait 1\nonly (H1)",  "Causal for Trait 2\nonly (H2)",  "Causal variant\ndistinct (H3)",  "Causal variant\nshared (H4)", "H4/(H3+H4)", "exposure")
cor_colocplot <- as.data.frame(pivot_longer(colocmerge, 
                                            cols = c("Causal for Trait 1\nonly (H1)",  "Causal for Trait 2\nonly (H2)",  "Causal variant\ndistinct (H3)",  "Causal variant\nshared (H4)", "H4/(H3+H4)")))
colnames(cor_colocplot) <- c('Protein', 'Posterior probability of hypothesis (%)', 'Posterior probability (%)')

cor <- ggplot(data = cor_colocplot, aes(x=`Posterior probability of hypothesis (%)`, y=Protein, fill=`Posterior probability (%)`)) +
  geom_tile() + scale_fill_gradient2(low="white",  high="darkseagreen3", na.value = "grey93") +  
  geom_text(aes(label = str_c(round(`Posterior probability (%)`*100, 1), '%')), color = "black", size = 3) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
pdf("~/Desktop/preec/resproteomic/heatmap_preec_coloc_decode.pdf", width = 5, height = 5.2)
cor 
dev.off()

#### HEATMAP - Main & Coloc - sunukb ####

head(mergedres)
coloc <- as.data.frame(fread('~/Desktop/preec/ressunukb_preec/sunukb_coloc_sighits.csv')[,-1])
colnames(coloc) <- c('exposure', "nsnps", "PP.H0.abf", "PP.H1.abf",  "PP.H2.abf",  "PP.H3.abf",  "PP.H4.abf")
coloc$exposure
sunukbcis_preec_sig$exposure
colocmerge <- merge(sunukbcis_preec_sig[,c('exposure', 'b', 'se', 'pval', 'padj')], coloc, by='exposure', all.x=TRUE, all.y=FALSE)

colocmerge$OR <- round(exp(colocmerge$b), 2)
colocmerge$LCI <- round(exp(colocmerge$b - 1.96*colocmerge$se), 2)
colocmerge$UCI <- round(exp(colocmerge$b + 1.96*colocmerge$se), 2)
colocmerge$AdjPround <- round(colocmerge$padj, 3)
colocmerge$`PP.H4.abf/(PP.H3+PP.H4)` <- colocmerge$`PP.H4.abf`/(colocmerge$`PP.H3.abf`+ colocmerge$`PP.H4.abf`)
write.csv(colocmerge, '~/Desktop/preec/resproteomic/mainres_significant_coloc_sunukb.csv')

# Plot colocalization results
colocmerge <- colocmerge[,c("PP.H1.abf",  "PP.H2.abf",  "PP.H3.abf",  "PP.H4.abf", "PP.H4.abf/(PP.H3+PP.H4)", "exposure")]
colnames(colocmerge) <- c("Causal for Trait 1\nonly (H1)",  "Causal for Trait 2\nonly (H2)",  "Causal variant\ndistinct (H3)",  "Causal variant\nshared (H4)", "H4/(H3+H4)", "exposure")
cor_colocplot <- as.data.frame(pivot_longer(colocmerge, 
                                            cols = c("Causal for Trait 1\nonly (H1)",  "Causal for Trait 2\nonly (H2)",  "Causal variant\ndistinct (H3)",  "Causal variant\nshared (H4)", "H4/(H3+H4)")))
colnames(cor_colocplot) <- c('Protein', 'Posterior probability of hypothesis (%)', 'Posterior probability (%)')

cor <- ggplot(data = cor_colocplot, aes(x=`Posterior probability of hypothesis (%)`, y=Protein, fill=`Posterior probability (%)`)) +
  geom_tile() + scale_fill_gradient2(low="white",  high="darkseagreen3", na.value = "grey93") +  
  geom_text(aes(label = str_c(round(`Posterior probability (%)`*100, 1), '%')), color = "black", size = 3) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
pdf("~/Desktop/preec/resproteomic/heatmap_preec_coloc_sunukb.pdf", width = 5, height = 3.2)
cor 
dev.off()

#### HEATMAP - Proteome and relevant transcriptome ####
transcriptome <- as.data.frame(fread('~/Desktop/preec/restranscriptomic/transcriptomic_mergedres_main_preec.csv'))[,-1]
head(transcriptome)

transcriptome_select <- filter(transcriptome, 
                               # 3MG
                               transcriptome$Gene == 'MPG' | 
                                 transcriptome$Gene == 'ALDH2'| 
                                 transcriptome$Gene == 'NPPA'| 
                                 transcriptome$Gene == 'ADAMTS13'| 
                                 transcriptome$Gene == 'FGL1'| 
                                 transcriptome$Gene == 'JUND'| 
                                 transcriptome$Gene == 'MANEA'| 
                                 transcriptome$Gene == 'METAP1'| 
                                 transcriptome$Gene == 'GDF15'| 
                                 transcriptome$Gene == 'NOTUM'| 
                                 transcriptome$Gene == 'PZP'|
                                 transcriptome$Gene == 'SERPINE2'| 
                                 transcriptome$Gene == 'RGS18' |
                                 transcriptome$Gene == 'SULT1A1'|
                                 transcriptome$Gene == 'SH2B3'|
                                 transcriptome$Gene == 'FGF5'|
                                 transcriptome$Gene == 'FES'|
                                 transcriptome$Gene == 'APOBR')

transcriptome_select<- transcriptome_select[,c('Gene', 'Tissue', 'b', 'se', 'pval', 'padj')]
transcriptome_select$Protein <- transcriptome_select$Gene

setwd("~/Desktop/preec/resdecode_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
decodecis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(decodecis_preec) <- NULL
decodecis_preec$analysis <- 'deCODE\n(SomaScan)'
decodecis_preec$Protein <- decodecis_preec$exposure
decodecis_preec$padj <- p.adjust(decodecis_preec$pval, method = 'fdr')
decodecis_preec$Protein <- gsub('3MG', 'MPG', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('ALDH-E2', 'ALDH2', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('NPPA', 'ANP', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('ATS13', 'ADAMTS13', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('jun-D', 'JUND', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('MIC-1', 'GDF15', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('Pregnancy zone protein', 'PZP', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('Protease nexin I', 'SERPINE2', decodecis_preec$Protein, fixed = TRUE)

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec', pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
sunukbcis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rm(files, alldat)
rownames(sunukbcis_preec) <- NULL
sunukbcis_preec$analysis <- 'UK Biobank\n(Olink)'
sunukbcis_preec$padj <- p.adjust(sunukbcis_preec$pval, method = 'fdr')
sunukbcis_preec$Protein <- gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec/', '', sunukbcis_preec$filename)
sunukbcis_preec$Protein <- gsub('_res.csv', '', sunukbcis_preec$Protein)
sunukbcis_preec$Protein <- gsub('3MG', 'MPG', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('ALDH-E2', 'ALDH2', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('NPPA', 'ANP', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('ATS13', 'ADAMTS13', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('jun-D', 'JUND', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('MIC-1', 'GDF15', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('Pregnancy zone protein', 'PZP', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('Protease nexin I', 'SERPINE2', sunukbcis_preec$Protein, fixed = TRUE)

decodecis_preec_sig <- filter(decodecis_preec, decodecis_preec$padj <0.05 & decodecis_preec$method == 'Inverse variance weighted' |
                                decodecis_preec$padj <0.05 & decodecis_preec$method == 'Wald ratio')
sunukbcis_preec_sig <- filter(sunukbcis_preec, sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Inverse variance weighted' |
                                sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Wald ratio')

genelist <- rbind(decodecis_preec_sig[,c('b', 'se', 'pval', 'padj', 'Protein')], sunukbcis_preec_sig[,c('b', 'se', 'pval', 'padj', 'Protein')])

proteome_res <- genelist[,c('Protein', 'b', 'se', 'pval', 'padj')]
proteome_res <- proteome_res[!duplicated(proteome_res$Protein),]
proteome_res$Tissue <- '1 pQTL result'
# genenames <- c('MPG','ALDH2','NPPA', 'ADAMTS13','FGL1','JUND','MANEA','METAP1', 'GDF15','NOTUM', 'PZP','SERPINE2','RGS18' )
# proteinnames <- c('3MG','ALDH-E2','ANP', 'ATS13','FGL1','jun-D','MANEA','METAP1', 'MIC-1','NOTUM', 'Pregnancy zone protein','Protease nexin I','RGS18' )
proteome_res$Gene <- proteome_res$Protein
rownames(proteome_res) <- NULL
proteome_res <- proteome_res[,c('Gene', 'Tissue', 'b', 'se', 'pval', 'padj', 'Protein')]

transcriptome_select <- filter(transcriptome_select, transcriptome_select$Tissue != 'Placenta')

transcriptome_select <- rbind(proteome_res, transcriptome_select)
proteome_res$outcome <- 'preec'
mergedres_transcriptome <- merge(proteome_res[,c('Protein', 'outcome')], transcriptome_select, by = 'Protein', all.x=TRUE, all.y=FALSE)
mergedres_transcriptome$Protein <- gsub('Pregnancy zone protein', 'PZP', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('MIC-1', 'GDF15', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('Protease nexin I', 'SERPINE2', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('ATS13', 'ADAMTS13', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('MPG', '3MG', mergedres_transcriptome$Protein, fixed = TRUE)

mergedres_transcriptome$Gene <- gsub('Pregnancy zone protein', 'PZP', mergedres_transcriptome$Gene, fixed = TRUE)
mergedres_transcriptome$Gene <- gsub('MIC-1', 'GDF15', mergedres_transcriptome$Gene, fixed = TRUE)
mergedres_transcriptome$Gene <- gsub('Protease nexin I', 'SERPINE2', mergedres_transcriptome$Gene, fixed = TRUE)
mergedres_transcriptome$Gene <- gsub('ATS13', 'ADAMTS13', mergedres_transcriptome$Gene, fixed = TRUE)
mergedres_transcriptome$Gene <- gsub('MPG', '3MG', mergedres_transcriptome$Protein, fixed = TRUE)


mergedres_transcriptome$signif <- ifelse(mergedres_transcriptome$pval <1,'ns',  '')
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$pval <0.05,'*', mergedres_transcriptome$signif)
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$pval <0.001,'**', mergedres_transcriptome$signif)
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$padj <0.05, '***', mergedres_transcriptome$signif)

mergedres_transcriptome <- mergedres_transcriptome %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres_transcriptome$signif!='ns'  ~ "1 Positive",
  b<0 & mergedres_transcriptome$signif!='ns'~ "2 Negative",
  b>0 & mergedres_transcriptome$signif=='ns' ~ "3 No significant association",
  b<0 & mergedres_transcriptome$signif=='ns'~ "3 No significant association",
  is.na(b) ~ as.character(NA)))
mergedres_transcriptome$signif <- gsub('ns', '', mergedres_transcriptome$signif)

mergedres_transcriptome$Tissue <- replace_na(mergedres_transcriptome$Tissue, 'Artery Aorta')

# Name protein and all genes in interaction network & pathway (via GeneCards / STRIING protein interaction network)

cor <- ggplot(data = mergedres_transcriptome, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Positive", "Negative", "No significant association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/resproteomic/heatmap_preec_main_withtranscript_strict.pdf", width = 8.5, height = 5.2)
cor
dev.off()


### HEATMAP - Proteome and relevant transcriptome

transcriptome <- as.data.frame(fread('~/Desktop/preec/restranscriptomic/transcriptomic_mergedres_main_preec.csv'))[,-1]
head(transcriptome)


setwd("~/Desktop/preec/resdecode_preec")
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/individual_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
decodecis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(decodecis_preec) <- NULL
decodecis_preec$analysis <- 'deCODE\n(SomaScan)'
decodecis_preec$Protein <- decodecis_preec$exposure
decodecis_preec$padj <- p.adjust(decodecis_preec$pval, method = 'fdr')
decodecis_preec$Protein <- gsub('3MG', 'MPG', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('ALDH-E2', 'ALDH2', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('NPPA', 'ANP', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('ATS13', 'ADAMTS13', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('jun-D', 'JUND', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('MIC-1', 'GDF15', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('Pregnancy zone protein', 'PZP', decodecis_preec$Protein, fixed = TRUE)
decodecis_preec$Protein <- gsub('Protease nexin I', 'SERPINE2', decodecis_preec$Protein, fixed = TRUE)

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec', pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
sunukbcis_preec <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rm(files, alldat)
rownames(sunukbcis_preec) <- NULL
sunukbcis_preec$analysis <- 'UK Biobank\n(Olink)'
sunukbcis_preec$padj <- p.adjust(sunukbcis_preec$pval, method = 'fdr')
sunukbcis_preec$Protein <- gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/individual_results_cis_preec/', '', sunukbcis_preec$filename)
sunukbcis_preec$Protein <- gsub('_res.csv', '', sunukbcis_preec$Protein)
sunukbcis_preec$Protein <- gsub('3MG', 'MPG', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('ALDH-E2', 'ALDH2', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('NPPA', 'ANP', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('ATS13', 'ADAMTS13', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('jun-D', 'JUND', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('MIC-1', 'GDF15', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('Pregnancy zone protein', 'PZP', sunukbcis_preec$Protein, fixed = TRUE)
sunukbcis_preec$Protein <- gsub('Protease nexin I', 'SERPINE2', sunukbcis_preec$Protein, fixed = TRUE)

decodecis_preec_sig <- filter(decodecis_preec, decodecis_preec$padj <0.05 & decodecis_preec$method == 'Inverse variance weighted' |
                                decodecis_preec$padj <0.05 & decodecis_preec$method == 'Wald ratio')
sunukbcis_preec_sig <- filter(sunukbcis_preec, sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Inverse variance weighted' |
                                sunukbcis_preec$padj <0.05 & sunukbcis_preec$method == 'Wald ratio')

genelist <- rbind(decodecis_preec_sig[,c('b', 'se', 'pval', 'padj', 'Protein')], sunukbcis_preec_sig[,c('b', 'se', 'pval', 'padj', 'Protein')])
proteome_res <- genelist[,c('Protein', 'b', 'se', 'pval', 'padj')]
proteome_res <- proteome_res[!duplicated(proteome_res$Protein),]
proteome_res$Tissue <- '1 pQTL result'
# genenames <- c('MPG','ALDH2','NPPA', 'ADAMTS13','FGL1','JUND','MANEA','METAP1', 'GDF15','NOTUM', 'PZP','SERPINE2','RGS18' )
# proteinnames <- c('3MG','ALDH-E2','ANP', 'ATS13','FGL1','jun-D','MANEA','METAP1', 'MIC-1','NOTUM', 'Pregnancy zone protein','Protease nexin I','RGS18' )
proteome_res$Gene <- proteome_res$Protein
rownames(proteome_res) <- NULL
proteome_res <- proteome_res[,c('Gene', 'Tissue', 'b', 'se', 'pval', 'padj', 'Protein')]

transcriptome_select <- filter(transcriptome, 
                               # 3MG
                               transcriptome$Gene == 'MPG' | 
                                 # ALDH-E2
                                 transcriptome$Gene == 'ALDH2'| transcriptome$Gene == 'AOC3'| transcriptome$Gene == 'ADH1C'| transcriptome$Gene == 'ADH1B'| transcriptome$Gene == 'MAOA'| transcriptome$Gene == 'ADH4'| 
                                 # ANP
                                 transcriptome$Gene == 'NPPA'| transcriptome$Gene == 'NPR1'| transcriptome$Gene == 'NPR2'| transcriptome$Gene == 'NPR3'| transcriptome$Gene == 'NPPB'| transcriptome$Gene == 'NPPC'| 
                                 # ADAMTS13
                                 transcriptome$Gene == 'ADAMTS13'| transcriptome$Gene == 'F8'| transcriptome$Gene == 'GP1BA'|  transcriptome$Gene == 'VWF'| 
                                 # FGL1
                                 transcriptome$Gene == 'FGL1'| transcriptome$Gene == 'FGG'| transcriptome$Gene == 'LAG3'| transcriptome$Gene == 'FGA'| 
                                 # jun-D
                                 transcriptome$Gene == 'JUND'| transcriptome$Gene == 'FOS'| transcriptome$Gene == 'FOSB'| transcriptome$Gene == 'FOSL1'| transcriptome$Gene == 'FOSL2'| transcriptome$Gene == 'JUNB'| 
                                 # MANEA
                                 transcriptome$Gene == 'MANEA'| transcriptome$Gene == 'MAN1B1'| 
                                 # METAP1
                                 transcriptome$Gene == 'METAP1'| transcriptome$Gene == 'METAP2'| 
                                 # GDF15
                                 transcriptome$Gene == 'GDF15'| transcriptome$Gene == 'RET'| transcriptome$Gene == 'GFRAL'| transcriptome$Gene == 'TGFBR2'| 
                                 # NOTUM
                                 transcriptome$Gene == 'NOTUM'| transcriptome$Gene == 'WNT7A'| transcriptome$Gene == 'WNT16'| transcriptome$Gene == 'WNT3A'| transcriptome$Gene == 'WNT1'| transcriptome$Gene == 'WNT11'| 
                                 # PZP
                                 transcriptome$Gene == 'PZP'| 
                                 # SERPINE2
                                 transcriptome$Gene == 'SERPINE2'| transcriptome$Gene == 'PLAU'| transcriptome$Gene == 'F2'| transcriptome$Gene == 'SNX1'| 
                                 # RGS18
                                 transcriptome$Gene == 'RGS18'|  transcriptome$Gene == 'GNAQ'| transcriptome$Gene == 'GNB5' | 
                                 # FGF5
                                 transcriptome$Gene == 'FGF5'| transcriptome$Gene == 'FGFR1'| transcriptome$Gene == 'FGFR4'| transcriptome$Gene == 'EGF'| transcriptome$Gene == 'KL'| transcriptome$Gene == 'FGFR2'| 
                                 # FES
                                 transcriptome$Gene == 'FES'| transcriptome$Gene == 'PLXNA1'| transcriptome$Gene == 'PLXNA2'| transcriptome$Gene == 'DPYSL5'| transcriptome$Gene == 'DPYSL2'| transcriptome$Gene == 'PLXNA3'| 
                                 # SULT1A1
                                 transcriptome$Gene == 'SULT1A1'| transcriptome$Gene == 'SULT1A3'| transcriptome$Gene == 'SULT1A4'| transcriptome$Gene == 'SULT1A2'| transcriptome$Gene == 'SULT1C2'|  
                                 # SH2B3
                                 transcriptome$Gene == 'SH2B3'| transcriptome$Gene == 'NTRK3'| transcriptome$Gene == 'SH2B1'| transcriptome$Gene == 'JAK2'| transcriptome$Gene == 'GRB2'| transcriptome$Gene == 'NTRK1'| 
                                 # APOBR
                                 transcriptome$Gene == 'APOBR'| transcriptome$Gene == 'APOB'
                               
)

transcriptome_select <- transcriptome_select[,c('Gene', 'Tissue', 'b', 'se', 'pval', 'padj')]
transcriptome_select$Protein <- recode_factor(transcriptome_select$Gene, 
                                              'MPG' = 'MPG', 
                                              'MANEA' = 'MANEA', 'MAN1B1' = 'MANEA', 
                                              'NPPA' = 'ANP', 'NPR1' = 'ANP', 'NPR2' = 'ANP', 'NPR3' = 'ANP', 'NPPB' = 'ANP', 'NPPC' = 'ANP', 
                                              'ALDH2' = 'ALDH2','AOC3' = 'ALDH2','ADH1C' = 'ALDH2','ADH1B' = 'ALDH2','MAOA' = 'ALDH2','ADH4' = 'ALDH2',
                                              'ADAMTS13' = 'ADAMTS13', 'F8' = 'ADAMTS13', 'VWF' = 'ADAMTS13', 'GP1BA' = 'ADAMTS13', 
                                              'JUND' = 'jun-D', 'FOS' = 'jun-D', 'FOSB' = 'jun-D', 'FOSL1' = 'jun-D', 'FOSL2' = 'jun-D', 'JUNB' = 'jun-D', 
                                              'FGL1' = 'FGL1','FGG' = 'FGL1','LAG3' = 'FGL1','FGA' = 'FGL1',
                                              'METAP1' = 'METAP1',  'METAP2' = 'METAP1',  
                                              'GDF15' = 'GDF15', 'RET' = 'GDF15', 'GFRAL' = 'GDF15', 'TGFBR2' = 'GDF15',
                                              'NOTUM' = 'NOTUM', 'WNT7A' = 'NOTUM', 'WNT16' = 'NOTUM', 'WNT3A' = 'NOTUM', 'WNT1' = 'NOTUM', 'WNT11' = 'NOTUM', 
                                              'PZP' = 'PZP', 
                                              'SERPINE2' = 'SERPINE2', 'PLAU' = 'SERPINE2', 'F2' = 'SERPINE2', 'SNX1' = 'SERPINE2', 
                                              'RGS18' = 'RGS18','GNAQ' = 'RGS18','GNB5' = 'RGS18',
                                              'APOBR' = 'APOBR','APOB' = 'APOBR',
                                              'FGF5' = 'FGF5','FGFR1' = 'FGF5','FGFR4' = 'FGF5','EGF' = 'FGF5','KL' = 'FGF5','FGFR2' = 'FGF5',
                                              'FES' = 'FES','PLXNA1' = 'FES','PLXNA2' = 'FES','DPYSL5' = 'FES','DPYSL2' = 'FES','PLXNA3' = 'FES',
                                              'SULT1A1' = 'SULT1A1','SULT1A3' = 'SULT1A1','SULT1A4' = 'SULT1A1','SULT1A2' = 'SULT1A1','SULT1C2' = 'SULT1A1',
                                              'SH2B3' = 'SH2B3','NTRK3' = 'SH2B3','SH2B1' = 'SH2B3','JAK2' = 'SH2B3','GRB2' = 'SH2B3','NTRK1' = 'SH2B3' )
transcriptome_select <- filter(transcriptome_select, transcriptome_select$Tissue != 'Placenta')

proteome_res$Tissue <- '1 pQTL result'
# genenames <- c('MPG','ALDH2','NPPA', 'ADAMTS13','FGL1','JUND','MANEA','METAP1', 'GDF15','NOTUM', 'PZP','SERPINE2','RGS18' )
# proteinnames <- c('3MG','ALDH-E2','ANP', 'ATS13','FGL1','jun-D','MANEA','METAP1', 'MIC-1','NOTUM', 'Pregnancy zone protein','Protease nexin I','RGS18' )
proteome_res$Protein <- gsub('3MG', 'MPG', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('ALDH-E2', 'ALDH2', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('NPPA', 'ANP', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('ATS13', 'ADAMTS13', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('jun-D', 'jun-D', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('MIC-1', 'GDF15', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('Pregnancy zone protein', 'PZP', proteome_res$Protein, fixed = TRUE)
proteome_res$Protein <- gsub('Protease nexin I', 'SERPINE2', proteome_res$Protein, fixed = TRUE)
proteome_res$Gene <- proteome_res$Protein
rownames(proteome_res) <- NULL
transcriptome_select <- rbind(proteome_res, transcriptome_select)
proteome_res$outcome <- 'preec'

mergedres_transcriptome <- merge(proteome_res[,c('Protein', 'outcome')], transcriptome_select, by = 'Protein', all.x=TRUE, all.y=FALSE)
mergedres_transcriptome$Protein <- gsub('Pregnancy zone protein', 'PZP', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('MIC-1', 'GDF15', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('LPH', 'LCT', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('Protease nexin I', 'SERPINE2', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('ATS13', 'ADAMTS13', mergedres_transcriptome$Protein, fixed = TRUE)
mergedres_transcriptome$Protein <- gsub('MPG', '3MG', mergedres_transcriptome$Protein, fixed = TRUE)
write.csv(mergedres_transcriptome, '~/Desktop/preec/resproteomic/mainres_significant_transcriptomics.csv')


mergedres_transcriptome$signif <- ifelse(mergedres_transcriptome$pval <1,'ns',  '')
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$pval <0.05,'*', mergedres_transcriptome$signif)
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$pval <0.001,'**', mergedres_transcriptome$signif)
mergedres_transcriptome$signif<- ifelse(mergedres_transcriptome$padj <0.05, '***', mergedres_transcriptome$signif)

mergedres_transcriptome <- mergedres_transcriptome %>% mutate(`Direction of association` = case_when(
  b>0 & mergedres_transcriptome$signif!='ns'  ~ "1 Positive",
  b<0 & mergedres_transcriptome$signif!='ns'~ "2 Negative",
  b>0 & mergedres_transcriptome$signif=='ns' ~ "3 No significant association",
  b<0 & mergedres_transcriptome$signif=='ns'~ "3 No significant association",
  is.na(b) ~ as.character(NA)))
mergedres_transcriptome$signif <- gsub('ns', '', mergedres_transcriptome$signif)

mergedres_transcriptome$Tissue <- replace_na(mergedres_transcriptome$Tissue, 'Artery Aorta')

# Name protein and all genes in interaction network & pathway (via GeneCards / STRIING protein interaction network)
mergedres_transcriptome$eQTL <- str_c(mergedres_transcriptome$Protein, ' -                      ', mergedres_transcriptome$Gene)
cor <- ggplot(data = mergedres_transcriptome, aes(x=Tissue, y=eQTL, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Positive", "Negative", "No significant association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

pdf("~/Desktop/preec/resproteomic/heatmap_preec_main_withtranscript_withinteractions.pdf", width = 10.5, height = 12.2)
cor
dev.off()


# Name protein and all genes in interaction network & pathway (via GeneCards / STRIING protein interaction network)
mergedres_transcriptome <- filter(mergedres_transcriptome, mergedres_transcriptome$Tissue != '1 pQTL result')
cor <- ggplot(data = mergedres_transcriptome, aes(x=Tissue, y=Gene, fill=`Direction of association`)) +
  geom_tile() + 
  scale_fill_brewer(palette = "Pastel2", labels = c("Positive", "Negative", "No significant association", 'Not available')) +
  geom_text(aes(label = signif), color = "black", size = 3) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 39, vjust = 1, hjust=1), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(size = 8, angle = 0)) +
  facet_grid(Protein~., scales = 'free', drop = TRUE,space = 'free')

pdf("~/Desktop/preec/resproteomic/heatmap_preec_main_withtranscript_withinteractions_facet.pdf", width = 11, height = 10.8)
cor
dev.off()




rm(list=ls())







#### Merge F-stats ####
fst1 <- as.data.frame(fread('~/desktop/preec/resdecode_preec/cis_rsqfstat.csv')[,-1])
fst1$Study <- 'deCODE'
fst2 <- as.data.frame(fread('~/desktop/preec/reszheng_preec/cis_rsqfstat.csv')[,-1])
fst2$Study <- 'Zheng et al'
fst3 <- as.data.frame(fread('~/desktop/preec/ressunukb_preec/cis_rsqfstat.csv')[,-1])
fst3$Study <- 'UK Biobank'

fst <- rbind(fst1, fst2)
fst <- rbind(fst, fst3)

fst$protein <- gsub('clp_', '', fst$protein)
rownames(fst) <- NULL
head(fst)
fst$fstat <- round(fst$fstat, 2)
fst$rsq <- round(fst$rsq, 5)

colnames(fst) <- c('Protein', 'R-squared', 'Number of SNPs', 'F-statistic', 'Data source')
fst <- fst[,c('Data source','Protein', 'R-squared', 'Number of SNPs', 'F-statistic')]


write.csv(fst, '~/desktop/preec/resproteomic/fstats.csv', row.names = FALSE)

rm(list=ls())

#### -----------------------------------   PHENOSCAN LOOKUP ---------------------------------------  ####
library(phenoscanner)

# Load decode
setwd("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","_decode",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
rm(files, data_list)

# load ukb
setwd("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","_sunukb", list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)

# Select
clplist <- c('clp_DNA_3_methyladenine_glycosylase_decode','clp_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_13_decode','clp_Aldehyde_dehydrogenase_mitochondrial_decode',
             'clp_Atrial_natriuretic_factor_decode','clp_Fibrinogen_like_protein_1_decode','clp_Growth_differentiation_factor_15_decode',
             'clp_Transcription_factor_jun_D_decode','clp_Glycoprotein_endo_alpha_1_2_mannosidase_decode','clp_Methionine_aminopeptidase_1_decode',
             'clp_Palmitoleoyl_protein_carboxylesterase_NOTUM_decode','clp_Pregnancy_zone_protein_decode','clp_Regulator_of_G_protein_signaling_18_decode','clp_Glia_derived_nexin_decode',
             'clp_ADAMTS13_sunukb','clp_APOBR_sunukb','clp_FES_sunukb','clp_FGF5_sunukb','clp_PZP_sunukb','clp_SERPINE2_sunukb','clp_SH2B3_sunukb','clp_SULT1A1_sunukb')

rm(list=ls()[!(ls() %in% clplist)]) 
clplist <- c('clp_DNA_3_methyladenine_glycosylase_decode','clp_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_13_decode','clp_Aldehyde_dehydrogenase_mitochondrial_decode',
             'clp_Atrial_natriuretic_factor_decode','clp_Fibrinogen_like_protein_1_decode','clp_Growth_differentiation_factor_15_decode',
             'clp_Transcription_factor_jun_D_decode','clp_Glycoprotein_endo_alpha_1_2_mannosidase_decode','clp_Methionine_aminopeptidase_1_decode',
             'clp_Palmitoleoyl_protein_carboxylesterase_NOTUM_decode','clp_Pregnancy_zone_protein_decode','clp_Regulator_of_G_protein_signaling_18_decode','clp_Glia_derived_nexin_decode',
             'clp_ADAMTS13_sunukb','clp_APOBR_sunukb','clp_FES_sunukb','clp_FGF5_sunukb','clp_PZP_sunukb','clp_SERPINE2_sunukb','clp_SH2B3_sunukb','clp_SULT1A1_sunukb')

# Extract SNPs from each 
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- dat$SNP
}
join_list<-lapply(clplist, formfunc)
names(join_list) <- clplist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Phenoscan results
phenoscanfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- phenoscanner(snpquery=dat, pvalue = 5e-8)[['results']]
}
join_list<-lapply(clplist, phenoscanfunc)
names(join_list) <- gsub('clp_', '', clplist)
list_drop_empty(join_list)
join_list <- join_list[sapply(join_list, function(x) dim(x)[1]) > 0]
join_list <- lapply(names(join_list),
                    function(current_name)
                      transform(join_list[[current_name]],
                                Protein = current_name))
res <- Reduce(rbind, join_list)
head(res)

res <- res[,c("Protein", "snp", "hg19_coordinates", "hg38_coordinates", "a1", 
              "a2", "trait", "efo", "study", "pmid", "ancestry", "year", "beta", 
              "se", "p")]
colnames(res) <- c("Protein", "SNP", "hg19 coordinates", "hg38 coordinates", "Allele 1", 
                   "Allele 2", "Phenotype", "EFO", "Study", "PubMed ID", "Ancestry", "Year", "Beta", 
                   "SE", "pvalue")

write.csv(res, '~/desktop/preec/resproteomic/phenoscan.csv')


rm(list=ls())





#### ---------------   Association w Cong malform - SNP not in FinnGen ----------------------------  ####
library(phenoscanner)
# Most complete SNP list (all others are subsets of this ) = on outcome of AF
# Therefore loading and harmonising on AF for phenoscanner

# #### MR OUTCOME = finngen_R9_Q17_VACTREL ####
# ##load exposures -  all CSVs from IV folder 
# # ANP
# anp_iv<- read_exposure_data(
#   filename="/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_unadj_cis_preec/p_Atrial_natriuretic_factor.csv", sep=",", snp_col = "SNP", 
#   beta_col = "beta.exposure", se_col = "se.exposure", ncase_col = "ncase.exposure", 
#   effect_allele_col = "effect_allele.exposure",other_allele_col = "other_allele.exposure",  eaf_col = "eaf.exposure",
#   pval_col = "pval.exposure", phenotype_col= "shortname")
# #  
# instruments <- anp_iv
# #
# #
# #extract outcome for all IVs - R9_Q17_VACTREL
# Q17_VACTREL <- as.data.frame(fread(("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_Q17_VACTREL.gz")))
# setnames(Q17_VACTREL, old=c("rsids", "alt", "beta", "sebeta", "pval", "#chrom", "pos","ref","af_alt"),
#          new=c("SNP", "Q17_VACTREL_ea", "Q17_VACTREL_beta", "Q17_VACTREL_se", "Q17_VACTREL_p", "chr", "pos","Q17_VACTREL_ref","Q17_VACTREL_maf"))
# Q17_VACTREL<-Q17_VACTREL[, c("Q17_VACTREL_ea", "Q17_VACTREL_beta", "Q17_VACTREL_se", "Q17_VACTREL_p","SNP","Q17_VACTREL_ref","Q17_VACTREL_maf")]
# Q17_VACTREL$Q17_VACTREL_ea<-toupper(Q17_VACTREL$Q17_VACTREL_ea)
# colnames(Q17_VACTREL)
# 
# # Split SNPs which contain two rsids. 
# # First we split dataframe in two parts, with SNPs that are duplicate and those that are not
# Q17_VACTRELNoDoubles<-Q17_VACTREL[-which(grepl(",", Q17_VACTREL$SNP)),]
# Q17_VACTRELDoubles<-Q17_VACTREL[which(grepl(",", Q17_VACTREL$SNP)),]
# # Then we duplicate the dataframe with doubles
# Q17_VACTRELDoublespart1<-Q17_VACTRELDoubles
# Q17_VACTRELDoublespart2<-Q17_VACTRELDoubles
# # In the first dataset we create values for the rsid *BEFORE* the comma (,)
# Q17_VACTRELDoublespart1$SNP<-trimws(Q17_VACTRELDoublespart1$SNP, whitespace = ",[^,]*")
# #In the first dataset we create values for the risd *AFTER* to the comma (,)
# Q17_VACTRELDoublespart2$SNP<-trimws(Q17_VACTRELDoublespart1$SNP,which="right", whitespace = ",[^,]*")
# # Then we remerge the data frame with (1) No doubles (2) First rsid (2) Second rsid
# Q17_VACTREL<-rbind(Q17_VACTRELNoDoubles,Q17_VACTRELDoublespart1,Q17_VACTRELDoublespart2)
# rm(Q17_VACTRELDoubles,Q17_VACTRELDoublespart1,Q17_VACTRELDoublespart2, Q17_VACTRELNoDoubles)
# 
# # Merge updated outcome data to instrument list 
# Q17_VACTREL <- Q17_VACTREL[!duplicated(Q17_VACTREL$SNP),]
# instruments_Q17_VACTREL<-dplyr::left_join(instruments,Q17_VACTREL, by=c("SNP"))
# 
# # Drop duplicates, and drop instruments with no estimate for Q17_VACTREL
# data<-instruments_Q17_VACTREL[-which(is.na(instruments_Q17_VACTREL$Q17_VACTREL_beta)),]
# head(data)
# 
# # Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
# exposuredata<-anp_iv
# 
# outcomedata<-data[,c("Q17_VACTREL_beta","Q17_VACTREL_ea", "Q17_VACTREL_se", "Q17_VACTREL_maf","SNP","Q17_VACTREL_ref","pos","chr","Q17_VACTREL_p")]
# outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
#                         beta_col = "Q17_VACTREL_beta",
#                         se_col = "Q17_VACTREL_se",
#                         eaf_col = "Q17_VACTREL_maf",
#                         effect_allele_col = "Q17_VACTREL_ea",
#                         other_allele_col = "Q17_VACTREL_ref",pos_col	="pos",chr_col="chr",pval_col = "Q17_VACTREL_p")
# 
# # Harmonize exposure and outcome data
# harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
# harmonizeddata<-harmonizeddata[which(harmonizeddata$mr_keep==TRUE),]
# rm(exposuredat,exposuredata,outcomedat,outcomedata)
# 
# setnames(harmonizeddata, old=c("chr.outcome"), new=c("chr"))
# 
# # Join harmonised effects to each IV set (by merging harmonised data with [SNP,chr] of each IV list)
# snplist<-harmonizeddata$SNP
# 
# select<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   filter(dat, SNP %in% snplist)
# }
# join_list<-sapply(genelist, select, simplify = FALSE)
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# rm(snplist)
# 
# joining<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   merge(dat[,c("SNP","chr")], harmonizeddata, by=c("SNP","chr"))
# }
# join_list<-sapply(genelist, joining, simplify = FALSE)
# invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
# 
# # # Delete rows with NAs frorm each IV list
# # nas<-function(dat){
# #   dat <- get(dat, envir = .GlobalEnv)
# #   dat<-na.omit(dat)
# # }
# # iv_list<- sapply(genelist, nas, simplify = FALSE)
# # invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
# 
# # Do MR
# mrfunc2<-function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   b<-0
#   se<-0
#   pval<-0
#   n.SNP<-0
#   
#   if(nrow(dat)==1){
#     dat$effect_allele_col<-"A"
#     dat$other_allele_col<-"C"
#     exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
#                              beta_col = "beta.exposure",
#                              se_col = "se.exposure",
#                              pval_col="pval.exposure",
#                              eaf_col = "eaf.exposure",
#                              effect_allele_col = "effect_allele_col",
#                              other_allele_col = "other_allele_col")
#     
#     
#     outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
#                             beta_col = "beta.outcome",
#                             se_col = "se.outcome",
#                             pval_col="pval.outcome",
#                             eaf_col = "eaf.exposure",
#                             effect_allele_col = "effect_allele_col",
#                             other_allele_col = "other_allele_col")
#     
#     harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
#     output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
#     b<-output[1,7]
#     se<-output[1,8]
#     pval<-output[1,9]
#     n.SNP<-1
#     
#   }else{
#     output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
#     b<-unlist(output)[1]
#     se<-unlist(output)[2]
#     pval<-unlist(output)[3]
#     n.SNP<-unlist(output)[4]
#   }
#   
#   outputobj<-list(b,se,pval,n.SNP)
#   return(unlist(outputobj))
# }
# mr_table2 <- as.data.frame(t(data.frame((sapply(genelist,mrfunc2)))))
# mr_table2$outcome <-  "Q17_VACTREL"
# write.csv(mr_table2,"~/Desktop/preec/phenomescan/Q17_VACTREL_ivw.csv")
# rm(list=ls())
# #
# 
# 

#### --------------------------------------------------------------------------------------------- ####
#### -------------------------------------- SENSITIVITY 1 ---------------------------------------- ####
#### -------------------------------------different thresholds------------------------------------ ####
#### ------------------------------ DECODE ------------------------------- ####
#### Format exposure data files - extract cis- (+- 500kb) p<1e-4 ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec")

# 3MG 16:77,007-85,851
gunzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt')
poslow <- 77007-500000
poshigh <- 85851+500000
chrom <- 'chr16'
p_3MG_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_3MG_ivs<-filter(p_3MG_ivs, p_3MG_ivs$Pos>poslow & p_3MG_ivs$Pos<poshigh & p_3MG_ivs$Chrom == chrom)
write.csv(p_3MG_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_3MG_ivs.csv')
gzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz')

# ADAMTS13 9:133414358 - 133459402
gunzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt')
poslow <- 133414358-500000
poshigh <- 133459402+500000
chrom <- 'chr9'
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$Pos>poslow & p_ADAMTS13_ivs$Pos<poshigh & p_ADAMTS13_ivs$Chrom == chrom)
write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_ADAMTS13_ivs.csv')
gzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz')


# ALDH2 12:111766887 - 111817532
gunzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt')
poslow <- 111766887-500000
poshigh <- 111817532+500000
chrom <- 'chr12'
p_ALDHE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ALDHE2_ivs<-filter(p_ALDHE2_ivs, p_ALDHE2_ivs$Pos>poslow & p_ALDHE2_ivs$Pos<poshigh & p_ALDHE2_ivs$Chrom == chrom)
write.csv(p_ALDHE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_ALDHE2_ivs.csv')
gzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz')

# ANP 1:11845709 - 11848345
gunzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt')
poslow <- 11845709-500000
poshigh <- 11848345+500000
chrom <- 'chr1'
p_ANP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ANP_ivs<-filter(p_ANP_ivs, p_ANP_ivs$Pos>poslow & p_ANP_ivs$Pos<poshigh & p_ANP_ivs$Chrom == chrom)
write.csv(p_ANP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_ANP_ivs.csv')
gzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz')

# FGL1 8:17864380 - 17910365
gunzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt')
poslow <- 17864380-500000
poshigh <- 17910365+500000
chrom <- 'chr8'
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$Pos>poslow & p_FGL1_ivs$Pos<poshigh & p_FGL1_ivs$Chrom == chrom)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_FGL1_ivs.csv')
gzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz')

# GDF15 19:18374731 - 18389176
gunzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt')
poslow <- 18374731-500000
poshigh <- 18389176+500000
chrom <- 'chr19'
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$Pos>poslow & p_GDF15_ivs$Pos<poshigh & p_GDF15_ivs$Chrom == chrom)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_GDF15_ivs.csv')
gzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz')

# JUND 19:18279694 - 18281622
gunzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt')
poslow <- 18279694-500000
poshigh <- 18281622+500000
chrom <- 'chr19'
p_JUND_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_JUND_ivs<-filter(p_JUND_ivs, p_JUND_ivs$Pos>poslow & p_JUND_ivs$Pos<poshigh & p_JUND_ivs$Chrom == chrom)
write.csv(p_JUND_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_JUND_ivs.csv')
gzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz')

# MANEA 6:95577485 - 95609470
gunzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt')
poslow <-  95577485-500000
poshigh <- 95609470+500000
chrom <- 'chr6'
p_MANEA_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_MANEA_ivs<-filter(p_MANEA_ivs, p_MANEA_ivs$Pos>poslow & p_MANEA_ivs$Pos<poshigh & p_MANEA_ivs$Chrom == chrom)
write.csv(p_MANEA_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_MANEA_ivs.csv')
gzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz')

# METAP1 4:98995659 - 99062809 
gunzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt')
poslow <- 98995659-500000
poshigh <- 99062809+500000
chrom <- 'chr4'
p_METAP1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_METAP1_ivs<-filter(p_METAP1_ivs, p_METAP1_ivs$Pos>poslow & p_METAP1_ivs$Pos<poshigh & p_METAP1_ivs$Chrom == chrom)
write.csv(p_METAP1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_METAP1_ivs.csv')
gzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz')

# NOTUM 17:81952507 - 81961840
gunzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt')
poslow <- 81952507-500000
poshigh <- 81961840+500000
chrom <- 'chr17'
f <- as.data.frame(fread('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt'))
p_NOTUM_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', separator = "\t")$
  filter(pl$col("Pval") < 5e-8)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_NOTUM_ivs<-filter(p_NOTUM_ivs, p_NOTUM_ivs$Pos>poslow & p_NOTUM_ivs$Pos<poshigh & p_NOTUM_ivs$Chrom == chrom)
write.csv(p_NOTUM_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_NOTUM_ivs.csv')
gzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz')

# PZP 12:9148840 - 9208395
gunzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt')
poslow <- 9148840-500000
poshigh <- 9208395+500000
chrom <- 'chr12'
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$Pos>poslow & p_PZP_ivs$Pos<poshigh & p_PZP_ivs$Chrom == chrom)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_PZP_ivs.csv')
gzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz')

# RGS18 1:192158462 - 192185815 
gunzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt')
poslow <- 192158462-500000
poshigh <- 192185815+500000
chrom <- 'chr1'
p_RGS18_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_RGS18_ivs<-filter(p_RGS18_ivs, p_RGS18_ivs$Pos>poslow & p_RGS18_ivs$Pos<poshigh & p_RGS18_ivs$Chrom == chrom)
write.csv(p_RGS18_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_RGS18_ivs.csv')
gzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz')

# SERPINE2 2:223975045 - 224039318
gunzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt')
poslow <- 223975045-500000
poshigh <- 224039318+500000
chrom <- 'chr2'
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', separator = "\t")$
  filter(pl$col("Pval") < 1e-4)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$Pos>poslow & p_SERPINE2_ivs$Pos<poshigh & p_SERPINE2_ivs$Chrom == chrom)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec/p_SERPINE2_ivs.csv')
gzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz')
rm(list=ls())

#### Format outcome data - select SNP 1e-4 ####
dir.create('~/desktop/preec/outdecode_1e4')
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
pqtls <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
pqtls$chr <- gsub('chr', '', pqtls$Chrom)
pqtls$chrpos <- str_c(pqtls$chr, ':', pqtls$Pos)
pqtls <- pqtls[,c('rsids', 'chrpos')]
colnames(pqtls) <- c('SNP', 'chrpos')
pqtls$SNP <- replace_na(pqtls$SNP, pqtls$chrpos)
pqtls$SNP <- ifelse(!is.na(pqtls$SNP), pqtls$SNP, pqtls$chrpos)


preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt', drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, '~/desktop/preec/outdecode_1e4/preec_out_unadj.csv')
rm(list=ls())

#### MR - cis - 5e6 r2 0.2 ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_1e4_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rsids",
                    beta_col = "Beta",
                    se_col = "SE",
                    pval_col="minus_log10_pval",
                    eaf_col = "ImpMAF",
                    effect_allele_col = "effectAllele",
                    other_allele_col = "otherAllele", 
                    phenotype_col =  "phenotype", log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Filter 5e-6
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- filter(dat, dat$pval.exposure < 5e-6)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outdecode_1e4/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
genelist<-ls(pattern = "har_", mget(ls()))
try_clp <- function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- clump_data(dat, clump_r2 = 0.2)
}
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- as.data.frame(mr.raps.overdispersed.robust(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,"tukey"))
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_raps')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
reslist<-ls(pattern = "_raps", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

# Other sensitivity analyses
# # re-clump 0.001 locally 
# genelist<-ls(pattern = "clp_", mget(ls()))
# try_clp <- function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <- clump_data(dat, clump_r2 = 0.001) 
# }
# iv_list <- sapply(genelist, try_clp, simplify = FALSE)
# names(iv_list) <- clplist 
# invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_sens')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
reslist<-ls(pattern = "_sens", mget(ls()))

setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}


# MR intercept 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr_pleiotropy_test(dat)
}

mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_intercept')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_intercept", mget(ls()))
setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

# save raps
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec', 
                    pattern = "_raps.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))

mergedres$or <- exp(mergedres$beta.hat)
mergedres$lci <- exp(mergedres$beta.hat - 1.96*mergedres$beta.se)
mergedres$uci <- exp(mergedres$beta.hat + 1.96*mergedres$beta.se)
mergedres <- mergedres[,c("beta.hat",  "beta.se",  "beta.p.value", 
                          "filename", "or", "lci", "uci")]
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec/', '', mergedres$filename))
mergedres$filename <- NULL
colnames(mergedres) <- c("beta", "beta", "pval", "or", "lci", "uci", 
                         "outcome", "exposure")
mergedres$method <- 'p<5e-6, r2<0.2, MR-RAPS'
write.csv(mergedres, "~/Desktop/preec/resdecode_preec/raps_cis_preec_5e-6.csv", row.names = FALSE)

rm(list=ls())

# save other sensitivity analyses
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec', 
                    pattern = "_sens.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL


files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec', 
                    pattern = "_intercept.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
intercepts <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(intercepts) <- NULL
setnames(intercepts, old = 'egger_intercept', new = 'b')
intercepts$nsnp <- NA
intercepts$method <- 'MR-Egger intercept'
intercepts <- intercepts[,colnames(mergedres)]

mergedres <- rbind(mergedres, intercepts)
mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/resdecode_5e-6_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres <- filter(mergedres, mergedres$method == 'MR Egger' | mergedres$method == 'Weighted median'| mergedres$method == 'MR-Egger intercept')
write.csv(mergedres, "~/Desktop/preec/resdecode_preec/sens_cis_preec_5e6.csv", row.names = FALSE)

mergedres
rm(list=ls())

#### -------------------------------- UKB ------------------------------- ####
#### Format exposure data files - extract cis- (+- 500kb) p<1e-4 ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec")
dir.create('/volumes/maddy2/sunukb/temptar')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ADAMTS13 9:133414358 - 133459402
untar('/volumes/maddy2/sunukb/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt')
poslow <- 133414358-500000
poshigh <- 133459402+500000
chrom <- 9
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(5e-6)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$GENPOS>poslow & p_ADAMTS13_ivs$GENPOS<poshigh & p_ADAMTS13_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_ADAMTS13_ivs <- cbind(p_ADAMTS13_ivs, data.frame(do.call('rbind', strsplit(as.character(p_ADAMTS13_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_ADAMTS13_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '9')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_ADAMTS13_ivs <- merge(p_ADAMTS13_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)

write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_ADAMTS13_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# APOBR 16:28494643 - 28498964
untar('/volumes/maddy2/sunukb/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt')
poslow <- 28494643-500000
poshigh <- 28498964+500000
chrom <- 16
p_APOBR_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_APOBR_ivs<-filter(p_APOBR_ivs, p_APOBR_ivs$GENPOS>poslow & p_APOBR_ivs$GENPOS<poshigh & p_APOBR_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_APOBR_ivs <- cbind(p_APOBR_ivs, data.frame(do.call('rbind', strsplit(as.character(p_APOBR_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_APOBR_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_APOBR_ivs <- merge(p_APOBR_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)

write.csv(p_APOBR_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_APOBR_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# IL27 16:28499362 - 28512051 
untar('/volumes/maddy2/sunukb/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EB13_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.txt')
poslow <- 28499362-500000
poshigh <- 28512051+500000
chrom <- 16
p_EB13_IL27_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EB13_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_EB13_IL27_ivs<-filter(p_EB13_IL27_ivs, p_EB13_IL27_ivs$GENPOS>poslow & p_EB13_IL27_ivs$GENPOS<poshigh & p_EB13_IL27_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_EB13_IL27_ivs <- cbind(p_EB13_IL27_ivs, data.frame(do.call('rbind', strsplit(as.character(p_EB13_IL27_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_EB13_IL27_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_EB13_IL27_ivs <- merge(p_EB13_IL27_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_EB13_IL27_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_EB13_IL27_ivs.csv')

unlink('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology', recursive = TRUE)
rm(list=ls())

# FES 15:90883695 - 90895776
untar('/volumes/maddy2/sunukb/FES_P07332_OID21207_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt')
poslow <- 90883695-500000
poshigh <- 90895776+500000
chrom <- 15
p_FES_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FES_ivs<-filter(p_FES_ivs, p_FES_ivs$GENPOS>poslow & p_FES_ivs$GENPOS<poshigh & p_FES_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FES_ivs <- cbind(p_FES_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FES_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FES_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '15')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FES_ivs <- merge(p_FES_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FES_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_FES_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology', recursive = TRUE)
rm(list=ls())

# FGF5 4:80266639 - 80336680
untar('/volumes/maddy2/sunukb/FGF5_P12034_OID20490_v1_Inflammation.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.gz', 
       '/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt')
poslow <- 80266639-500000
poshigh <- 80336680+500000
chrom <- 4
p_FGF5_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGF5_ivs<-filter(p_FGF5_ivs, p_FGF5_ivs$GENPOS>poslow & p_FGF5_ivs$GENPOS<poshigh & p_FGF5_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGF5_ivs <- cbind(p_FGF5_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGF5_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGF5_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '4')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGF5_ivs <- merge(p_FGF5_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGF5_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_FGF5_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation', recursive = TRUE)
rm(list=ls())

# FGL1 8:17864380 - 17910365
untar('/volumes/maddy2/sunukb/FGL1_Q08830_OID30702_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt')
poslow <- 17864380-500000
poshigh <- 17910365+500000
chrom <- 8
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$GENPOS>poslow & p_FGL1_ivs$GENPOS<poshigh & p_FGL1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGL1_ivs <- cbind(p_FGL1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGL1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGL1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '8')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGL1_ivs <- merge(p_FGL1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_FGL1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())
x
# GDF15 19:18374731 - 18389176
untar('/volumes/maddy2/sunukb/GDF15_Q99988_OID20251_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt')
poslow <- 18374731-500000
poshigh <- 18389176+500000
chrom <- 19
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$GENPOS>poslow & p_GDF15_ivs$GENPOS<poshigh & p_GDF15_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_GDF15_ivs <- cbind(p_GDF15_ivs, data.frame(do.call('rbind', strsplit(as.character(p_GDF15_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_GDF15_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '19')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_GDF15_ivs <- merge(p_GDF15_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_GDF15_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# PZP 12:9148840 - 9208395
untar('/volumes/maddy2/sunukb/PZP_P20742_OID30730_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt')
poslow <- 9148840-500000
poshigh <- 9208395+500000
chrom <- 12
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$GENPOS>poslow & p_PZP_ivs$GENPOS<poshigh & p_PZP_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_PZP_ivs <- cbind(p_PZP_ivs, data.frame(do.call('rbind', strsplit(as.character(p_PZP_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_PZP_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_PZP_ivs <- merge(p_PZP_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_PZP_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())

# SERPINE2 2:223975045 - 224039318
untar('/volumes/maddy2/sunukb/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt')
poslow <- 223975045-500000
poshigh <- 224039318+500000
chrom <- 2
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$GENPOS>poslow & p_SERPINE2_ivs$GENPOS<poshigh & p_SERPINE2_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SERPINE2_ivs <- cbind(p_SERPINE2_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SERPINE2_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SERPINE2_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '2')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SERPINE2_ivs <- merge(p_SERPINE2_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_SERPINE2_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# SH2B3 12:111405923 - 111451623
untar('/volumes/maddy2/sunukb/SH2B3_Q9UQQ2_OID21222_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt')
poslow <- 111405923-500000
poshigh <- 111451623+500000
chrom <- 12
p_SH2B3_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SH2B3_ivs<-filter(p_SH2B3_ivs, p_SH2B3_ivs$GENPOS>poslow & p_SH2B3_ivs$GENPOS<poshigh & p_SH2B3_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SH2B3_ivs <- cbind(p_SH2B3_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SH2B3_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SH2B3_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SH2B3_ivs <- merge(p_SH2B3_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SH2B3_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_SH2B3_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology', recursive = TRUE)
rm(list=ls())

# SULT1A1 16:28605196 - 28614279
untar('/volumes/maddy2/sunukb/SULT1A1_P50225_OID21031_v1_Neurology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.gz', 
       '/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt')
poslow <- 28605196-500000
poshigh <- 28614279+500000
chrom <- 16
p_SULT1A1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 5.30103)$ # -log10(1e-4)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SULT1A1_ivs<-filter(p_SULT1A1_ivs, p_SULT1A1_ivs$GENPOS>poslow & p_SULT1A1_ivs$GENPOS<poshigh & p_SULT1A1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SULT1A1_ivs <- cbind(p_SULT1A1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SULT1A1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SULT1A1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SULT1A1_ivs <- merge(p_SULT1A1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SULT1A1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec/p_SULT1A1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology', recursive = TRUE)
rm(list=ls())



#### Format outcome data - select SNP 1e-4 ####
dir.create('~/desktop/preec/outsunukb_1e4')
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
pqtls <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
pqtls$chrpos <- str_c(pqtls$CHR, ':', pqtls$GENPOS)
pqtls <- pqtls[,c('SNP', 'chrpos')]
rownames(pqtls) <- NULL

preec <- as.data.frame(fread('/Volumes/MADDY2/datasets/preeclampsia/metal_preec_European_allBiobanks_omitNone_1.txt', drop = c('FreqSE', 'MinFreq', 'MaxFreq', 'Direction')))
preec$chrpos <- str_c(preec$Chromosome, ':', preec$Position)
preec_out <- merge(pqtls, preec, by='chrpos', all.x=FALSE, all.y = FALSE)
preec_out <- preec_out[,c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'chrpos')]
preec_out$Allele1 <- toupper(preec_out$Allele1)
preec_out$Allele2 <- toupper(preec_out$Allele2)
setnames(preec_out, old=c('SNP', 'Chromosome', 'Position', 'Allele1', 'Allele2', 'Freq1', 'Effect','StdErr', 'P-value', 'chrpos'), 
         new = c('SNP', 'chr.outcome', 'pos.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'beta.outcome','se.outcome', 'pval.outcome', 'chrpos'))
preec_out$phenotype <- 'preec'
preec_out <- preec_out[!duplicated(preec_out$SNP),]
write.csv(preec_out, '~/desktop/preec/outsunukb_1e4/preec_out_unadj.csv')
rm(list=ls())


#### MR - cis - 5e6 r2 0.2 ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_1e4_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
                    beta_col = "BETA",
                    se_col = "SE",
                    pval_col="LOG10P",
                    eaf_col = "A1FREQ",
                    effect_allele_col = "ALLELE1",
                    other_allele_col = "ALLELE0", 
                    phenotype_col =  "new_column",  # CHECK THIS IS THE EXPOSURE NAME
                    log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Filter 5e-6
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- filter(dat, dat$pval.exposure < 5e-6)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outsunukb_1e4/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
genelist<-ls(pattern = "har_", mget(ls()))
try_clp <- function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- clump_data(dat, clump_r2 = 0.2)
}
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- as.data.frame(mr.raps.overdispersed.robust(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,"tukey"))
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_raps')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
reslist<-ls(pattern = "_raps", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

# Other sensitivity analyses
# re-clump 0.001 locally 
# genelist<-ls(pattern = "clp_", mget(ls()))
# try_clp <- function(dat) {
#   dat <- get(dat, envir = .GlobalEnv)
#   dat <- clump_data(dat, clump_r2 = 0.001) 
# }
# iv_list <- sapply(genelist, try_clp, simplify = FALSE)
# names(iv_list) <- clplist 
# invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
# to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
# rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
# rm(to.rm)
# rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_sens')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)
reslist<-ls(pattern = "_sens", mget(ls()))

setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}


# MR intercept 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr_pleiotropy_test(dat)
}

mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_intercept')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_intercept", mget(ls()))
setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

# save raps
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec', 
                    pattern = "_raps.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))

mergedres$or <- exp(mergedres$beta.hat)
mergedres$lci <- exp(mergedres$beta.hat - 1.96*mergedres$beta.se)
mergedres$uci <- exp(mergedres$beta.hat + 1.96*mergedres$beta.se)
mergedres <- mergedres[,c("beta.hat",  "beta.se",  "beta.p.value", 
                          "filename", "or", "lci", "uci")]
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec/', '', mergedres$filename))
mergedres$filename <- NULL
colnames(mergedres) <- c("beta", "beta", "pval", "or", "lci", "uci", 
                         "outcome", "exposure")
mergedres$method <- 'p<5e-6, r2<0.2, MR-RAPS'
write.csv(mergedres, "~/Desktop/preec/ressunukb_preec/raps_cis_preec_5e-6.csv", row.names = FALSE)

rm(list=ls())

# save other sensitivity analyses
files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec', 
                    pattern = "_sens.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(mergedres) <- NULL


files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec', 
                    pattern = "_intercept.csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
intercepts <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))
rownames(intercepts) <- NULL
setnames(intercepts, old = 'egger_intercept', new = 'b')
intercepts$nsnp <- NA
intercepts$method <- 'MR-Egger intercept'
intercepts <- intercepts[,colnames(mergedres)]

mergedres <- rbind(mergedres, intercepts)
mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_5e-6_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres <- filter(mergedres, mergedres$method == 'MR Egger' | mergedres$method == 'Weighted median'| mergedres$method == 'MR-Egger intercept')
write.csv(mergedres, "~/Desktop/preec/ressunukb_preec/sens_cis_preec_5e6.csv", row.names = FALSE)

mergedres
rm(list=ls())


#### --------------------------------------------------------------------------------------------- ####
#### -------------------------------------- SENSITIVITY 2 ---------------------------------------- ####
#### -----------------------------------outside gene region--------------------------------------- ####
#### ------------------------------ DECODE ------------------------------- ####
#### Format exposure data files - extract cis- (+- 500kb) outside ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec")

# 3MG 16:77,007-85,851
gunzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt')
poslow <- 77007
poshigh <- 85851
chrom <- 'chr16'
p_3MG_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_3MG_ivs<-filter(p_3MG_ivs, p_3MG_ivs$Pos<poslow & p_3MG_ivs$Chrom == chrom | p_3MG_ivs$Pos>poshigh & p_3MG_ivs$Chrom == chrom)
write.csv(p_3MG_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_3MG_ivs.csv')
gzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz')

# ADAMTS13 9:133414358 - 133459402
gunzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt')
poslow <- 133414358
poshigh <- 133459402
chrom <- 'chr9'
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$Pos<poslow & p_ADAMTS13_ivs$Chrom == chrom| p_ADAMTS13_ivs$Pos>poshigh & p_ADAMTS13_ivs$Chrom == chrom)
write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_ADAMTS13_ivs.csv')
gzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz')


# ALDH2 12:111766887 - 111817532
gunzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt')
poslow <- 111766887
poshigh <- 111817532
chrom <- 'chr12'
p_ALDHE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ALDHE2_ivs<-filter(p_ALDHE2_ivs, p_ALDHE2_ivs$Pos<poslow & p_ALDHE2_ivs$Chrom == chrom | p_ALDHE2_ivs$Pos>poshigh & p_ALDHE2_ivs$Chrom == chrom)
write.csv(p_ALDHE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_ALDHE2_ivs.csv')
gzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz')

# ANP 1:11845709 - 11848345
gunzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt')
poslow <- 11845709
poshigh <- 11848345
chrom <- 'chr1'
p_ANP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ANP_ivs<-filter(p_ANP_ivs, p_ANP_ivs$Pos<poslow & p_ANP_ivs$Chrom == chrom | p_ANP_ivs$Pos>poshigh & p_ANP_ivs$Chrom == chrom)
write.csv(p_ANP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_ANP_ivs.csv')
gzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz')

# FGL1 8:17864380 - 17910365
gunzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt')
poslow <- 17864380
poshigh <- 17910365
chrom <- 'chr8'
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$Pos<poslow & p_FGL1_ivs$Chrom == chrom| p_FGL1_ivs$Pos>poshigh & p_FGL1_ivs$Chrom == chrom)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_FGL1_ivs.csv')
gzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz')

# GDF15 19:18374731 - 18389176
gunzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt')
poslow <- 18374731
poshigh <- 18389176
chrom <- 'chr19'
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$Pos<poslow & p_GDF15_ivs$Chrom == chrom| p_GDF15_ivs$Pos>poshigh & p_GDF15_ivs$Chrom == chrom)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_GDF15_ivs.csv')
gzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz')

# JUND 19:18279694 - 18281622
gunzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt')
poslow <- 18279694
poshigh <- 18281622
chrom <- 'chr19'
p_JUND_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_JUND_ivs<-filter(p_JUND_ivs, p_JUND_ivs$Pos<poslow & p_JUND_ivs$Chrom == chrom| p_JUND_ivs$Pos>poshigh & p_JUND_ivs$Chrom == chrom)
write.csv(p_JUND_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_JUND_ivs.csv')
gzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz')

# MANEA 6:95577485 - 95609470
gunzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt')
poslow <-  95577485
poshigh <- 95609470
chrom <- 'chr6'
p_MANEA_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_MANEA_ivs<-filter(p_MANEA_ivs, p_MANEA_ivs$Pos<poslow & p_MANEA_ivs$Chrom == chrom| p_MANEA_ivs$Pos>poshigh & p_MANEA_ivs$Chrom == chrom)
write.csv(p_MANEA_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_MANEA_ivs.csv')
gzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz')

# METAP1 4:98995659 - 99062809 
gunzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt')
poslow <- 98995659
poshigh <- 99062809
chrom <- 'chr4'
p_METAP1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_METAP1_ivs<-filter(p_METAP1_ivs, p_METAP1_ivs$Pos<poslow & p_METAP1_ivs$Chrom == chrom| p_METAP1_ivs$Pos>poshigh & p_METAP1_ivs$Chrom == chrom)
write.csv(p_METAP1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_METAP1_ivs.csv')
gzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz')

# NOTUM 17:81952507 - 81961840
gunzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt')
poslow <- 81952507
poshigh <- 81961840
chrom <- 'chr17'
p_NOTUM_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_NOTUM_ivs<-filter(p_NOTUM_ivs, p_NOTUM_ivs$Pos<poslow & p_NOTUM_ivs$Chrom == chrom| p_NOTUM_ivs$Pos>poshigh & p_NOTUM_ivs$Chrom == chrom)
write.csv(p_NOTUM_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_NOTUM_ivs.csv')
gzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz')

# PZP 12:9148840 - 9208395
gunzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt')
poslow <- 9148840
poshigh <- 9208395
chrom <- 'chr12'
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$Pos<poslow & p_PZP_ivs$Chrom == chrom| p_PZP_ivs$Pos>poshigh & p_PZP_ivs$Chrom == chrom)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_PZP_ivs.csv')
gzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz')

# RGS18 1:192158462 - 192185815 
gunzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt')
poslow <- 192158462
poshigh <- 192185815
chrom <- 'chr1'
p_RGS18_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_RGS18_ivs<-filter(p_RGS18_ivs, p_RGS18_ivs$Pos<poslow  & p_RGS18_ivs$Chrom == chrom| p_RGS18_ivs$Pos>poshigh & p_RGS18_ivs$Chrom == chrom)
write.csv(p_RGS18_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_RGS18_ivs.csv')
gzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz')

# SERPINE2 2:223975045 - 224039318
gunzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt')
poslow <- 223975045
poshigh <- 224039318
chrom <- 'chr2'
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$Pos<poslow & p_SERPINE2_ivs$Chrom == chrom| p_SERPINE2_ivs$Pos>poshigh & p_SERPINE2_ivs$Chrom == chrom)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec/p_SERPINE2_ivs.csv')
gzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz')
rm(list=ls())

beep(2)

#### MR - cis - outside ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_outside_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rsids",
                    beta_col = "Beta",
                    se_col = "SE",
                    pval_col="minus_log10_pval",
                    eaf_col = "ImpMAF",
                    effect_allele_col = "effectAllele",
                    other_allele_col = "otherAllele", 
                    phenotype_col =  "phenotype", log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Select only exposures with SNPs in outcome 
out <- as.data.frame(fread("~/Desktop/preec/outdecode/preec_out_unadj.csv")[,-1])
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- dat[which(dat$SNP %in% out$SNP),]
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% genelist)]) # keep only clumped files
genelist<-ls(pattern = "_ivs", mget(ls())) # Make list with names of all cluped sets


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outdecode/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/resdecode_outside_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/resdecode_outside_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_outside_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_outside_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/resdecode_outside_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres$analysis <- 'Outside cis-region'
mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted' | mergedres$method == 'Wald ratio')
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', length(mergedres$pval))
write.csv(mergedres, "~/Desktop/preec/resdecode_preec/sens_cis_preec_outside.csv", row.names = FALSE)

rm(list=ls())



#### -------------------------------- UKB ------------------------------- ####
#### Format exposure data files - extract cis- (+- 500kb) outside ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec")
dir.create('/volumes/maddy2/sunukb/temptar')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ADAMTS13 9:133414358 - 133459402
untar('/volumes/maddy2/sunukb/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt')
poslow <- 133414358
poshigh <- 133459402
chrom <- 9
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(P<1.7×10−11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$GENPOS<poslow | p_ADAMTS13_ivs$GENPOS>poshigh & p_ADAMTS13_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_ADAMTS13_ivs <- cbind(p_ADAMTS13_ivs, data.frame(do.call('rbind', strsplit(as.character(p_ADAMTS13_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_ADAMTS13_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '9')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_ADAMTS13_ivs <- merge(p_ADAMTS13_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)

write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_ADAMTS13_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# APOBR 16:28494643 - 28498964
untar('/volumes/maddy2/sunukb/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt')
poslow <- 28494643
poshigh <- 28498964
chrom <- 16
p_APOBR_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_APOBR_ivs<-filter(p_APOBR_ivs, p_APOBR_ivs$GENPOS<poslow | p_APOBR_ivs$GENPOS>poshigh & p_APOBR_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_APOBR_ivs <- cbind(p_APOBR_ivs, data.frame(do.call('rbind', strsplit(as.character(p_APOBR_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_APOBR_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_APOBR_ivs <- merge(p_APOBR_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)

write.csv(p_APOBR_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_APOBR_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# IL27 16:28499362 - 28512051 
untar('/volumes/maddy2/sunukb/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EB13_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.txt')
poslow <- 28499362
poshigh <- 28512051
chrom <- 16
p_EB13_IL27_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EB13_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_EB13_IL27_ivs<-filter(p_EB13_IL27_ivs, p_EB13_IL27_ivs$GENPOS<poslow | p_EB13_IL27_ivs$GENPOS>poshigh & p_EB13_IL27_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_EB13_IL27_ivs <- cbind(p_EB13_IL27_ivs, data.frame(do.call('rbind', strsplit(as.character(p_EB13_IL27_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_EB13_IL27_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_EB13_IL27_ivs <- merge(p_EB13_IL27_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_EB13_IL27_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_EB13_IL27_ivs.csv')

unlink('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology', recursive = TRUE)
rm(list=ls())

# FES 15:90883695 - 90895776
untar('/volumes/maddy2/sunukb/FES_P07332_OID21207_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt')
poslow <- 90883695
poshigh <- 90895776
chrom <- 15
p_FES_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FES_ivs<-filter(p_FES_ivs, p_FES_ivs$GENPOS<poslow | p_FES_ivs$GENPOS>poshigh & p_FES_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FES_ivs <- cbind(p_FES_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FES_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FES_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '15')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FES_ivs <- merge(p_FES_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FES_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_FES_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology', recursive = TRUE)
rm(list=ls())

# FGF5 4:80266639 - 80336680
untar('/volumes/maddy2/sunukb/FGF5_P12034_OID20490_v1_Inflammation.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.gz', 
       '/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt')
poslow <- 80266639
poshigh <- 80336680
chrom <- 4
p_FGF5_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGF5_ivs<-filter(p_FGF5_ivs, p_FGF5_ivs$GENPOS<poslow | p_FGF5_ivs$GENPOS>poshigh & p_FGF5_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGF5_ivs <- cbind(p_FGF5_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGF5_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGF5_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '4')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGF5_ivs <- merge(p_FGF5_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGF5_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_FGF5_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation', recursive = TRUE)
rm(list=ls())

# FGL1 8:17864380 - 17910365
untar('/volumes/maddy2/sunukb/FGL1_Q08830_OID30702_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt')
poslow <- 17864380
poshigh <- 17910365
chrom <- 8
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$GENPOS<poslow | p_FGL1_ivs$GENPOS>poshigh & p_FGL1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGL1_ivs <- cbind(p_FGL1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGL1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGL1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '8')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGL1_ivs <- merge(p_FGL1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_FGL1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())
x
# GDF15 19:18374731 - 18389176
untar('/volumes/maddy2/sunukb/GDF15_Q99988_OID20251_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt')
poslow <- 18374731
poshigh <- 18389176
chrom <- 19
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$GENPOS<poslow | p_GDF15_ivs$GENPOS>poshigh & p_GDF15_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_GDF15_ivs <- cbind(p_GDF15_ivs, data.frame(do.call('rbind', strsplit(as.character(p_GDF15_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_GDF15_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '19')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_GDF15_ivs <- merge(p_GDF15_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_GDF15_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# PZP 12:9148840 - 9208395
untar('/volumes/maddy2/sunukb/PZP_P20742_OID30730_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt')
poslow <- 9148840
poshigh <- 9208395
chrom <- 12
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$GENPOS<poslow | p_PZP_ivs$GENPOS>poshigh & p_PZP_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_PZP_ivs <- cbind(p_PZP_ivs, data.frame(do.call('rbind', strsplit(as.character(p_PZP_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_PZP_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_PZP_ivs <- merge(p_PZP_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_PZP_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())

# SERPINE2 2:223975045 - 224039318
untar('/volumes/maddy2/sunukb/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt')
poslow <- 223975045
poshigh <- 224039318
chrom <- 2
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$GENPOS<poslow | p_SERPINE2_ivs$GENPOS>poshigh & p_SERPINE2_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SERPINE2_ivs <- cbind(p_SERPINE2_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SERPINE2_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SERPINE2_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '2')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SERPINE2_ivs <- merge(p_SERPINE2_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_SERPINE2_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# SH2B3 12:111405923 - 111451623
untar('/volumes/maddy2/sunukb/SH2B3_Q9UQQ2_OID21222_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt')
poslow <- 111405923
poshigh <- 111451623
chrom <- 12
p_SH2B3_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SH2B3_ivs<-filter(p_SH2B3_ivs, p_SH2B3_ivs$GENPOS<poslow | p_SH2B3_ivs$GENPOS>poshigh & p_SH2B3_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SH2B3_ivs <- cbind(p_SH2B3_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SH2B3_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SH2B3_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SH2B3_ivs <- merge(p_SH2B3_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SH2B3_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_SH2B3_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology', recursive = TRUE)
rm(list=ls())

# SULT1A1 16:28605196 - 28614279
untar('/volumes/maddy2/sunukb/SULT1A1_P50225_OID21031_v1_Neurology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.gz', 
       '/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt')
poslow <- 28605196
poshigh <- 28614279
chrom <- 16
p_SULT1A1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.76955)$ # -log10(outside)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SULT1A1_ivs<-filter(p_SULT1A1_ivs, p_SULT1A1_ivs$GENPOS<poslow | p_SULT1A1_ivs$GENPOS>poshigh & p_SULT1A1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SULT1A1_ivs <- cbind(p_SULT1A1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SULT1A1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SULT1A1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SULT1A1_ivs <- merge(p_SULT1A1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SULT1A1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec/p_SULT1A1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology', recursive = TRUE)
rm(list=ls())




#### MR - cis - outside ####
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_outside_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
                    beta_col = "BETA",
                    se_col = "SE",
                    pval_col="LOG10P",
                    eaf_col = "A1FREQ",
                    effect_allele_col = "ALLELE1",
                    other_allele_col = "ALLELE0", 
                    phenotype_col =  "new_column",  # CHECK THIS IS THE EXPOSURE NAME
                    log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)


# Select only exposures with SNPs in outcome 
out <- as.data.frame(fread("~/Desktop/preec/outdecode/preec_out_unadj.csv")[,-1])
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- dat[which(dat$SNP %in% out$SNP),]
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% genelist)]) # keep only clumped files
genelist<-ls(pattern = "_ivs", mget(ls())) # Make list with names of all  sets


# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- data.frame()
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outsunukb/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)


# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% genelist)]) # keep only clumped files
genelist<-ls(pattern = "_ivs", mget(ls())) # Make list with names of all  sets

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all clumped sets

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub('_ivs', '', gsub('p_', '', gsub("clp_","",clplist)))
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/ressunukb_outside_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_outside_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_outside_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_outside_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))


mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)
mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_outside_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres$analysis <- 'Outside cis-region'
mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted' | mergedres$method == 'Wald ratio')
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', length(mergedres$pval))
write.csv(mergedres, "~/Desktop/preec/ressunukb_preec/sens_cis_preec_outside.csv", row.names = FALSE)

rm(list=ls())






#### --------------------------------------------------------------------------------------------- ####
#### -------------------------------------- REVERSE MR ------------------------------------------- ####
#### --------------------------------------------------------------------------------------------- ####
#### decode ####
# Import all IVs
preec_ivs <- read.csv('~/desktop/preec/preec_ivs/poret_iv.csv')
genelist <- 'preec_ivs'
preec_ivs$phenotype <- 'Pre-eclampsia'

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- format_data(dat,type="exposure",snp_col = "RefSNP_id",
                     beta_col = "poret_beta",
                     se_col = "poret_se",
                     pval_col="poret_p",
                     eaf_col = "poret_eaf",
                     effect_allele_col = "poret_ea",
                     other_allele_col = "poret_nea", 
                     phenotype_col =  "phenotype")
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates

#3MG
out_3MG <- read_outcome_data(snps = preec_ivs$SNP, 
                             filename = "/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz", 
                             sep = "\t", 
                             snp_col = "rsids",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col="minus_log10_pval",
                             eaf_col = "ImpMAF",
                             effect_allele_col = "effectAllele",
                             other_allele_col = "otherAllele", log_pval=TRUE)
out_3MG$phenotype <- '3MG'

# PZP
out_PZP <- read_outcome_data(snps = preec_ivs$SNP, 
                             filename = "/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz", 
                             sep = "\t", 
                             snp_col = "rsids",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col="minus_log10_pval",
                             eaf_col = "ImpMAF",
                             effect_allele_col = "effectAllele",
                             other_allele_col = "otherAllele", log_pval=TRUE)
out_PZP$phenotype <- 'PZP'

#ANP
out_ANP <- read_outcome_data(snps = preec_ivs$SNP, 
                             filename = "/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz", 
                             sep = "\t", 
                             snp_col = "rsids",
                             beta_col = "Beta",
                             se_col = "SE",
                             pval_col="minus_log10_pval",
                             eaf_col = "ImpMAF",
                             effect_allele_col = "effectAllele",
                             other_allele_col = "otherAllele", log_pval=TRUE)
out_ANP$phenotype <- 'ANP'


#METAP1
out_METAP1 <- read_outcome_data(snps = preec_ivs$SNP, 
                                filename = "/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz", 
                                sep = "\t", 
                                snp_col = "rsids",
                                beta_col = "Beta",
                                se_col = "SE",
                                pval_col="minus_log10_pval",
                                eaf_col = "ImpMAF",
                                effect_allele_col = "effectAllele",
                                other_allele_col = "otherAllele", log_pval=TRUE)
out_METAP1$phenotype <- 'METAP1'


#NOTUM
out_NOTUM <- read_outcome_data(snps = preec_ivs$SNP, 
                               filename = "/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz", 
                               sep = "\t", 
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col="minus_log10_pval",
                               eaf_col = "ImpMAF",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele", log_pval=TRUE)
out_NOTUM$phenotype <- 'NOTUM'


#JUND
out_JUND <- read_outcome_data(snps = preec_ivs$SNP, 
                              filename = "/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz", 
                              sep = "\t", 
                              snp_col = "rsids",
                              beta_col = "Beta",
                              se_col = "SE",
                              pval_col="minus_log10_pval",
                              eaf_col = "ImpMAF",
                              effect_allele_col = "effectAllele",
                              other_allele_col = "otherAllele", log_pval=TRUE)
out_JUND$phenotype <- 'JUND'


#ADAMTS13
out_ADAMTS13 <- read_outcome_data(snps = preec_ivs$SNP, 
                                  filename = "/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz", 
                                  sep = "\t", 
                                  snp_col = "rsids",
                                  beta_col = "Beta",
                                  se_col = "SE",
                                  pval_col="minus_log10_pval",
                                  eaf_col = "ImpMAF",
                                  effect_allele_col = "effectAllele",
                                  other_allele_col = "otherAllele", log_pval=TRUE)
out_ADAMTS13$phenotype <- 'ADAMTS13'


#GDF15
out_GDF15 <- read_outcome_data(snps = preec_ivs$SNP, 
                               filename = "/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz", 
                               sep = "\t", 
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col="minus_log10_pval",
                               eaf_col = "ImpMAF",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele", log_pval=TRUE)
out_GDF15$phenotype <- 'GDF15'


#FGL1
out_FGL1 <- read_outcome_data(snps = preec_ivs$SNP, 
                              filename = "/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz", 
                              sep = "\t", 
                              snp_col = "rsids",
                              beta_col = "Beta",
                              se_col = "SE",
                              pval_col="minus_log10_pval",
                              eaf_col = "ImpMAF",
                              effect_allele_col = "effectAllele",
                              other_allele_col = "otherAllele", log_pval=TRUE)
out_FGL1$phenotype <- 'FGL1'


#SERPINE2
out_SERPINE2 <- read_outcome_data(snps = preec_ivs$SNP, 
                                  filename = "/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz", 
                                  sep = "\t", 
                                  snp_col = "rsids",
                                  beta_col = "Beta",
                                  se_col = "SE",
                                  pval_col="minus_log10_pval",
                                  eaf_col = "ImpMAF",
                                  effect_allele_col = "effectAllele",
                                  other_allele_col = "otherAllele", log_pval=TRUE)
out_SERPINE2$phenotype <- 'SERPINE2'


#RGS18
out_RGS18 <- read_outcome_data(snps = preec_ivs$SNP, 
                               filename = "/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz", 
                               sep = "\t", 
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col="minus_log10_pval",
                               eaf_col = "ImpMAF",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele", log_pval=TRUE)
out_RGS18$phenotype <- 'RGS18'


#ALDH2
out_ALDH2 <- read_outcome_data(snps = preec_ivs$SNP, 
                               filename = "/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz", 
                               sep = "\t", 
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col="minus_log10_pval",
                               eaf_col = "ImpMAF",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele", log_pval=TRUE)
out_ALDH2$phenotype <- 'ALDH2'


#MANEA
out_MANEA <- read_outcome_data(snps = preec_ivs$SNP, 
                               filename = "/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz", 
                               sep = "\t", 
                               snp_col = "rsids",
                               beta_col = "Beta",
                               se_col = "SE",
                               pval_col="minus_log10_pval",
                               eaf_col = "ImpMAF",
                               effect_allele_col = "effectAllele",
                               other_allele_col = "otherAllele", log_pval=TRUE)
out_MANEA$phenotype <- 'MANEA'

beep(2)

outlist <- ls(pattern = "out_", mget(ls())) 

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- gsub('out_', 'har_', outlist)
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

unlink("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_reverse_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_reverse_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_reverse_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/decode_harm_clump_reverse_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
clplist<-ls(pattern = "clp_", mget(ls()))

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/reverse_results_cis_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/reverse_results_cis_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/reverse_results_cis_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/reverse_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)


mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted')
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', n = nrow(mergedres))


mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/resdecode_preec/reverse_results_cis_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres[,1:3] <- NULL
mergedres$analysis <- 'Reverse direction'
write.csv(mergedres, "~/Desktop/preec/resdecode_preec/sens_reverse_cis_preec.csv", row.names = FALSE)


rm(list=ls())





#### sunukb - format exposure data ####
dir.create("/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse")
preec_ivs <- read.csv('~/desktop/preec/preec_ivs/poret_iv.csv')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

table(preec_ivs$Chromosome)


chr1_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '1'))[,c('pos', 'RefSNP_id')]
colnames(chr1_ref) <- c('pos37', 'SNP')
chr1_ref <- chr1_ref[which(chr1_ref$SNP %in% preec_ivs$RefSNP_id),]

chr2_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '2'))[,c('pos', 'RefSNP_id')]
colnames(chr2_ref) <- c('pos37', 'SNP')
chr2_ref <- chr2_ref[which(chr2_ref$SNP %in% preec_ivs$RefSNP_id),]

chr3_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '3'))[,c('pos', 'RefSNP_id')]
colnames(chr3_ref) <- c('pos37', 'SNP')
chr3_ref <- chr3_ref[which(chr3_ref$SNP %in% preec_ivs$RefSNP_id),]

chr4_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '4'))[,c('pos', 'RefSNP_id')]
colnames(chr4_ref) <- c('pos37', 'SNP')
chr4_ref <- chr4_ref[which(chr4_ref$SNP %in% preec_ivs$RefSNP_id),]

chr6_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '6'))[,c('pos', 'RefSNP_id')]
colnames(chr6_ref) <- c('pos37', 'SNP')
chr6_ref <- chr6_ref[which(chr6_ref$SNP %in% preec_ivs$RefSNP_id),]

chr7_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '7'))[,c('pos', 'RefSNP_id')]
colnames(chr7_ref) <- c('pos37', 'SNP')
chr7_ref <- chr7_ref[which(chr7_ref$SNP %in% preec_ivs$RefSNP_id),]

chr8_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '8'))[,c('pos', 'RefSNP_id')]
colnames(chr8_ref) <- c('pos37', 'SNP')
chr8_ref <- chr8_ref[which(chr8_ref$SNP %in% preec_ivs$RefSNP_id),]

chr10_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '10'))[,c('pos', 'RefSNP_id')]
colnames(chr10_ref) <- c('pos37', 'SNP')
chr10_ref <- chr10_ref[which(chr10_ref$SNP %in% preec_ivs$RefSNP_id),]

chr12_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12'))[,c('pos', 'RefSNP_id')]
colnames(chr12_ref) <- c('pos37', 'SNP')
chr12_ref <- chr12_ref[which(chr12_ref$SNP %in% preec_ivs$RefSNP_id),]

chr13_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '13'))[,c('pos', 'RefSNP_id')]
colnames(chr13_ref) <- c('pos37', 'SNP')
chr13_ref <- chr13_ref[which(chr13_ref$SNP %in% preec_ivs$RefSNP_id),]

chr15_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '15'))[,c('pos', 'RefSNP_id')]
colnames(chr15_ref) <- c('pos37', 'SNP')
chr15_ref <- chr15_ref[which(chr15_ref$SNP %in% preec_ivs$RefSNP_id),]

chr16_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16'))[,c('pos', 'RefSNP_id')]
colnames(chr16_ref) <- c('pos37', 'SNP')
chr16_ref <- chr16_ref[which(chr16_ref$SNP %in% preec_ivs$RefSNP_id),]

chr19_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '19'))[,c('pos', 'RefSNP_id')]
colnames(chr19_ref) <- c('pos37', 'SNP')
chr19_ref <- chr19_ref[which(chr19_ref$SNP %in% preec_ivs$RefSNP_id),]

chr20_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '20'))[,c('pos', 'RefSNP_id')]
colnames(chr20_ref) <- c('pos37', 'SNP')
chr20_ref <- chr20_ref[which(chr20_ref$SNP %in% preec_ivs$RefSNP_id),]

chr21_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '21'))[,c('pos', 'RefSNP_id')]
colnames(chr21_ref) <- c('pos37', 'SNP')
chr21_ref <- chr21_ref[which(chr21_ref$SNP %in% preec_ivs$RefSNP_id),]

chr22_ref = data.frame(snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '22'))[,c('pos', 'RefSNP_id')]
colnames(chr22_ref) <- c('pos37', 'SNP')
chr22_ref <- chr22_ref[which(chr22_ref$SNP %in% preec_ivs$RefSNP_id),]

#-## ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic #-##

# untar
untar('/volumes/maddy2/sunukb/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
ADAMTS13_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr1_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr1 <- cbind(ADAMTS13_chr1, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr1$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr1)[15] <- 'pos37'
ADAMTS13_chr1 <- merge(ADAMTS13_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
ADAMTS13_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr2_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr2 <- cbind(ADAMTS13_chr2, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr2$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr2)[15] <- 'pos37'
ADAMTS13_chr2 <- merge(ADAMTS13_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
ADAMTS13_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr3_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr3 <- cbind(ADAMTS13_chr3, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr3$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr3)[15] <- 'pos37'
ADAMTS13_chr3 <- merge(ADAMTS13_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
ADAMTS13_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr4_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr4 <- cbind(ADAMTS13_chr4, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr4$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr4)[15] <- 'pos37'
ADAMTS13_chr4 <- merge(ADAMTS13_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
ADAMTS13_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr6_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr6 <- cbind(ADAMTS13_chr6, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr6$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr6)[15] <- 'pos37'
ADAMTS13_chr6 <- merge(ADAMTS13_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
ADAMTS13_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr7_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr7 <- cbind(ADAMTS13_chr7, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr7$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr7)[15] <- 'pos37'
ADAMTS13_chr7 <- merge(ADAMTS13_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
ADAMTS13_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr8_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr8 <- cbind(ADAMTS13_chr8, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr8$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr8)[15] <- 'pos37'
ADAMTS13_chr8 <- merge(ADAMTS13_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
ADAMTS13_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr10_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr10 <- cbind(ADAMTS13_chr10, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr10$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr10)[15] <- 'pos37'
ADAMTS13_chr10 <- merge(ADAMTS13_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
ADAMTS13_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr12_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr12 <- cbind(ADAMTS13_chr12, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr12$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr12)[15] <- 'pos37'
ADAMTS13_chr12 <- merge(ADAMTS13_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
ADAMTS13_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr13_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr13 <- cbind(ADAMTS13_chr13, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr13$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr13)[15] <- 'pos37'
ADAMTS13_chr13 <- merge(ADAMTS13_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
ADAMTS13_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr15_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr15 <- cbind(ADAMTS13_chr15, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr15$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr15)[15] <- 'pos37'
ADAMTS13_chr15 <- merge(ADAMTS13_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
ADAMTS13_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr16_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr16 <- cbind(ADAMTS13_chr16, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr16$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr16)[15] <- 'pos37'
ADAMTS13_chr16 <- merge(ADAMTS13_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
ADAMTS13_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr19_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr19 <- cbind(ADAMTS13_chr19, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr19$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr19)[15] <- 'pos37'
ADAMTS13_chr19 <- merge(ADAMTS13_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
ADAMTS13_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr20_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr20 <- cbind(ADAMTS13_chr20, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr20$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr20)[15] <- 'pos37'
ADAMTS13_chr20 <- merge(ADAMTS13_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
ADAMTS13_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr21_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr21 <- cbind(ADAMTS13_chr21, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr21$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr21)[15] <- 'pos37'
ADAMTS13_chr21 <- merge(ADAMTS13_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
ADAMTS13_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr22_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz'))
ADAMTS13_chr22 <- cbind(ADAMTS13_chr22, data.frame(do.call('rbind', strsplit(as.character(ADAMTS13_chr22$ID),':',fixed=TRUE)))[,2])
colnames(ADAMTS13_chr22)[15] <- 'pos37'
ADAMTS13_chr22 <- merge(ADAMTS13_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


ADAMTS13_out <- rbind(ADAMTS13_chr1, ADAMTS13_chr2)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr3)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr4)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr6)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr7)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr8)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr10)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr12)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr13)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr15)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr16)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr19)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr20)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr21)
ADAMTS13_out <- rbind(ADAMTS13_out, ADAMTS13_chr22)
write.csv(ADAMTS13_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/ADAMTS13.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic', recursive = TRUE)



#-## APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II #-##

# untar
untar('/volumes/maddy2/sunukb/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
APOBR_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr1_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr1 <- cbind(APOBR_chr1, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr1$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr1)[15] <- 'pos37'
APOBR_chr1 <- merge(APOBR_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
APOBR_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr2_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr2 <- cbind(APOBR_chr2, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr2$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr2)[15] <- 'pos37'
APOBR_chr2 <- merge(APOBR_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
APOBR_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr3_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr3 <- cbind(APOBR_chr3, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr3$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr3)[15] <- 'pos37'
APOBR_chr3 <- merge(APOBR_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
APOBR_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr4_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr4 <- cbind(APOBR_chr4, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr4$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr4)[15] <- 'pos37'
APOBR_chr4 <- merge(APOBR_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
APOBR_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr6_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr6 <- cbind(APOBR_chr6, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr6$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr6)[15] <- 'pos37'
APOBR_chr6 <- merge(APOBR_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
APOBR_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr7_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr7 <- cbind(APOBR_chr7, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr7$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr7)[15] <- 'pos37'
APOBR_chr7 <- merge(APOBR_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
APOBR_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr8_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr8 <- cbind(APOBR_chr8, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr8$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr8)[15] <- 'pos37'
APOBR_chr8 <- merge(APOBR_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
APOBR_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr10_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr10 <- cbind(APOBR_chr10, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr10$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr10)[15] <- 'pos37'
APOBR_chr10 <- merge(APOBR_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
APOBR_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr12_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr12 <- cbind(APOBR_chr12, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr12$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr12)[15] <- 'pos37'
APOBR_chr12 <- merge(APOBR_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
APOBR_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr13_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr13 <- cbind(APOBR_chr13, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr13$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr13)[15] <- 'pos37'
APOBR_chr13 <- merge(APOBR_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
APOBR_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr15_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr15 <- cbind(APOBR_chr15, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr15$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr15)[15] <- 'pos37'
APOBR_chr15 <- merge(APOBR_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
APOBR_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr16 <- cbind(APOBR_chr16, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr16$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr16)[15] <- 'pos37'
APOBR_chr16 <- merge(APOBR_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
APOBR_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr19_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr19 <- cbind(APOBR_chr19, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr19$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr19)[15] <- 'pos37'
APOBR_chr19 <- merge(APOBR_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
APOBR_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr20_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr20 <- cbind(APOBR_chr20, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr20$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr20)[15] <- 'pos37'
APOBR_chr20 <- merge(APOBR_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
APOBR_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr21_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr21 <- cbind(APOBR_chr21, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr21$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr21)[15] <- 'pos37'
APOBR_chr21 <- merge(APOBR_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
APOBR_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr22_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz'))
APOBR_chr22 <- cbind(APOBR_chr22, data.frame(do.call('rbind', strsplit(as.character(APOBR_chr22$ID),':',fixed=TRUE)))[,2])
colnames(APOBR_chr22)[15] <- 'pos37'
APOBR_chr22 <- merge(APOBR_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


APOBR_out <- rbind(APOBR_chr1, APOBR_chr2)
APOBR_out <- rbind(APOBR_out, APOBR_chr3)
APOBR_out <- rbind(APOBR_out, APOBR_chr4)
APOBR_out <- rbind(APOBR_out, APOBR_chr6)
APOBR_out <- rbind(APOBR_out, APOBR_chr7)
APOBR_out <- rbind(APOBR_out, APOBR_chr8)
APOBR_out <- rbind(APOBR_out, APOBR_chr10)
APOBR_out <- rbind(APOBR_out, APOBR_chr12)
APOBR_out <- rbind(APOBR_out, APOBR_chr13)
APOBR_out <- rbind(APOBR_out, APOBR_chr15)
APOBR_out <- rbind(APOBR_out, APOBR_chr16)
APOBR_out <- rbind(APOBR_out, APOBR_chr19)
APOBR_out <- rbind(APOBR_out, APOBR_chr20)
APOBR_out <- rbind(APOBR_out, APOBR_chr21)
APOBR_out <- rbind(APOBR_out, APOBR_chr22)
write.csv(APOBR_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/APOBR.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II', recursive = TRUE)



#-## EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology.tar #-##

# untar
untar('/volumes/maddy2/sunukb/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
EBI3_IL27_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr1_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr1 <- cbind(EBI3_IL27_chr1, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr1$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr1)[15] <- 'pos37'
EBI3_IL27_chr1 <- merge(EBI3_IL27_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
EBI3_IL27_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr2_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr2 <- cbind(EBI3_IL27_chr2, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr2$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr2)[15] <- 'pos37'
EBI3_IL27_chr2 <- merge(EBI3_IL27_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
EBI3_IL27_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr3_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr3 <- cbind(EBI3_IL27_chr3, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr3$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr3)[15] <- 'pos37'
EBI3_IL27_chr3 <- merge(EBI3_IL27_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
EBI3_IL27_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr4_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr4 <- cbind(EBI3_IL27_chr4, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr4$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr4)[15] <- 'pos37'
EBI3_IL27_chr4 <- merge(EBI3_IL27_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
EBI3_IL27_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr6_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr6 <- cbind(EBI3_IL27_chr6, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr6$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr6)[15] <- 'pos37'
EBI3_IL27_chr6 <- merge(EBI3_IL27_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
EBI3_IL27_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr7_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr7 <- cbind(EBI3_IL27_chr7, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr7$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr7)[15] <- 'pos37'
EBI3_IL27_chr7 <- merge(EBI3_IL27_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
EBI3_IL27_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr8_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr8 <- cbind(EBI3_IL27_chr8, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr8$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr8)[15] <- 'pos37'
EBI3_IL27_chr8 <- merge(EBI3_IL27_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
EBI3_IL27_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr10_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr10 <- cbind(EBI3_IL27_chr10, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr10$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr10)[15] <- 'pos37'
EBI3_IL27_chr10 <- merge(EBI3_IL27_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
EBI3_IL27_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr12_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr12 <- cbind(EBI3_IL27_chr12, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr12$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr12)[15] <- 'pos37'
EBI3_IL27_chr12 <- merge(EBI3_IL27_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
EBI3_IL27_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr13_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr13 <- cbind(EBI3_IL27_chr13, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr13$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr13)[15] <- 'pos37'
EBI3_IL27_chr13 <- merge(EBI3_IL27_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
EBI3_IL27_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr15_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr15 <- cbind(EBI3_IL27_chr15, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr15$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr15)[15] <- 'pos37'
EBI3_IL27_chr15 <- merge(EBI3_IL27_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
EBI3_IL27_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr16_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr16 <- cbind(EBI3_IL27_chr16, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr16$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr16)[15] <- 'pos37'
EBI3_IL27_chr16 <- merge(EBI3_IL27_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
EBI3_IL27_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr19_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr19 <- cbind(EBI3_IL27_chr19, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr19$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr19)[15] <- 'pos37'
EBI3_IL27_chr19 <- merge(EBI3_IL27_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
EBI3_IL27_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr20_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr20 <- cbind(EBI3_IL27_chr20, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr20$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr20)[15] <- 'pos37'
EBI3_IL27_chr20 <- merge(EBI3_IL27_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
EBI3_IL27_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr21_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr21 <- cbind(EBI3_IL27_chr21, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr21$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr21)[15] <- 'pos37'
EBI3_IL27_chr21 <- merge(EBI3_IL27_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
EBI3_IL27_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology/discovery_chr22_EBI3_IL27:Q14213_Q8NEV9:OID21389:v1:Oncology.gz'))
EBI3_IL27_chr22 <- cbind(EBI3_IL27_chr22, data.frame(do.call('rbind', strsplit(as.character(EBI3_IL27_chr22$ID),':',fixed=TRUE)))[,2])
colnames(EBI3_IL27_chr22)[15] <- 'pos37'
EBI3_IL27_chr22 <- merge(EBI3_IL27_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


EBI3_IL27_out <- rbind(EBI3_IL27_chr1, EBI3_IL27_chr2)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr3)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr4)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr6)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr7)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr8)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr10)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr12)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr13)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr15)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr16)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr19)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr20)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr21)
EBI3_IL27_out <- rbind(EBI3_IL27_out, EBI3_IL27_chr22)
write.csv(EBI3_IL27_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/EBI3_IL27.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology', recursive = TRUE)



#-## FES_P07332_OID21207_v1_Oncology.tar #-##

# untar
untar('/volumes/maddy2/sunukb/FES_P07332_OID21207_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
FES_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr1_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr1 <- cbind(FES_chr1, data.frame(do.call('rbind', strsplit(as.character(FES_chr1$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr1)[15] <- 'pos37'
FES_chr1 <- merge(FES_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
FES_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr2_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr2 <- cbind(FES_chr2, data.frame(do.call('rbind', strsplit(as.character(FES_chr2$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr2)[15] <- 'pos37'
FES_chr2 <- merge(FES_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
FES_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr3_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr3 <- cbind(FES_chr3, data.frame(do.call('rbind', strsplit(as.character(FES_chr3$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr3)[15] <- 'pos37'
FES_chr3 <- merge(FES_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
FES_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr4_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr4 <- cbind(FES_chr4, data.frame(do.call('rbind', strsplit(as.character(FES_chr4$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr4)[15] <- 'pos37'
FES_chr4 <- merge(FES_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
FES_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr6_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr6 <- cbind(FES_chr6, data.frame(do.call('rbind', strsplit(as.character(FES_chr6$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr6)[15] <- 'pos37'
FES_chr6 <- merge(FES_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
FES_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr7_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr7 <- cbind(FES_chr7, data.frame(do.call('rbind', strsplit(as.character(FES_chr7$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr7)[15] <- 'pos37'
FES_chr7 <- merge(FES_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
FES_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr8_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr8 <- cbind(FES_chr8, data.frame(do.call('rbind', strsplit(as.character(FES_chr8$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr8)[15] <- 'pos37'
FES_chr8 <- merge(FES_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
FES_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr10_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr10 <- cbind(FES_chr10, data.frame(do.call('rbind', strsplit(as.character(FES_chr10$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr10)[15] <- 'pos37'
FES_chr10 <- merge(FES_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
FES_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr12_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr12 <- cbind(FES_chr12, data.frame(do.call('rbind', strsplit(as.character(FES_chr12$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr12)[15] <- 'pos37'
FES_chr12 <- merge(FES_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
FES_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr13_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr13 <- cbind(FES_chr13, data.frame(do.call('rbind', strsplit(as.character(FES_chr13$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr13)[15] <- 'pos37'
FES_chr13 <- merge(FES_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
FES_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr15 <- cbind(FES_chr15, data.frame(do.call('rbind', strsplit(as.character(FES_chr15$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr15)[15] <- 'pos37'
FES_chr15 <- merge(FES_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
FES_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr16_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr16 <- cbind(FES_chr16, data.frame(do.call('rbind', strsplit(as.character(FES_chr16$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr16)[15] <- 'pos37'
FES_chr16 <- merge(FES_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
FES_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr19_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr19 <- cbind(FES_chr19, data.frame(do.call('rbind', strsplit(as.character(FES_chr19$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr19)[15] <- 'pos37'
FES_chr19 <- merge(FES_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
FES_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr20_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr20 <- cbind(FES_chr20, data.frame(do.call('rbind', strsplit(as.character(FES_chr20$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr20)[15] <- 'pos37'
FES_chr20 <- merge(FES_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
FES_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr21_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr21 <- cbind(FES_chr21, data.frame(do.call('rbind', strsplit(as.character(FES_chr21$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr21)[15] <- 'pos37'
FES_chr21 <- merge(FES_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
FES_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr22_FES:P07332:OID21207:v1:Oncology.gz'))
FES_chr22 <- cbind(FES_chr22, data.frame(do.call('rbind', strsplit(as.character(FES_chr22$ID),':',fixed=TRUE)))[,2])
colnames(FES_chr22)[15] <- 'pos37'
FES_chr22 <- merge(FES_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


FES_out <- rbind(FES_chr1, FES_chr2)
FES_out <- rbind(FES_out, FES_chr3)
FES_out <- rbind(FES_out, FES_chr4)
FES_out <- rbind(FES_out, FES_chr6)
FES_out <- rbind(FES_out, FES_chr7)
FES_out <- rbind(FES_out, FES_chr8)
FES_out <- rbind(FES_out, FES_chr10)
FES_out <- rbind(FES_out, FES_chr12)
FES_out <- rbind(FES_out, FES_chr13)
FES_out <- rbind(FES_out, FES_chr15)
FES_out <- rbind(FES_out, FES_chr16)
FES_out <- rbind(FES_out, FES_chr19)
FES_out <- rbind(FES_out, FES_chr20)
FES_out <- rbind(FES_out, FES_chr21)
FES_out <- rbind(FES_out, FES_chr22)
write.csv(FES_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FES.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology', recursive = TRUE)



#-## FGF5_P12034_OID20490_v1_Inflammation.tar #-##

# untar
untar('/volumes/maddy2/sunukb/FGF5_P12034_OID20490_v1_Inflammation.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
FGF5_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr1_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr1 <- cbind(FGF5_chr1, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr1$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr1)[15] <- 'pos37'
FGF5_chr1 <- merge(FGF5_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
FGF5_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr2_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr2 <- cbind(FGF5_chr2, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr2$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr2)[15] <- 'pos37'
FGF5_chr2 <- merge(FGF5_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
FGF5_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr3_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr3 <- cbind(FGF5_chr3, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr3$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr3)[15] <- 'pos37'
FGF5_chr3 <- merge(FGF5_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
FGF5_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr4 <- cbind(FGF5_chr4, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr4$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr4)[15] <- 'pos37'
FGF5_chr4 <- merge(FGF5_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
FGF5_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr6_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr6 <- cbind(FGF5_chr6, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr6$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr6)[15] <- 'pos37'
FGF5_chr6 <- merge(FGF5_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
FGF5_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr7_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr7 <- cbind(FGF5_chr7, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr7$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr7)[15] <- 'pos37'
FGF5_chr7 <- merge(FGF5_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
FGF5_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr8_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr8 <- cbind(FGF5_chr8, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr8$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr8)[15] <- 'pos37'
FGF5_chr8 <- merge(FGF5_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
FGF5_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr10_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr10 <- cbind(FGF5_chr10, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr10$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr10)[15] <- 'pos37'
FGF5_chr10 <- merge(FGF5_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
FGF5_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr12_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr12 <- cbind(FGF5_chr12, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr12$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr12)[15] <- 'pos37'
FGF5_chr12 <- merge(FGF5_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
FGF5_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr13_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr13 <- cbind(FGF5_chr13, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr13$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr13)[15] <- 'pos37'
FGF5_chr13 <- merge(FGF5_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
FGF5_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr15_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr15 <- cbind(FGF5_chr15, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr15$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr15)[15] <- 'pos37'
FGF5_chr15 <- merge(FGF5_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
FGF5_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr16_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr16 <- cbind(FGF5_chr16, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr16$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr16)[15] <- 'pos37'
FGF5_chr16 <- merge(FGF5_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
FGF5_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr19_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr19 <- cbind(FGF5_chr19, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr19$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr19)[15] <- 'pos37'
FGF5_chr19 <- merge(FGF5_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
FGF5_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr20_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr20 <- cbind(FGF5_chr20, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr20$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr20)[15] <- 'pos37'
FGF5_chr20 <- merge(FGF5_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
FGF5_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr21_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr21 <- cbind(FGF5_chr21, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr21$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr21)[15] <- 'pos37'
FGF5_chr21 <- merge(FGF5_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
FGF5_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr22_FGF5:P12034:OID20490:v1:Inflammation.gz'))
FGF5_chr22 <- cbind(FGF5_chr22, data.frame(do.call('rbind', strsplit(as.character(FGF5_chr22$ID),':',fixed=TRUE)))[,2])
colnames(FGF5_chr22)[15] <- 'pos37'
FGF5_chr22 <- merge(FGF5_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


FGF5_out <- rbind(FGF5_chr1, FGF5_chr2)
FGF5_out <- rbind(FGF5_out, FGF5_chr3)
FGF5_out <- rbind(FGF5_out, FGF5_chr4)
FGF5_out <- rbind(FGF5_out, FGF5_chr6)
FGF5_out <- rbind(FGF5_out, FGF5_chr7)
FGF5_out <- rbind(FGF5_out, FGF5_chr8)
FGF5_out <- rbind(FGF5_out, FGF5_chr10)
FGF5_out <- rbind(FGF5_out, FGF5_chr12)
FGF5_out <- rbind(FGF5_out, FGF5_chr13)
FGF5_out <- rbind(FGF5_out, FGF5_chr15)
FGF5_out <- rbind(FGF5_out, FGF5_chr16)
FGF5_out <- rbind(FGF5_out, FGF5_chr19)
FGF5_out <- rbind(FGF5_out, FGF5_chr20)
FGF5_out <- rbind(FGF5_out, FGF5_chr21)
FGF5_out <- rbind(FGF5_out, FGF5_chr22)
write.csv(FGF5_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FGF5.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation', recursive = TRUE)



#-## FGL1_Q08830_OID30702_v1_Inflammation_II.tar #-##

# untar
untar('/volumes/maddy2/sunukb/FGL1_Q08830_OID30702_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
FGL1_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr1_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr1 <- cbind(FGL1_chr1, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr1$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr1)[15] <- 'pos37'
FGL1_chr1 <- merge(FGL1_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
FGL1_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr2_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr2 <- cbind(FGL1_chr2, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr2$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr2)[15] <- 'pos37'
FGL1_chr2 <- merge(FGL1_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
FGL1_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr3_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr3 <- cbind(FGL1_chr3, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr3$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr3)[15] <- 'pos37'
FGL1_chr3 <- merge(FGL1_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
FGL1_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr4_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr4 <- cbind(FGL1_chr4, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr4$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr4)[15] <- 'pos37'
FGL1_chr4 <- merge(FGL1_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
FGL1_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr6_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr6 <- cbind(FGL1_chr6, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr6$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr6)[15] <- 'pos37'
FGL1_chr6 <- merge(FGL1_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
FGL1_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr7_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr7 <- cbind(FGL1_chr7, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr7$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr7)[15] <- 'pos37'
FGL1_chr7 <- merge(FGL1_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
FGL1_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr8 <- cbind(FGL1_chr8, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr8$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr8)[15] <- 'pos37'
FGL1_chr8 <- merge(FGL1_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
FGL1_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr10_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr10 <- cbind(FGL1_chr10, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr10$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr10)[15] <- 'pos37'
FGL1_chr10 <- merge(FGL1_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
FGL1_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr12_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr12 <- cbind(FGL1_chr12, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr12$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr12)[15] <- 'pos37'
FGL1_chr12 <- merge(FGL1_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
FGL1_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr13_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr13 <- cbind(FGL1_chr13, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr13$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr13)[15] <- 'pos37'
FGL1_chr13 <- merge(FGL1_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
FGL1_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr15_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr15 <- cbind(FGL1_chr15, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr15$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr15)[15] <- 'pos37'
FGL1_chr15 <- merge(FGL1_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
FGL1_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr16_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr16 <- cbind(FGL1_chr16, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr16$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr16)[15] <- 'pos37'
FGL1_chr16 <- merge(FGL1_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
FGL1_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr19_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr19 <- cbind(FGL1_chr19, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr19$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr19)[15] <- 'pos37'
FGL1_chr19 <- merge(FGL1_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
FGL1_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr20_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr20 <- cbind(FGL1_chr20, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr20$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr20)[15] <- 'pos37'
FGL1_chr20 <- merge(FGL1_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
FGL1_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr21_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr21 <- cbind(FGL1_chr21, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr21$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr21)[15] <- 'pos37'
FGL1_chr21 <- merge(FGL1_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
FGL1_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr22_FGL1:Q08830:OID30702:v1:Inflammation_II.gz'))
FGL1_chr22 <- cbind(FGL1_chr22, data.frame(do.call('rbind', strsplit(as.character(FGL1_chr22$ID),':',fixed=TRUE)))[,2])
colnames(FGL1_chr22)[15] <- 'pos37'
FGL1_chr22 <- merge(FGL1_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


FGL1_out <- rbind(FGL1_chr1, FGL1_chr2)
FGL1_out <- rbind(FGL1_out, FGL1_chr3)
FGL1_out <- rbind(FGL1_out, FGL1_chr4)
FGL1_out <- rbind(FGL1_out, FGL1_chr6)
FGL1_out <- rbind(FGL1_out, FGL1_chr7)
FGL1_out <- rbind(FGL1_out, FGL1_chr8)
FGL1_out <- rbind(FGL1_out, FGL1_chr10)
FGL1_out <- rbind(FGL1_out, FGL1_chr12)
FGL1_out <- rbind(FGL1_out, FGL1_chr13)
FGL1_out <- rbind(FGL1_out, FGL1_chr15)
FGL1_out <- rbind(FGL1_out, FGL1_chr16)
FGL1_out <- rbind(FGL1_out, FGL1_chr19)
FGL1_out <- rbind(FGL1_out, FGL1_chr20)
FGL1_out <- rbind(FGL1_out, FGL1_chr21)
FGL1_out <- rbind(FGL1_out, FGL1_chr22)
write.csv(FGL1_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FGL1.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II', recursive = TRUE)



#-## GDF15_Q99988_OID20251_v1_Cardiometabolic.tar #-##

# untar
untar('/volumes/maddy2/sunukb/GDF15_Q99988_OID20251_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
GDF15_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr1_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr1 <- cbind(GDF15_chr1, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr1$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr1)[15] <- 'pos37'
GDF15_chr1 <- merge(GDF15_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
GDF15_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr2_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr2 <- cbind(GDF15_chr2, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr2$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr2)[15] <- 'pos37'
GDF15_chr2 <- merge(GDF15_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
GDF15_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr3_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr3 <- cbind(GDF15_chr3, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr3$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr3)[15] <- 'pos37'
GDF15_chr3 <- merge(GDF15_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
GDF15_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr4_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr4 <- cbind(GDF15_chr4, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr4$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr4)[15] <- 'pos37'
GDF15_chr4 <- merge(GDF15_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
GDF15_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr6_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr6 <- cbind(GDF15_chr6, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr6$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr6)[15] <- 'pos37'
GDF15_chr6 <- merge(GDF15_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
GDF15_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr7_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr7 <- cbind(GDF15_chr7, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr7$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr7)[15] <- 'pos37'
GDF15_chr7 <- merge(GDF15_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
GDF15_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr8_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr8 <- cbind(GDF15_chr8, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr8$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr8)[15] <- 'pos37'
GDF15_chr8 <- merge(GDF15_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
GDF15_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr10_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr10 <- cbind(GDF15_chr10, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr10$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr10)[15] <- 'pos37'
GDF15_chr10 <- merge(GDF15_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
GDF15_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr12_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr12 <- cbind(GDF15_chr12, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr12$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr12)[15] <- 'pos37'
GDF15_chr12 <- merge(GDF15_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
GDF15_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr13_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr13 <- cbind(GDF15_chr13, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr13$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr13)[15] <- 'pos37'
GDF15_chr13 <- merge(GDF15_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
GDF15_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr15_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr15 <- cbind(GDF15_chr15, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr15$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr15)[15] <- 'pos37'
GDF15_chr15 <- merge(GDF15_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
GDF15_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr16_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr16 <- cbind(GDF15_chr16, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr16$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr16)[15] <- 'pos37'
GDF15_chr16 <- merge(GDF15_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
GDF15_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr19 <- cbind(GDF15_chr19, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr19$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr19)[15] <- 'pos37'
GDF15_chr19 <- merge(GDF15_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
GDF15_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr20_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr20 <- cbind(GDF15_chr20, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr20$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr20)[15] <- 'pos37'
GDF15_chr20 <- merge(GDF15_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
GDF15_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr21_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr21 <- cbind(GDF15_chr21, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr21$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr21)[15] <- 'pos37'
GDF15_chr21 <- merge(GDF15_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
GDF15_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr22_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz'))
GDF15_chr22 <- cbind(GDF15_chr22, data.frame(do.call('rbind', strsplit(as.character(GDF15_chr22$ID),':',fixed=TRUE)))[,2])
colnames(GDF15_chr22)[15] <- 'pos37'
GDF15_chr22 <- merge(GDF15_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


GDF15_out <- rbind(GDF15_chr1, GDF15_chr2)
GDF15_out <- rbind(GDF15_out, GDF15_chr3)
GDF15_out <- rbind(GDF15_out, GDF15_chr4)
GDF15_out <- rbind(GDF15_out, GDF15_chr6)
GDF15_out <- rbind(GDF15_out, GDF15_chr7)
GDF15_out <- rbind(GDF15_out, GDF15_chr8)
GDF15_out <- rbind(GDF15_out, GDF15_chr10)
GDF15_out <- rbind(GDF15_out, GDF15_chr12)
GDF15_out <- rbind(GDF15_out, GDF15_chr13)
GDF15_out <- rbind(GDF15_out, GDF15_chr15)
GDF15_out <- rbind(GDF15_out, GDF15_chr16)
GDF15_out <- rbind(GDF15_out, GDF15_chr19)
GDF15_out <- rbind(GDF15_out, GDF15_chr20)
GDF15_out <- rbind(GDF15_out, GDF15_chr21)
GDF15_out <- rbind(GDF15_out, GDF15_chr22)
write.csv(GDF15_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/GDF15.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic', recursive = TRUE)



#-## PZP_P20742_OID30730_v1_Inflammation_II.tar #-##

# untar
untar('/volumes/maddy2/sunukb/PZP_P20742_OID30730_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
PZP_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr1_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr1 <- cbind(PZP_chr1, data.frame(do.call('rbind', strsplit(as.character(PZP_chr1$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr1)[15] <- 'pos37'
PZP_chr1 <- merge(PZP_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
PZP_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr2_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr2 <- cbind(PZP_chr2, data.frame(do.call('rbind', strsplit(as.character(PZP_chr2$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr2)[15] <- 'pos37'
PZP_chr2 <- merge(PZP_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
PZP_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr3_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr3 <- cbind(PZP_chr3, data.frame(do.call('rbind', strsplit(as.character(PZP_chr3$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr3)[15] <- 'pos37'
PZP_chr3 <- merge(PZP_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
PZP_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr4_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr4 <- cbind(PZP_chr4, data.frame(do.call('rbind', strsplit(as.character(PZP_chr4$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr4)[15] <- 'pos37'
PZP_chr4 <- merge(PZP_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
PZP_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr6_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr6 <- cbind(PZP_chr6, data.frame(do.call('rbind', strsplit(as.character(PZP_chr6$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr6)[15] <- 'pos37'
PZP_chr6 <- merge(PZP_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
PZP_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr7_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr7 <- cbind(PZP_chr7, data.frame(do.call('rbind', strsplit(as.character(PZP_chr7$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr7)[15] <- 'pos37'
PZP_chr7 <- merge(PZP_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
PZP_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr8_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr8 <- cbind(PZP_chr8, data.frame(do.call('rbind', strsplit(as.character(PZP_chr8$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr8)[15] <- 'pos37'
PZP_chr8 <- merge(PZP_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
PZP_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr10_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr10 <- cbind(PZP_chr10, data.frame(do.call('rbind', strsplit(as.character(PZP_chr10$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr10)[15] <- 'pos37'
PZP_chr10 <- merge(PZP_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
PZP_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr12 <- cbind(PZP_chr12, data.frame(do.call('rbind', strsplit(as.character(PZP_chr12$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr12)[15] <- 'pos37'
PZP_chr12 <- merge(PZP_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
PZP_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr13_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr13 <- cbind(PZP_chr13, data.frame(do.call('rbind', strsplit(as.character(PZP_chr13$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr13)[15] <- 'pos37'
PZP_chr13 <- merge(PZP_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
PZP_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr15_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr15 <- cbind(PZP_chr15, data.frame(do.call('rbind', strsplit(as.character(PZP_chr15$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr15)[15] <- 'pos37'
PZP_chr15 <- merge(PZP_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
PZP_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr16_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr16 <- cbind(PZP_chr16, data.frame(do.call('rbind', strsplit(as.character(PZP_chr16$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr16)[15] <- 'pos37'
PZP_chr16 <- merge(PZP_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
PZP_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr19_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr19 <- cbind(PZP_chr19, data.frame(do.call('rbind', strsplit(as.character(PZP_chr19$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr19)[15] <- 'pos37'
PZP_chr19 <- merge(PZP_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
PZP_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr20_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr20 <- cbind(PZP_chr20, data.frame(do.call('rbind', strsplit(as.character(PZP_chr20$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr20)[15] <- 'pos37'
PZP_chr20 <- merge(PZP_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
PZP_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr21_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr21 <- cbind(PZP_chr21, data.frame(do.call('rbind', strsplit(as.character(PZP_chr21$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr21)[15] <- 'pos37'
PZP_chr21 <- merge(PZP_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
PZP_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr22_PZP:P20742:OID30730:v1:Inflammation_II.gz'))
PZP_chr22 <- cbind(PZP_chr22, data.frame(do.call('rbind', strsplit(as.character(PZP_chr22$ID),':',fixed=TRUE)))[,2])
colnames(PZP_chr22)[15] <- 'pos37'
PZP_chr22 <- merge(PZP_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


PZP_out <- rbind(PZP_chr1, PZP_chr2)
PZP_out <- rbind(PZP_out, PZP_chr3)
PZP_out <- rbind(PZP_out, PZP_chr4)
PZP_out <- rbind(PZP_out, PZP_chr6)
PZP_out <- rbind(PZP_out, PZP_chr7)
PZP_out <- rbind(PZP_out, PZP_chr8)
PZP_out <- rbind(PZP_out, PZP_chr10)
PZP_out <- rbind(PZP_out, PZP_chr12)
PZP_out <- rbind(PZP_out, PZP_chr13)
PZP_out <- rbind(PZP_out, PZP_chr15)
PZP_out <- rbind(PZP_out, PZP_chr16)
PZP_out <- rbind(PZP_out, PZP_chr19)
PZP_out <- rbind(PZP_out, PZP_chr20)
PZP_out <- rbind(PZP_out, PZP_chr21)
PZP_out <- rbind(PZP_out, PZP_chr22)
write.csv(PZP_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/PZP.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II', recursive = TRUE)



#-## SERPINE2_P07093_OID30359_v1_Cardiometabolic_II.tar #-##

# untar
untar('/volumes/maddy2/sunukb/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
SERPINE2_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr1_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr1 <- cbind(SERPINE2_chr1, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr1$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr1)[15] <- 'pos37'
SERPINE2_chr1 <- merge(SERPINE2_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
SERPINE2_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr2 <- cbind(SERPINE2_chr2, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr2$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr2)[15] <- 'pos37'
SERPINE2_chr2 <- merge(SERPINE2_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
SERPINE2_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr3_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr3 <- cbind(SERPINE2_chr3, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr3$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr3)[15] <- 'pos37'
SERPINE2_chr3 <- merge(SERPINE2_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
SERPINE2_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr4_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr4 <- cbind(SERPINE2_chr4, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr4$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr4)[15] <- 'pos37'
SERPINE2_chr4 <- merge(SERPINE2_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
SERPINE2_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr6_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr6 <- cbind(SERPINE2_chr6, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr6$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr6)[15] <- 'pos37'
SERPINE2_chr6 <- merge(SERPINE2_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
SERPINE2_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr7_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr7 <- cbind(SERPINE2_chr7, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr7$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr7)[15] <- 'pos37'
SERPINE2_chr7 <- merge(SERPINE2_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
SERPINE2_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr8_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr8 <- cbind(SERPINE2_chr8, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr8$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr8)[15] <- 'pos37'
SERPINE2_chr8 <- merge(SERPINE2_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
SERPINE2_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr10_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr10 <- cbind(SERPINE2_chr10, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr10$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr10)[15] <- 'pos37'
SERPINE2_chr10 <- merge(SERPINE2_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
SERPINE2_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr12_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr12 <- cbind(SERPINE2_chr12, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr12$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr12)[15] <- 'pos37'
SERPINE2_chr12 <- merge(SERPINE2_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
SERPINE2_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr13_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr13 <- cbind(SERPINE2_chr13, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr13$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr13)[15] <- 'pos37'
SERPINE2_chr13 <- merge(SERPINE2_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
SERPINE2_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr15_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr15 <- cbind(SERPINE2_chr15, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr15$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr15)[15] <- 'pos37'
SERPINE2_chr15 <- merge(SERPINE2_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
SERPINE2_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr16_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr16 <- cbind(SERPINE2_chr16, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr16$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr16)[15] <- 'pos37'
SERPINE2_chr16 <- merge(SERPINE2_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
SERPINE2_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr19_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr19 <- cbind(SERPINE2_chr19, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr19$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr19)[15] <- 'pos37'
SERPINE2_chr19 <- merge(SERPINE2_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
SERPINE2_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr20_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr20 <- cbind(SERPINE2_chr20, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr20$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr20)[15] <- 'pos37'
SERPINE2_chr20 <- merge(SERPINE2_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
SERPINE2_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr21_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr21 <- cbind(SERPINE2_chr21, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr21$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr21)[15] <- 'pos37'
SERPINE2_chr21 <- merge(SERPINE2_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
SERPINE2_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr22_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz'))
SERPINE2_chr22 <- cbind(SERPINE2_chr22, data.frame(do.call('rbind', strsplit(as.character(SERPINE2_chr22$ID),':',fixed=TRUE)))[,2])
colnames(SERPINE2_chr22)[15] <- 'pos37'
SERPINE2_chr22 <- merge(SERPINE2_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


SERPINE2_out <- rbind(SERPINE2_chr1, SERPINE2_chr2)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr3)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr4)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr6)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr7)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr8)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr10)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr12)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr13)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr15)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr16)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr19)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr20)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr21)
SERPINE2_out <- rbind(SERPINE2_out, SERPINE2_chr22)
write.csv(SERPINE2_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SERPINE2.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II', recursive = TRUE)



#-## SH2B3_Q9UQQ2_OID21222_v1_Oncology.tar #-##

# untar
untar('/volumes/maddy2/sunukb/SH2B3_Q9UQQ2_OID21222_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
SH2B3_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr1_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr1 <- cbind(SH2B3_chr1, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr1$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr1)[15] <- 'pos37'
SH2B3_chr1 <- merge(SH2B3_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
SH2B3_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr2_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr2 <- cbind(SH2B3_chr2, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr2$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr2)[15] <- 'pos37'
SH2B3_chr2 <- merge(SH2B3_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
SH2B3_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr3_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr3 <- cbind(SH2B3_chr3, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr3$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr3)[15] <- 'pos37'
SH2B3_chr3 <- merge(SH2B3_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
SH2B3_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr4_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr4 <- cbind(SH2B3_chr4, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr4$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr4)[15] <- 'pos37'
SH2B3_chr4 <- merge(SH2B3_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
SH2B3_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr6_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr6 <- cbind(SH2B3_chr6, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr6$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr6)[15] <- 'pos37'
SH2B3_chr6 <- merge(SH2B3_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
SH2B3_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr7_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr7 <- cbind(SH2B3_chr7, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr7$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr7)[15] <- 'pos37'
SH2B3_chr7 <- merge(SH2B3_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
SH2B3_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr8_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr8 <- cbind(SH2B3_chr8, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr8$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr8)[15] <- 'pos37'
SH2B3_chr8 <- merge(SH2B3_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
SH2B3_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr10_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr10 <- cbind(SH2B3_chr10, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr10$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr10)[15] <- 'pos37'
SH2B3_chr10 <- merge(SH2B3_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
SH2B3_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr12 <- cbind(SH2B3_chr12, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr12$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr12)[15] <- 'pos37'
SH2B3_chr12 <- merge(SH2B3_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
SH2B3_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr13_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr13 <- cbind(SH2B3_chr13, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr13$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr13)[15] <- 'pos37'
SH2B3_chr13 <- merge(SH2B3_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
SH2B3_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr15_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr15 <- cbind(SH2B3_chr15, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr15$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr15)[15] <- 'pos37'
SH2B3_chr15 <- merge(SH2B3_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
SH2B3_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr16_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr16 <- cbind(SH2B3_chr16, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr16$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr16)[15] <- 'pos37'
SH2B3_chr16 <- merge(SH2B3_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
SH2B3_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr19_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr19 <- cbind(SH2B3_chr19, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr19$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr19)[15] <- 'pos37'
SH2B3_chr19 <- merge(SH2B3_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
SH2B3_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr20_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr20 <- cbind(SH2B3_chr20, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr20$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr20)[15] <- 'pos37'
SH2B3_chr20 <- merge(SH2B3_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
SH2B3_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr21_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr21 <- cbind(SH2B3_chr21, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr21$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr21)[15] <- 'pos37'
SH2B3_chr21 <- merge(SH2B3_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
SH2B3_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr22_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz'))
SH2B3_chr22 <- cbind(SH2B3_chr22, data.frame(do.call('rbind', strsplit(as.character(SH2B3_chr22$ID),':',fixed=TRUE)))[,2])
colnames(SH2B3_chr22)[15] <- 'pos37'
SH2B3_chr22 <- merge(SH2B3_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


SH2B3_out <- rbind(SH2B3_chr1, SH2B3_chr2)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr3)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr4)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr6)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr7)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr8)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr10)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr12)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr13)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr15)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr16)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr19)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr20)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr21)
SH2B3_out <- rbind(SH2B3_out, SH2B3_chr22)
write.csv(SH2B3_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SH2B3.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology', recursive = TRUE)



#-## SULT1A1_P50225_OID21031_v1_Neurology.tar #-##

# untar
untar('/volumes/maddy2/sunukb/SULT1A1_P50225_OID21031_v1_Neurology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
#chr1
SULT1A1_chr1 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr1_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr1 <- cbind(SULT1A1_chr1, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr1$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr1)[15] <- 'pos37'
SULT1A1_chr1 <- merge(SULT1A1_chr1, chr1_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr2
SULT1A1_chr2 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr2_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr2 <- cbind(SULT1A1_chr2, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr2$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr2)[15] <- 'pos37'
SULT1A1_chr2 <- merge(SULT1A1_chr2, chr2_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr3
SULT1A1_chr3 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr3_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr3 <- cbind(SULT1A1_chr3, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr3$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr3)[15] <- 'pos37'
SULT1A1_chr3 <- merge(SULT1A1_chr3, chr3_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr4
SULT1A1_chr4 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr4_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr4 <- cbind(SULT1A1_chr4, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr4$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr4)[15] <- 'pos37'
SULT1A1_chr4 <- merge(SULT1A1_chr4, chr4_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr6
SULT1A1_chr6 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr6_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr6 <- cbind(SULT1A1_chr6, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr6$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr6)[15] <- 'pos37'
SULT1A1_chr6 <- merge(SULT1A1_chr6, chr6_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr7
SULT1A1_chr7 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr7_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr7 <- cbind(SULT1A1_chr7, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr7$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr7)[15] <- 'pos37'
SULT1A1_chr7 <- merge(SULT1A1_chr7, chr7_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr8
SULT1A1_chr8 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr8_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr8 <- cbind(SULT1A1_chr8, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr8$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr8)[15] <- 'pos37'
SULT1A1_chr8 <- merge(SULT1A1_chr8, chr8_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr10
SULT1A1_chr10 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr10_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr10 <- cbind(SULT1A1_chr10, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr10$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr10)[15] <- 'pos37'
SULT1A1_chr10 <- merge(SULT1A1_chr10, chr10_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr12
SULT1A1_chr12 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr12_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr12 <- cbind(SULT1A1_chr12, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr12$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr12)[15] <- 'pos37'
SULT1A1_chr12 <- merge(SULT1A1_chr12, chr12_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr13
SULT1A1_chr13 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr13_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr13 <- cbind(SULT1A1_chr13, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr13$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr13)[15] <- 'pos37'
SULT1A1_chr13 <- merge(SULT1A1_chr13, chr13_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr15
SULT1A1_chr15 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr15_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr15 <- cbind(SULT1A1_chr15, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr15$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr15)[15] <- 'pos37'
SULT1A1_chr15 <- merge(SULT1A1_chr15, chr15_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr16
SULT1A1_chr16 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr16 <- cbind(SULT1A1_chr16, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr16$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr16)[15] <- 'pos37'
SULT1A1_chr16 <- merge(SULT1A1_chr16, chr16_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr19
SULT1A1_chr19 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr19_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr19 <- cbind(SULT1A1_chr19, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr19$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr19)[15] <- 'pos37'
SULT1A1_chr19 <- merge(SULT1A1_chr19, chr19_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr20
SULT1A1_chr20 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr20_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr20 <- cbind(SULT1A1_chr20, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr20$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr20)[15] <- 'pos37'
SULT1A1_chr20 <- merge(SULT1A1_chr20, chr20_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr21
SULT1A1_chr21 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr21_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr21 <- cbind(SULT1A1_chr21, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr21$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr21)[15] <- 'pos37'
SULT1A1_chr21 <- merge(SULT1A1_chr21, chr21_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)

#chr22
SULT1A1_chr22 <- as.data.frame(fread('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr22_SULT1A1:P50225:OID21031:v1:Neurology.gz'))
SULT1A1_chr22 <- cbind(SULT1A1_chr22, data.frame(do.call('rbind', strsplit(as.character(SULT1A1_chr22$ID),':',fixed=TRUE)))[,2])
colnames(SULT1A1_chr22)[15] <- 'pos37'
SULT1A1_chr22 <- merge(SULT1A1_chr22, chr22_ref, by = 'pos37', all.x=FALSE, all.y=FALSE)


SULT1A1_out <- rbind(SULT1A1_chr1, SULT1A1_chr2)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr3)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr4)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr6)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr7)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr8)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr10)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr12)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr13)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr15)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr16)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr19)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr20)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr21)
SULT1A1_out <- rbind(SULT1A1_out, SULT1A1_chr22)
write.csv(SULT1A1_out, '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SULT1A1.csv')

# delete untarred file
unlink('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology', recursive = TRUE)





#-## CLEAR WE #-##
rm(list=ls())

#### sunukb - reverse mr ####
# Import all IVs
preec_ivs <- read.csv('~/desktop/preec/preec_ivs/poret_iv.csv')
genelist <- 'preec_ivs'
preec_ivs$phenotype <- 'Pre-eclampsia'

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- format_data(dat,type="exposure",snp_col = "RefSNP_id",
                     beta_col = "poret_beta",
                     se_col = "poret_se",
                     pval_col="poret_p",
                     eaf_col = "poret_eaf",
                     effect_allele_col = "poret_ea",
                     other_allele_col = "poret_nea", 
                     phenotype_col =  "phenotype")
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates
out_ADAMTS13 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/ADAMTS13.csv', 
                                  sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_ADAMTS13$phenotype <- 'ADAMTS13'
out_APOBR <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/APOBR.csv', 
                               sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_APOBR$phenotype <- 'APOBR'
out_EBI3_IL27 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/EBI3_IL27.csv', 
                                   sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_EBI3_IL27$phenotype <- 'EBI3_IL27'
out_FES <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FES.csv', 
                             sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_FES$phenotype <- 'FES'
out_FGF5 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FGF5.csv', 
                              sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_FGF5$phenotype <- 'FGF5'
out_FGL1 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/FGL1.csv', 
                              sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_FGL1$phenotype <- 'FGL1'
out_GDF15 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/GDF15.csv', 
                               sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_GDF15$phenotype <- 'GDF15'
out_PZP <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/PZP.csv', 
                             sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_PZP$phenotype <- 'PZP'
out_SERPINE2 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SERPINE2.csv', 
                                  sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_SERPINE2$phenotype <- 'SERPINE2'
out_SH2B3 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SH2B3.csv', 
                               sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_SH2B3$phenotype <- 'SH2B3'
out_SULT1A1 <- read_outcome_data(snps = preec_ivs$SNP,  filename = '/Volumes/MADDY2/datasets/preeclampsia/outsunukb_reverse/SULT1A1.csv', 
                                 sep = ",",  snp_col = "SNP", beta_col = "BETA", se_col = "SE", pval_col="LOG10P", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0", log_pval=TRUE)
out_SULT1A1$phenotype <- 'SULT1A1'

outlist <- ls(pattern = "out_", mget(ls())) 

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(outlist)], SIMPLIFY = FALSE)
names(join_list) <- gsub('out_', 'har_', outlist)
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 

# Clump locally 
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
} # UPDATED local clump - includes error handler to return null for unclumpables and continue running
iv_list <- sapply(genelist, try_clp, simplify = FALSE)
names(iv_list) <- gsub("^har_","clp_",names(iv_list)) 
clplist <- gsub("^har_","clp_",names(iv_list)) 
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # remove empties
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
rm(list=ls()[!(ls() %in% clplist)]) # keep only clumped files
clplist<-ls(pattern = "clp_", mget(ls())) # Make list with names of all cluped sets
length(clplist) # XXX????? After clump

unlink("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_reverse_cis_preec")
dir.create("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_reverse_cis_preec")
setwd("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_reverse_cis_preec")

files <- mget(ls(pattern = '^clp_')) 
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

#mr & save results (use mr(dat) as this keeps only mr.keep==T)
setwd("/Volumes/MADDY2/datasets/preeclampsia/sunukb_harm_clump_reverse_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]  #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(pattern="*.csv", full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)],
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
clplist<-ls(pattern = "clp_", mget(ls()))

mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  # dat$mr_keep <- TRUE
  res<- mr(dat)
}
mr_table2 <- lapply(clplist, mrfunc2)
names(mr_table2) <- gsub("clp_","",clplist) 
names(mr_table2) <- str_c(names(mr_table2), '_res')
invisible(lapply(names(mr_table2), function(x) assign(x,mr_table2[[x]],envir=.GlobalEnv)))

reslist <- names(mr_table2)

rm(list=ls()[!(ls() %in% reslist)]) 
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<1)) # or ncol<3?
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)
rm(to.rm)
reslist<-ls(pattern = "_res", mget(ls()))

unlink("/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/reverse_results_cis_preec")
dir.create('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/reverse_results_cis_preec')
setwd('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/reverse_results_cis_preec')

files <- mget(reslist)
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}

rm(list=ls())

files <- list.files('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/reverse_results_cis_preec', 
                    pattern = ".csv", recursive = TRUE, full.names = TRUE)
alldat <- lapply(setNames(nm = files), read.csv)
mergedres <- do.call(rbind, Map(function(x, nm) transform(x, filename = nm), alldat, names(alldat)))

mergedres$or <- exp(mergedres$b)
mergedres$lci <- exp(mergedres$b - 1.96*mergedres$se)
mergedres$uci <- exp(mergedres$b + 1.96*mergedres$se)

mergedres <- filter(mergedres, mergedres$method == 'Inverse variance weighted')
mergedres$padj <- p.adjust(mergedres$pval, method = 'fdr', n = nrow(mergedres))

mergedres$outcome <- 'Pre-eclampsia'
mergedres$exposure <- gsub('_res.csv', '', gsub('/Volumes/MADDY2/datasets/preeclampsia/ressunukb_preec/reverse_results_cis_preec/', '', mergedres$filename))
mergedres$filename <- NULL
mergedres[,1:3] <- NULL
mergedres$analysis <- 'Reverse direction'
write.csv(mergedres, "~/Desktop/preec/ressunukb_preec/sens_reverse_cis_preec.csv", row.names = FALSE)


rm(list=ls())





#### --------------------------------------------------------------------------------------------- ####
#### -------------------------------------- COLOC ------------------------------------------------ ####
#### --------------------------------------------------------------------------------------------- ####
#### ------------ DECODE --------------- ####
#### Format exposure data files - extract cis- (+- 500kb) ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec")

# 3MG 16:77,007-85,851
gunzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt')
poslow <- 77007-500000
poshigh <- 85851+500000
chrom <- 'chr16'
p_3MG_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_3MG_ivs<-filter(p_3MG_ivs, p_3MG_ivs$Pos>poslow & p_3MG_ivs$Pos<poshigh & p_3MG_ivs$Chrom == chrom)
write.csv(p_3MG_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_3MG_ivs.csv')
gzip('/volumes/maddy2/decode/12438_127_MPG_3MG.txt', '/volumes/maddy2/decode/12438_127_MPG_3MG.txt.gz')

# ADAMTS13 9:133414358 - 133459402
gunzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt')
poslow <- 133414358-500000
poshigh <- 133459402+500000
chrom <- 'chr9'
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$Pos>poslow & p_ADAMTS13_ivs$Pos<poshigh & p_ADAMTS13_ivs$Chrom == chrom)
write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_ADAMTS13_ivs.csv')
gzip('/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt', '/volumes/maddy2/decode/3175_51_ADAMTS13_ATS13.txt.gz')


# ALDH2 12:111766887 - 111817532
gunzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt')
poslow <- 111766887-500000
poshigh <- 111817532+500000
chrom <- 'chr12'
p_ALDHE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ALDHE2_ivs<-filter(p_ALDHE2_ivs, p_ALDHE2_ivs$Pos>poslow & p_ALDHE2_ivs$Pos<poshigh & p_ALDHE2_ivs$Chrom == chrom)
write.csv(p_ALDHE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_ALDHE2_ivs.csv')
gzip('/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt', '/volumes/maddy2/decode/18381_16_ALDH2_ALDH_E2.txt.gz')

# ANP 1:11845709 - 11848345
gunzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt')
poslow <- 11845709-500000
poshigh <- 11848345+500000
chrom <- 'chr1'
p_ANP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ANP_ivs<-filter(p_ANP_ivs, p_ANP_ivs$Pos>poslow & p_ANP_ivs$Pos<poshigh & p_ANP_ivs$Chrom == chrom)
write.csv(p_ANP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_ANP_ivs.csv')
gzip('/volumes/maddy2/decode/5443_62_NPPA_ANP.txt', '/volumes/maddy2/decode/5443_62_NPPA_ANP.txt.gz')

# FGL1 8:17864380 - 17910365
gunzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt')
poslow <- 17864380-500000
poshigh <- 17910365+500000
chrom <- 'chr8'
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$Pos>poslow & p_FGL1_ivs$Pos<poshigh & p_FGL1_ivs$Chrom == chrom)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_FGL1_ivs.csv')
gzip('/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt', '/volumes/maddy2/decode/5581_28_FGL1_FGL1.txt.gz')

# GDF15 19:18374731 - 18389176
gunzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt')
poslow <- 18374731-500000
poshigh <- 18389176+500000
chrom <- 'chr19'
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$Pos>poslow & p_GDF15_ivs$Pos<poshigh & p_GDF15_ivs$Chrom == chrom)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_GDF15_ivs.csv')
gzip('/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt', '/volumes/maddy2/decode/4374_45_GDF15_MIC_1.txt.gz')

# JUND 19:18279694 - 18281622
gunzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt')
poslow <- 18279694-500000
poshigh <- 18281622+500000
chrom <- 'chr19'
p_JUND_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_JUND_ivs<-filter(p_JUND_ivs, p_JUND_ivs$Pos>poslow & p_JUND_ivs$Pos<poshigh & p_JUND_ivs$Chrom == chrom)
write.csv(p_JUND_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_JUND_ivs.csv')
gzip('/volumes/maddy2/decode/19602_36_JUND_jun_D.txt', '/volumes/maddy2/decode/19602_36_JUND_jun_D.txt.gz')

# MANEA 6:95577485 - 95609470
gunzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt')
poslow <-  95577485-500000
poshigh <- 95609470+500000
chrom <- 'chr6'
p_MANEA_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_MANEA_ivs<-filter(p_MANEA_ivs, p_MANEA_ivs$Pos>poslow & p_MANEA_ivs$Pos<poshigh & p_MANEA_ivs$Chrom == chrom)
write.csv(p_MANEA_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_MANEA_ivs.csv')
gzip('/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt', '/volumes/maddy2/decode/8014_359_MANEA_MANEA.txt.gz')

# METAP1 4:98995659 - 99062809 
gunzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt')
poslow <- 98995659-500000
poshigh <- 99062809+500000
chrom <- 'chr4'
p_METAP1_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_METAP1_ivs<-filter(p_METAP1_ivs, p_METAP1_ivs$Pos>poslow & p_METAP1_ivs$Pos<poshigh & p_METAP1_ivs$Chrom == chrom)
write.csv(p_METAP1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_METAP1_ivs.csv')
gzip('/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt', '/volumes/maddy2/decode/3210_1_METAP1_METAP1.txt.gz')

# NOTUM 17:81952507 - 81961840
gunzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt')
poslow <- 81952507-500000
poshigh <- 81961840+500000
chrom <- 'chr17'
p_NOTUM_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_NOTUM_ivs<-filter(p_NOTUM_ivs, p_NOTUM_ivs$Pos>poslow & p_NOTUM_ivs$Pos<poshigh & p_NOTUM_ivs$Chrom == chrom)
write.csv(p_NOTUM_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_NOTUM_ivs.csv')
gzip('/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt', '/volumes/maddy2/decode/8252_2_NOTUM_NOTUM.txt.gz')

# PZP 12:9148840 - 9208395
gunzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt')
poslow <- 9148840-500000
poshigh <- 9208395+500000
chrom <- 'chr12'
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$Pos>poslow & p_PZP_ivs$Pos<poshigh & p_PZP_ivs$Chrom == chrom)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_PZP_ivs.csv')
gzip('/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt', '/volumes/maddy2/decode/6580_29_PZP_Pregnancy_zone_protein.txt.gz')

# RGS18 1:192158462 - 192185815 
gunzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt')
poslow <- 192158462-500000
poshigh <- 192185815+500000
chrom <- 'chr1'
p_RGS18_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_RGS18_ivs<-filter(p_RGS18_ivs, p_RGS18_ivs$Pos>poslow & p_RGS18_ivs$Pos<poshigh & p_RGS18_ivs$Chrom == chrom)
write.csv(p_RGS18_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_RGS18_ivs.csv')
gzip('/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt', '/volumes/maddy2/decode/13982_33_RGS18_RGS18.txt.gz')

# SERPINE2 2:223975045 - 224039318
gunzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt')
poslow <- 223975045-500000
poshigh <- 224039318+500000
chrom <- 'chr2'
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', separator = "\t")$
  filter(pl$col("Pval") < 1.8e-9)$
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$Pos>poslow & p_SERPINE2_ivs$Pos<poshigh & p_SERPINE2_ivs$Chrom == chrom)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec/p_SERPINE2_ivs.csv')
gzip('/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt', '/volumes/maddy2/decode/19154_41_SERPINE2_Protease_nexin_I.txt.gz')
rm(list=ls())

#### coloc ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivsdecode_coloc_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "rsids",
                    beta_col = "Beta",
                    se_col = "SE",
                    pval_col="minus_log10_pval",
                    eaf_col = "ImpMAF",
                    effect_allele_col = "effectAllele",
                    other_allele_col = "otherAllele", 
                    phenotype_col =  "phenotype", log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outdecode_1e4/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) 

#do coloc
colfunc<-function(d) {
  d <- get(d, envir = .GlobalEnv)
  (coloc.abf((list(type = "quant", beta = d$beta.exposure, varbeta = d$se.exposure^2, 
                   N = 35559, sdY = 1,  snp = d$SNP )), 
             (list( type = "cc", beta = d$beta.outcome, varbeta = d$se.outcome^2, 
                    N = 611484,  s = 0.03076104, snp = d$SNP )),
             MAF = NULL,  p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)$summary)
}
table1 <- data.frame(t(sapply(genelist,colfunc, simplify = TRUE)))
table1<-tibble::rownames_to_column(table1, var = "exposure")
table1$exposure <- gsub('_ivs', '', gsub('har_p_', '', table1$exposure))
table1 <- apply(table1,2,as.character)
write.csv(table1, '~/Desktop/preec/resdecode_preec/decode_coloc_sighits.csv')

rm(list=ls())




#### --------------- UKB --------------- ####
#### Format exposure data files - extract cis- (+- 500kb) ####

dir.create("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec")
dir.create('/volumes/maddy2/sunukb/temptar')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ADAMTS13 9:133414358 - 133459402
untar('/volumes/maddy2/sunukb/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt')
poslow <- 133414358-500000
poshigh <- 133459402+500000
chrom <- 9
p_ADAMTS13_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic/discovery_chr9_ADAMTS13:Q76LX8:OID20249:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_ADAMTS13_ivs<-filter(p_ADAMTS13_ivs, p_ADAMTS13_ivs$GENPOS>poslow & p_ADAMTS13_ivs$GENPOS<poshigh & p_ADAMTS13_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_ADAMTS13_ivs <- cbind(p_ADAMTS13_ivs, data.frame(do.call('rbind', strsplit(as.character(p_ADAMTS13_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_ADAMTS13_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '9')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_ADAMTS13_ivs <- merge(p_ADAMTS13_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_ADAMTS13_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_ADAMTS13_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/ADAMTS13_Q76LX8_OID20249_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# APOBR 16:28494643 - 28498964
untar('/volumes/maddy2/sunukb/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt')
poslow <- 28494643-500000
poshigh <- 28498964+500000
chrom <- 16
p_APOBR_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II/discovery_chr16_APOBR:Q0VD83:OID30246:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_APOBR_ivs<-filter(p_APOBR_ivs, p_APOBR_ivs$GENPOS>poslow & p_APOBR_ivs$GENPOS<poshigh & p_APOBR_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_APOBR_ivs <- cbind(p_APOBR_ivs, data.frame(do.call('rbind', strsplit(as.character(p_APOBR_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_APOBR_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_APOBR_ivs <- merge(p_APOBR_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_APOBR_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_APOBR_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/APOBR_Q0VD83_OID30246_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# FES 15:90883695 - 90895776
untar('/volumes/maddy2/sunukb/FES_P07332_OID21207_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt')
poslow <- 90883695-500000
poshigh <- 90895776+500000
chrom <- 15
p_FES_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology/discovery_chr15_FES:P07332:OID21207:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FES_ivs<-filter(p_FES_ivs, p_FES_ivs$GENPOS>poslow & p_FES_ivs$GENPOS<poshigh & p_FES_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FES_ivs <- cbind(p_FES_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FES_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FES_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '15')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FES_ivs <- merge(p_FES_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FES_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_FES_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FES_P07332_OID21207_v1_Oncology', recursive = TRUE)
rm(list=ls())

# FGF5 4:80266639 - 80336680
untar('/volumes/maddy2/sunukb/FGF5_P12034_OID20490_v1_Inflammation.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.gz', 
       '/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt')
poslow <- 80266639-500000
poshigh <- 80336680+500000
chrom <- 4
p_FGF5_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation/discovery_chr4_FGF5:P12034:OID20490:v1:Inflammation.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGF5_ivs<-filter(p_FGF5_ivs, p_FGF5_ivs$GENPOS>poslow & p_FGF5_ivs$GENPOS<poshigh & p_FGF5_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGF5_ivs <- cbind(p_FGF5_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGF5_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGF5_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '4')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGF5_ivs <- merge(p_FGF5_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGF5_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_FGF5_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGF5_P12034_OID20490_v1_Inflammation', recursive = TRUE)
rm(list=ls())

# FGL1 8:17864380 - 17910365
untar('/volumes/maddy2/sunukb/FGL1_Q08830_OID30702_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt')
poslow <- 17864380-500000
poshigh <- 17910365+500000
chrom <- 8
p_FGL1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II/discovery_chr8_FGL1:Q08830:OID30702:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_FGL1_ivs<-filter(p_FGL1_ivs, p_FGL1_ivs$GENPOS>poslow & p_FGL1_ivs$GENPOS<poshigh & p_FGL1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_FGL1_ivs <- cbind(p_FGL1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_FGL1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_FGL1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '8')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_FGL1_ivs <- merge(p_FGL1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_FGL1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_FGL1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/FGL1_Q08830_OID30702_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())
x
# GDF15 19:18374731 - 18389176
untar('/volumes/maddy2/sunukb/GDF15_Q99988_OID20251_v1_Cardiometabolic.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.gz', 
       '/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt')
poslow <- 18374731-500000
poshigh <- 18389176+500000
chrom <- 19
p_GDF15_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic/discovery_chr19_GDF15:Q99988:OID20251:v1:Cardiometabolic.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_GDF15_ivs<-filter(p_GDF15_ivs, p_GDF15_ivs$GENPOS>poslow & p_GDF15_ivs$GENPOS<poshigh & p_GDF15_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_GDF15_ivs <- cbind(p_GDF15_ivs, data.frame(do.call('rbind', strsplit(as.character(p_GDF15_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_GDF15_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '19')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_GDF15_ivs <- merge(p_GDF15_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_GDF15_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_GDF15_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/GDF15_Q99988_OID20251_v1_Cardiometabolic', recursive = TRUE)
rm(list=ls())

# PZP 12:9148840 - 9208395
untar('/volumes/maddy2/sunukb/PZP_P20742_OID30730_v1_Inflammation_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.gz', 
       '/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt')
poslow <- 9148840-500000
poshigh <- 9208395+500000
chrom <- 12
p_PZP_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II/discovery_chr12_PZP:P20742:OID30730:v1:Inflammation_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_PZP_ivs<-filter(p_PZP_ivs, p_PZP_ivs$GENPOS>poslow & p_PZP_ivs$GENPOS<poshigh & p_PZP_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_PZP_ivs <- cbind(p_PZP_ivs, data.frame(do.call('rbind', strsplit(as.character(p_PZP_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_PZP_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_PZP_ivs <- merge(p_PZP_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_PZP_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_PZP_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/PZP_P20742_OID30730_v1_Inflammation_II', recursive = TRUE)
rm(list=ls())

# SERPINE2 2:223975045 - 224039318
untar('/volumes/maddy2/sunukb/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.gz', 
       '/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt')
poslow <- 223975045-500000
poshigh <- 224039318+500000
chrom <- 2
p_SERPINE2_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II/discovery_chr2_SERPINE2:P07093:OID30359:v1:Cardiometabolic_II.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SERPINE2_ivs<-filter(p_SERPINE2_ivs, p_SERPINE2_ivs$GENPOS>poslow & p_SERPINE2_ivs$GENPOS<poshigh & p_SERPINE2_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SERPINE2_ivs <- cbind(p_SERPINE2_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SERPINE2_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SERPINE2_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '2')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SERPINE2_ivs <- merge(p_SERPINE2_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SERPINE2_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_SERPINE2_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SERPINE2_P07093_OID30359_v1_Cardiometabolic_II', recursive = TRUE)
rm(list=ls())

# SH2B3 12:111405923 - 111451623
untar('/volumes/maddy2/sunukb/SH2B3_Q9UQQ2_OID21222_v1_Oncology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.gz', 
       '/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt')
poslow <- 111405923-500000
poshigh <- 111451623+500000
chrom <- 12
p_SH2B3_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology/discovery_chr12_SH2B3:Q9UQQ2:OID21222:v1:Oncology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SH2B3_ivs<-filter(p_SH2B3_ivs, p_SH2B3_ivs$GENPOS>poslow & p_SH2B3_ivs$GENPOS<poshigh & p_SH2B3_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SH2B3_ivs <- cbind(p_SH2B3_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SH2B3_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SH2B3_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '12')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SH2B3_ivs <- merge(p_SH2B3_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SH2B3_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_SH2B3_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SH2B3_Q9UQQ2_OID21222_v1_Oncology', recursive = TRUE)
rm(list=ls())

# SULT1A1 16:28605196 - 28614279
untar('/volumes/maddy2/sunukb/SULT1A1_P50225_OID21031_v1_Neurology.tar', 
      exdir = '/volumes/maddy2/sunukb/temptar')
gunzip('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.gz', 
       '/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt')
poslow <- 28605196-500000
poshigh <- 28614279+500000
chrom <- 16
p_SULT1A1_ivs <- pl$
  scan_csv('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology/discovery_chr16_SULT1A1:P50225:OID21031:v1:Neurology.txt', separator = " ")$
  filter(pl$col("LOG10P") > 10.7695511)$ # -log10(1.7e-11)
  collect()$
  to_data_frame() %>%
  as_tibble()
p_SULT1A1_ivs<-filter(p_SULT1A1_ivs, p_SULT1A1_ivs$GENPOS>poslow & p_SULT1A1_ivs$GENPOS<poshigh & p_SULT1A1_ivs$CHROM == chrom)

# Use biomart to get SNP names
p_SULT1A1_ivs <- cbind(p_SULT1A1_ivs, data.frame(do.call('rbind', strsplit(as.character(p_SULT1A1_ivs$ID),':',fixed=TRUE)))[,2])
colnames(p_SULT1A1_ivs)[15] <- 'pos37'
all_snps_info = snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, '16')
all_snps_info_DF <- data.frame(all_snps_info)[,c('pos', 'RefSNP_id')]
colnames(all_snps_info_DF) <- c('pos37', 'SNP')
p_SULT1A1_ivs <- merge(p_SULT1A1_ivs, all_snps_info_DF, by = 'pos37', all.x=FALSE, all.y=FALSE)
write.csv(p_SULT1A1_ivs, '/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec/p_SULT1A1_ivs.csv')
unlink('/volumes/maddy2/sunukb/temptar/SULT1A1_P50225_OID21031_v1_Neurology', recursive = TRUE)
rm(list=ls())




#### coloc ####
# Import all IVs
setwd("/Volumes/MADDY2/datasets/preeclampsia/ivssunukb_coloc_unadj_cis_preec")
files = list.files(pattern="*.csv")[which(file.info(list.files(pattern="*.csv"))$size>3)]   #make list of all csv names
files <- files[which(file.info(files)$size>3)]
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
data_list <- lapply(names(data_list),
                    function(current_name)
                      transform(data_list[[current_name]],
                                new_column = current_name))
names(data_list) <- gsub(".csv","", list.files(pattern="*.csv", 
                                               full.names=FALSE)[which(file.info(list.files(pattern="*.csv", full.names=FALSE))$size>3)], fixed = TRUE)
list2env(data_list, envir = .GlobalEnv)
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
rm(files, data_list)
setwd("~/Desktop/preec")

# Format IVs as exposure data
formfunc<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <-format_data(dat,type="exposure",snp_col = "SNP",
                    beta_col = "BETA",
                    se_col = "SE",
                    pval_col="LOG10P",
                    eaf_col = "A1FREQ",
                    effect_allele_col = "ALLELE1",
                    other_allele_col = "ALLELE0", 
                    phenotype_col =  "new_column",  # CHECK THIS IS THE EXPOSURE NAME
                    log_pval = TRUE)
}
join_list<-lapply(genelist, formfunc)
names(join_list) <- genelist
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(formfunc, join_list)

# Import outcome association estimates
outex<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  dat <- read_outcome_data(snps = dat$SNP, 
                           filename = "~/Desktop/preec/outsunukb_1e4/preec_out_unadj.csv", 
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
outlist <- str_c("out_",genelist) 
rm(outex, join_list)

# Harmonise data
hrm<-function(exp, out){
  exp <- get(exp, envir = .GlobalEnv)
  out <- get(out, envir = .GlobalEnv)
  dat <- harmonise_data(exp, out, action=2)
}
join_list<- mapply(hrm, genelist[1:length(genelist)], outlist[1:length(genelist)], SIMPLIFY = FALSE)
names(join_list) <- str_c("har_",names(join_list) ) 
invisible(lapply(names(join_list), function(x) assign(x,join_list[[x]],envir=.GlobalEnv)))
rm(hrm, join_list)
genelist<-ls(pattern = "har_", mget(ls())) 
rm(list=ls()[!(ls() %in% genelist)]) 
genelist<-ls(pattern = "har_", mget(ls())) 

#do coloc
colfunc<-function(d) {
  d <- get(d, envir = .GlobalEnv)
  (coloc.abf((list(type = "quant", beta = d$beta.exposure, varbeta = d$se.exposure^2, 
                   N = 34557, sdY = 1,  snp = d$SNP )), 
             (list( type = "cc", beta = d$beta.outcome, varbeta = d$se.outcome^2, 
                    N = 611484,  s =  0.03076104, snp = d$SNP )),
             MAF = NULL,  p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)$summary)
}
table1 <- data.frame(t(sapply(genelist,colfunc, simplify = TRUE)))
table1<-tibble::rownames_to_column(table1, var = "exposure")
table1$exposure <- gsub('_ivs', '', gsub('har_p_', '', table1$exposure))
table1 <- apply(table1,2,as.character)
write.csv(table1, '~/Desktop/preec/ressunukb_preec/sunukb_coloc_sighits.csv')

rm(list=ls())





#### -------------------------------------- protein annotation ------------------------------------------------ ####
# Load required libraries
library(readxl)
library(dplyr)

# Read the Excel file, skipping the first row and using the second row as column names
data <- readxl::read_excel("~/desktop/preec/manuscript/Tables3.xlsx", sheet = "SuppTable11", skip = 1)

# Tally counts of values in the Description column overall
overall_counts <- data %>%
  group_by(Description) %>%
  summarize(count = n())

# Tally counts of values in the Description column split by Category
category_counts <- data %>%
  group_by(Category, Description) %>%
  summarize(count = n())

# Print overall counts
print("Overall Counts:")
print(overall_counts)

# Print counts split by Category
print("Counts by Category:")
print(category_counts)


# Select and save top 10 counts by category
top_10_by_category <- category_counts %>%
  group_by(Category) %>%
  top_n(10, wt = count) %>%
  arrange(Category, desc(count)) %>%
  ungroup()

# Save the top 10 counts by category to a new data frame
top_10_counts <- data.frame(top_10_by_category)

# Print the top 10 counts by category
print("Top 10 Counts by Category:")
print(top_10_counts)