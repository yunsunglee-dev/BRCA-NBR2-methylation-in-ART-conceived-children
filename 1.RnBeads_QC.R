#####################################################
## What  : MoBa met013 DNA methylation QC; attempt 2
## Who   : Christian M Page (CMP)
## Where : NIPH, Oslo, Norway
## When  : Summer, 2024
##
## Description: 
## Using RnBeads (instead of Minfi), since this has
## support for EPICv2. 
##
##
##
##
####################################################


## Windows only script!!! 
## QC using RnBeads
library(RnBeads)
library(RnBeads.hg38)

library(sesame)
library(RPMM)

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#debug(read.idat.files)
#debug(rnb.detect.infinium.platform)

#debug(rnb.plot.control.barplot)

## Data loader: 
setwd("N:/durable/Projects/MoBa_genomics/MoBa13/src")


idat.dir 	<- file.path("../iDat/")

batch <- 1

sample.annotation <- file.path(paste("../doc/2024_05_07_PDB2374_Met013_sample_sheet_batch_",batch, ".csv", sep="")) 
analysis.dir 	<- file.path("../res")#Sys.getenv("resFolder") # Defines a directory where the data will be saved
report.dir 	<- file.path(analysis.dir,paste(format(Sys.Date(), "%Y%m%d"),"RnBeads", "batch", batch,sep="_")) # Difines a folder within analysis.dir where the results are saved


#Setting up analysis
# ?rnb.options() gives quite alot of information regarding the variables to be defined.
rnb.options(
  import.idat.platform     = "auto", ## probesEPICv2
  analysis.name            = "MoBa13",          ## How can we annotate if there is multiple data (ie. some cordBlood and som PBMC?
  #inference.reference.methylome.column	= "CellType",
  identifiers.column       = "SAMPLE_NAME", #"barcode",        ## Depending on the sample sheet
  filtering.blacklist      = "../dat/blacklist_crossHybridizing_probes_UpdatedCpGandNonCpH.txt", 
  filtering.cross.reactive = TRUE, ## Should take out any from the balcklist
  filtering.snp            = "no",
  import.sex.prediction    = TRUE, 
  qc.snp.heatmap	         = TRUE, 
  qc.snp.distances	       = TRUE,
  qc.snp.barplot	         = TRUE,
  filtering.sex.chromosomes.removal     = FALSE,  ## Should we remove the X & Y chomosome? Might want to look into X-Chromsome? 
  filtering.greedycut.pvalue.threshold  = 0.01,   ## Specific value is probably ok, but sould be discussed
  qc.sample.batch.size      = 61,                 ## changing batch size for plots from 500; should not be devidable by the number of samples!!
  filtering.context.removal = NULL,               ## Do not remove any CpH based no context alone
  normalization.background.method       = "sesame.noob",  # default method: methylumi.noob
  normalization.method      = "bmiq", 
  normalization             = TRUE,
  filtering.greedycut       = TRUE,
  export.to.trackhub        = NULL,
  filtering.missing.value.quantile  = 1,     # retains all CpG methylation, does not remove based on NAs across samples
  filtering.coverage.threshold      = 5,      # minimum coverage for each CpG
  export.to.bed             = FALSE,
  qc.boxplots               = TRUE, ## Set to false if QC keeps failing
  qc.barplots               = TRUE, ## Set to false if QC keeps failing
  qc.negative.boxplot       = TRUE, ## Set to false if QC keeps failing
  disk.dump.big.matrices    = FALSE,
  disk.dump.bigff           = FALSE,
  exploratory               = FALSE,
  exploratory.beta.distribution    = FALSE, 
  exploratory.clustering    = "none",
  logging.exit.on.error	    = TRUE,
  differential              = FALSE)



rnb.run.analysis(
  dir.reports = report.dir, 
  sample.sheet= sample.annotation, 
  data.dir    = idat.dir, 
  data.type   = "idat.dir")


#as.matrix(meth(load.rnb.set("../res/20240715_RnBeads_batch_1/rnbSet_preprocessed/"), row.names = TRUE))
bmat_1 <- meth(load.rnb.set(file.path(report.dir,"rnbSet_preprocessed")))
save(bmat_1, file = "../bmat/bmat_1.RData")


batch <- 2
gc()
sample.annotation <- file.path(paste("../doc/2024_05_07_PDB2374_Met013_sample_sheet_batch_",batch, ".csv", sep="")) 
report.dir 	<- file.path(analysis.dir,paste(format(Sys.Date(), "%Y%m%d"),"RnBeads", "batch", batch,sep="_")) # Difines a folder within analysis.dir where the results are saved

rnb.run.analysis(
  dir.reports = report.dir, 
  sample.sheet= sample.annotation, 
  data.dir    = idat.dir, 
  data.type   = "idat.dir")


bmat_2 <- meth(load.rnb.set(file.path(report.dir,"rnbSet_preprocessed")))
save(bmat_2, file = "../bmat/bmat_2.RData")


## Merge the two data sets
load("../annot/annot_EPIC_v2-A2.RData")
cg   <- Reduce(intersect, list(rownames(bmat_1), rownames(bmat_2)))
cg.auto <- intersect(cg, subset(annot.epic8.v2_A2$IlmnID, annot.epic8.v2_A2$CHR %in% 1:22))

rownames(annot.epic8.v2_A2) <- annot.epic8.v2_A2$IlmnID

## May need fix for collum names! 
bmat <- cbind(bmat_1[cg,],bmat_2[cg,])


key  <- read.table(file.path("../doc/2024_05_07_PDB2374_Met013_sample_sheet.csv", sep=""), sep = ",", header=TRUE, stringsAsFactors = FALSE)

key_ubc <- subset(key, key$AGE_DATESAMPLE < 10 & key$SAMPLE_NAME %in% colnames(bmat))
key_wbc <- subset(key, key$AGE_DATESAMPLE > 10 & key$SAMPLE_NAME %in% colnames(bmat))


bmat.child.umbilical.autosomal.BMIQ.NewStart <- bmat[cg.auto, key_ubc$SAMPLE_NAME]
bmat.child.wb.autosomal.BMIQ.NewStart <- bmat[cg.auto, key_wbc$SAMPLE_NAME]


rownames(bmat.child.umbilical.autosomal.BMIQ.NewStart) <- annot.epic8.v2_A2[rownames(bmat.child.umbilical.autosomal.BMIQ.NewStart),"Name"]
rownames(bmat.child.wb.autosomal.BMIQ.NewStart) <- annot.epic8.v2_A2[rownames(bmat.child.wb.autosomal.BMIQ.NewStart),"Name"]


cc.ubc <- setNames(paste("PREG_ID_2374",key_ubc$PREG_ID_2374, key_ubc$BARN_NR, sep="_"), key_ubc$SAMPLE_NAME)
cc.wb <- setNames(paste("PREG_ID_2374",key_wbc$PREG_ID_2374,  key_wbc$BARN_NR, sep="_"), key_wbc$SAMPLE_NAME)

bmat.child.umbilical.autosomal.BMIQ.NewStart <- bmat.child.umbilical.autosomal.BMIQ.NewStart[,names(cc.ubc)]
colnames(bmat.child.umbilical.autosomal.BMIQ.NewStart) <- cc.ubc

bmat.child.wb.autosomal.BMIQ.NewStart <- bmat.child.wb.autosomal.BMIQ.NewStart[,names(cc.wb)]
colnames(bmat.child.wb.autosomal.BMIQ.NewStart) <- cc.wb

bmat.child.umbilical.autosomal.BMIQ.NewStart[1:10,1:5]
bmat.child.wb.autosomal.BMIQ.NewStart[1:5,1:5]

save(bmat.child.umbilical.autosomal.BMIQ.NewStart, file = "../met_13/Clean_after_QC/Autosomal/noRep/bmat_child_umbilical_autosomal_BMIQ_NewStart.RData")
save(bmat.child.wb.autosomal.BMIQ.NewStart, file = "../met_13/Clean_after_QC/Autosomal/noRep/bmat_child_WB_autosomal_BMIQ_NewStart.RData")


## Key
key.NewStart <- rbind(rnb_b1@pheno, rnb_b2@pheno)
key.umbilical.NewStart <- subset(key.NewStart, key.NewStart$AGE_DATESAMPLE < 10)
key.wb.NewStart <- subset(key.NewStart, key.NewStart$AGE_DATESAMPLE > 10)

rownames(key.umbilical.NewStart) <- paste("PREG_ID_2374", key.umbilical.NewStart$PREG_ID_2374, key.umbilical.NewStart$BARN_NR, sep="_")
rownames(key.wb.NewStart) <- paste("PREG_ID_2374", key.wb.NewStart$PREG_ID_2374, key.wb.NewStart$BARN_NR, sep="_")
save(key.umbilical.NewStart, key.wb.NewStart, file = "../met_13/Clean_after_QC/key/PDB2374_PregID_Key_Technical_sheet_NewStart.RData")


## Annot
annot.epic8.v2_A2.orig <- annot.epic8.v2_A2
annot.epic8.v2_A2 <- annot.epic8.v2_A2[cg,c("IlmnID", "Name", "Infinium_Design_Type", "CHR", "MAPINFO",
                                            "Color_Channel", "Probe_Type", "UCSC_RefGene_Group",
                                            "UCSC_RefGene_Name", "UCSC_RefGene_Accession","MAPINFO_37",
                                            "Strand_FR","Strand_TB","Strand_CO")]
rownames(annot.epic8.v2_A2) <- annot.epic8.v2_A2$Name
save(annot.epic8.v2_A2, file = "../met_13/Clean_after_QC/key/annot_EPIC_v2_A2.RData")


cg.X <- intersect(cg, subset(annot.epic8.v2_A2$IlmnID, annot.epic8.v2_A2$CHR == "X"))

bmat.child.umbilical.chrX.BMIQ.NewStart <- bmat[cg.X, names(cc.ubc)]
colnames(bmat.child.umbilical.chrX.BMIQ.NewStart) <- cc.ubc

bmat.child.wb.chrX.BMIQ.NewStart <- bmat[cg.X, names(cc.wb)]
colnames(bmat.child.wb.chrX.BMIQ.NewStart) <- cc.wb

save(bmat.child.umbilical.chrX.BMIQ.NewStart, file = "../met_13/Clean_after_QC/ChrX/bmat_child_umbilical_chrX_BMIQ_NewStart.RData")
save(bmat.child.wb.chrX.BMIQ.NewStart, file = "../met_13/Clean_after_QC/ChrX/bmat_child_wb_chrX_BMIQ_NewStart.RData")

## Y chromosome
cg.Y <- intersect(cg, subset(annot.epic8.v2_A2$IlmnID, annot.epic8.v2_A2$CHR == "Y"))
cc.ubc.y <- subset(setNames(paste("PREG_ID_2374",key_ubc$PREG_ID_2374, key_ubc$BARN_NR, sep="_"), key_ubc$SAMPLE_NAME), key_ubc$SEX == "M")
cc.wb.y  <- subset(setNames(paste("PREG_ID_2374",key_wbc$PREG_ID_2374, key_wbc$BARN_NR, sep="_"), key_wbc$SAMPLE_NAME), key_wbc$SEX == "M")

  
bmat.child.umbilical.chrY.BMIQ.NewStart <- bmat[cg.Y, names(cc.ubc.y)]
colnames(bmat.child.umbilical.chrY.BMIQ.NewStart) <- cc.ubc.y

bmat.child.wb.chrY.BMIQ.NewStart <- bmat[cg.Y, names(cc.wb.y)]
colnames(bmat.child.wb.chrY.BMIQ.NewStart) <- cc.wb.y

save(bmat.child.umbilical.chrY.BMIQ.NewStart, file = "../met_13/Clean_after_QC/ChrY/bmat_child_umbilical_chrY_BMIQ_NewStart.RData")
save(bmat.child.wb.chrY.BMIQ.NewStart, file = "../met_13/Clean_after_QC/ChrY/bmat_child_wb_chrY_BMIQ_NewStart.RData")







