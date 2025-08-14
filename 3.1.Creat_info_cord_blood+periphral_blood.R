source("setwd.R")

### Load bmat
load("/mnt/scratch/yunsung.lee/5.met013/Clean_after_QC/Autosomal/noRep//bmat_child_umbilical_autosomal_BMIQ_NewStart.RData")
load("/mnt/scratch/yunsung.lee/5.met013/Clean_after_QC/Autosomal/noRep//bmat_child_WB_autosomal_BMIQ_NewStart.RData")

### Load the UniLIFE reference panel and apply
require(EpiDISH)
centUniLIFE.m=readRDS("data/centUniLIFE.m.RData")

cell_wb <- epidish(beta.m = bmat.child.wb.autosomal.BMIQ.NewStart, ref.m = centUniLIFE.m[,c(8:19)], method = "RPC")$estF
cell_um <- epidish(beta.m = bmat.child.umbilical.autosomal.BMIQ.NewStart, ref.m = centUniLIFE.m[,c(1:7)], method = "RPC")$estF

table(rownames(bmat.child.umbilical.autosomal.BMIQ.NewStart) == rownames(bmat.child.wb.autosomal.BMIQ.NewStart))
colnames(bmat.child.umbilical.autosomal.BMIQ.NewStart)=paste0(as.character(colnames(bmat.child.umbilical.autosomal.BMIQ.NewStart)),"_um")
colnames(bmat.child.wb.autosomal.BMIQ.NewStart)=paste0(as.character(colnames(bmat.child.wb.autosomal.BMIQ.NewStart)),"_wb")
bmat=cbind(bmat.child.umbilical.autosomal.BMIQ.NewStart,bmat.child.wb.autosomal.BMIQ.NewStart)

### Transform beta values to M values
bmat=log2(bmat/(1-bmat));gc(TRUE)
print(head(bmat[,1:4]));print(dim(bmat))

### EPIC annotation.
anno=data.table::fread("/mnt/work/yunsung.lee/data/MethylationEPIC v2.0 Files/EPIC-8v2-0_A2.csv",data.table=F,skip=7)
anno$UCSC_RefGene_Name2=unlist(lapply(1:nrow(anno),FUN=function(i){
paste0(unique(unlist(strsplit(anno$UCSC_RefGene_Name[i],";"))),collapse = ";")    
}))
anno_sub=subset(anno,anno$CHR_37=="17"& anno$MAPINFO_37 > 41277059 & anno$MAPINFO_37 < 41278712)
# ### Select BRCA1 and NBR2 region
# anno_sub1=anno[grep(pattern = "BRCA1",anno$UCSC_RefGene_Name),]
# anno_sub2=anno[grep(pattern = "NBR2",anno$UCSC_RefGene_Name),]
# anno_sub=unique(rbind(anno_sub1,anno_sub2))
# anno_sub=anno_sub[order(anno_sub$MAPINFO_37),]

common=intersect(anno_sub$Name,rownames(bmat))
anno_sub2=anno_sub[match(common,anno_sub$Name),c("Name","CHR","MAPINFO_37")]
anno_sub2=anno_sub2[order(anno_sub2$MAPINFO_37),]
common=anno_sub2$Name
bmat_sub=t(bmat[match(common, rownames(bmat)),]);gc(TRUE)

### Load info
info_um=readRDS("data/key_umbilical_NewStart.RData");info_um$type="um"
info_wb=readRDS("data/key_wb_NewStart.RData");info_wb$type="wb"

### Add Cell type composition to info
cell_um=cell_um[match(paste0("PREG_ID_2374_",info_um$PREG_ID_2374_r),rownames(cell_um)),]
cell_wb=cell_wb[match(paste0("PREG_ID_2374_",info_wb$PREG_ID_2374_r),rownames(cell_wb)),]
table(info_um$PREG_ID_2374_r == gsub("PREG_ID_2374_","",rownames(cell_um)))
table(info_wb$PREG_ID_2374_r == gsub("PREG_ID_2374_","",rownames(cell_wb)))
info_um=cbind(info_um,cell_um);info_wb=cbind(info_wb,cell_wb)

for(i in setdiff(colnames(info_um),colnames(info_wb))){info_wb[,i]=NA}
for(i in setdiff(colnames(info_wb),colnames(info_um))){info_um[,i]=NA}
table(colnames(info_wb) %in% colnames(info_um))

info=rbind(info_um,info_wb)
info$PREG_ID_2374_r2=paste0("PREG_ID_2374_",info$PREG_ID_2374_r,"_",info$type)
table(info$PREG_ID_2374_r2 %in% colnames(bmat))

### Align info and bmat
common_ID=intersect(info$PREG_ID_2374_r2,rownames(bmat_sub))
info=info[match(common_ID,info$PREG_ID_2374_r2),]
bmat_sub=bmat_sub[match(common_ID,rownames(bmat_sub)),]
print(dim(bmat_sub));gc(TRUE)
print(table(info$PREG_ID_2374_r2==rownames(bmat_sub)))

### Create info
info=cbind(info,bmat_sub)
table(info$ANY_ART,info$type)
table(table(info$PREG_ID_2374))

saveRDS(info,file="data/info_bmat/info_BRCA1.RData")
saveRDS(common,file="data/info_bmat/common_BRCA1.RData")

info=readRDS("data/info_bmat/info_BRCA1.RData")

### Select one twin from twin pairs
tmp=unique(info[,c("PREG_ID_2374","PREG_ID_2374_r","multiple","ANY_ART")])
tmp=data.frame(table(tmp$PREG_ID_2374));table(tmp$Freq)
twin_ids=as.character(subset(tmp,tmp$Freq==2)$Var1)
selector=function(i){
    one=info[info$PREG_ID_2374 == twin_ids[i],c("PREG_ID_2374","PREG_ID_2374_r","PREG_ID_2374_r2","multiple")]
    tmp=data.frame(table(one$PREG_ID_2374_r))
    if(length(unique(tmp$Freq))==2){
        return(as.character(subset(tmp,tmp$Freq==max(tmp$Freq))$Var1))
    }else{
        return(sample(as.character(tmp$Var1),1))
    }
}
set.seed(1);twin_ids_s=unlist(lapply(1:length(twin_ids),selector))
info_sub=subset(info,info$PREG_ID_2374_r %in% c(subset(info,!info$PREG_ID_2374 %in% twin_ids)$PREG_ID_2374_r,twin_ids_s))
tmp=data.frame(table(info_sub$PREG_ID_2374));table(tmp$Freq)
table(info_sub$type,info_sub$ANY_ART)

table(unique(info_sub[,c("PREG_ID_2374_r","PREG_ID_2374","ANY_ART")])$ANY_ART)

saveRDS(info_sub,file="data/info_bmat/info_BRCA1_notwinpairs.RData")