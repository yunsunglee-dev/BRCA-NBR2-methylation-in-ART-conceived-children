source("setwd.R")

### Load info without twin "pairs" (individuals subject to multiple births are still present)
info=readRDS("data/info_bmat/info_BRCA1_notwinpairs.RData")

### Set the ART variable (could be non_fresh, non_frozen, non_ivf..)
exposure="ANY_ART"#commandArgs(TRUE)[1];
cpgs=readRDS("data/info_bmat/common_BRCA1.RData")
length(cpgs)

table(info$type)

### Read cell types available
cells_um=colnames(readRDS("data/centUniLIFE.m.RData"))[1:7]
cells_wb=colnames(readRDS("data/centUniLIFE.m.RData"))[8:19]

### Define the work horse.
work_horse=function(i,um_wb,cell=FALSE, twin_exc=FALSE){
    
    ### Define different adjusting variables
    if(um_wb=="um"){
        adj_vars=c('mat.age','mat.smk','mat.bmi','mat.parity','child.sex','multiple')
    }else{
        adj_vars=c('mat.age','mat.smk','mat.bmi','mat.parity','child.sex','multiple','child.age')
    }      
    ### Cell type adjustment
    if(cell==TRUE & um_wb=="um"){adj_vars=c(adj_vars,cells_um[-1])}else if(cell==TRUE & um_wb=="wb"){adj_vars=c(adj_vars,cells_wb[-1])}
    # print(adj_vars)
    
    ### Load info
    if(twin_exc==TRUE){info=subset(info,info$multiple==0);adj_vars=setdiff(adj_vars,"multiple")}
    info_sub=subset(info,info$type==um_wb);#print(dim(info_sub))
    ### Generate dat
    dat=na.omit(data.frame(y=info_sub[,cpgs[i]],info_sub[,c(exposure,adj_vars,"plate")]))
    mad=median(abs(dat$y-median(dat$y)))
    num_exc=sum(dat$y > (median(dat$y) + mad*5) | dat$y < (median(dat$y) - mad*5))
    #hist(dat$y);
    #abline(v=c(median(dat$y),median(dat$y) + mad*5,median(dat$y) - mad*5),col="red",lty=2)
    dat=subset(dat,dat$y <= (median(dat$y) + mad*5) & dat$y >= (median(dat$y) - mad*5))
    model=try(nlme::lme(as.formula(paste0("y ~ ",paste0(c(exposure,adj_vars),collapse = "+"))),random = ~1|plate,data=dat),silent=T)
    #model=try(MASS::rlm(as.formula(paste0("y ~ ",paste0(c(exposure,adj_vars,"plate"),"*child.sex",collapse = "+"))),
    #                     data=dat),silent=T)
    if(class(model)[1]=="try-error"){
        tab=c(rep(NA,5),as.numeric(table(dat[,exposure])))
    }else{
        #tab=summary(model)$coef
        tab=summary(model)$tTable
        tab=c(as.numeric(tab[which(rownames(tab)==exposure),]), as.numeric(table(dat[,exposure])))
    }
    return(c(i,num_exc,tab))
}


for(twin_exc in c("","_twin")){
    for(cell_adj in c("","_w_cell")){
        for(one in c("um","wb")){    
            ### Run it.
            system.time({ewas=parallel::mclapply(1:length(cpgs),work_horse,um_wb=one,
                                                 cell=ifelse(cell_adj=="",FALSE,TRUE),                                             
                                                 twin_exc=ifelse(twin_exc=="",FALSE,TRUE),   
                                                 mc.cores=1)})
            ### Transform the list to a dataframe.
            print(length(ewas))
            ewas=as.data.frame(do.call("rbind",ewas))
            ewas$V1=cpgs
            colnames(ewas)=c("Name","num_outliers","b","se","df","t","p","n0","n1")
            print(head(ewas))
            saveRDS(ewas,file=paste0("result_data/assoc_",one,cell_adj,twin_exc,".RData"))
            
            #print(table(info$ANY_ART,info$art_iiff,useNA = "always"))
        }
    }
}
