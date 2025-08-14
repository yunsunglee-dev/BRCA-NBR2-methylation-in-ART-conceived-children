source("setwd.R")

info=readRDS("data/info_bmat/info_BRCA1_notwinpairs.RData");dim(info)
common=readRDS("data/info_bmat/common_BRCA1.RData")
length(common)

### Read cell types available
cells_um=colnames(readRDS("data/centUniLIFE.m.RData"))[1:7]
cells_wb=colnames(readRDS("data/centUniLIFE.m.RData"))[8:19]

### When running original analysis
info_sub=info;dim(info_sub)

### When running sensitivity analysis without twin individuals
# print(paste0("# of twin individuals:",length(unique(subset(info,info$multiple==1)$PREG_ID_2374_r))))
# info_sub=subset(info,info$multiple==0);dim(info_sub)

horse=function(cpg,cell=FALSE){
    ### Remove and count outliers and NAs.
    bundle=lapply(c("um","wb"),FUN=function(i){
        info_sub=subset(info_sub,info_sub$type==i)    
        #print(nrow(info_sub))
        ### Remove outliers
        mad=median(abs(info_sub[,cpg]-median(info_sub[,cpg],na.rm=T)),na.rm=T)
        num_exc=sum(info_sub[,cpg] > (median(info_sub[,cpg],na.rm=T) + mad*5) | info_sub[,cpg] < (median(info_sub[,cpg],na.rm=T) - mad*5),na.rm=T)
        num_na=sum(is.na(info_sub[,cpg]))
        info_sub=subset(info_sub,info_sub[,cpg] <= (median(info_sub[,cpg],na.rm=T) + mad*5) & info_sub[,cpg] >= (median(info_sub[,cpg],na.rm=T) - mad*5))
        #print(nrow(info_sub))
        #cat("\n")
        return(list(num_exc,num_na,info_sub))
    })
    num_exc=bundle[[1]][[1]] + bundle[[2]][[1]]
    num_na=bundle[[1]][[2]] + bundle[[2]][[2]]
    info_sub=rbind(bundle[[1]][[3]], bundle[[2]][[3]])
    info_sub[,setdiff(common,cpg)]=c()
    
    ### Split um and wb
    info_um=subset(info_sub,info_sub$type=="um");dim(info_um)
    info_wb=subset(info_sub,info_sub$type=="wb");dim(info_wb)
    
    ### Regress out child.age, child.sex, and plate
    if(cell==TRUE){
        model=nlme::lme(as.formula(paste0(cpg,"~child.sex+",paste0(cells_um[-1],collapse = "+"))),random = ~1|plate, data=info_um,na.action = "na.exclude")
        info_um[,paste0(cpg)]=as.numeric(resid(model))
        model=nlme::lme(as.formula(paste0(cpg,"~child.age + child.sex+",paste0(cells_wb[-1],collapse = "+"))),random = ~1|plate, data=info_wb,na.action = "na.exclude")
        info_wb[,paste0(cpg)]=as.numeric(resid(model))
    }else{
        model=nlme::lme(as.formula(paste0(cpg,"~child.sex")),random = ~1|plate, data=info_um,na.action = "na.exclude")
        info_um[,paste0(cpg)]=as.numeric(resid(model))
        model=nlme::lme(as.formula(paste0(cpg,"~child.age + child.sex")),random = ~1|plate, data=info_wb,na.action = "na.exclude")
        info_wb[,paste0(cpg)]=as.numeric(resid(model))
    }
    ### Flatten data
    info=merge(x=info_um, y=info_wb[,c("PREG_ID_2374_r","child.age",cpg)],by="PREG_ID_2374_r")
    dim(info)
    ### Calculate the age gap
    #info$child.age=info$child.age.y-info$child.age.x
    ### Create interaction term
    #info$ANY_ART_child.age=info$ANY_ART*info$child.age
    
    
    ### Calculate the difference
    info[,paste0(cpg)] = info[,paste0(cpg,".y")] - info[,paste0(cpg,".x")]
    info[,paste0(cpg,".x")]=c();info[,paste0(cpg,".y")]=c()
    
    adj_vars=c('mat.age','mat.smk','mat.bmi','mat.parity','multiple')
    interest=c("ANY_ART")
    
    ### Ordinary linear regression
    model=glm(as.formula(paste0(cpg,"~ ",paste0(interest,collapse = "+"),"+",paste0(adj_vars,collapse = "+"))),                    
                    data=info,na.action="na.exclude")
    tab=summary(model)$coef;#print(tab)
    tab=as.numeric(tab[which(rownames(tab) %in% interest),])
    return(tab)
}

for(twin_inc in c("","_twin")){
    if(twin_inc==""){
        info_sub=info;print(dim(info_sub))
    }else{
        print(paste0("# of twin individuals:",length(unique(subset(info,info$multiple==1)$PREG_ID_2374_r))))
        info_sub=subset(info,info$multiple==0);dim(info_sub)
    }
    for(cell_adj in c("","_w_cell")){    
        result=data.frame(Name=common,do.call("rbind",lapply(common,horse,cell=ifelse(cell_adj=="_w_cell",TRUE,FALSE))))
        colnames(result)=c("Name","b","se","z","p")
        result$q=p.adjust(result$p,method = "BH")
        result$lower=result$b -qnorm(.975)*result$se;result$upper=result$b +qnorm(.975)*result$se
        print(head(result[order(result$p),]))
        saveRDS(result,file=paste0("result_data/assoc_paired_resid",cell_adj,twin_inc,".RData"))
    }
}
