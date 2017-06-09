#get.END.from.vcf = function(x){ 
#    a = as.character(x) 
#    b = strsplit(a,split = ";")
#    c = strsplit(b[[1]][2],split = '=')
#    d = as.integer(c[[1]][2])
#}

get.END.from.vcf = function(x){ 
    a = as.character(x) 
    b = strsplit(a,split = ";")
    b = unlist(b)
    b = b[grep("^END=",b)]
    c = unlist(sapply(b,strsplit,split = '='))
    d = as.integer(c[complete.cases(as.integer(c))],echo = F)
}

addRow.id = function(df){ if(colnames(df)[1]!="row.id")
    colnames(df)[1:3] =c("V1","V2","V3")
df$row.id = as.character(paste(df$V1,df$V2,df$V3,sep = "."))
df = as.data.frame(cbind(df$row.id,df[,1:(ncol(df)-1)]))
colnames(df)[1:4] =c("row.id","chr","start","stop")
return(df)
}

setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data/")
#The files we read in here were produced by a shell command to subset based on the SVTYPE
#of the phase3 structural variants VCF file

ALL.wgs.integrated_sv_map_v2.20130502.CNV_ONLY <- read.delim("./ALL.wgs.integrated_sv_map_v2.20130502.CNV_ONLY.vcf", header=FALSE)
endCol = sapply(ALL.wgs.integrated_sv_map_v2.20130502.CNV_ONLY$V8,get.END.from.vcf)
m = cbind(ALL.wgs.integrated_sv_map_v2.20130502.CNV_ONLY[,c(1,2)],endCol)
colnames(m) = c('chr','start','stop')
m$chr = gsub(x = m[,1],pattern = '^',replacement = 'chr')
m = addRow.id(m)
m = m[complete.cases(m),]
write.table(x = m[,2:4],file = "Phase3_1KG_CNV.bed",
            col.names = F,row.names = F, sep = "\t",quote = F)


ALL.wgs.integrated_sv_map_v2.20130502.DUP_ONLY <- read.delim("./ALL.wgs.integrated_sv_map_v2.20130502.DUP_ONLY.vcf", header=FALSE)
endCol = sapply(ALL.wgs.integrated_sv_map_v2.20130502.DUP_ONLY$V8,get.END.from.vcf)
m = cbind(ALL.wgs.integrated_sv_map_v2.20130502.DUP_ONLY[,c(1,2)],endCol)
colnames(m) = c('chr','start','stop')
m$chr = gsub(x = m[,1],pattern = '^',replacement = 'chr')
m = addRow.id(m)
m = m[complete.cases(m),]
write.table(x = m[,2:4],file = "Phase3_1KG_DUP.bed",
            col.names = F,row.names = F, sep = "\t",quote = F)

ALL.wgs.integrated_sv_map_v2.20130502.DEL_ONLY <- read.delim("./ALL.wgs.integrated_sv_map_v2.20130502.DEL_ONLY.vcf", header=FALSE)
endCol = sapply(ALL.wgs.integrated_sv_map_v2.20130502.DEL_ONLY$V8,get.END.from.vcf)
m = cbind(ALL.wgs.integrated_sv_map_v2.20130502.DEL_ONLY[,c(1,2)],endCol)
colnames(m) = c('chr','start','stop')
m$chr = gsub(x = m[,1],pattern = '^',replacement = 'chr')
m = addRow.id(m)
m = m[complete.cases(m),]
write.table(x = m[,2:4],file = "Phase3_1KG_DEL.bed",
            col.names = F,row.names = F, sep = "\t",quote = F)

ALL.wgs.integrated_sv_map_v2.20130502.INV_ONLY <- read.delim("./ALL.wgs.integrated_sv_map_v2.20130502.INV_ONLY.vcf", header=FALSE)
endCol = sapply(ALL.wgs.integrated_sv_map_v2.20130502.INV_ONLY$V8,get.END.from.vcf)
m = cbind(ALL.wgs.integrated_sv_map_v2.20130502.INV_ONLY[,c(1,2)],endCol)
colnames(m) = c('chr','start','stop')
m$chr = gsub(x = m[,1],pattern = '^',replacement = 'chr')
m = addRow.id(m)
m = m[complete.cases(m),]
write.table(x = m[,2:4],file = "Phase3_1KG_INV.bed",
            col.names = F,row.names = F, sep = "\t",quote = F)


#Now, dos2unix ... sort -Vk1.4,2 ... these in unix and use intersect beds to bin in windows...


###########
