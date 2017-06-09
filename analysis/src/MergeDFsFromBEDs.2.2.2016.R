#MergeDFsFromBEDs.R
MergeDFsFromBEDs = function(){
source('H:/Dropbox/BRL/Methylation_Struct_Mut/MyScripts/getTileEnrichment.R', echo=F)
setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data")
##Import the rest of the Structural Variant Data Sets from  Li/Harris and 1KG
##Automate the import
PGC_7.19W_emb1_CpGMe <- read.delim("./GSE63818_PGC_7-19W_embryo1_M_100kb_hg19_methylation.txt", header=FALSE)
colnames(PGC_7.19W_emb1_CpGMe)[4:8] =  c('mPGC.WK7.me','mPGC.WK10.me','mPGC.WK11.me','mPGC.WK13.me','mPGC.WK19.me')
sperm_CpGMe <- read.delim("./GSE30340_human_sperm_CpG_100kbTile_hg19_molaro.sorted.wig", header=FALSE, stringsAsFactors=FALSE)
colnames(sperm_CpGMe)[4] = "sperm_CpGMe"
#Somatic Control Methylomes
Somatic_100kb <- read.delim("./GSE63818_Somatic_100kb_hg19_methylation.txt", header=FALSE)
colnames(Somatic_100kb)[4:ncol(Somatic_100kb)] =  c('Heart_5W_CpGMe','Soma_7w' ,'Soma_10w','Soma_11w','Soma_17w','Soma_19w')


#Following represents count data output from intersectBeds with -c switch
#and in this case compared to the set of 100kb tiling windows across hg19:
Rearrangement_Hum_Spec <- read.delim("./Rearrangement_Human_Specific.0.1Mb.hg19.bed", header=FALSE)
colnames(Rearrangement_Hum_Spec)[4] = "Rearrangement_Hum_Spec"
Phase3_all_bkpts <- read.delim("./1KG_phase3_all_bkpts.v5_0.1Mb.hg19.bed", header=FALSE)
colnames(Phase3_all_bkpts)[4] = "Phase3_all_bkpts"
CNV_Bipolar_Case <- read.delim("./singletonCNV_Bipolar_Case.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Bipolar_Case)[4] = "CNV_Bipolar_Case"
CNV_Bipolar_Control <- read.delim("./singletonCNV_Bipolar_Control.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Bipolar_Control)[4] = "CNV_Bipolar_Control"
CNV_Autism_Case <- read.delim("./rareCNV_Autism_Case.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Autism_Case)[4] = "CNV_Autism_Case"
CNV_Autism_Control <- read.delim("./rareCNV_Autism_Control.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Autism_Control)[4] = "CNV_Autism_Control"
CNV_DevDelay_Case <- read.delim("./rareCNV_DevDelay_Case.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_DevDelay_Case)[4] = "CNV_DevDelay_Case"
CNV_DevDelay_Control <- read.delim("./rareCNV_DevDelay_Control.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_DevDelay_Control)[4] = "CNV_DevDelay_Control"
CNV_Schiz_Case  <- read.delim("./rareCNV_Schizophrenia_Case.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Schiz_Case)[4] = "CNV_Schiz_Case"
CNV_Schiz_Control <- read.delim("./rareCNV_Schizophrenia_Control.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_Schiz_Control)[4] = "CNV_Schiz_Control"
CNV_270_HapMap <- read.delim("./CNV_270HapMap.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_270_HapMap)[4] = "CNV_270_HapMap"
CNV_450_HapMap <- read.delim("./CNV_450HapMap.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_450_HapMap)[4] = "CNV_450_HapMap"
CNV_400_MGL <- read.delim("./CNV_400MGL.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_400_MGL)[4] = "CNV_400_MGL"
CNV_WTCCC <- read.delim("./CNV_WTCCC.0.1Mb.hg19.bed", header=FALSE)
colnames(CNV_WTCCC)[4] = "CNV_WTCCC"
Phase3_mCNV <- read.delim("./Phase3_1KG_CNV_hg19_0.1Mb.bed", header=FALSE)
colnames(Phase3_mCNV)[4] = "Phase3_mCNV"
Phase3_DUP <- read.delim("./Phase3_1KG_DUP_hg19_0.1Mb.bed", header=FALSE)
colnames(Phase3_DUP)[4] = "Phase3_DUP"
Phase3_DEL <- read.delim("./Phase3_1KG_DEL_hg19_0.1Mb.bed", header=FALSE)
colnames(Phase3_DEL)[4] = "Phase3_DEL"
Phase3_INV <- read.delim("./Phase3_1KG_INV_hg19_0.1Mb.bed", header=FALSE)
colnames(Phase3_INV)[4] = "Phase3_INV"

########
#Import ENCODE and Chia TF DATA:
Cmyc.ENCODE <- read.delim("./Cmyc.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Cmyc.ENCODE)[4] = "Cmyc.ENCODE"
Ctcf.ENCODE <- read.delim("./Ctcf.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Ctcf.ENCODE)[4] = "Ctcf.ENCODE"
Foxa1.ENCODE <- read.delim("./Foxa1.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Foxa1.ENCODE)[4] = "Foxa1.ENCODE"
Jund.ENCODE <- read.delim("./Jund.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Jund.ENCODE)[4] = "Jund.ENCODE"
Mef2a.ENCODE <- read.delim("./Mef2a.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Mef2a.ENCODE)[4] = "Mef2a.ENCODE"
Mef2c.ENCODE <- read.delim("./Mef2c.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Mef2c.ENCODE)[4] = "Mef2c.ENCODE"
Pou5f1.ENCODE <- read.delim("./Pou5f1.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Pou5f1.ENCODE)[4] = "Pou5f1.ENCODE"
Prdm1.ENCODE <- read.delim("./Prdm1.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Prdm1.ENCODE)[4] = "Prdm1.ENCODE"
PRDM14_ES_PEAKS <- read.delim("./PRDM14_ES_PEAKS_hg19_0.1Mb.bed", header=FALSE)
colnames(PRDM14_ES_PEAKS)[4] = "PRDM14_ES_PEAKS"
Rad21.ENCODE <- read.delim("./Rad21.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Rad21.ENCODE)[4] = "Rad21.ENCODE"
Suz12.ENCODE <- read.delim("./Suz12.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Suz12.ENCODE)[4] = "Suz12.ENCODE"
Yy1.ENCODE <- read.delim("./Yy1.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Yy1.ENCODE)[4] = "Yy1.ENCODE"
Znf143.ENCODE <- read.delim("./Znf143.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Znf143.ENCODE)[4] = "Znf143.ENCODE"
Znf263.ENCODE <- read.delim("./Znf263.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Znf263.ENCODE)[4] = "Znf263.ENCODE"
Znf274.ENCODE <- read.delim("./Znf274.ENCODE_hg19_0.1Mb.bed", header=FALSE)
colnames(Znf274.ENCODE)[4] = "Znf274.ENCODE"


##Import Annotations
Havana.gene.gencode <- read.delim("./gencode.v19.hg19_0.1Mb_HAVANA.gene.bed", header=FALSE)
colnames(Havana.gene.gencode)[4] = "Havana.gene.gencode"
#####I need a way to reliably and appropriately name the columns....
#####

#IncorporatePreprocessing steps on the data:
# add a list of case control pairs for the CNV data
# for the sets in this list make a function that removes the CNVs that
# Case and control have in common
# 









MethylationSets = list(Somatic_100kb = Somatic_100kb,
                    PGC_7.19W_emb1_CpGMe = PGC_7.19W_emb1_CpGMe,
                    sperm_CpGMe = sperm_CpGMe)

CNV_sets = list( Rearrangement_Hum_Spec = Rearrangement_Hum_Spec,
                 Phase3_all_bkpts = Phase3_all_bkpts,
                 Phase3_mCNV = Phase3_mCNV,
                 Phase3_DUP = Phase3_DUP,
                 Phase3_DEL = Phase3_DEL,
                 Phase3_INV = Phase3_INV,
                 CNV_Bipolar_Case = CNV_Bipolar_Case,
                 CNV_Bipolar_Control = CNV_Bipolar_Control,
                 CNV_Autism_Case = CNV_Autism_Case,
                 CNV_Autism_Control = CNV_Autism_Control,
                 CNV_DevDelay_Case = CNV_DevDelay_Case,
                 CNV_DevDelay_Control = CNV_DevDelay_Control,
                 CNV_Schiz_Case = CNV_Schiz_Case,
                 CNV_Schiz_Control = CNV_Schiz_Control,
                 CNV_WTCCC = CNV_WTCCC,
                 CNV_270_HapMap = CNV_270_HapMap,
                 CNV_450_HapMap = CNV_450_HapMap,
                 CNV_400_MGL = CNV_400_MGL)

Peak_Sets = list (Cmyc.ENCODE = Cmyc.ENCODE,
                  Ctcf.ENCODE = Ctcf.ENCODE,
                  Foxa1.ENCODE = Foxa1.ENCODE,
                  Jund.ENCODE = Jund.ENCODE,
                  Mef2a.ENCODE = Mef2a.ENCODE,
                  Mef2c.ENCODE = Mef2c.ENCODE,
                  Pou5f1.ENCODE = Pou5f1.ENCODE,
                  Prdm1.ENCODE =Prdm1.ENCODE,
                  PRDM14_ES_PEAKS = PRDM14_ES_PEAKS,
                  Rad21.ENCODE =Rad21.ENCODE,
                  Suz12.ENCODE = Suz12.ENCODE,
                  Yy1.ENCODE = Yy1.ENCODE,
                  Znf143.ENCODE = Znf143.ENCODE,
                  Znf263.ENCODE = Znf263.ENCODE,
                  Znf274.ENCODE = Znf274.ENCODE
                  )

Annotation_Sets = list(Havana.gene.gencode = Havana.gene.gencode)

allFeatures = c(MethylationSets,CNV_sets,Peak_Sets,Annotation_Sets)


rm(PGC_7.19W_emb1_CpGMe, sperm_CpGMe)
rm( Rearrangement_Hum_Spec, Phase3_all_bkpts,
    CNV_Bipolar_Case,CNV_Bipolar_Control,CNV_Autism_Case,CNV_Autism_Control,
    CNV_DevDelay_Case,CNV_DevDelay_Control,CNV_Schiz_Case,CNV_Schiz_Control,
    CNV_WTCCC,CNV_270_HapMap,CNV_450_HapMap,CNV_400_MGL)
rm(Ctcf.ENCODE,Foxa1.ENCODE,Pou5f1.ENCODE,
    Prdm1.ENCODE,PRDM14_ES_PEAKS,Suz12.ENCODE)

#This works for the  SV counts dfs but not the methylation average dfs
#because the number of colums is inconsistent
###
#####..Fix that...


addRow.id = function(df){ if(colnames(df)[1]!="row.id")
    df$row.id = as.character(paste(df$V1,df$V2,df$V3,sep = "."))
    df = as.data.frame(cbind(df$row.id,df[,1:(ncol(df)-1)]))
    colnames(df)[1:4] =c("row.id","chr","start","stop")
    return(df)
}


#After adding Row.id column, merge the data to a single data frame for
#processing by getTileEnrichment()...

merge2.dfs = function(dfList,x,y){
    df.out = merge.data.frame(x = dfList[[x]],y =  dfList[[y]],
                                all.x = T,by = c("row.id","chr","start","stop"))#'row.id')
    return(df.out[order(df.out[,2],df.out[,3]),])
}

#Must run addRow.id first
mergeAll.dfs = function(dfList){
    df.out =dfList[[1]]
    for(i in 1:(length(allFeatures)-1)){
    df.out = merge.data.frame(x = df.out,y =  dfList[[i+1]],
                              all.x = T,by = c("row.id","chr","start","stop"))
    }
    return(df.out[order(df.out[,2],df.out[,3]),])
}


allFeatures = lapply(allFeatures, addRow.id)
z = mergeAll.dfs(dfList = allFeatures)
}
