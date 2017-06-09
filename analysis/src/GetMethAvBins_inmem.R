#Title GetMethAvBins_inmem.R
#Hard Coded for now, The script should take as input argument the list of 
#chomosome sizes and a parameter describing the window size
#!/bin/bash
#args = (commandArgs(trailingOnly=TRUE))
cat("starting script!")

#setwd("./MethylomeWindowFiles")
setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data")
options(scipen = 999)
WinCoordsFile = "hg19_0.1Mb_windowCoords.sorted"
GenomeWindows = as.data.frame(read.table(paste('./',WinCoordsFile,sep = ''),sep = "\t", header=F,as.is = T))
#currentMethylome = "./testMethylome.bed.sorted"
currentMethylome = "./GSE63818_PGC_7W_embryo1_M_methylation_calling.CpG.bed.sorted"
#This script should take the genome windows as input, then read sorted 
#methylome (CpG) files,line by line and take either the average of the
#methylation percentage at each CpG in the window or
#the methylation percent in the window as the ration of all methylated versus
#total CpGs (coverage is >=5 at each CpG) from a preprocessing step

MethCalls = read.table(currentMethylome,header = F,sep = "\t",stringsAsFactors = F)
ByCpGAvs = rep(9, nrow(GenomeWindows))
ByCoverageAvs = rep(9, nrow(GenomeWindows))

GetChromInt = function(ChromStrng){ #Takes string,returns 4th character as int
    tmp = ChromStrng[[1]]
    try(as.numeric(strsplit(tmp,"")[1][[1]][4]))
}
winlen = as.integer(GenomeWindows[1,3])- as.integer(GenomeWindows[1,2]) +1

#Initialization
i = 1 #represents current CpGRead position
for(genWindow in 1:nrow(GenomeWindows[,])){
    #update window coordinates
    WinChr = GetChromInt(GenomeWindows[genWindow,1])
    LoLim = as.numeric(GenomeWindows[genWindow,2])
    HiLim = as.numeric(GenomeWindows[genWindow,3])
    
    CpGRead = MethCalls[i,]
    CpGPos = c(GetChromInt(CpGRead[1,1]),as.numeric(CpGRead[1,2]))
    
    #Ensure the two tables are "synchronized"
    while(WinChr != CpGPos[1] | !(CpGPos[2] >= LoLim && CpGPos[2] <= HiLim)){
        #if current line in methylome file is before winstart, read more lines
        #this assumes that the chromosomes will match
        while(WinChr > CpGPos[1]){
            i = i+1
            CpGRead = MethCalls[i,]
            CpGPos = c(GetChromInt(CpGRead[1,1]),as.numeric(CpGRead[1,2]))
        }
        if (WinChr == CpGPos[[1]] & CpGPos[[2]] <= LoLim){  
            i = i+1
            CpGRead = MethCalls[i,]
            CpGPos = c(GetChromInt(CpGRead[1,1]),as.numeric(CpGRead[1,2]))
        }
        if (WinChr == CpGPos[[1]] & CpGPos[[2]] >= HiLim){  
            #Need to catch genomic window coordinates up to methylome file
            break
        }
        if (i >=nrow(MethCalls)){break}
    }
    
    AvByCpG = as.integer(rep(NA,winlen))
    MethylatedCnt = as.integer(rep(NA,winlen))
    UnMethylatedCnt = as.integer(rep(NA,winlen))
    CpGPos = c(GetChromInt(CpGRead[1,1]),as.numeric(CpGRead[1,2]))
    j = 1 # represents the index of the CpG within the current window
    while(WinChr == CpGPos[1] && (CpGPos[2] >= LoLim && CpGPos[2] <= HiLim)){
        avmet = as.numeric(CpGRead[8])
        metreads = as.integer(CpGRead[6])
        unmetreads = as.integer(CpGRead[7])
        AvByCpG[j] = avmet
        MethylatedCnt[j] = metreads
        UnMethylatedCnt[j] =  unmetreads
        if (i >=nrow(MethCalls)){break}
        j = j+1
        i = i+1
        CpGRead = MethCalls[i,]
        CpGPos = c(GetChromInt(CpGRead[1,1]),as.numeric(CpGRead[1,2]))
        
    }
    
    winCpgMean =mean(AvByCpG,na.rm = T)
    ByCpGAvs[genWindow] = winCpgMean
    
    metSum = sum(MethylatedCnt,na.rm = T)
    CountsAv = metSum/(metSum + sum(UnMethylatedCnt,na.rm = T))
    ByCoverageAvs[genWindow] = CountsAv
    if (i >=nrow(MethCalls)){break}
}
    

colnames(GenomeWindows) = c('Chromosome','start','stop')
outframe = cbind(GenomeWindows,ByCpGAvs,ByCoverageAvs)
colnames(outframe)[4:5] = c('Methylation by CpG','Methylation by combined reads') 
write.table(x = outframe,file = paste((winlen/1000000),'Mb_window_methylation_avg.txt'),
            quote = F,row.names = F,col.names = T,sep = "\t")

