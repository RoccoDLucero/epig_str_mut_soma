#Title GetMethylationAvinGenomicBins.R
#Hard Coded for now, The script should take as input argument the list of 
#chomosome sizes and a parameter describing the window size
#!/bin/bash
#args = (commandArgs(trailingOnly=TRUE))

setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Processed_Data")
try(close.connection(myFileCon),silent = T)
options(scipen = 999)
winlen = 0.1 #Window Length in Mb
WinCoordsFile = "hg19_0.1Mb_windowCoords.sorted"
GenomeWindows = as.data.frame(read.table(paste('./',WinCoordsFile,sep = ''),sep = "\t", header=T,as.is = T))
currentMethylome = "./GSE63818_PGC_7W_embryo1_M_methylation_calling.CpG.bed.sorted"
#This script should take the genome windows as input, then read sorted 
#methylome (CpG) files,line by line and take either the average of the
#methylation percentage at each CpG in the window or
#the methylation percent in the window as the ration of all methylated versus
#total CpGs (coverage is >=5 at each CpG) from a preprocessing step


myFileCon = file(currentMethylome, open = "r")
#Store window methylation averages in two ways
ByCpGAvs = rep(9, nrow(GenomeWindows))
ByCoverageAvs = rep(9, nrow(GenomeWindows))
#ByCpGAvs = vector(mode = "numeric",length = nrow(GenomeWindows)) 
#ByCoverageAvs = vector(mode = "numeric",length = nrow(GenomeWindows)) 
#I need to write code that will identify the lines in the sorted methylome
#file that demarcate the beginning and end of each chromosome.
#I can use these min and max indices to ensure that the current line matches
# the current chromosome window
GetChromInt = function(ChromStrng){
    tmp = ChromStrng[[1]]
    as.numeric(strsplit(tmp,"")[1][[1]][4])
}

#Initialization
CpGRead = unlist(strsplit(x =readLines(con = myFileCon,n = 1),split = "\t")) 
CpGPos = list(GetChromInt(CpGRead[1]),as.numeric(CpGRead[2]))
for(genWindow in 1:nrow(GenomeWindows[,])){
    Chr = as.numeric(strsplit(GenomeWindows[genWindow,1],"")[1][[1]][4])
    #update window borders
    LoLim = as.numeric(GenomeWindows[genWindow,2])
    HiLim = as.numeric(GenomeWindows[genWindow,3])
    
    #Ensure the two files are "synchronized"
    while(Chr != CpGPos[[1]] | !(CpGPos[[2]] >= LoLim && CpGPos[[2]] <= HiLim)){
        #print(paste("checking position",CpGPos[[2]]))
        #if current line in methylome file is before winstart, read more lines
        #this assumes that the chromosomes will match
        if(Chr > CpGPos[[1]]){
            CpGRead = unlist(strsplit(x =readLines(con = myFileCon,n = 1),split = "\t"))
            CpGPos = list(GetChromInt(CpGRead[1]),as.numeric(CpGRead[2]))
        }
        if (Chr == CpGPos[[1]] & CpGPos[[2]] <= LoLim){  
            CpGRead = unlist(strsplit(x =readLines(con = myFileCon,n = 1),split = "\t"))
            CpGPos = list(GetChromInt(CpGRead[1]),as.numeric(CpGRead[2]))
            
        }
        if (Chr == CpGPos[[1]] & CpGPos[[2]] >= HiLim){  
            #Need to catch genomic window coordinates up to methylome file
            break
            #print("increment window coordinates")
        }
        
    }
    #print('Escaped_while_loop1')
    AvByCpG = as.numeric(c())
    MethylatedCnt = as.integer(c())
    UnMethylatedCnt = as.integer(c())
    CpGPos = list(GetChromInt(CpGRead[1]),as.numeric(CpGRead[2]))
    while(Chr == CpGPos[[1]] && (CpGPos[[2]] >= LoLim && CpGPos[[2]] <= HiLim)){
        avmet = as.numeric(CpGRead[8])
        metreads = as.integer(CpGRead[6])
        unmetreads = as.integer(CpGRead[7])
        AvByCpG = c(AvByCpG, avmet)
        MethylatedCnt = c(MethylatedCnt, metreads)
        UnMethylatedCnt = c(UnMethylatedCnt, unmetreads)
        CpGRead = unlist(strsplit(x =readLines(con = myFileCon,n = 1),split = "\t"))
        CpGPos = list(GetChromInt(CpGRead[1]),as.numeric(CpGRead[2]))
        #print(c(avmet,metreads,unmetreads))#
    }
    #print('escaped_while_loop2')#
    winCpgMean =mean(AvByCpG)
    #print(paste('mean by means =', winCpgMean))#
    ByCpGAvs[genWindow] = winCpgMean
    CountsAv = sum(MethylatedCnt)/(sum(MethylatedCnt)+sum(UnMethylatedCnt))
    #print(paste("countsAv =",CountsAv))#
    ByCoverageAvs[genWindow] = CountsAv
    #print(paste("window",genWindow))#
}
close.connection(myFileCon)
colnames(GenomeWindows) = c('Chromosome','start','stop')
outframe = cbind(GenomeWindows,ByCpGAvs,ByCoverageAvs)
colnames(outframe)[4:5] = c('Methylation by CpG','Methylation by combined reads') 
write.table(x = outframe,file = paste(winlen,'Mb_window_methylation_avg.txt'),
            quote = F,row.names = F,col.names = T,sep = "\t",)

