#Hard Coded for now, The script should take as input argument the list of chomosome sizes
#and a parameter describing the window size
setwd("H:/Dropbox/BRL/Methylation_Struct_Mut/Raw_Data")
options(scipen = 999)
ChromSizes = read.delim("./hg19.chrom.sizes", header=F)
ChromSizes = ChromSizes[1:24,]
GenomeBuild = "hg19"
windowLen =  (10^6) #c(10^3,10^4,10^6,10^7)
winMbs = windowLen/10^6
NumWindows = ceiling(sum(as.numeric(ChromSizes[,2]))/windowLen)

#Write this as a function that intakes a chromsize file and a vector of desired window lengths
#into which to bin the given genome.
##This should probably be a perl or awk script....
ChromWindows = matrix(nrow = NumWindows, ncol = 3,dimnames = list(c(),c("Chr","winStart", "winStop")))
i = 1
for (chrom in ChromSizes[,1]){
    winstart = 1
    winstop = windowLen
    foo = 0
    upperLim = ChromSizes$V2[ChromSizes$V1 == chrom]
    while(winstop <= upperLim){
        ChromWindows[i,] = c(chrom,winstart,winstop)
        winstart = winstart + windowLen
        winstop = winstop + windowLen
        i = i+1
    }
    if (winstop > upperLim){
        ChromWindows[i,] = c(chrom,winstart,upperLim)
    }
    
}

ChromWindows = ChromWindows[complete.cases(ChromWindows),]
write.table(x = ChromWindows,sep = '\t',
            file = paste("../Processed_Data/",GenomeBuild,'_',winMbs,"Mb_windowCoords",sep =""),quote = F,row.names = F)

