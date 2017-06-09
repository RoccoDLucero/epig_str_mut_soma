#Analyze Transcription Factor Binding Site Motif Enrichment
#at the mm9 PRDM14 locus and surrounding genomic regions
#for Monica Justice lab

library(grid)
library(parallel)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade",suppressUpdates = T)
biocLite('biomaRt',suppressUpdates = T)
library(biomaRt)
dep.pckgs <- c('BiocGenerics', 'Biostrings','MotifDb','testthat',
               'gtools','BiocStyle', 'knitr')
biocLite(dep.pckgs,suppressUpdates = T)

pckgs <- c('PWMEnrich','PWMEnrich.Hsapiens.background')
biocLite(pckgs,type = 'source',suppressUpdates = T)
biocLite('regioneR',suppressUpdates = T)
library(regioneR)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
data("PWMCutoff4.hg19.MotifDb.Hsap")
data("PWMLogn.hg19.MotifDb.Hsap")
library(MotifDb)
data("MotifDb.Hsap")

biocLite("BSgenome",suppressUpdates = T)
library(BSgenome)
#available.genomes()
biocLite("BSgenome.Hsapiens.UCSC.hg19",type ='source',suppressUpdates = T) #installs the hg19 genome
library("BSgenome.Hsapiens.UCSC.hg19") 
Hg19 <- BSgenome.Hsapiens.UCSC.hg19

biocLite("BSgenome.Hsapiens.UCSC.hg19.masked",type ='source',suppressUpdates = T) #installs the masked hg19 genome
library("BSgenome.Hsapiens.UCSC.hg19.masked")
Hg19.mask <- BSgenome.Hsapiens.UCSC.hg19.masked

biocLite("GenomicFeatures",suppressUpdates = T)
library("GenomicFeatures")
#supportedUCSCtables()[1:4, ]
#Hg19KG <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
#saveFeatures(Hg19KG, file= "./input/Hg19.UCSC,KnownGenes.sqlite")
Hg19KG <- loadFeatures("./input/Hg19.UCSC,KnownGenes.sqlite")
Hg19Transcripts <- transcripts(Hg19KG)

#Example of data frame format for genomic regions:
#   Gene    Gene.Start  Gene.Stop   Chromosome  Strand  
#1  Lbx1    45308217    45309726    chr19       -      
#2  Adam6b  114727775   114730046   chr12       +

#PRDM14 coordinates from UCSC hg19: chr8:70,963,886-70,983,562
#PRDM14 coordinates from UCSC hg38: chr8:70,051,651-70,071,327


Prdm14.locus <- data.frame(t(c('Prdm14',70963886,70983562,'chr8', '-')),stringsAsFactors = F)
colnames(Prdm14.locus) <- c('Gene','Gene.Start','Gene.Stop','Chromosome','Strand')
Prdm14.locus[,2:3] <- sapply(Prdm14.locus[2:3],as.integer)
#sapply(Prdm14.locus,class)
Prdm14.length <- -(70963886 - 70983562)

if(T){Hg19 <- Hg19.mask}

################################################################################
# Get sequences for ROIs and random regions
################################################################################
###GET THE GENOMIC COORDINATES AND SEQUENCES FOR:
###(PRDM14:GENE BODY,100 KB UPSTEAM,100 KB DOWNSTREAM,
###RANDOM REGIONS OF SAME SIZE AS PRDM14:FROM WHOLE GENOME, BASED ON GENE STARTS 
###RANDOM REGIONS OF 100KB SIZE: FROM WHOLE GENOME, BASED ON GENE STARTS, BASED ON GENE ENDS)
###

#for random regions based on UCSC mm9 transcripts
set.seed(20160905)
geneBodSampSize <- 4000
minTxLen <- 2*1000
maxTxLen <- 5*Prdm14.length
bgdepth <- 1000

my.hg19.genebody.ranges <- sample(Hg19Transcripts,geneBodSampSize)
my.hg19.genebody.ranges <- my.hg19.genebody.ranges[my.hg19.genebody.ranges@ranges@width < maxTxLen]
my.hg19.genebody.ranges <- my.hg19.genebody.ranges[my.hg19.genebody.ranges@ranges@width > minTxLen]
#my.mm9.genebody.ranges <- my.mm9.genebody.ranges[1:(bgdepth)]
hist(my.hg19.genebody.ranges@ranges@width,breaks = 100)
 
my.hg19.whole.genebody.seqs <- getSeq(Hg19, my.hg19.genebody.ranges@seqnames,
                                      my.hg19.genebody.ranges@ranges@start,
                                     (my.hg19.genebody.ranges@ranges@start + my.hg19.genebody.ranges@ranges@width))

my.hg19.whole.genebody.seqs <- my.hg19.whole.genebody.seqs[-grep('N',my.hg19.whole.genebody.seqs)]
my.hg19.whole.genebody.seqs <- my.hg19.whole.genebody.seqs[1:bgdepth]

my.hg19.whole.genebody.logn.bg <- readRDS('./output/my.hg19.whole.genebody.logn.bg')
if(!exists('my.hg19.whole.genebody.logn.bg')){
    my.hg19.whole.genebody.logn.bg <- makePWMLognBackground(motifs = MotifDb.Hsap,
                                                       bg.seq = my.hg19.whole.genebody.seqs,
                                                       bg.len = 500,
                                                       bg.source = "hg19 UCSC transcripts of lengths between 2-70kb split into 500bp regions",
                                                       algorithm = 'human',verbose = T)

saveRDS(my.hg19.whole.genebody.logn.bg,'./output/my.hg19.whole.genebody.logn.bg')
} 
if(F){
my.mm9.genebody.P14len.seqs <- getSeq(Mmus, my.mm9.genebody.ranges@seqnames,
                                      my.mm9.genebody.ranges@ranges@start,
                                      (my.mm9.genebody.ranges@ranges@start + Prdm14.length))

my.mm9.genebody.P14len.logn.bg <- makeBackground(motifs = PWMLogn.mm9.MotifDb.Mmus,
                                     bg.seq = my.mm9.genebody.P14len.seqs,
                                     type = 'logn')
}
################################################################################
################################################################################

if(F){
    #FOR PRDM14 GENE BODY:
    Prdm14.length <- (Prdm14.locus$Gene.Stop - Prdm14.locus$Gene.Start)
    Prdm14.seq.gene  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                               Prdm14.locus$Gene.Start,
                               Prdm14.locus$Gene.Stop)
        
        #for random regions from entire mm9
        my.randomRegs <- createRandomRegions(nregions = geneBodSampSize,
                                             length.mean = Prdm14.length,
                                             mask = NULL,
                                             genome = 'mm9')
        
        my.randomSeqs <- getSeq(Mmus, my.randomRegs@seqnames,
                                my.randomRegs@ranges@start,
                                (my.randomRegs@ranges@start + my.randomRegs@ranges@width -1))
        
        
        #FOR UPSTREAM AND DOWNSTREAM FLANKING REGIONS:
        flank <- 10*1000
        flankSampSize <- 1000
        
        #For the PRDM14 gene flanking regions
        Prdm14.seq.upstream10k  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                                           (if(Prdm14.locus$Gene.Start-flank > flank){Prdm14.locus$Gene.Start-flank}else NA),
                                           Prdm14.locus$Gene.Start)
        
        
        Prdm14.seq.downstream10k  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                                             Prdm14.locus$Gene.Stop,
                                             (Prdm14.locus$Gene.Stop+flank))
        
        #for random regions based on UCSC mm9 transcripts
        my.randomTx.ranges <- sample(mm9Transcripts,flankSampSize)
        
        my.randomTx.10kb.up.seqs <- getSeq(Mmus, my.randomTx.ranges@seqnames,
                                            my.randomTx.ranges@ranges@start - flank,
                                            (my.randomTx.ranges@ranges@start))
        
        my.randomTx.10kb.dn.seqs <- getSeq(Mmus, my.randomTx.ranges@seqnames,
                                            my.randomTx.ranges@ranges@start,
                                            (my.randomTx.ranges@ranges@start + flank))
        
        
        #for random regions from entire mm9
        my.randomRegs10k <- createRandomRegions(nregions = flankSampSize,
                                                 length.mean = (flank),
                                                 mask = NULL,
                                                 genome = 'mm9')
        
        my.randomSeqs10k <- getSeq(Mmus, my.randomRegs10k@seqnames,
                                    my.randomRegs10k@ranges@start,
                                    (my.randomRegs10k@ranges@start + my.randomRegs10k@ranges@width -1))

#This needs to account for occasional sampling from unsequenced regions
#my.randomSeqs100k <- my.randomSeqs100k[-c(5,15,28)]
################################################################################
## We now have sequences on which to perform motif enrichment
## ROI enrichments will be contrasted agains the random region enrichments
################################################################################
}


################################################################################
# perform motif enrichment analysis, relative to Prdm14 gene body as background
################################################################################
"Prdm14.gene.enr <- motifEnrichment(sequences = Prdm14.seq.gene, pwms = MotifDb.Mmus)
Prdm14.gene.enr <- motifEnrichment(sequences = Prdm14.seq.gene, pwms = my.mm9.whole.genebody.logn.bg)
prdm14.report.gn <- groupReport(Prdm14.gene.enr)
prdm14.report.gn[1:20]

Prdm14.upstream100k.enr <- motifEnrichment(sequences = Prdm14.seq.upstream100k,
                                           pwms = PWMLogn.mm9.MotifDb.Mmus,
                                           bg = my.mm9.whole.genebody.logn.bg)

Prdm14.downstream100k.enr <- motifEnrichment(sequences = Prdm14.seq.downstream100k,
                                             pwms = PWMLogn.mm9.MotifDb.Mmus,
                                             bg = my.mm9.whole.genebody.logn.bg)

Prdm14.upstream10k.enr <- motifEnrichment(sequences = Prdm14.seq.upstream10k ,pwms = PWMLogn.mm9.MotifDb.Mmus)
Prdm14.downstream10k.enr <- motifEnrichment(sequences = Prdm14.seq.downstream10k ,pwms = PWMLogn.mm9.MotifDb.Mmus)

################################################################################
##Generate enrichment reports
prdm14.report.gn <- groupReport(Prdm14.gene.enr)
my.randomTx.report <-groupReport(my.randomTx.enr)
my.randoms.report <- groupReport(my.randoms.enr)

prdm14.report.up100k <- groupReport(Prdm14.upstream100k.enr)
prdm14.report.dn100k <- groupReport(Prdm14.downstream100k.enr)
my.randomTx.100kb.up.report <- groupReport(my.randomTx.100kb.up.enr)
my.randomTx.100kb.dn.report <- groupReport(my.randomTx.100kb.dn.enr)
my.randoms.100k.report <- groupReport(my.randoms100k.enr)

prdm14.report.up10k <- groupReport(Prdm14.upstream10k.enr)
prdm14.report.dn10k <- groupReport(Prdm14.downstream10k.enr)
my.randomTx.10kb.up.report <- groupReport(my.randomTx.10kb.up.enr)
my.randomTx.10kb.dn.report <- groupReport(my.randomTx.10kb.dn.enr)
my.randoms.10k.report <- groupReport(my.randoms10k.enr)
################################################################################
## We now have the enrichment reports for the ROIs and random regions
## We need to contrast these to clarify the relevant factors
################################################################################


################################################################################
#make a plot that shows when the overlaps begin to occur and also draw a pval threshold
################################################################################
plotOverlapsByRankThresh <- function(rpt1,rpt2,resolution =50,title){
    rankThresh <- seq(1,length(intersect(rpt1@d$target,rpt2@d$target)),resolution)
    OLs <-c()
    for(i in rankThresh){
        overlaps <- length(intersect(rpt1@d$target[1:i],rpt2@d$target[1:i]))
        OLs <- c(OLs,overlaps)
    }
    plot(rankThresh ,OLs, main = title, xlab = 'Rank Threshold', ylab = 'Shared Motifs')
}

#This shows that random regions enrich for different factors than prdm14 gene body
res <- 10
plotOverlapsByRankThresh(prdm14.report.gn,my.randomTx.report,
                         resolution = res, title = "randomTX vs. PRDM14 gene body")

plotOverlapsByRankThresh(prdm14.report.gn,my.randoms.report,
                         resolution = res, title = "randomRegion vs. PRDM14 gene body")

##Upstream
plotOverlapsByRankThresh(prdm14.report.up100k,my.randomTx.100kb.up.report,
                         resolution = res, title = "randomTX vs. 100k upstream")

plotOverlapsByRankThresh(prdm14.report.up100k,my.randoms.100k.report,
                         resolution = res, title = "random vs. 100k upstream")

plotOverlapsByRankThresh(prdm14.report.up10k,my.randomTx.10kb.up.report,
                         resolution = res, title = "randomTX vs. 10k upstream")

plotOverlapsByRankThresh(prdm14.report.up10k,my.randoms.10k.report,
                         resolution = res, title = "random vs. 10k upstream")

#Downstream
plotOverlapsByRankThresh(prdm14.report.dn100k,my.randomTx.100kb.dn.report,
                         resolution = res, title = "randomTX vs. 100k dnstream")

plotOverlapsByRankThresh(prdm14.report.dn100k,my.randoms.100k.report,
                         resolution = res, title = "random vs. 100k dnstream")

plotOverlapsByRankThresh(prdm14.report.dn10k,my.randomTx.100kb.dn.report,
                         resolution = res, title = "randomTX vs. 10k dnstream")

plotOverlapsByRankThresh(prdm14.report.dn10k,my.randoms.100k.report,
                         resolution = res, title = "random vs. 10k dnstream")

#Random v Random
plotOverlapsByRankThresh(my.randomTx.100kb.up.report,my.randoms.100k.report,
                         resolution = res, title = "randomTX 100kbup vs. random 100kb")

plotOverlapsByRankThresh(my.randomTx.100kb.dn.report,my.randoms.100k.report,
                         resolution = res, title = "randomTX 100kbdn vs. random 100kb")

################################################################################
## These plots show that the top enriched candidates in the PRDM14 gene body
## are different from those in random regions
################################################################################

################################################################################
#Report the unique factors within the top 50 for our ROI
################################################################################


getTopUnique <- function(roi,bg,rankLim=50){
    uniqueTFs <- setdiff(roi@d$target[1:rankLim],bg@d$target[1:rankLim])
    roi@d[which(uniqueTFs %in% roi@d$target[1:rankLim]),]
}

write.table(getTopUnique(prdm14.report.gn,my.randomTx.report),
            file = './output/UniqueMotifsTop50_mm9_PRDM14gene.tab',
            quote = F,sep = '\t',row.names = F)
write.table(getTopUnique(prdm14.report.up100k,my.randomTx.100kb.up.report),
            file = './output/UniqueMotifsTop50_mm9_PRDM14100kbUp.tab',
            quote = F,sep = '\t',row.names = F)
write.table(getTopUnique(prdm14.report.dn100k,my.randomTx.100kb.dn.report),
            file = './output/UniqueMotifsTop50_mm9_PRDM14100kbDn.tab',
            quote = F,sep = '\t',row.names = F)
write.table(getTopUnique(prdm14.report.up10k,my.randomTx.10kb.up.report),
            file = './output/UniqueMotifsTop50_mm9_PRDM1410kbUp.tab',
            quote = F,sep = '\t',row.names = F)
write.table(getTopUnique(prdm14.report.dn10k,my.randomTx.10kb.dn.report),
            file = './output/UniqueMotifsTop50_mm9_PRDM1410kbDn.tab',
            quote = F,sep = '\t',row.names = F)

"