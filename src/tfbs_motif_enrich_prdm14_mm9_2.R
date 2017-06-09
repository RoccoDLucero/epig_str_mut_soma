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

pckgs <- c('PWMEnrich','PWMEnrich.Mmusculus.background')
biocLite(pckgs,type = 'source',suppressUpdates = T)
biocLite('regioneR',suppressUpdates = T)
library(regioneR)
library(PWMEnrich)
library(PWMEnrich.Mmusculus.background)

data(MotifDb.Mmus.PFM)
data('MotifDb.Mmus')
data('PWMLogn.mm9.MotifDb.Mmus')
data('PWMCutoff4.mm9.MotifDb.Mmus')
data('PWMCutoff5.mm9.MotifDb.Mmus')
data('PWMPvalueCutoff1e2.mm9.MotifDb.Mmus')
data('PWMPvalueCutoff1e3.mm9.MotifDb.Mmus')
data('PWMPvalueCutoff1e4.mm9.MotifDb.Mmus')
data('mm9.upstream2000')

biocLite("BSgenome",suppressUpdates = T) 
biocLite("BSgenome.Mmusculus.UCSC.mm9",type ='source',suppressUpdates = T) #installs the mouse mm9 genome
library("BSgenome.Mmusculus.UCSC.mm9") 
Mmus <- BSgenome.Mmusculus.UCSC.mm9

biocLite("BSgenome.Mmusculus.UCSC.mm9.masked",type ='source',suppressUpdates = T) #installs the mouse mm9 genome
library("BSgenome.Mmusculus.UCSC.mm9.masked")
Mmus.mask <- BSgenome.Mmusculus.UCSC.mm9.masked

biocLite("GenomicFeatures",suppressUpdates = T)
library("GenomicFeatures")
#supportedUCSCtables()[1:4, ]
mm9KG <- makeTxDbFromUCSC(genome = "mm9", tablename = "knownGene")
#saveFeatures(mm9KG, file= "./input//Mmus.UCSC,KnownGenes.sqlite")
#mm9KG <- loadFeatures("./input//Mmus.UCSC,KnownGenes.sqlite")
mm9Transcripts <- transcripts(mm9KG)

#Example of data frame format for genomic regions:
#   Gene    Gene.Start  Gene.Stop   Chromosome  Strand  
#1  Lbx1    45308217    45309726    chr19       -      
#2  Adam6b  114727775   114730046   chr12       +

#PRDM14 coordinates from UCSC mm9: Chr1:13113457-13127163 bp, - strand
#Position: mm9 chr1:13,103,509-13,117,244 Size: 13,736 Total Exon Count: 8 Strand: -


#Prdm14.locus <- data.frame(t(c('Prdm14',13113457,13127163,'chr1', '-')),stringsAsFactors = F)
Prdm14.locus <- data.frame(t(c('Prdm14',13103509,13117244,'chr1', '-')),stringsAsFactors = F)
colnames(Prdm14.locus) <- c('Gene','Gene.Start','Gene.Stop','Chromosome','Strand')
Prdm14.locus[,2:3] <- sapply(Prdm14.locus[2:3],as.integer)
#sapply(Prdm14.locus,class)

if(T){Mmus <- Mmus.mask}

################################################################################
# Get sequences for ROIs and random regions
################################################################################
###GET THE GENOMIC COORDINATES AND SEQUENCES FOR:
###(PRDM14:GENE BODY,100 KB UPSTEAM,100 KB DOWNSTREAM,
###RANDOM REGIONS OF SAME SIZE AS PRDM14:FROM WHOLE GENOME, BASED ON GENE STARTS 
###RANDOM REGIONS OF 100KB SIZE: FROM WHOLE GENOME, BASED ON GENE STARTS, BASED ON GENE ENDS)
###
#FOR PRDM14 GENE BODY:
Prdm14.length <- (Prdm14.locus$Gene.Stop - Prdm14.locus$Gene.Start)
Prdm14.seq.gene  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                        Prdm14.locus$Gene.Start,
                        Prdm14.locus$Gene.Stop)

#for random regions based on UCSC mm9 transcripts
set.seed(20160901)
geneBodSampSize <- 100
my.randomTx.ranges <- sample(mm9Transcripts,geneBodSampSize)

my.randomTx.seqs <- getSeq(Mmus, my.randomTx.ranges@seqnames,
                           my.randomTx.ranges@ranges@start,
                           (my.randomTx.ranges@ranges@start + Prdm14.length))

#for random regions from entire mm9
my.randomRegs <- createRandomRegions(nregions = geneBodSampSize,
                                     length.mean = Prdm14.length,
                                     mask = NULL,
                                     genome = 'mm9')

my.randomSeqs <- getSeq(Mmus, my.randomRegs@seqnames,
                        my.randomRegs@ranges@start,
                        (my.randomRegs@ranges@start + my.randomRegs@ranges@width -1))

#FOR UPSTREAM AND DOWNSTREAM FLANKING REGIONS:
flank <- 100*1000
flankSampSize <- 100

Prdm14.seq.upstream100k  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                                (if(Prdm14.locus$Gene.Start-flank > flank){Prdm14.locus$Gene.Start-flank}else NA),
                                Prdm14.locus$Gene.Start)


Prdm14.seq.downstream100k  <- getSeq(Mmus, Prdm14.locus$Chromosome,
                                  Prdm14.locus$Gene.Stop,
                                  (Prdm14.locus$Gene.Stop+flank))

#for random regions based on UCSC mm9 transcripts
my.randomTx.ranges <- sample(mm9Transcripts,flankSampSize)

my.randomTx.100kb.up.seqs <- getSeq(Mmus, my.randomTx.ranges@seqnames,
                           my.randomTx.ranges@ranges@start - flank,
                           (my.randomTx.ranges@ranges@start))

my.randomTx.100kb.dn.seqs <- getSeq(Mmus, my.randomTx.ranges@seqnames,
                                    my.randomTx.ranges@ranges@start,
                                    (my.randomTx.ranges@ranges@start + flank))


#for random regions from entire mm9
my.randomRegs100k <- createRandomRegions(nregions = flankSampSize,
                                     length.mean = (flank),
                                     mask = NULL,
                                     genome = 'mm9')

my.randomSeqs100k <- getSeq(Mmus, my.randomRegs100k@seqnames,
                        my.randomRegs100k@ranges@start,
                        (my.randomRegs100k@ranges@start + my.randomRegs100k@ranges@width -1))

#FOR UPSTREAM AND DOWNSTREAM FLANKING REGIONS:
flank <- 10*1000
flankSampSize <- 100

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



################################################################################
# perform motif enrichments on the selected sequences from the previous step
################################################################################
Prdm14.gene.enr <- motifEnrichment(sequences = Prdm14.seq.gene ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randomTx.enr <- motifEnrichment(sequences = my.randomTx.seqs ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randoms.enr <-  motifEnrichment(sequences = my.randomSeqs,pwms = PWMLogn.mm9.MotifDb.Mmus)

Prdm14.upstream100k.enr <- motifEnrichment(sequences = Prdm14.seq.upstream100k ,pwms = PWMLogn.mm9.MotifDb.Mmus)
Prdm14.downstream100k.enr <- motifEnrichment(sequences = Prdm14.seq.downstream100k ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randomTx.100kb.up.enr <- motifEnrichment(sequences = my.randomTx.100kb.up.seqs  ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randomTx.100kb.dn.enr <- motifEnrichment(sequences = my.randomTx.100kb.dn.seqs  ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randoms100k.enr <-  motifEnrichment(sequences = my.randomSeqs100k,pwms = PWMLogn.mm9.MotifDb.Mmus)

Prdm14.upstream10k.enr <- motifEnrichment(sequences = Prdm14.seq.upstream10k ,pwms = PWMLogn.mm9.MotifDb.Mmus)
Prdm14.downstream10k.enr <- motifEnrichment(sequences = Prdm14.seq.downstream10k ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randomTx.10kb.up.enr <- motifEnrichment(sequences = my.randomTx.10kb.up.seqs  ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randomTx.10kb.dn.enr <- motifEnrichment(sequences = my.randomTx.10kb.dn.seqs  ,pwms = PWMLogn.mm9.MotifDb.Mmus)
my.randoms10k.enr <-  motifEnrichment(sequences = my.randomSeqs10k,pwms = PWMLogn.mm9.MotifDb.Mmus)
#my.randomSeqs10k <- my.randomSeqs10k[-c(31,90,91,95)]

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