################################################################################
# Get background sequences for ROIs and random regions and make lognormal backgrounds
################################################################################
##GET THE GENOMIC COORDINATES AND SEQUENCES FOR:
##GENE BODIES,10 KB UPSTEAM,10 KB DOWNSTREAM,

##REQUIRES##
source("https://bioconductor.org/biocLite.R")
library(BSgenome)#, lib.loc = './Rlibs')
library(PWMEnrich)#, lib.loc = './Rlibs')
library(PWMEnrich.Hsapiens.background)
data("MotifDb.Hsap")
data("MotifDb.Hsap.PFM")
library("BSgenome.Hsapiens.UCSC.hg19")#, lib.loc = './Rlibs')
library("BSgenome.Hsapiens.UCSC.hg19.masked")#,lib.loc = './Rlibs')
library("GenomicFeatures")#,lib.loc = './Rlibs')
library(MotifDb)#, lib.loc = './Rlibs')
library(seqLogo)


#To make the gene-centric backgrounds we sample UCSC Hg19 known genes
#for random transcripts coordinates then
#get sequences for gene body and flanking regions
#Then compute the lognormal background for the PWMs in the PWMEnrich package.

#Import From Environment
ref.Tx.len <- Prdm14.length

Hg19 <- BSgenome.Hsapiens.UCSC.hg19
Hg19.mask <- BSgenome.Hsapiens.UCSC.hg19.masked
Hg19KG <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
#Hg19KG <- makeTranscriptDbFromUCSC(genome = "hg19", tablename = "knownGene")
Hg19Transcripts <- transcripts(Hg19KG)

h.sap.motifs <- query(MotifDb, 'sapiens',ignore.case = T)


my.genome <- Hg19.mask
TxSet <- Hg19Transcripts
minTxLen <- 2*1000
maxTxLen <- 3*ref.Tx.len
bgdepth <- 10000
#FOR UPSTREAM AND DOWNSTREAM FLANKING REGIONS:
up.flank <- 10*1000
dn.flank <- 10*1000



#Randomly sample transcripts of set size range:
set.seed(20160905)
geneBodSampSize <- 1000
genebody.ranges <- sample(TxSet,geneBodSampSize)
genebody.ranges <- genebody.ranges[genebody.ranges@ranges@width < maxTxLen]
genebody.ranges <- genebody.ranges[genebody.ranges@ranges@width > minTxLen]

#Get Sequences 
genebody.seqs <- getSeq(my.genome, genebody.ranges@seqnames,
                              genebody.ranges@ranges@start,
                              (genebody.ranges@ranges@start + genebody.ranges@ranges@width))

up.flank.seqs <- getSeq(my.genome, genebody.ranges@seqnames,
                              (genebody.ranges@ranges@start - up.flank),
                              genebody.ranges@ranges@start)

dn.flank.seqs <- getSeq(my.genome, genebody.ranges@seqnames,
                        (genebody.ranges@ranges@start + genebody.ranges@ranges@width),
                        (genebody.ranges@ranges@start + genebody.ranges@ranges@width + dn.flank))

#Remove transcripts that lie in masked regions and select desired number of
#sequences that will be included in the background computation
pick.good.seqs <- function(dna.set,depth=bgdepth){
    if(length(grep('N',dna.set)) != 0){dna.set <- dna.set[-grep('N',dna.set)]}
    dna.set <- dna.set[1:depth]
}

genebody.seqs <- pick.good.seqs(genebody.seqs)
up.flank.seqs <- pick.good.seqs(up.flank.seqs)
dn.flank.seqs <- pick.good.seqs(dn.flank.seqs)

#Create backgrounds
bg.source.string <- "hg19 UCSC transcripts of lengths between 2-70kb split into 500bp regions"
my.hg19.whole.genebody.logn.bg1 <- makePWMLognBackground(motifs = MotifDb.Hsap.PFM,
                                                        bg.seq = genebody.seqs,
                                                        bg.len = 500,
                                                        bg.source = bg.source.string,
                                                        algorithm = 'human',verbose = T)

saveRDS(my.hg19.whole.genebody.logn.bg,'./output/my.hg19.whole.genebody.logn.bg')
rm(my.hg19.whole.genebody.logn.bg)

bg.source.string <- "hg19 UCSC TSS minus 10kb split into 500bp regions"
my.hg19.TSS.minus10bk.logn.bg<- makePWMLognBackground(motifs = MotifDb.Hsap,
                                                        bg.seq = up.flank.seqs,
                                                        bg.len = 500,
                                                        bg.source = bg.source.string,
                                                        algorithm = 'human',verbose = T)

saveRDS(my.hg19.TSS.minus10bk.logn.bg,'./output/my.hg19.TSS.minus10bk.logn.bg')
rm(my.hg19.TSS.minus10bk.logn.bg)

bg.source.string <- "hg19 UCSC TTS plus 10kb split into 500bp regions"
my.hg19.TSS.plus10bk.logn.bg <- makePWMLognBackground(motifs = MotifDb.Hsap,
                                                        bg.seq = dn.flank.seqs,
                                                        bg.len = 500,
                                                        bg.source = bg.source.string,
                                                        algorithm = 'human',verbose = T)

saveRDS(my.hg19.TSS.plus10bk.logn.bg,'./output/my.hg19.TSS.plus10bk.logn.bg')
rm(my.hg19.TSS.plus10bk.logn.bg)
################################################################################
################################################################################
## We now have backgrounds against which to perform motif enrichment
##
################################################################################