#Analyze Transcription Factor Binding Site Motif Enrichment
#In the 2012 Simko et. al. PRDM14 paper from Monica Justice lab
#

library(biomaRt)
#library(org.Hs.eg.db)
#library(annotate)
library(grid)
library(parallel)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite('biomaRt')
library(biomaRt)
dep.pckgs <- c('BiocGenerics', 'Biostrings','MotifDb','testthat',
               'gtools','BiocStyle', 'knitr')
biocLite(dep.pckgs,suppressUpdates = T)
pckgs <- c('PWMEnrich','PWMEnrich.Mmusculus.background')
biocLite(pckgs[2],type = 'source',suppressUpdates = T)
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


CN_GT_3_FilteredGenes <- read.delim("./input/CN_GT_3_FilteredGenesTableVIEW.txt",
                                    stringsAsFactors = F,check.names = T)
CN_LT_0.5_FilteredGenes <- read.delim("./input/CN_LT_0.5_FilteredGenesTableVIEW.txt",
                                      stringsAsFactors = F, check.names = T)


biocLite("BSgenome",suppressUpdates = T) 
biocLite("BSgenome.Mmusculus.UCSC.mm9",type ='source',suppressUpdates = T) #installs the mouse mm9 genome) 
library("BSgenome.Mmusculus.UCSC.mm9") 
Mmus <- BSgenome.Mmusculus.UCSC.mm9

geneset <- CN_GT_3_FilteredGenes
geneset$Gene.Start <- as.integer(gsub(pattern = ',',replacement = '',x = geneset$Gene.Start))
geneset$Gene.Stop <-  as.integer(gsub(pattern = ',',replacement = '',x = geneset$Gene.Stop))


#Possibly need to take strand information into account here
my.seqs.gene  <- getSeq(Mmus, geneset$Chromosome,
                        geneset$Gene.Start,
                        geneset$Gene.Stop)

my.seqs.upstream2k  <- getSeq(Mmus, geneset$Chromosome,
                              (geneset$Gene.Start-2000),
                              geneset$Gene.Start)

my.seqs.downstream2k  <- getSeq(Mmus, geneset$Chromosome,
                                geneset$Gene.Stop,
                                (geneset$Gene.Stop+2000))

my.seqs.upstream100k  <- getSeq(Mmus, geneset$Chromosome,
                                (if(geneset$Gene.Start-100000 > 10000){geneset$Gene.Start-100000}else NA),
                                geneset$Gene.Start)
                                

my.seqs.downstream100k  <- getSeq(Mmus, geneset$Chromosome,
                                geneset$Gene.Stop,
                                (geneset$Gene.Stop+100000))

upstream2kEnr <- motifEnrichment(sequences = my.seqs.upstream2k,pwms = PWMLogn.mm9.MotifDb.Mmus)
downstream2kEnr <- motifEnrichment(sequences = my.seqs.downstream2k,pwms = PWMLogn.mm9.MotifDb.Mmus)
upstream100kEnr <- motifEnrichment(sequences = my.seqs.upstream100k,pwms = PWMLogn.mm9.MotifDb.Mmus)
downstream100kEnr <- motifEnrichment(sequences = my.seqs.downstream100k,pwms = PWMLogn.mm9.MotifDb.Mmus)
genebodyEnr <- motifEnrichment(sequences = my.seqs.gene,pwms = PWMLogn.mm9.MotifDb.Mmus)



enrichment.report.up2k <- groupReport(upstream2kEnr)
enrichment.report.up100k <- groupReport(upstream100kEnr)
enrichment.report.dn2k <- groupReport(downstream2kEnr)
enrichment.report.dn100k <- groupReport(downstream100kEnr)
enrichment.report.gn <- groupReport(genebodyEnr)


write.table(x = enrichment.report.up2k@d, file = './output/2kupstream.motifs.txt',quote = F,sep = '\t',row.names = F)
write.table(x = enrichment.report.dn2k@d, file = './output/2kdnstream.motifs.txt',quote = F,sep = '\t',row.names = F)
write.table(x = enrichment.report.up100k@d, file = './output/100kupstream.motifs.txt',quote = F,sep = '\t',row.names = F)
write.table(x = enrichment.report.dn100k@d, file = './output/100kdnstream.motifs.txt',quote = F,sep = '\t',row.names = F)
write.table(x = enrichment.report.gn@d, file = './output/genebody.motifs.txt',quote = F,sep = '\t',row.names = F)
