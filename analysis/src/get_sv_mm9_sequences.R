my.get.SV.mm9.sequences <- function(my.SVs){
    biocLite("BSgenome",suppressUpdates = T) 
    biocLite("BSgenome.Mmusculus.UCSC.mm9",type ='source',suppressUpdates = T) #installs the human genome (~850 MB download). 
    library("BSgenome.Mmusculus.UCSC.mm9") 
    
    geneset <- my.SVs
    #??DO I NEED TO LIFTOVER from MM9 to MM10 in order to get sequencers??
    my.seqs.gene  <- getSeq(Mmus, geneset$Chromosome,
                       geneset$Gene.Start,
                       geneset$Gene.Stop)
    
    #Possibly need to take strand information into account here
    my.seqs.upstream2k  <- getSeq(Mmus, geneset$Chromosome,
                            (geneset$Gene.Start-2000),
                            geneset$Gene.Start)
    
    my.seqs.downstream2k  <- getSeq(Mmus, geneset$Chromosome,
                            geneset$Gene.Stop,
                            (geneset$Gene.Stop+2000))

}
