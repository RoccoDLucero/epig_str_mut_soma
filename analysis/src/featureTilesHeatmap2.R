featureTilesHeatmap2 = function(featureTilesDF,firstSet,secondSet,
                                pct= 1,rev = F,
                                prec =1000,p.cutoff = 0.05,reportTiles = T,
                                title = ""){
    
    #featureTilesDF =z;firstSet = 5;secondSet = 26;pct=1;rev =F;prec = 100;reportTiles = T
    targTiles = list()
    HILO = NULL
    if(rev == T){HILO = c("highest","most positive")}
    if(rev== F){HILO = c("lowest","most negative")}
    All.Enrichments1 = data.frame(matrix(data = NA,nrow =500,ncol = 3))
    #compute enrichment for each tile
    k =1 #initialize counter
    for(meth in firstSet){
        for(sv in secondSet){
            my.foi = c(meth,sv)
            tmp.enr = getTileEnrichment2(featureTilesDF = featureTilesDF,my.foi,pct = pct,precision = prec,rev = rev)
            rownames(All.Enrichments1)[k] = tmp.enr[[1]]
            All.Enrichments1[k,] = c(tmp.enr[[1]],tmp.enr[[2]][[1]],tmp.enr[[3]])
            if(reportTiles == T){targTiles[[k]] = tmp.enr[[2]][[2]]
                                    names(targTiles)[[k]] = tmp.enr$comparison
            }
            k = k+1
        }
    }
    colnames(All.Enrichments1) = names(tmp.enr)
    All.Enrichments1 = All.Enrichments1[complete.cases(All.Enrichments1),]
    #Bonf = length(secondSet)
    #All.Enrichments1[as.numeric(All.Enrichments1[,3])>(p.cutoff/Bonf),2] = NA
    #write.table(x = All.Enrichments1,file = "EnrichmentComparisons1.txt",
    #            col.names = T,row.names = F, sep = "\t",quote = F)
    dimsr = list(SV_names,Methylome_names)
    r = matrix((as.numeric(All.Enrichments1[,2])),nrow = length(SV_names), dimnames = dimsr)
    r[r==0] =NA
    notes = r
    notes[abs(log2(as.numeric(notes)))>=1] = '*'
    notes[abs(log2(as.numeric(notes)))<=1] = ''
    #notes
    #max(abs(log2(r)))
    dev.new()
    par(bg = "transparent")
    #par(bg = "grey85")
    numBreaks = 256
    pairs.breaks = seq(from=-3.2,to=3.2,length.out=257)
    heatmap.2(log(r,base = 2),na.rm = T,trace = 'none',margins = c(10,18),
              breaks = pairs.breaks,
              symbreaks = any(breaks<0,na.rm = T),
              density.info = 'none', key.xlab = 'Log2 Fold Enrichment', key.title = NA,
              cexCol = 1.9,cexRow = 1.5,cellnote = notes,notecol = "black",
              col = bluered(numBreaks),Colv = F,Rowv = F,
              main = paste("\nEnrichment of Structural Variants in\n 100kb Tiles with",
                           HILO[1],pct,"% CpG Methylation\n",HILO[2],"first and second deltas\n", 
                           "in human male developmental germline",sep = " "))
    targTiles
}