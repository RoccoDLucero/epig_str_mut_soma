######  Rename this or generalize it to allow enrichment by
######  modes: (percent thresholds v count, count v count, pres abs v count
######                    pres abs v pres abs)
######  
######  getTileEnrichment.R
######  By: Rocco Lucero
######  Created: February 2, 2016
######  Last Modified February 10, 2016 by Rocco Lucero
######
######  getTileEnrichment. R is a script that takes as input a data frame of
######  genome feature  counts (or levels) in uniformly-sized genomic
######  ranges (tiles) and computes the enrichment of the overlap
######  between a pair of the available features within the set of tiles
######  defined by a fixed percentage supplied by the user. Fold enrichment
######  is reported.
######
######  FOI IS A VECTOR OF COL INDICES OR COLNAMES
        
getTileEnrichment = function(featureTilesDF,foi,pct = 1,rev = F,precision = 100){
    set.seed(20160202)
    #library(coin)
    df = featureTilesDF[order(featureTilesDF[,foi[1]],decreasing = rev),foi]
    df = df[complete.cases(df),]
    #print(head(df))
    
    comp.Enrichment = function(df=df,pct=pct){
        totalTiles = nrow(df)
        testTiles = round(pct/100*totalTiles,digits = 0)
        testEnrich = sum(df[1:testTiles,2])/testTiles 
        otherEnrich = sum(df[(testTiles+1):totalTiles,2])/(totalTiles - testTiles) 
        Target.enrichment = testEnrich/otherEnrich
        
    }
    
    perm.P.Val = function(df=df,pct=pct,trgEnr,prec = precision){
        #sample "testTiles" number of tiles randomly 1000 times
        #and compute enrichmnet p value
        perm.dist = c(rep(NULL,prec))
        for(i in 1:prec){
            testTiles = round(pct/100*nrow(df),digits = 0)
            my.samp = sample(x = rownames(df),size = testTiles,replace = F)
            df = df[my.samp,]
            sampledEnr = comp.Enrichment(df,pct)
            perm.dist[i] = sampledEnr
        }
        ####
        if(trgEnr>=1){length(perm.dist[which(perm.dist>trgEnr)])/prec}
        else{length(perm.dist[which(perm.dist<trgEnr)])/prec}
    }
    
    enr = comp.Enrichment(df = df,pct = pct)
    p.val = perm.P.Val(df = df,pct = 5,trgEnr = enr,prec = precision)
    if (p.val == 0){p.val = paste("<",(1/precision),sep = " ")}
    result = list(paste(colnames(df)[1],"vs.",colnames(df)[2],sep=" "),enr,p.val)
    names(result) = c("comparison","fold enrichment","p-value")
    result
}



getTileEnrichment.pres = function(featureTilesDF,featurePair,rev = T,precision = 100,mode = "chi"){
    set.seed(20160210)
    primaryFeature = featurePair[1]
    df = featureTilesDF[order(featureTilesDF[,primaryFeature],decreasing = rev),featurePair]
    df = df[complete.cases(df),]
    
    #in this function df will now be a 2column data frame

    comp.Enrichment = function(df=df){
        testTiles = df[df[,1]>0,]
        otherTiles = df[df[,1]==0,]
        a = testTiles[testTiles[,2]>0,] #Tiles that Share the primary feature
        b = testTiles[testTiles[,2]==0,] #Tiles with only the Primary feature
        c = otherTiles[otherTiles[,2]>0,]# Tiles with only Secondary Feature
        d = otherTiles[otherTiles[,2]==0,] #Neither Feature
        chiTab = matrix(c(nrow(a),nrow(b),nrow(c),nrow(d)),ncol = 2,byrow = T)
        if(mode == "tab"){r = chiTab}
        if(mode == "pct"){r = (nrow(a)/(nrow(a)+nrow(b)))}
        if(mode == "chi"){r = chisq.test(x = chiTab)}
        ##if(mode == "enr"){r = (nrow(a)/(nrow(a)+nrow(b)))/(nrow(a)/(nrow(a)+nrow(c)))}
        #Needs Fix This is intended to give a crude enrichment value
        r
        
    }
    comp.Enrichment(df)
  
}





