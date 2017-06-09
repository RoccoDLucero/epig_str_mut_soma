getTile.KS = function(featureTilesDF,foi,rev = T, plots =F){
    set.seed(20160213)
    #featureTilesDF = z
    #foi = c(5,13)
    #print(colnames(z)[foi])
    df = featureTilesDF[order(featureTilesDF[,foi[2]],decreasing = rev),foi]
    df = df[complete.cases(df),]
    SV_T = df[df[,2]>0,1]
    SV_F = df[df[,2]==0,1]
    
    if(plots == T){
        dev.new()
        #b = lowess(x = SV_T,f=5)
        plot(ecdf(SV_F),pch = '.',cex = 4,col = 'black', lwd =5)
        lines(ecdf(SV_T),pch = '.',cex = 2, col = 'dodgerblue3',lwd = 5)
    }
    #qqplot(x = SV_T,y = SV_F)
    #print(ks.test(SV_T,SV_F))
    ks.test(SV_T,SV_F)
}