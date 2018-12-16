library(fdapace)
library(data.table)

#fgf4 (mm10)
chr=7
start=144861386  # + strand
end=144865243
EXT=500

files <-list.files(path = ".", pattern = "acc_processed.tsv", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) #[1:5]

length(files)

# read scNMT-seq data in the folder
Ly <- list()
Lt <- list()
counter <- 0
for(i in 1:length(files) ){
    print(i)
    c1 <- fread(files[i], head=FALSE)
    c1f <- c1[  which (c1$V1==chr & c1$V2>= start-EXT & c1$V2 <= start+EXT)  , ]
    if (length(c1f$V3)>0  ){
        counter <-  counter + 1
        Ly[[counter]] <- c1f$V3
        Lt[[counter]] <- c1f$V2
    }
    rm(c1,c1f)
}

# Number of cells with at least one GpC value
length(Ly)
# Proportion of cells with data
100* ( length(Ly)/length(files) )


#  Each vector in t should be in ascending order in fdapace
Ly_sorted <- list()
Lt_sorted <- list()
ID <- list()
for (j in 1:length(Lt) ){
    temp <- sort(Lt[[j]], index.return=TRUE, decreasing=FALSE)$ix
    Lt_sorted[[j]] <- Lt[[j]][temp]
    Ly_sorted[[j]] <- Ly[[j]][temp]
    #ID[[j]] <- rep(j, length(unlist(Ly[[j]][temp]) ) )  # for face R package
     rm(temp)

}

#--------------------------------------------------------------------------------------
# fdapace
#--------------------------------------------------------------------------------------
pace <- FPCA(Ly=Ly_sorted, Lt=Lt_sorted, optns = list(maxK=30, nRegGrid=100, plot=TRUE) ) #methodSelectK='BIC', userBwMu = 10,,  nRegGrid=100,kFoldMuCov=10   dataType='Sparse'
# , methodBwCov='GMeanAndGCV' , methodBwMu ='GMeanAndGCV', 
pdf('fig1.pdf',height=6.5, width=6.25)
CreatePathPlot( pace, K=10, subset = 2:5, main = "", pch = 1, showMean=TRUE,xlab='chr7', ylab='GpC accessibility')
dev.off()


library(wesanderson)
cellcolor <- wes_palette("BottleRocket2",n=5, type="discrete")

pdf('fig1.pdf',height=5, width=5.2)
par(mfrow=c(2,3))
par(mar=c(4,4,1,1))
for (i in 2:5){
    if (i==4 | i==5) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], xlab='chr7', ylab='GpC accessibility',main=paste('cell ', i-1) )
    if (i==2 | i==3) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], ylab='GpC accessibility',xlab='',main=paste('cell ', i-1) )
#Sys.sleep(5)
}
dev.off()

pdf('fig2.pdf',height=5.5, width=5.25)
CreateScreePlot(pace)
dev.off()

par(mfrow=c(1,2))
    library(ks)
  CreateOutliersPlot(pace, optns = list(K = 3, variant = 'KDE'))
  CreateFuncBoxPlot(pace, xlab = 'Days', ylab = '# of eggs laid', optns = list(K =3, variant='bagplot'))
#plot(pace)

 CreateFuncBoxPlot(pace, xlab = 'Days', ylab = '# of eggs laid', optns = list(K =3, variant='bagplot'))
#plot(pace)

pdf('fig3.pdf',height=4.6, width=4.25)
 CreateFuncBoxPlot(pace, xlab = 'chr7', ylab = 'GpC accessibility',main='Functional box-plot',col='red')
dev.off()


pdf('FIGURE1.pdf',height=5, width=7)
par(mfrow=c(2,3))
par(mar=c(5,4,2,1))
for (i in 2:5){
    if (i==4 | i==5) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], xlab='chr7', ylab='GpC accessibility',main=paste('cell ', i-1) )
    if (i==2 | i==3) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], ylab='GpC accessibility',xlab='chr7',main=paste('cell ', i-1) )
#Sys.sleep(5)
}
CreateScreePlot(pace)
CreateFuncBoxPlot(pace, xlab = 'chr7', ylab = 'GpC accessibility',main='Functional box-plot',col='red')
dev.off()

