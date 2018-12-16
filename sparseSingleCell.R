library(fdapace)
library(data.table)

#fgf4 (mm10), + strand
chr=7
start=144861386  
end=144865243
EXT=500

# read scNMT-seq data 
files <-list.files(path = ".", pattern = "acc_processed.tsv", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

length(files)

#---------------------------------
# Input for FPCA
#---------------------------------
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
    rm(temp)
}

#---------------------------------
# fdapace
#---------------------------------
pace <- FPCA(Ly=Ly_sorted, Lt=Lt_sorted, optns = list(maxK=30, nRegGrid=100, plot=TRUE) )

library(wesanderson)
cellcolor <- wes_palette("BottleRocket2",n=5, type="discrete")

pdf('FIGURE1.pdf',height=5, width=7)
par(mfrow=c(2,3))
par(mar=c(5,4,2,1))
for (i in 2:5){
    if (i==4 | i==5) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], xlab='chr7', ylab='GpC accessibility',main=paste('cell ', i-1) )
    if (i==2 | i==3) CreatePathPlot( pace, K=9, subset = i, main = "", pch = 16, showMean=FALSE,col=cellcolor[i-1], ylab='GpC accessibility',xlab='chr7',main=paste('cell ', i-1) )
}
CreateScreePlot(pace)
CreateFuncBoxPlot(pace, xlab = 'chr7', ylab = 'GpC accessibility',main='Functional box-plot',col='red')
dev.off()

