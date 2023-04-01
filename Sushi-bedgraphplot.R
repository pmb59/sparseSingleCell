library("Sushi")

bw <- read.table("ENCFF974HWL.bedGraph", head = FALSE)
colnames(bw) <- c("chrom", "start", "end", "value")

# Fgf4 (mm10)
chrom <- "chr7"
chromstart <- 144861386 - 20000
chromend <- 144861386 + 20000

pdf("fgf4.pdf", width = 8, height = 5)
layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
par(mar = c(3, 4, 1, 1))

plotBedgraph(bw, chrom, chromstart, chromend, colorbycol = SushiColors(5), linecolor = "black", addscale = FALSE)
labelgenome(chrom, chromstart, chromend, n = 4, scale = "Mb")
mtext("Read-depth normalized signal", side = 2, line = 1.75, cex = 0.8, font = 2)
axis(side = 2, las = 2, tcl = .2)
zoomregion1 <- c(144861386 - 500, 144861386 + 500)
zoomsregion(zoomregion1, extend = c(0.01, 0.13), wideextend = 0.05, offsets = c(0, 0))

plotBedgraph(bw, chrom, chromstart = zoomregion1[1], chromend = zoomregion1[2], colorbycol = SushiColors(5), linecolor = "black")
labelgenome(chrom, chromstart = zoomregion1[1], chromend = zoomregion1[2], n = 4, scale = "bp", edgeblankfraction = 0.2, cex.axis = .75)
zoombox()
mtext("", side = 2, line = 1.75, cex = 0.5, font = 2)
axis(side = 2, las = 2, tcl = .2)
dev.off()
