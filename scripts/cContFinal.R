# crosscont 2015/12/16 with Rcircos
# modif 2016/05/01 to change colour at lower values than 1 %

library(reshape2)
#
library(RCircos)
library("circlize")

#setwd("/home/cristina/R/crosscont/")
# table<- read.table(file="/home/cristina/Work/Projects/crossCont/SARCSU_SARCS_contam/CcEnd.csv", header=TRUE,sep="\t")

# Run using args.=> Rscript cContFinal.R "/home/cristina/Work/Projects/crossCont/SARCSU_SARCS/CcEnd.csv" "./SARCSU_plot.pdf"
#                                                                                             table         path plot's name

args<-commandArgs(TRUE)
end.table <- args[1]
proj<-args[2]     #

table<- read.table(file=end.table, header=TRUE,sep="\t")

###

table$X <- NULL    
colnames(table)<-gsub("\\.","-",colnames(table))      
tt<-t(table)        
colnames(tt)<-tt[1,]                                     
tt<-tt[-1,]                                           
class(tt) <- "numeric"                               

# col_mat

col_mat<-tt
row.names(col_mat) <- NULL
color <- rand_color(nrow(col_mat), transparency = 0)      
num.col<-ncol(col_mat)                                 
                                               
for (x in 1:num.col) {                                
   col_mat[,x] <- color
}
for (x in 1:num.col) { 
  for (y in 1:num.col) { 
    if(x==y){col_mat[x,y]<-"#e0ebeb"}
  }
}
col_mat[tt < 1] = "#e0ebeb"                               # change colour if value is < 1, to highlight the other values.

grid.c<- NULL
samp <-table$Samples
samp<-as.vector(samp)
grid.c<- cbind(samp,color)
row.names(grid.c)<-samp
grid.c<-grid.c[,-1]
# exportarlo como pdf con tamaño: 13, 13 inches(default), landscape

pdf(file=proj,width=13, height=13)
#pdf(file='/home/cristina/Work/Projects/crossCont/SARCSU_SARCS/x2.pdf',width=13, height=13)         

#chordDiagram
chordDiagram(tt, grid.col = grid.c, col = col_mat, annotationTrack = "grid", directional = 1, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.arr.length = 0.02, diffHeight = 0.04, preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.7)
  circos.axis(h = "top", labels.cex = 0.03, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
circos.clear()

dev.off()
