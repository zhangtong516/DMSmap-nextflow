#!/usr/bin/env Rscript

args = commandArgs(T);
if (length(args) < 2) {
    stop(paste( "\n", sub("--file=","",commandArgs(F)[4]),"mapCsv sampleTsv", sep=" "))
}
mapCsv = args[1];
sampleTsv = args[2];

library('data.table')
library('stringr')
library('ggplot2')
library('doMC')


sampleInfo = fread(sampleTsv)
sampleInfo[, UniqueID:= str_c(SeqBatch, "_", LibID)]
mapRate = fread(mapCsv)
colnames(mapRate) = c("libID", "totalReads","mappedReads", "uniqMapped","unMapped")
mapRate[, mapRate := (mappedReads / totalReads) * 100]
mapRate[,Lib :=  lapply(.SD, function(x) {sampleInfo[UniqueID==x,str_c(Sample,"__",Replicates)]}),
        .SDcols="libID", by="libID" ]
mapRate[,LibType :=  lapply(.SD, function(x) {sampleInfo[UniqueID==x,LibType]}), .SDcols="libID", by="libID" ]
mapRate[,Sample :=  lapply(.SD, function(x) {sampleInfo[UniqueID==x,Sample]}), .SDcols="libID", by="libID" ]
mapRateL = melt(mapRate, id = c("libID","Lib","LibType","Sample"), measure=c("totalReads","mappedReads","mapRate"));
mapRateL = mapRateL[order(Lib)]

pdfFile = gsub("txt", "pdf", mapCsv)
pdf(pdfFile,width=9,height=5,bg="white")
p<-ggplot(data=mapRateL[variable != "mapRate",], aes(x=Lib, y=value,fill=variable, color=Sample)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(data=mapRateL[variable == "mapRate",],
              aes(x=Lib, y=4e+7, label=format(value,digits=3), vjust=0.9, fill=NULL),,
              position = position_dodge(width=0.1), angle=50, size=6) +
    facet_grid(LibType~.,  scales = "free") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p = p + theme(axis.text.x = element_text(angle = 50, hjust = 1))
p = p + theme(axis.text.x= element_text(size = 15 ), axis.text.y=element_text(size=15),
              axis.title.x= element_text(size = 15 ), axis.title.y=element_text(size=15),
              legend.text=element_text(size=15))
print(p)
dev.off()
