#!/usr/bin/env Rscript

# easyRNA
#
# Copyright (C) Ping Zhu
# Contact: Ping Zhu <pingzhu.work@gmail.com>
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


# Library installation check
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is.installed("RColorBrewer")){
  warning("Detect package \"RColorBrewer\" is not installed in your R enviroment.")
  warning("Trying to install the \"RColorBrewer\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("RColorBrewer")
}
if(!is.installed("optparse")){
  warning("Detect package \"optparse\" is not installed in your R enviroment.")
  warning("Trying to install the \"optparse\" package.")
  warning("If failed, please try to install it mannually.")
  install.packages("optparse")
}

## Load libraries
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="Input file necessary for each mode"),
  make_option(c("-m", "--mode"), dest = "mode", default = "",
              help="Data analysis or visualization. Required arguments are listed below each mode.\
              \
              deg\tDifferentially expressed genes. \
                   - two groups of cells\
              *barplot\tBarplot of specified cells and genes \
                   - rData image file. \
                   - one or more groups of cells. \
                   - a list of genes\
              *heatmap\tHeatmap of specified cells and genes \
                   - rData image file. \
                   - one or more groups of cells. \
                   - a list of genes\
              *volcano\tVolcano plot of differentially expressed genes\
                   - output from deg"),
  make_option(c("-c", "--cells"), dest = "cells", default = "", type = "character",
              help = "Groups of cells"),
  make_option(c("-G", "--genes"), dest = "genes", default = "", type = "character",
              help = "List of genes separated by comma or a file containing genes by rows"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] Output file name. [default: modes.SysDate]"),
  make_option(c("-l", "--labels"), dest = "labels", default = T,
              help = "[opt] Show labels of DEGs in volcano plot. TRUE (T) or FALSE (F). [default: T]"),
  make_option(c("-W", "--width"), dest = "figure.width", default = 7,
              help = "[opt] Width of figure (inch). [default: 7]"),
  make_option(c("-H", "--height"), dest = "figure.height", default = 7,
              help = "[opt] Height of figure (inch). [default: 7]")
  #make_option(c("-f","--format"), dest = "figure.format", default = "pdf",
  #            help = "[opt] Format of output figure. Alternative: png. [default: pdf]"),
  #make_option(c("-R","--resolution"), dest = "figure.resolution", default = 300,
  #            help = "[opt] Resolution in ppi. Only available for png format. [default: 300]")
)

parser <- OptionParser(usage = "Rscript %prog [options]", option_list=option_list)
## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
#if(infile == ""){ # default, STDIN
#  infile <- file("stdin")
#} else { # user specified
  if( infile == "" | file.access(infile) == -1){ # file not exists
    #print_help(parser)
    stop("Please specify the input file using '-i'. Use '-h' to see detailed usage.")
  }
#}

modes <- arguments$mode
if( ! modes %in% c("deg", "barplot", "heatmap", "volcano")){
    stop("Please specify the mode using '-m'. Use '-h' to see detailed usage.")
}

labels <- arguments$labels
if( ! labels %in% c(T, TRUE, F, FALSE)){
    stop("Only TRUE (T) or FALSE (F) could be assigned to '-l'. Use '-h' to see detailed usage.")
}


outfile <- arguments$outfile
if(outfile == ""){ # default, "modes.Date"
  outfile <- modes
} else { # user specified
  outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
#figure.format <- arguments$figure.format
#if(! figure.format %in% c("pdf", "png")){ # format not support
#  print_help(parser)
#} else {
#  outfile <- paste(outfile, figure.format, sep = ".")
#}
#figure.resolution <- arguments$figure.resolution

# date marker
sysDate <- suppressWarnings(Sys.Date())

# subrutines

cellIndex <- function(){
    # cells
    cells <- arguments$cells
    # format
    # Group1:Type=='A' & Class=='B';Group2:Type=='C' & Class=='D'
    cellGroups <- unlist(strsplit(cells, ":|;", perl =T))
    cell.idx <- list()
    for( i in 1:(length(cellGroups)/2)){
        j <- i * 2
        cell.idx[[cellGroups[j-1]]] <- with( meta,  which(eval(parse(text=cellGroups[j]))) )
    }
    return(cell.idx)
}

geneIndex <- function(){
    # genes
    genes <- unique(arguments$genes)
    if (! file.exists(genes)){
        genes <- unlist(strsplit(genes, split=",|,\\s|\\s", perl = T))
    } else {
        genes <- read.table(genes, stringsAsFactors = F)[,1]
    }
    gene.idx <- match(genes, row.names(exp))
    gene.idx <- gene.idx[ !is.na(gene.idx) ]
    return(gene.idx)
}

Barplot <- function(cell, gene, unit=8){
    gene.num <- length(gene)
    this.data <- exp[gene, as.numeric(unlist(cell))]
    pages <- ceiling(gene.num / unit)
    cell.colors <- colorPick(cell)
    for(i in 1:pages){
        pdf(file=paste(outfile, "barplot", sysDate, i, "pdf", sep="."), width=figure.width, height=figure.height,
            onefile=FALSE)
        par(mar=c(1,1,1,2), oma=c(2,1,2,3))
        layout(matrix(c(1, rep(2, unit), 3, seq(4, length.out=unit)), ncol = 2), width=c(1.5,8.5))
        
        # panel numbers for multiple genes
        start <- (i-1)*unit + 1
        end <- min(nrow(this.data), i * unit)
        
        # panel 1, null
        plot(x=0,y=0, xlim=c(0,1), type = "n", axes = F, xlab = "", ylab = "")

        # panel 2, ylab
        plot(x=0,y=0, xlim=c(0,1), ylim=c(0,1), type = "n", axes = F, xlab = "", ylab = "", yaxs = "i")
        text(x=1, y=1 - (end-start-1)/(unit*2), adj=c(1,1), srt=90, "Relative expression", xpd=T, cex=1.5)
        #text(x=1, y=1 - (end-start-1)/(unit*2), srt=90, "Relative expression", xpd=T, cex=1.5)

        # panel 3, figure legend
        plot(x=0,y=0, xlim=c(0,1), type = "n", axes = F, xlab = "", ylab = "")
        #legend(x=0.588, y=0, xjust=0.5, yjust = 1, legend = names(cell), col = unique(cell.colors), ncol = 6, bty = "n", lwd=2, xpd = T)
        legend.ncol <- length(cell)
        if (legend.ncol > 6){ legend.ncol <- 6}
        legend("bottom", legend = names(cell), col = unique(cell.colors), bty = "n", lwd=2, xpd = T, ncol=legend.ncol)

        # main figure
        if(is.null(nrow(this.data))){
            barplot(height = this.data, names.arg = F, col = cell.colors, border = cell.colors, lwd=2, las=1)
            mtext(side=3, text = gene)
        } else {
            for(j in start:end){
                barplot(height = this.data[j,], names.arg = F, col = cell.colors, border = cell.colors, lwd=2, las=1)
                mtext(side=3, text = row.names(this.data)[j])
            }
        }
        
        invisible(dev.off())
    }
}

colorPick <- function(cells, type="discrete"){
    # colors
    library(RColorBrewer)
    cols <- colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(cells))
    col.cells <- c()
    for( i in 1:length(cells)){
        col.cells <- c(col.cells, rep(cols[i], length(cells[[i]])))
    }
    return(col.cells)
}

findDEG <- function(cell, fc=2, fdr=0.01, expCut=2){
    if (length(cell[[1]]) < 2 || length(cell[[2]]) < 2){
        warning("Error: at least two samples needed to be sepcified in each group.")
    } else {
        gene.test <- apply(exp, 1, function(x){ # students' t-test
                               g1.mean <- mean(x[cell[[1]]], na.rm=T)
                               g2.mean <- mean(x[cell[[2]]], na.rm=T)
                               this.fc <- g1.mean/g2.mean
                               this.test <- try(t.test(x[cell[[1]]], x[cell[[2]]]), silent = TRUE)
                               if(class(this.test) == "try-error"){
                                   this.pvalue <- NA
                               } else {
                                   this.pvalue <- this.test$p.value
                               }
                               c(g1.mean, g2.mean, this.fc, this.pvalue)
            })

        test.matrix <- matrix(rep(NA, nrow(exp)*5), ncol = 5)
        for(i in 1:nrow(exp)){
            test.matrix[i,] <- c(row.names(exp)[i], gene.test[, i])
        }
        test.matrix <- as.data.frame(test.matrix, stringsAsFactors=F)
        colnames(test.matrix) <- c("Gene", paste("Mean", names(cell)[1], sep="."), 
                                   paste("Mean", names(cell)[2], sep="."), "FC.Group1_vs_Group2", "Pvalue")
        test.matrix$CorrectedPvalue <- p.adjust(as.numeric(test.matrix$Pvalue), method = "fdr")
        # DEG
        #fc.idx <- which(as.numeric(test.matrix$FC.Group1_vs_Group2) >= fc | as.numeric(test.matrix$FC.Group1_vs_Group2) <= 1/fc)
        exp.idx <- which(as.numeric(test.matrix$Mean.Group1) >= expCut | as.numeric(test.matrix$Mean.Group2) >= expCut)
        fdr.idx <- which(as.numeric(test.matrix$CorrectedPvalue) <= fdr | is.na(test.matrix$CorrectedPvalue))
        test.matrix$UpDown <- "Stable"
        test.matrix$UpDown[ intersect(intersect(exp.idx, fdr.idx), which(as.numeric(test.matrix$FC.Group1_vs_Group2) >= fc)) ] <- "Up" 
        test.matrix$UpDown[  intersect(intersect(exp.idx, fdr.idx), which(as.numeric(test.matrix$FC.Group1_vs_Group2) <= 1/fc)) ] <- "Down" 
        write.table(test.matrix, file=paste(outfile, "DEG", sysDate, "xls", sep="."), quote=F, sep="\t", col.names=T, row.names=F)
    }
    # else { # compare multiple groups
    #   clusters <- levels(meta[[groupBy]])
    #   gene.test <- sapply(clusters, function(x){
    #     group1.idx <- which(meta[[groupBy]] == x)
    #     group2.idx <- which(meta[[groupBy]] != x)
    #     apply(data, 1, function(y){
    #       g1.mean <- mean(y[group1.idx], na.rm=T)
    #       g2.mean <- mean(y[group2.idx], na.rm=T)
    #       this.fc <- g1.mean/g2.mean
    #       this.test <- try(t.test(y[group1.idx], y[group2.idx]), silent = TRUE)
    #       if(class(this.test) == "try-error"){
    #         this.pvalue <- NA
    #       } else {
    #         this.pvalue <- this.test$p.value
    #       }
    #       c(x, g1.mean, g2.mean, this.fc, this.pvalue)
    #     })
    #   })
    #   test.matrix <- matrix(rep(NA, nrow(data)*length(clusters)*6), ncol = 6)
    #   for(i in 1:ncol(gene.test)){
    #     for(j in 1:nrow(data)){
    #       k <- (j-1) * 5 + 1
    #       test.matrix[(i-1)*nrow(data)+j,] <- c(row.names(data)[j], gene.test[k:(k+4), i])
    #     }
    #   }
    #   test.matrix <- as.data.frame(test.matrix, stringsAsFactors=F)
    #   colnames(test.matrix) <- c("Gene", "Cell.group", "Mean.thisGroup", "Mean.otherGroups", "FC.this_vs_other", "Pvalue")
    #   test.matrix$CorrectedPvalue <- p.adjust(as.numeric(test.matrix$Pvalue), method = "fdr")
    #   # DEG
    #   test.matrix <- test.matrix[ which(as.numeric(test.matrix$FC.this_vs_other) >= fc & as.numeric(test.matrix$Mean.thisGroup) >= exp), ]
    #   test.matrix <- test.matrix[ which(as.numeric(test.matrix$CorrectedPvalue) <= fdr | is.na(test.matrix$CorrectedPvalue)), ]
    #   test.matrix$Cell.group <- factor(test.matrix$Cell.group, levels = clusters, ordered = T)
    #   return(test.matrix)
    #   #return(gene.test)
    # }
}

VolcanoPlot <- function(fdrCut=0.01, fcCut=log2(2), ymax=5, labelsShow=labels){
    if(!is.installed("ggrepel")){
      warning("Detect package \"ggrepel\" is not installed in your R enviroment.")
      warning("Trying to install the \"ggrepel\" package.")
      warning("If failed, please try to install it mannually.")
      install.packages("ggrepel")
    }

    ## Load libraries
    suppressPackageStartupMessages(library(ggrepel))

    this.data <- read.table(infile, head=T, stringsAsFactors=F)
    up.num <- sum(this.data$UpDown == "Up")
    down.num <- sum(this.data$UpDown == "Down")
    this.data <- this.data[ ! (is.na(this.data$CorrectedPvalue) | is.infinite(log(this.data$FC.Group1_vs_Group2))), ]
    this.data$CorrectedPvalue <- -log10(this.data$CorrectedPvalue)
    this.data$FC.Group1_vs_Group2 <- log2(this.data$FC.Group1_vs_Group2)
    this.data$CorrectedPvalue[this.data$CorrectedPvalue >= ymax] <- ymax 
    
    theme.publication <- theme(panel.background = element_rect(fill = "white", colour = "black"),
                               axis.title=element_text(size=rel(1.5)), 
                               plot.title=element_text(hjust=0.5, size=rel(2)))

    p <- ggplot( this.data, aes(x=FC.Group1_vs_Group2, y=CorrectedPvalue)) +
        geom_point(aes(color=UpDown)) + 
        scale_color_manual( values=c("Up" = "red", "Down" = "blue", "Stable" = "grey"), 
                           breaks=c("Up", "Down", "Stable"),
                           labels=c(paste(up.num, "Up"), paste(down.num, "Down"), "Not Sig")) +
        xlab("log2 fold change") + ylab("-log10 corrected p value") +
        labs(title=outfile, color=paste(sum(up.num, down.num), "significant DEGs")) +
        geom_vline(xintercept=c(-fcCut, fcCut), color="grey", linetype="dotted") + 
        geom_hline(yintercept=-log10(fdrCut), color="grey", linetype="dotted")

    if(isTRUE(labelsShow)){
        p <- p + geom_text_repel(data = subset(this.data, UpDown != "Stable"), 
                        aes(label=Gene), size=5, box.padding = 0.25, point.padding = 0.3)
    }
    p + theme.publication
    ggsave(paste(outfile, "Volcano", sysDate, "pdf", sep="."), width=figure.width+2, height=figure.height)
    #pdf(file=paste(outfile, "Volcano", sysDate, "pdf", sep="."), width=figure.width, height=figure.height,
    #    onefile=FALSE)
    #plot(x=fc, y=fdr, pch=20, col=cols, ylim=c(0, ymax), xlab="log2 fold change", ylab="-log10 corrected p value", las=1,
    #     main=outfile) 
    #legend("topleft", legend=c(paste( sum(cols=="red"), "Up", sep=" "), paste( sum(cols=="blue"), "Down", sep=" ")), 
    #       pch=20, col=c("red", "blue"), 
    #       title=paste("Significant", sum(cols!="grey"), "DEGs", sep=" "), bty="n")
    #abline(v=c(-fcCut, fcCut), col="grey", lty="dotted")
    #abline(h=-log10(fdrCut), col="grey", lty="dotted")
    #invisible(dev.off())
}

# content for each mode
# Barplot
if ( modes == "barplot"){
    load(file=infile)
    cell.idx <- cellIndex()
    gene.idx <- geneIndex()
    Barplot(cell=cell.idx, gene=gene.idx)
}

# Heatmap
if (modes == "heatmap"){

}

# DEG
if (modes == "deg"){
    load(file=infile)
    cell.idx <- cellIndex()
    findDEG(cell.idx)
}

# Volcano 
if (modes == "volcano"){
    VolcanoPlot()   
}
