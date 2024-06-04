# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values for nperm, nproc, and input file
nperm <- 100000
nproc <- 40
input_file <- 'data/signed_log_pv_BioID.tsv'  # Default input file name

# Parse command line arguments
for (i in seq_along(args)) {
  arg <- strsplit(args[i], "=")[[1]]
  if (length(arg) == 2) {
    if (arg[1] == "nperm") {
      nperm <- as.numeric(arg[2])
    } else if (arg[1] == "nproc") {
      nproc <- as.numeric(arg[2])
    } else if (arg[1] == "input") {
      input_file <- arg[2]
    } else if (arg[1] == "condapath") {
      .libPaths( paste(arg[2],"/lib/R/library/",sep="") )
    }
  }
}

library(rlist)
library(fgsea)
library(gplots)
library(colorspace)
library(tidyverse)

cos.dissim <- function(x)
{
    outp <- matrix(nrow = length(x[,1]), ncol = length(x[,1]))
    row.names(outp) <- row.names(x)
    colnames(outp) <- row.names(x)
    for(i in seq(from=1, to=length(x[,1]), by=1)){
        for(j in seq(from=1, to=length(x[,1]), by=1)){
            A = x[i,]
            B = x[j,]
            outp[i,j] <- 1 - sum(A*B) / sqrt(sum(A^2) * sum(B^2))
        }
    }
    return(outp)
}

z_norm <- function(x)
{
    (x - mean(x)) / sd(x)
}

#Obtaining go sets, for reproducibility we will use preset data 
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#library(gage)
#go_gs <- go.gsets(species = "human", pkg.name='org.Hs.eg.db', id.type = "eg", keep.evidence=FALSE)
#go_bp=go_gs$go.sets[go_gs$go.subs$BP]
#go_cc=go_gs$go.sets[go_gs$go.subs$CC]
#go_mf=go_gs$go.sets[go_gs$go.subs$MF]
go_bp <- list.load('bin/go_bp.yml')
go_cc <- list.load('bin/go_cc.yml')
go_mf <- list.load('bin/go_mf.yml')

#Obtaining z-scores table with ENTREZID
signed_pv_scores <- read.table(input_file, sep='\t', header=TRUE, row.names = 1)
signed_pv_scores$Prots <- row.names(signed_pv_scores)
signed_pv_scores <- data.frame(separate_rows(signed_pv_scores, Prots, sep = ';'))

#Annotating genes by entrez ids, for reproducibility we will use preset data
#symb_to_ids <- data.frame(mapIds(org.Hs.eg.db, signed_pv_scores$Prots, "ENTREZID", "SYMBOL", multiVals = "filter"))
#symb_to_ids$SYMB <- names(mapIds(org.Hs.eg.db, signed_pv_scores$Prots, "ENTREZID", "SYMBOL", multiVals = "filter"))
#colnames(symb_to_ids) <- c('ENTREZ','SYMB')
symb_to_ids <- read.table('bin/Gene_Symbol_To_Entrez.tsv', sep='\t', header = T, stringsAsFactors = F)
signed_pv_scores <- merge(signed_pv_scores, symb_to_ids, by.x='Prots', by.y='SYMB', all.x=TRUE)
signed_pv_scores <- as.matrix(data.frame(signed_pv_scores[!is.na(signed_pv_scores$ENTREZ),2:(ncol(signed_pv_scores)-1)], row.names=signed_pv_scores$ENTREZ[!is.na(signed_pv_scores$ENTREZ)]))

z_scores_signed_pv <- signed_pv_scores
for(i in seq(from=1,to=nrow(z_scores_signed_pv), by=1)){
    z_scores_signed_pv[i,] <- z_norm(z_scores_signed_pv[i,])
}
z_scores_signed_pv <- data.frame(z_scores_signed_pv)
z_scores_signed_pv$ID <- row.names(z_scores_signed_pv)
z_scores_signed_pv <- z_scores_signed_pv[,c((ncol(z_scores_signed_pv)),1:(ncol(z_scores_signed_pv)-1))]

#Running fgsea
set.seed(13)

#GO BP
fgseaResList_signed_pv_BP <- list()
for(i in seq(from=2,to=ncol(z_scores_signed_pv),by=1)){
    print(i)
    fgseaRes_BP <- fgsea(pathways=go_bp, stats=deframe(z_scores_signed_pv[,c(1,i)])[order(deframe(z_scores_signed_pv[,c(1,i)]), decreasing = FALSE)], nperm=nperm, nproc=nproc)
    fgseaResList_signed_pv_BP[[(i-1)]] <- fgseaRes_BP
}

#GO MF
fgseaResList_signed_pv_MF <- list()
for(i in seq(from=2,to=ncol(z_scores_signed_pv),by=1)){
    print(i)
    fgseaRes_MF <- fgseaSimple(pathways=go_mf, stats=deframe(z_scores_signed_pv[,c(1,i)])[order(deframe(z_scores_signed_pv[,c(1,i)]), decreasing = FALSE)], nperm=nperm, nproc=nproc)
    fgseaResList_signed_pv_MF[[(i-1)]] <- fgseaRes_MF
}

#GO CC
fgseaResList_signed_pv_CC <- list()
for(i in seq(from=2,to=ncol(z_scores_signed_pv),by=1)){
    print(i)
    fgseaRes_CC <- fgseaSimple(pathways=go_cc, stats=deframe(z_scores_signed_pv[,c(1,i)])[order(deframe(z_scores_signed_pv[,c(1,i)]), decreasing = FALSE)], nperm=nperm, nproc=nproc)
    fgseaResList_signed_pv_CC[[(i-1)]] <- fgseaRes_CC
}

#Obtaining tables for the BioID-based annotation heatmap (BP)
#Table containing NES
fgseaResList_signed_pv_BP_NES <- data.frame('Path'=fgseaResList_signed_pv_BP[[1]][,1], 'Prot_name'=fgseaResList_signed_pv_BP[[1]][,5])
colnames(fgseaResList_signed_pv_BP_NES)[2] <- colnames(z_scores_signed_pv)[2]
for(i in seq(from=2, to=length(fgseaResList_signed_pv_BP), by=1)){
    fgseaResList_signed_pv_BP_NES <- merge(fgseaResList_signed_pv_BP_NES, data.frame('Path'=fgseaResList_signed_pv_BP[[i]][,1], 'Prot_name'=fgseaResList_signed_pv_BP[[i]][,5]), by = 'pathway')
    colnames(fgseaResList_signed_pv_BP_NES)[i+1] <- colnames(z_scores_signed_pv)[i+1]
}

fgseaResList_signed_pv_BP_NES_filtered <- fgseaResList_signed_pv_BP_NES[rowSums(fgseaResList_signed_pv_BP_NES[,2:ncol(fgseaResList_signed_pv_BP_NES)]>=2, na.rm = TRUE)>0,]
for(i in seq(from=2,to=ncol(fgseaResList_signed_pv_BP_NES_filtered), by=1)){
    fgseaResList_signed_pv_BP_NES_filtered[is.na(fgseaResList_signed_pv_BP_NES_filtered[,i])|fgseaResList_signed_pv_BP_NES_filtered[,i]<0,i] <- 0
}

#Table containing adjusted p-values
fgseaResList_signed_pv_BP_Padj <- data.frame('Path'=fgseaResList_signed_pv_BP[[1]][,1], 'Prot_name'=fgseaResList_signed_pv_BP[[1]][,3])
colnames(fgseaResList_signed_pv_BP_Padj)[2] <- colnames(z_scores_signed_pv)[2]
for(i in seq(from=2, to=length(fgseaResList_signed_pv_BP), by=1)){
    fgseaResList_signed_pv_BP_Padj <- merge(fgseaResList_signed_pv_BP_Padj, data.frame('Path'=fgseaResList_signed_pv_BP[[i]][,1], 'Prot_name'=fgseaResList_signed_pv_BP[[i]][,3]), by = 'pathway')
    colnames(fgseaResList_signed_pv_BP_Padj)[i+1] <- colnames(z_scores_signed_pv)[i+1]
}

fgseaResList_signed_pv_BP_Padj_filtered <- fgseaResList_signed_pv_BP_Padj[fgseaResList_signed_pv_BP_Padj$pathway%in%fgseaResList_signed_pv_BP_NES_filtered$pathway, ]
for(i in seq(from=2,to=ncol(fgseaResList_signed_pv_BP_Padj_filtered), by=1)){
    fgseaResList_signed_pv_BP_Padj_filtered[as.numeric(fgseaResList_signed_pv_BP_Padj_filtered[,i])<=0.05&as.numeric(fgseaResList_signed_pv_BP_NES_filtered[,i])>0&!is.na(fgseaResList_signed_pv_BP_Padj_filtered[,i]),i] <- '*'
    fgseaResList_signed_pv_BP_Padj_filtered[fgseaResList_signed_pv_BP_Padj_filtered[,i]!='*'|is.na(fgseaResList_signed_pv_BP_Padj_filtered[,i]),i] <- ''
}

#Biclusterization and drawing
fgseaResList_signed_pv_BP_d <- as.dist(cos.dissim(t(as.matrix(data.frame(fgseaResList_signed_pv_BP_NES_filtered[,2:ncol(fgseaResList_signed_pv_BP_NES_filtered)], row.names = fgseaResList_signed_pv_BP_NES_filtered[,1])))))
fgseaResList_signed_pv_BP_hclust <- hclust(fgseaResList_signed_pv_BP_d, method = 'ward.D2')
fgseaResList_signed_pv_BP_hclust_d <- as.dendrogram(fgseaResList_signed_pv_BP_hclust)

pdf('BioID_based_annotation_BP.pdf', width = 75, height=40)
heatmap.2(t(as.matrix(data.frame(fgseaResList_signed_pv_BP_NES_filtered[,2:ncol(fgseaResList_signed_pv_BP_NES_filtered)], row.names = fgseaResList_signed_pv_BP_NES_filtered[,1]))), trace = "none", col=sequential_hcl(n=299,palette='Blues 3', rev = TRUE), density.info="none", dist = function(x) as.dist(cos.dissim(x)), hclustfun = function(x) hclust(x, method = 'ward.D2'), key.title = "", key.xlab = 'NES', keysize=0.75, key.par = list(cex=2), srtCol = 75, margins = c(90,15), cexRow = 2, cexCol = 2, breaks = seq(0,2,length.out=300), cellnote = t(as.matrix(data.frame(fgseaResList_signed_pv_BP_Padj_filtered[,2:ncol(fgseaResList_signed_pv_BP_Padj_filtered)], row.names = fgseaResList_signed_pv_BP_Padj_filtered[,1]))), notecol = "red", notecex = 1.5, lhei = c(1,7), Rowv = fgseaResList_signed_pv_BP_hclust_d)
dev.off()
