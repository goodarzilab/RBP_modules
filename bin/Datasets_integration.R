library(metap)
library(stringr)
library(colorspace)
library(gplots)
library(extrafont)
library(extrafontdb)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables for filenames with default values if necessary
bioID_filename <- "signed_log_pv_BioID.tsv"  # Default filename for BioID
eCLIP_filename <- "signed_log_pv_eCLIP.tsv"  # Default filename for eCLIP
perturbSeq_filename <- "signed_log_pv_Perturbseq.tsv"  # Default filename for PerturbSeq

# Parse command line arguments
for (arg in args) {
  key_value <- strsplit(arg, "=")[[1]]
  if (length(key_value) == 2) {
    if (key_value[1] == "BioID") {
      bioID_filename <- key_value[2]
    } else if (key_value[1] == "eCLIP") {
      eCLIP_filename <- key_value[2]
    } else if (key_value[1] == "PerturbSeq") {
      perturbSeq_filename <- key_value[2]
    }
  }
}


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

obs_pvs <- function(x){
    res = mean(x<x[1])
    for(i in x[2:length(x)]){
        res = c(res,mean(x<i))
    }
    res
}

##Obtaining signed logPvs
#BioID
bioid_z <- read.table(bioID_filename, sep='\t', header=TRUE, row.names = 1)
for(i in seq(from=1,to=nrow(bioid_z), by=1)){
    bioid_z[i,] <- z_norm(as.numeric(bioid_z[i,]))
    bioid_z[i,bioid_z[i,]<0] <- 0
}
bioid_z <- bioid_z[rowSums(bioid_z!=0)>0,]
bioid_cosine <- cos.dissim(t(bioid_z))

#eCLIP
eclip_z <- read.table(eCLIP_filename, sep='\t', header=TRUE, row.names = 1)
for(i in seq(from=1,to=nrow(eclip_z), by=1)){
    eclip_z[i,][!is.na(eclip_z[i,])] <- z_norm(as.numeric(eclip_z[i,][!is.na(eclip_z[i,])]))
    eclip_z[i,is.na(eclip_z[i,])] <- 0
    eclip_z[i,eclip_z[i,]<0] <- 0
}
eclip_z <- eclip_z[rowSums(eclip_z!=0)>0,]
eclip_cosine <- cos.dissim(t(eclip_z))

#Perturb-seq
perseq_z <- read.table(perturbSeq_filename, sep='\t', header=TRUE, row.names = 1)
for(i in seq(from=1,to=nrow(perseq_z), by=1)){
    perseq_z[i,] <- z_norm(as.numeric(perseq_z[i,]))
}
perseq_cosine <- cos.dissim(t(perseq_z))

#Making table for logitp
Rank_products_df <- data.frame(matrix(ncol=4, nrow=1))
colnames(Rank_products_df) <- c('Pair','BioID','eCLIP','PerturbSeq')
curr_names <- colnames(eclip_cosine)[order(colnames(eclip_cosine))]
for(i in curr_names){
    curr_partners = curr_names[curr_names<i]
    for(j in curr_partners){
        curr_name = paste(i,j,sep = '_')
        Rank_products_df <- rbind(Rank_products_df,data.frame('Pair'=curr_name, 'BioID'=NA, 'eCLIP'=NA, 'PerturbSeq'=NA))
        if(i%in%colnames(bioid_cosine) & j%in%colnames(bioid_cosine)){
            Rank_products_df[nrow(Rank_products_df),2] <- as.numeric(bioid_cosine[row.names(bioid_cosine)==j,colnames(bioid_cosine)==i])}
        if(i%in%colnames(eclip_cosine) & j%in%colnames(eclip_cosine)){
            Rank_products_df[nrow(Rank_products_df),3] <- as.numeric(eclip_cosine[row.names(eclip_cosine)==j,colnames(eclip_cosine)==i])}
        if(i%in%colnames(perseq_cosine) & j%in%colnames(perseq_cosine)){
            Rank_products_df[nrow(Rank_products_df),4] <- as.numeric(perseq_cosine[row.names(perseq_cosine)==j,colnames(perseq_cosine)==i])}
    }
}
Rank_products_df <- Rank_products_df[2:nrow(Rank_products_df),]
pseudo <- 0.000001
Rank_products_df$BioID[!is.na(Rank_products_df$BioID)] <- obs_pvs(Rank_products_df$BioID[!is.na(Rank_products_df$BioID)]) + pseudo
Rank_products_df$eCLIP[!is.na(Rank_products_df$eCLIP)] <- obs_pvs(Rank_products_df$eCLIP[!is.na(Rank_products_df$eCLIP)]) + pseudo
Rank_products_df$PerturbSeq[!is.na(Rank_products_df$PerturbSeq)] <- obs_pvs(Rank_products_df$PerturbSeq[!is.na(Rank_products_df$PerturbSeq)]) + pseudo

#Calculating logitp for all
RP_res_df_all <- Rank_products_df[,c(1,2)]
colnames(RP_res_df_all)[2] <- 'P-value'
for(i in RP_res_df_all$Pair){
    curr_pair = as.numeric(Rank_products_df[Rank_products_df$Pair==i,c(2:4)])
    if(sum(!is.na(curr_pair))==1){
        RP_res_df_all[RP_res_df_all$Pair==i,2] <- curr_pair[!is.na(curr_pair)]
    }else{
        RP_res_df_all[RP_res_df_all$Pair==i,2] <- logitp(p = curr_pair[!is.na(curr_pair)])$p
    }
}
logitp_res_matrix <- matrix(nrow=length(curr_names), ncol=length(curr_names))
row.names(logitp_res_matrix) <- curr_names
colnames(logitp_res_matrix) <- curr_names
for(i in colnames(logitp_res_matrix)){
    for(j in colnames(logitp_res_matrix)[colnames(logitp_res_matrix)<i]){
        logitp_res_matrix[colnames(logitp_res_matrix)==i,row.names(logitp_res_matrix)==j] <- RP_res_df_all[str_split_fixed(RP_res_df_all$Pair, pattern = '_', n=2)[,1]==i&str_split_fixed(RP_res_df_all$Pair, pattern = '_', n=2)[,2]==j,2]
        logitp_res_matrix[colnames(logitp_res_matrix)==j,row.names(logitp_res_matrix)==i] <- RP_res_df_all[str_split_fixed(RP_res_df_all$Pair, pattern = '_', n=2)[,1]==i&str_split_fixed(RP_res_df_all$Pair, pattern = '_', n=2)[,2]==j,2]
        logitp_res_matrix[colnames(logitp_res_matrix)==i,row.names(logitp_res_matrix)==i] <- 0
        logitp_res_matrix[colnames(logitp_res_matrix)==j,row.names(logitp_res_matrix)==j] <- 0
    }
}

#Clusterization and drawing
logitp_res_d <- as.dist(cos.dissim(logitp_res_matrix))
logitp_res_hclust <- hclust(logitp_res_d, method = 'ward.D2')
logitp_res_hclust_d <- as.dendrogram(logitp_res_hclust)

cairo_pdf("Logitp_heatmap.pdf", family="Lato Semibold", width = 75, height = 75)
heatmap.2(logitp_res_matrix, trace = "none", col=sequential_hcl(n=319,palette='Purples 3', rev = FALSE)[1:300], Rowv=logitp_res_hclust_d, Colv=logitp_res_hclust_d, density.info="none", dist = function(x) as.dist(cos.dissim(x)), hclustfun = function(x) hclust(x, method = 'ward.D2'), key.title = "none", key.xlab = 'Logitp distance', keysize=0.75, srtCol = 75, notecex = 1.5, lhei = c(1,6), margins = c(15,15), cexRow = 3, cexCol = 3, key.par = list(cex=5))
dev.off()
