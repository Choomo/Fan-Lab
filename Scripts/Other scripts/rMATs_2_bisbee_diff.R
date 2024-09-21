rm(list = ls())
library(data.table)

input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]
output_file_2 <- commandArgs(trailingOnly = TRUE)[3]
control <- commandArgs(trailingOnly = TRUE)[4]
test <- commandArgs(trailingOnly = TRUE)[5]

in_data <- as.data.frame(fread(input_file))
colnames(in_data)[which(colnames(in_data) == "strand")] <- "Strand"

################output diff######################
out_data <- in_data
out_data$gene <- in_data$GeneID
out_data$strand <- in_data$Strand
out_data$contig <- in_data$chr
out_data$event_id <- paste0("exon_skip_",c(1:nrow(in_data)))
out_data$event_jid <- paste0(in_data$chr, "s", in_data$strand, ":g.", 
                             in_data$upstreamEE, "j", in_data$exonStart_0base+1,"_",
                             in_data$exonEnd, "j", in_data$downstreamES+1,">",
                             in_data$upstreamEE, "j",in_data$downstreamES+1,"[splES]")
out_data$PSI_higher <- NA
out_data$PSI_higher[which(out_data$IncLevelDifference < 0)] <- test
out_data$PSI_higher[which(out_data$IncLevelDifference > 0)] <- control
out_data$PSI_lower <- NA
out_data$PSI_lower[which(out_data$IncLevelDifference < 0)] <- control
out_data$PSI_lower[which(out_data$IncLevelDifference > 0)] <- test
out_data$ll_ratio <- 8

psi_df_control <- as.data.frame(t(apply(X = in_data, MARGIN = 1, FUN = function(x){
  tmp <- as.character(x)[21]
  tmp <- strsplit(x = tmp, split = ",", fixed = T)[[1]]
  return(tmp)
})))
psi_df_test <- as.data.frame(t(apply(X = in_data, MARGIN = 1, FUN = function(x){
  tmp <- as.character(x)[22]
  tmp <- strsplit(x = tmp, split = ",", fixed = T)[[1]]
  return(tmp)
})))
psi_df <- cbind.data.frame(psi_df_control, psi_df_test)

colnames(psi_df) <- c(paste0("psi_control", c(1:ncol(psi_df_control))), 
                      paste0("psi_test", c(1:ncol(psi_df_test))))
out_data <- cbind.data.frame(out_data, psi_df)
depth_df <- as.data.frame(matrix(100, nrow = nrow(in_data), ncol = ncol(psi_df)))
colnames(depth_df) <- c(paste0("depth_control", c(1:ncol(psi_df_control))), 
                      paste0("depth_test", c(1:ncol(psi_df_test))))
out_data <- cbind.data.frame(out_data, depth_df)
out_data$fitPSI_control <- 0.1
out_data$fitPSI_test <- 0.1
out_data$fitPSI_all <- 0.1
out_data$fitW_2group <- 0.1
out_data$fitW_all <- 0.1
out_data$event_type <- "exon_skip"

index <- seq(which(colnames(out_data) == "gene"), ncol(out_data))
write.csv(x = out_data[,index], file = output_file, quote = F, sep = ",", row.names = T, col.names = T)

################output ######################

out_data_2 <- in_data
out_data_2$gene <- in_data$GeneID
out_data_2$strand <- in_data$Strand
out_data_2$contig <- in_data$chr
out_data_2$event_id <- paste0("exon_skip_",c(1:nrow(in_data)))
out_data_2$confirmed <- 6
out_data_2$exon_pre_end <- in_data$upstreamEE
out_data_2$exon_start <- in_data$exonStart_0base+1
out_data_2$exon_end <- in_data$exonEnd
out_data_2$exon_aft_start <- in_data$downstreamES+1
out_data_2$event_jid <- out_data$event_jid
out_data_2$ll_ratio <- out_data$ll_ratio
tmp_psi <- out_data[,which(startsWith(x = colnames(out_data), prefix = "psi")),]
tmp_depth <- out_data[,which(startsWith(x = colnames(out_data), prefix = "depth")),]
tmp_fit <- out_data[,which(startsWith(x = colnames(out_data), prefix = "fit")),]
out_data_2 <- cbind.data.frame(out_data_2, tmp_psi, tmp_depth, tmp_fit)
index_n <- seq(which(colnames(out_data_2) == "gene"), ncol(out_data_2))
write.csv(x = out_data_2[,index], file = output_file_2, quote = F, sep = ",", row.names = T, col.names = T)















