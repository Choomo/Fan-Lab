rm(list = ls())
gc()

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

permutation_test <- function(A, B, bg, ntimes){
  A_length <- A@strand@lengths
  bg_length <- bg@strand@lengths
  observed_counts <- length(countOverlaps(A, B)[countOverlaps(A, B) == 1])
  
  permu_test <- data.frame("index" = 1:ntimes)
  permu_test$random <- apply(X = permu_test, MARGIN = 1, FUN = function(x){
    bg_tmp <- bg[sample(1:bg_length, A_length),]
    tmp_counts <- length(countOverlaps(bg_tmp, B)[countOverlaps(bg_tmp, B) == 1])
    return(tmp_counts)
  })
  
  ttest_tmp <- t.test(x = permu_test$random, mu = observed_counts)
  out_list <- list("observed_counts" = observed_counts, "test_df" = permu_test, "t_test" = ttest_tmp)
  return(out_list)
}


HeLa_dms <- import.bw(con = "D:/Fan_Lab/Research-data/D1_KD_HeLa/rrbs/HeLa_D1_KD_dmc_filted.bw")
D1_binding <- import.bed(con = "D:/Fan_Lab/Research-data/Wang Jing data/DNMT1-eCLiP/GSM6733468_HeLa_Sample1.narrowPeak.bed")
CG_bg <- import.bw(con = "D:/Fan_Lab/Research-data/Genome_rf/hg38_UCSC/hg38_CG_strand.bw")
HeLa_events <- import.bed("D:/Fan_Lab/Research-data/Genome_rf/hg38_UCSC/hg38_HeLa_events_region.bed")

dms_D1_binding <- permutation_test(A = HeLa_dms, B = D1_binding, bg = CG_bg, ntimes = 100)
dms_D1_events <- permutation_test(A = HeLa_dms, B = HeLa_events, bg = CG_bg, ntimes = 100)

ggplot(data = dms_D1_binding$test_df, mapping = aes(x = random)) + 
  geom_histogram() + 
  geom_vline(xintercept = dms_D1_binding$observed_counts)







