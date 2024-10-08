library("sleuth")
library("tidyr")
library("plyr")
library("dplyr")

# read in sample ids 
sample_id1 <- dir(file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20230815/con"))
sample_id2 <- dir(file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20230815/tub"))

# kallisto directories 
kal_dirs1 <- file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20230815/con",sample_id1)
kal_dirs2 <- file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20230815/tub",sample_id2)
  
kal_dirs <- c(kal_dirs1, kal_dirs2)

# sample to condition table
s2c <- read.table(file.path("/Volumes/albright_postdoc/tubulinRNAseq/bin", "conditions_ap.txt"), header = TRUE, sep="")
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# prep sleuth object, keep features with 5+ estimated counts in at least 14% of the samples
# this comes out to just 1 sample, which is quite permissive; however, 
# additional filters come later 

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, filter_fun=function(x){basic_filter(x, 5, 0.14)})

# if you want to look at the shiny
sleuth_live(so)

# create sleuth dataframe containing normalized estimated tpm
sleuth_matrix <- as.data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
head(sleuth_matrix)

write.csv(sleuth_matrix, file = "/Volumes/albright_postdoc/2022_AADAWM_v3/tub_norm_tpm.csv", row.names = TRUE)


# filter matrix for anterior and posterior samples
ap_df <- dplyr::select(sleuth_matrix, matches("^[a-z][a-z][a-z][1-5][A|P]"))

# log(1+x) transform
log_ap_df <- as.data.frame(lapply(ap_df, FUN = log1p))
rownames(log_ap_df) <- rownames(ap_df)

# funneled this into the jupyter notebook skew analysis
write.csv(log_ap_df, file = "/Volumes/albright_postdoc/tubulinRNAseq/bin/20230831_log_ap_df.csv")
