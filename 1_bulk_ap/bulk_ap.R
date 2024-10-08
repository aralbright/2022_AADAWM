library("sleuth")
library("tidyr")
library("plyr")
library("dplyr")
library("sva")

# sample IDs from path name
sample_id <- dir(file.path("/Volumes/albright_postdoc/ap_rna/kallisto_output/"))

# full paths
kal_dirs <- file.path("/Volumes/albright_postdoc/ap_rna/kallisto_output",sample_id)

# data frame (sample, condition, batch, path) from pre-filled table 
s2c <- read.table(file.path("/Volumes/albright_postdoc/ap_rna/bin", "conditions.csv"), header = TRUE, sep=",")
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# plotting specific colors needs batch to be a factor, not a numeric
s2c$batch <- factor(s2c$batch)

# prep sleuth object  
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# check PCA by batch 
p1 <- plot_pca(so, color_by = 'batch')
print(p1)

# remove outlier 
s2c <- s2c[s2c$batch != 5,]

# re-prep sleuth object  
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# PCA again
p2 <- plot_pca(so, color_by = 'batch')
print(p2)

# get estimated counts for input in ComBat
estimated_counts <- sleuth_to_matrix(so, "obs_raw", "est_counts")

write.csv(estimated_counts, file = "/Volumes/albright_postdoc/2022_AADAWM_v3/raw_counts.csv", row.names = TRUE)


# Extract batch information and batch correct 
batch <- so$sample_to_covariates$batch

batch_corrected_counts <- ComBat(dat = estimated_counts, batch = batch)

# save batch corrected counts
write.csv(batch_corrected_counts, file = "/Volumes/albright_postdoc/2022_AADAWM_v3/bc_counts.csv", row.names = TRUE)



