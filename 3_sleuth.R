library("sleuth")
library("dplyr")

sample_id <- dir(file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20221118"))

kal_dirs <- file.path("/Volumes/albright_postdoc/tubulinRNAseq/output_20221118",sample_id)

# s2c = samples to conditions
s2c <- read.table(file.path("/Volumes/albright_postdoc/tubulinRNAseq/bin/", "conditions_all.txt"), header = TRUE, sep="")
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# so = sleuth object
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, filter_fun=function(x){basic_filter(x, 5, 0.14)})

# if you want to look at the shiny
sleuth_live(so)

# create sleuth dataframe containing normalized estimated tpm
sleuth_matrix <- as.data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
head(sleuth_matrix)

# filter matrix for anterior and posterior samples
ap_df <- dplyr::select(sleuth_matrix, matches("^[a-z][a-z][a-z][1-5][A|P]"))

# log(1+x) transform
log_ap_df <- as.data.frame(lapply(ap_df, FUN = log1p))

# funneled this into the jupyter notebook skew analysis
write.csv(log_ap_df, file = "/Volumes/albright_postdoc/tubulinRNAseq/bin/log_ap_df.csv")