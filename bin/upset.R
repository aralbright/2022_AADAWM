library("tidyr")
library("ggplot2")
library("plyr")
library("dplyr")
library("UpSetR")
library("readxl")

# load all excel files in folder as dataframe
setwd("/Volumes/albright_postdoc/tubulinRNAseq/bin/upset/tables")
temp <- list.files(pattern = ".xlsx")
df.list <- lapply(temp, read_excel)

# list of lists of genes in these dataframes, assuming gene name in first column
set.list <- sapply(df.list, "[", 0:1)

G04 <- read.table(file.path("/Volumes/albright_postdoc/2022_AADAWM_v3/3_dynein/df_G04_sig.csv"), header = TRUE, sep=",")
tub <- read.table(file.path("/Volumes/albright_postdoc/2022_AADAWM_v3/2_tubulin/df_tub_sig.csv"), header = TRUE, sep=",")
ap <- read.table(file.path("/Volumes/albright_postdoc/2022_AADAWM_v3/1_bulk_ap/ap_sig_df_qval0pt2.csv"), header = TRUE, sep=",")

G04_list <- as.list((G04$gene))
tub_list <- as.list((tub$gene))
ap_list <- as.list((ap$gene))

set.list <- append(set.list, list(G04_list, tub_list, ap_list))

temp <- gsub("\\..*", "", temp) 
names(set.list) <- temp

clusters <- names(set.list)

formatted_data <- fromList(set.list)
str(formatted_data)
any(is.na(formatted_data))

colnames(formatted_data)[9:11] <- c("G04", "tub", "ap")

# Basic upset plot to check if the function works
upset(formatted_data, sets = colnames(formatted_data))

upset(formatted_data, sets = colnames(formatted_data), order.by = "freq", nintersects = NA, 
      mainbar.y.label='Intersection Size',
      sets.x.label = 'Set Size',
      #text.scale = c(1.75, 1.3, 1.75, 1, 2, 1.3),
      point.size = 2.8,
      keep.order = TRUE, 
      line.size = 0.25,
      mb.ratio= c(0.5, 0.5))




# Create an expanded queries list
# Add more queries incrementally
queries_list <- list(
  list(query = intersects, params = list("Set_1"), color = "skyblue", active = TRUE),
  list(query = intersects, params = list("Set_2"), color = "purple", active = TRUE)
)


# Generate the UpSet plot with queries
upset(upset_data, 
      sets = names(filtered_wald_test_results_list),  # Set names
      order.by = "freq",                            # Order intersections by frequency
      mainbar.y.label = 'Intersection Size',
      sets.x.label = 'Set Size',
      keep.order = TRUE, 
      queries = queries_list                    # Apply queries for custom colors
)

# Add more queries incrementally
queries_list_step1 <- list(
  list(query = intersects, params = list("Set_1"), color = "red", active = TRUE),
  list(query = intersects, params = list("Set_2"), color = "green", active = TRUE)
)

# Plot with incremental queries
upset(upset_data, 
      sets = names(filtered_wald_test_results_list),  
      order.by = "freq",                            
      mainbar.y.label = 'Intersection Size',
      sets.x.label = 'Set Size',
      text.scale = c(1.75, 1.3, 1.75, 1, 2, 1.3),
      show.numbers = "no",
      point.size = 2.8,
      keep.order = TRUE,
      line.size = 0.25,
      mb.ratio = c(0.5, 0.5),                       
      queries = queries_list_step1
)