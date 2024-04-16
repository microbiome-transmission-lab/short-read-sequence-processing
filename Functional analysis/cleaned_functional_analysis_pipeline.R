#Cleaned functional metagenomic anlysis pipeline/sop
#from funprofiler (sequences in box) and existing metadata file (same as metaphlan)

library(tidyverse)
library(readxl)
library(janitor)
library(gt)
library(ggplot2)
library(viridis)

#### Read in files ####
# set directory to folder containing all raw profiles

# prevent scientific notation
options(scipen = 999)

# identify sequencing files
file_list <-  list.files(pattern = "profiles")

# create empty df to hold results
df_clean <- data.frame()

# loop through sequencing files

for (file in file_list) {
  df <- read_csv(file) %>%
    mutate(file_id = substr(file, 1, 32)) 
  df_clean <- rbind(df_clean, df)
  print(file)
}

df_clean |> # data is too large right now to put into repo so store in box
  write_csv("C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\funprofiler_ko11_all_raw.csv")

###########################
####Clean data for use ####
###########################

df_clean<-read_csv("C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\funprofiler_ko11_all_raw.csv")  #stored in box due to size


##### Pull out sample ID's and attach metadata #####

# Create string of sample IDs from existing metadata to pull them out
metadata_dat<-read_csv("data/metaphlan_w_parent.csv")  #Sync box folder for this data
  

#List of id's to pull from full filenames
metadata_dat |>
  mutate(sample_id = gsub("_", "-", sample_id)) |>
  pull(sample_id) |>
  paste0(collapse = "|") |>
  identity() -> fmt_id_search_string

df_clean |>
  mutate(sample_id = file_id) |>
  mutate(sample_id = gsub("_","-", sample_id)) |>
  mutate(sample_id = str_extract(string = sample_id,fmt_id_search_string )) |>
  mutate(read_pair = case_when(str_detect(file_id,"R1") ~ "R1",
                               str_detect(file_id,"R2") ~ "R2")) |> #This is relevant for read-paired data
  left_join(metadata_dat, by = "sample_id") |> 
  identity() -> clean_k11_ko_codes
unique(clean_k11_ko_codes$sample_id) #All subjects accounted for!

write_csv(clean_k11_ko_codes, "C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\ko11_codes_with_metadata.csv") #Stored in box due to size

#################################################
##### Average abundances between read pairs #####
#################################################

clean_k11_ko_codes<-read_csv("C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\ko11_codes_with_metadata.csv")

clean_k11_ko_codes |>
  group_by(external_participant_id,collection_date,sample_id,sample_type,sample_cat,read_pair) |>
  mutate(total_abundance = sum(abundance))|> #Adds up to 1 per read pair
select(total_abundance, read_pair)

clean_k11_ko_codes |> #Average read pair abundances
  group_by(ko_id, sample_id, immediate_parent, external_participant_id, collection_date, sample_cat, sample_type, original_sample_type, parent) |>
  mutate(avg_abundance = mean(abundance)) |>
  select(external_participant_id, collection_date, sample_id, ko_id,avg_abundance, immediate_parent, sample_cat, sample_type, original_sample_type, parent) |>
 # select(external_participant_id, collection_date,ko_id,avg_abundance,sample_id,sample_type, sample_cat) |>
  distinct() |>
  ungroup() |>
  identity() -> clean_k11_ko_avg_read

clean_k11_ko_avg_read |>
  group_by(sample_id, immediate_parent, external_participant_id, collection_date, sample_cat, sample_type, original_sample_type, parent) |>
  mutate(total_abundance = sum(avg_abundance)) #Adds up to 1.11 each ? *****

clean_k11_ko_avg_read |> #Convert new abundances into relative abundances
  group_by(sample_id, immediate_parent, external_participant_id, collection_date, sample_cat, sample_type, original_sample_type, parent) |>
  mutate(total_sum = sum(avg_abundance)) |>
  ungroup() |>
  group_by(sample_id,ko_id, immediate_parent, external_participant_id, collection_date, sample_cat, sample_type, original_sample_type, parent) |>
  mutate(avg_rel_abundance = (avg_abundance/total_sum)) |>
  ungroup() |>
  identity() -> clean_k11_ko_relative_abs

clean_k11_ko_relative_abs |>
  group_by(external_participant_id, collection_date,sample_id,sample_type, sample_cat) |>
  mutate(totals = sum(avg_rel_abundance)) #adds up to 1

###############################################
##### Add in metadata to use for analysis #####
###############################################

clean_k11_ko_relative_abs |>
  filter(sample_cat == "donor") |>
  mutate(sample_source = case_when(original_sample_type == "Stool Swab" | sample_type == "Stool Swab" ~ "Stool Source",
                                   original_sample_type != "Stool Swab" & sample_type != "Stool Swab" ~ "Product Source")) |> 
  identity()-> clean_k11_ko_codes_donor

clean_k11_ko_relative_abs |>
  filter(sample_cat == "recipient") |>
  mutate(sample_source = "recipient") |> 
  identity()-> clean_k11_ko_codes_recipient

clean_k11_ko_metadata_all<-rbind(clean_k11_ko_codes_donor, clean_k11_ko_codes_recipient)

unique(clean_k11_ko_metadata_all$sample_source)

write_csv(clean_k11_ko_metadata_all, "C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\clean_k11_ko_metadata_all.csv")

clean_k11_ko_metadata_all |>
  filter(sample_source != "Product Source") |>
  write_csv("C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\clean_k11_ko_metadata_fmt.csv")

clean_k11_ko_metadata_all |>
  filter(sample_cat == "donor") |>
  write_csv("C:\\Users\\emiliacw\\Box Sync\\Funprofiler large datasets\\clean_k11_ko_metadata_donor.csv")

clean_k11_ko_metadata_all |>
  filter(sample_cat == "donor") |>
  count(external_participant_id, sample_source) #all product and stool source accounted for!
