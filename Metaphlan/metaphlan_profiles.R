###This script is for cleaning metaphlan outputted profiles into a row with 
#each classification and its relative abundance (kingdom, phylum, etc.)

#This follows the process of saving the output from box into github under
# data/metaphlan/profiles

#Set up
library(tidyverse)
library(janitor)
library(gt)
library(ggplot2)
library(viridis)

#Read in data

#If all files are from the same batch, read in the text files using a loop
#
list.files("data/metaphlan/profiles_unknown") |>
  tibble(file_name = _) |>
  mutate(specimen_name = file_name) |>
  mutate(specimen_name = gsub("_R1_001.fastq.gz_profiled.txt", "", specimen_name)) |>
  mutate(long_file_name = list.files("data/metaphlan/profiles_unknown", full.names = TRUE)) |>
  mutate(mpa_dat = map(.x = long_file_name, .f = ~ read_tsv(.x, skip = 4) |> mutate_all(.funs = ~ as.character(.x)))) |>
  select(specimen_name, mpa_dat) |>
  unnest(cols = c(mpa_dat)) |>
  mutate(relative_abundance = as.numeric(relative_abundance) / 100) |>
  identity() -> all_mpa_files
all_mpa_files

all_mpa_files |>
  rename_all(.funs = ~ gsub("#","", .x)) |>
  janitor::clean_names() |>
  mutate(clade_name = replace(clade_name, clade_name == "UNCLASSIFIED", "t__UNCLASSIFIED"))|>
  filter(grepl("t__",clade_name))|>
  select(-additional_species) |>
  separate_wider_delim(clade_name,"|", names = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"), too_few = "align_start") |>
  select(-ncbi_tax_id) |>
  rename(specimen_id = specimen_name) |>
  mutate(Kingdom = gsub("k__|t__", "", Kingdom)) |>
  mutate(Phylum = gsub("p__", "", Phylum)) |>
  mutate(Class = gsub("c__", "", Class)) |>
  mutate(Order = gsub("o__", "", Order)) |>
  mutate(Family = gsub("f__", "", Family)) |>
  mutate(Genus = gsub("g__", "", Genus)) |>
  mutate(Species = gsub("s__", "", Species)) |>
  mutate(Strain = gsub("t__", "", Strain)) |>
  identity() -> all_mpa_cleaned
all_mpa_cleaned

write_csv(all_mpa_cleaned, "data/metaphlan_profiles_unknown_cleaned.csv")

#Can add this code to above block to classify alphanumeric outputs as 'Unknown'
# mutate(Class = case_when(str_detect(Class, "CFGB") ~ "Unknown", .default = as.character(Class))) |>
# mutate(Order = case_when(str_detect(Order, "OFGB") ~ "Unknown", .default = as.character(Order))) |>
# mutate(Family = case_when(str_detect(Family, "FGB") ~ "Unknown", .default = as.character(Family))) |>
# mutate(Genus = case_when(str_detect(Genus, "GGB") ~ "Unknown", .default = as.character(Genus))) |>
# mutate(Species = case_when(str_detect(Species, "SGB") ~ "Unknown", .default = as.character(Species))) |>


#######################
#####Reads Processed#####
#######################

list.files("data/metaphlan/profiles_unknown") |>
  tibble(file_name = _) |>
  mutate(specimen_name = file_name) |>
  mutate(specimen_name = gsub("_R1_001.fastq.gz_profiled.txt", "", specimen_name)) |>
  mutate(long_file_name = list.files("data/metaphlan/profiles_unknown", full.names = TRUE)) |> 
  mutate(mpa_dat = map(.x = long_file_name, .f = ~ read_lines(.x, skip = 2, n_max = 1))) |>
  select(specimen_name, mpa_dat) |>
  unnest(cols = c(mpa_dat)) |>
  mutate(reads_processed = as.numeric(gsub("#| reads processed","", mpa_dat))) |> 
  select(-mpa_dat) |> 
  identity() -> mpa_reads_processed
mpa_reads_processed



#######################
#####Visualization#####
#######################

#stacked barplot
all_mpa_cleaned |>
  select(Family, relative_abundance,specimen_id) |>
  mutate(Family = gsub("f__", "", Family)) |>
  mutate(Family = case_when(str_detect(Family, "FGB") ~ "Unknown", .default = as.character(Family))) |>
  ggplot(aes(fill = Family, y = relative_abundance, x = specimen_id)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(legend.position = "bottom") +
  scale_fill_viridis(discrete = T) +
  theme(axis.text.x = ggtext::element_markdown(color = "black", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.text = element_text(size = 5)) +
#  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  labs(x = "Specimen ID", fill = "Bacterial Family", y = "Relative Abundance") |>
  identity() -> music_metaphlan_family_composition

ggsave(music_metaphlan_family_composition, filename = "figs/music_metaphlan_family_comp.png", width = 23, height = 30, units = "cm")

