
#### This is a test script ####
library(tidyverse)
#Upload data from GitHub

#Set session in github file
demo_data<-read.csv("Demo/Data/demo_data.csv")

#Change dataset
demo_data |>
  mutate(start_date = as.Date(start_date, "%m/%d/%Y")) |>
  mutate(end_date = as.Date(end_date, "%m/%d/%Y")) |>
  mutate(days_of_cdiff = as.numeric(end_date - start_date)) |>
  identity() -> demo_data_2

write_csv(demo_data_2,"Demo/Data/demo_data_2.csv")

#Saving to other folders

demo_data_2 |>
  ggplot(data = _, aes(x = days_of_cdiff)) +
  geom_bar() +
  theme_bw() |>
  identity() -> demo_figure

ggsave(demo_figure, filename= "Demo/Figs/demo_figure.png")
