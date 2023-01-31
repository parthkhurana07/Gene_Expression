#install packages
install.packages("dplyr")
install.packages("tidyverse")
BiocManager::install("GEOquery")
BiocManager::install("edgeR")

#load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

#read data
dat<-read_excel("GSE71562_E14R012_raw_counts (1).xlsx")
dat

#get metadata
gse<-getGEO(GEO = "GSE71562", GSEMatrix = TRUE)
gse
metadata<-pData(phenoData(gse[[1]]))
head(metadata)
metadata.modified<-metadata %>%
  select(1,9,12,19) %>%
  rename(Organism = organism_ch1) %>%
  rename(Time = characteristics_ch1.2) %>%
  mutate(Time = gsub("time: ","",Time))
  

head(dat)

#####MANIPULATION#####

#reshaping data
data.long<-dat %>%
  rename(Gene = 1) %>%
  gather(key = "Samples", value = "FPKM", -Gene) 
   

#join data.long and metadata.modified
data.long<-data.long %>%
  left_join(.,metadata.modified, by = c("Samples" = "description")) 
  
  
#explore data
data.long %>%
  filter(Gene == "aaeA" | Gene == "abgA") %>%
  group_by(Gene, Organism) %>%
  summarise(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM),
            sd_FPKM = sd(FPKM)) %>%
  head()


#####VISUALIZATION#####

#barPlot
data.long %>%
  filter(Gene == 'aaeA') %>%
  ggplot(.,aes(x = Samples, y = FPKM, fill = Time)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2")

#densityPlot
data.long %>%
  filter(Gene == "aaeA") %>%
  ggplot(.,aes(FPKM, fill = Time)) +
  geom_density(alpha = 0.3) +
  scale_fill_brewer(palette = "Dark2")

#scatterPlot
data.long %>%
  filter(Gene == "aaeA" | Gene == "abgA") %>%
  spread(key = Gene, value = FPKM) %>%
  ggplot(.,aes(x = aaeA, y = abgA, color = Time)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

#heatmap
genes.of.interest<-c("aaeA", "abgA", "aaeB", "abgB")
data.long %>%
  filter(Gene %in% genes.of.interest) %>%
  ggplot(.,aes(x = Samples, y = Gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')
  #barPlot
b<-data.long %>%
  filter(Gene == 'aaeA') %>%
  ggplot(.,aes(x = Samples, y = FPKM, fill = Time)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2")

#densityPlot
d<-data.long %>%
  filter(Gene == "aaeA") %>%
  ggplot(.,aes(FPKM, fill = Time)) +
  geom_density(alpha = 0.3) +
  scale_fill_brewer(palette = "Dark2")

#scatterPlot
s<-data.long %>%
  filter(Gene == "aaeA" | Gene == "abgA") %>%
  spread(key = Gene, value = FPKM) %>%
  ggplot(.,aes(x = aaeA, y = abgA, color = Time)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

#heatmap
genes.of.interest<-c("aaeA", "abgA", "aaeB", "abgB")
h<-data.long %>%
  filter(Gene %in% genes.of.interest) %>%
  ggplot(.,aes(x = Samples, y = Gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

#saving the plots
ggsave(b, filename = 'wullfen_barplot.pdf', width = 10,height = 8)
ggsave(d, filename = 'wullfen_densityplot.pdf', width = 10,height = 8)
ggsave(s, filename = 'wullfen_scatterplot.pdf', width = 10,height = 8)
ggsave(h, filename = 'wullfen_heatmapt.pdf', width = 10,height = 8)
