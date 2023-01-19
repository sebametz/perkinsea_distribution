##############################################################################################################
# Script: script_analysis.R
# Created: Dec 2022
# Author: Metz Sebastian
# Description: R script with all the analysis from "Perkinsea: a global perspective of environmental distribution
# and diversity explored by a meta-analysis of eDNA surveys"  
# NOTE: for more information refere to manuscript. 
##############################################################################################################

rm(list = ls(all.names = TRUE)) # Clean environment
gc() # free up memory

## Libreries ----
library(tidyverse) # data manipulation
library(vroom) # open files quickly
library(nVennR) # venn diagram
library(vegan) # ecology analysis
library(maps) # for world map
library(scatterpie) # pie chart for global map
library(circlize) # chord diagram
library(ComplexHeatmap) # to create labels for chord diagram
library(picante) # phylogenetic index analysis
library(phyloseq) # data manipulation
library("pheatmap")


## First part: Preprocessing of raw data for analysis of Perkinsea distribution ----
# NOTE 1: All the data is available in GitHub, for execution control the directories

### 1) get abundance table = ASVs as row and samples as columns
# We extracted the ASVs from two different tables from EukBank, one with Perkinsea ASVs and the other off 
# unknown Alveolata therefore we need to check both tables for the ASVs that we have and create only one table

# uploading ASVs from the Phylogenetic classification
taxonomy <- vroom("data/PERKINSEA_TABLE/taxonomy_per_query.tsv", delim = "\t")

# uploading Perkinsea table
perkinsea <- vroom("data/EUKBANK_TABLES/EukBank_2020_10_27_Perkinsea/eukbank_18SV4_asv.table", num_threads = 8)

# uploading unknown Alveolata table
alveolata <- vroom("data/EUKBANK_TABLES/EukBank_2020_10_27_unknown_alveolata/eukbank_18SV4_asv.table")

# creating a new table
# amplicon format of the tables don't have the ;size=XX as we have in the taxonomy file in name columns
perkinsea_new_table <- tibble(amplicon = str_extract(taxonomy$name, "\\w+")) %>%
  left_join(add_case(perkinsea, alveolata), by = "amplicon") %>% # merge Perkinsea and alveolata and filter ASVs
  select(amplicon, where( ~ is.numeric(.x) && sum(.x) != 0)) # select the samples with Perkinsea reads

### 2) get the metadata file 
eukbank_meta <- read_tsv("data/EUKBANK_TABLES/eukbank_18SV4_asv.metadata_28062021")
perkinsea_meta <- filter(eukbank_meta, sample %in% colnames(perkinsea_new_table))

# there are 758 samples that don't have metadata, most of them are repeated samples, 
# therefore we will delete them

perkinsea_new_table_filtered <- perkinsea_new_table[,colnames(perkinsea_new_table) %in% c("amplicon", perkinsea_meta$sample)]

# save tables
write_tsv(perkinsea_new_table_filtered, "data/PERKINSEA_TABLES/Perkinsea_abundance.table")
write_tsv(perkinsea_meta, "data/PERKINSEA_TABLES/Perkinsea_metadata.table")

### 3) prepare taxonomy table
ncol <- max(str_count(taxonomy$taxopath, ";")) + 1
colnm <- str_c("TAX", 1:ncol)
  
taxonomyF <- taxonomy %>%
  separate(col = taxopath, sep = ";", into = colnm, remove = F)

colnames(taxonomyF)[1] <- "amplicon"
write_tsv(taxonomyF, "data/PERKINSEA_TABLES/Perkinsea_taxonomy.table")

# The table of samples were then manually revised and the samples were 
# classified according to their environmental type

### 4) get rarefied tables <- from the eukbank database
table <- read_tsv("data/EUKBANK_TABLES/Perkinsea_18SV4_rareffied_10K.table")


rm(list = ls(all.names = TRUE)) # Clean environment
gc() # free up memory
## Second part: Analysis of Perkinsea data ----
# NOTE 2: All the figures were modified for a final version of the paper using Inkscape (https://inkscape.org/)

# Read tables:
abundance <- vroom("data/PERKINSEA_TABLES/Perkinsea_abundance.table", delim = "\t")
taxonomy <- vroom("data/PERKINSEA_TABLES/Perkinsea_taxonomy.table", delim = "\t")
samples <- vroom("data/PERKINSEA_TABLES/Perkinsea_metadata_reduced.table", delim = "\t")
rarefy <- vroom("data/PERKINSEA_TABLES/Perkinsea_filtered_rarrefy_10k.table")

## Some corrections, colors definition and theme ----
# Those ASVs that were not assigned to any cluster or to Perkinsea_X were classified as 'unclassified Perkinsea'
taxo <- taxonomy %>% 
  mutate(amplicon = str_extract(amplicon, "[a-z0-9]+"), 
         tax = ifelse(is.na(TAX5) | TAX5 == "Perkinsea_X", "unclassified Perkinsea", TAX5))

# Environments colors
env_colours <- c("#66CCEE", "#228833", "#CCBB44") # Environment colors
names(env_colours) <- c("Marine", "Land water", "Soil")

# Taxonomy colors
tax_colours <- c("#77AADD", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", 
                 "#864086","#EEDD88", "#EE8866", "#FFAABB","#DDDDDD")

names(tax_colours) <- c("NAG01", "PERK01", "PERK02", "PERK03", "PERK04",
                        "Pararosarium","Parviluciferaceae", "Perkinsidae", 
                        "Xcellidae","unclassified Perkinsea")


# Environments sub-categories colours

sub_colours <- c("#a0f4f5", "#7BAFDE", "#5289C7", "#1965B0", "#882E72", "#AE76A3",
                 "#D1BBD7","gray70","#90C987","#CAE0AB","#F7F056", 
                 "#e769f8","#F6C141","#0dc923")

names(sub_colours) <- c("Coastal zone", "Epipelagic zone", "Mesopelagic zone",
                        "Bathypelagic zone","Abyssal zone", "Pelagic zone", 
                        "Marine sediment","Other marine environment", 
                        "Lake", "River", "Freshwater sediment","Other land water",
                        "Temperate soil", "Tropical soil")

# I defined this theme just to homogenise the output of all the plots

tm <- theme(panel.border = element_blank(), 
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(color = "black", fill = "white"),
          text = element_text(family = "serif", size = 12, colour = "black"), 
          title = element_text(family = "serif", size = 14, colour = "black"))

### Figure 1A: Venn diagram - Shared ASVs between environments (Main and subcategories) ----

venn <- list(Marine = abundance$amplicon[ifelse(rowSums(abundance[colnames(abundance) %in% samples$sample[samples$environment == "Marine"]])>0, 1,0)>0],
             `Land water`= abundance$amplicon[ifelse(rowSums(abundance[colnames(abundance) %in% samples$sample[samples$environment == "Land water"]])>0, 1,0)>0],
             Soil = abundance$amplicon[ifelse(rowSums(abundance[colnames(abundance) %in% samples$sample[samples$environment == "Soil"]])>0, 1,0)>0])

# Venn Diagram plot using the library nVennR
venn_plot <- plotVenn(venn, setColors = env_colours, borderWidth = 2,
                      opacity = 0.4, outFile = "FIGURES/RAW/Figure1A.svg")

#Shared ASVss between the three environments
shared <- tibble(
  amplicon =intersect(intersect(venn$Marine,venn$`Land water`), venn$Soil)) %>%
  left_join(select(taxo, amplicon, tax))

shared %>%
  group_by(tax) %>%
  summarise(value = n())

# # A tibble: 4 × 2
# tax                    value
# <chr>                  <int>
# 1 NAG01                     16
# 2 PERK01                     2
# 3 PERK02                     8
# 4 unclassified Perkinsea    13

### Figure 1B: NMDS of selected samples representatives from each environment ----

# Extract ASVs with at least one read after rarefaction
aux_r <- rarefy %>%
  select(where( ~ is.numeric(.x) && sum(.x) != 0))

aux_r <- t(aux_r)

div <- tibble(sample = rownames(aux_r), 
              shannon = diversity(aux_r, index = "shannon"), 
              simpson = diversity(aux_r, index = "simpson"),
              inv = diversity(aux_r, index = "inv"),
              richness = specnumber(aux_r))

div_f <- left_join(div, select(samples, sample, environment, subcategory))

div_fn <- div_f 
samp <- NULL
for (i in unique(div_fn$subcategory)) {
  aux_s <- div_fn %>% 
    filter(subcategory == i)
  aux_s <- aux_s[order(aux_s$shannon, decreasing = T),]
  if(nrow(aux_s)>=15){aux_s <- aux_s[1:15,]}
  if(is.null(samp)){
    samp <- aux_s
  }else{samp <- add_case(samp, aux_s)}
}

# Table with selected samples
sampF <- read_tsv("data/PERKINSEA_TABLES/samples_NMDS.txt")

aux_rf <- aux_r[rownames(aux_r) %in% sampF$sample,]
aux_rf <- aux_rf[, -grep(T, colSums(aux_rf) == 0)]

table.h <- decostand(aux_rf,  method = "hellinger")
nmds <- metaMDS(table.h, distance = "bray", trymax = 500, parallel = 8, k = 2)
nmds$stress

aux_plot <- tibble(sample = names(nmds$points[,1]), NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2])
aux_plot <- aux_plot %>%
  left_join(select(samples, sample, environment, subcategory))

p <- ggplot(aux_plot, aes(x = NMDS1, y = NMDS2, fill = environment)) +
  geom_point(size = 4, colour = "black", pch = 21) + 
  scale_fill_manual(values = env_colours)+
  labs(fill = "Environment")+
  annotate("text", x = -4.3, y = 2.4, 
           label = str_c("Stress=",round(nmds$stress, digits = 2))) +
  theme_bw()+
  tm
p

pdf("FIGURES/RAW/Figure1B.pdf", height = 6, width = 8)
 p
dev.off()

### Figure 2: Map of distribution ----

# Create a grill of 10X10
xy <- select(samples, c("sample", "latitude", "longitude", "environment"))
xy$grid <- NA

for (i in unique(xy$environment)) {
  grid <- 0  
  for(x in seq(from = -180.0, to = 180.0, by = 10)){
    for (y in seq(from = -90.0, to = 90.0, by = 10)) {
      grid <- grid + 1
      xmin <- x
      xmax <- x + 9.9999999 # longitude max
      ymin <- y
      ymax <- y + 9.9999999 # latitude max
      xy$grid[between(xy$latitude, ymin, ymax) & between(xy$longitude, xmin, xmax) & xy$environment == i] <- str_c(i, "_", grid)
    }
  }
}
xy <- filter(xy, !is.na(grid))

# Now we collapse all the samples that are in the same square in XY
aux <- xy %>%
  filter(!is.na(grid)) %>%
  group_by(grid) %>%
  summarise(medianLat = median(latitude), medianLong = median(longitude), samples = n())

table_aux_wider <- NULL

# for echa squard in aux we search all the samples in the squard and collapse the amplicons
for (i in aux$grid) {
  table_aux <- select(abundance, "amplicon", xy$sample[xy$grid == i])
  
  table_aux$abundance <- rowSums(table_aux[,2:ncol(table_aux)])
  
  aux_t <- left_join(select(table_aux, c("amplicon", "abundance")), 
                     select(taxo, c("amplicon", "tax"))) %>%
    group_by(tax) %>% 
    summarise(value = sum(abundance>0)) %>% # number of amplicons
    select(tax, value) %>%
    pivot_wider(names_from = tax, values_from = value) %>%
    mutate(grid = i)
  
  if(is.null(table_aux_wider)){
    table_aux_wider <- aux_t
  }else{
    table_aux_wider <- table_aux_wider %>%
      add_case(aux_t)
  }
}

aux_f <- left_join(aux, table_aux_wider) %>% 
  mutate(environment = str_extract(grid, "([A-Za-z ]+)"))

#### World map plotting ----
world <- map_data("world")

for_map <- aux_f
for_map$N <- rowSums(for_map[,5:14])
for_map$radius <- log(for_map$N+1)
for_map$size <- for_map$radius/3
for_map <- for_map[order(for_map$samples, decreasing = T),]

p <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="gray90", 
               color = "gray75", size = 0.1, alpha=0.5) +
  geom_scatterpie(data = for_map, aes(y=medianLat, x=medianLong, r = radius, group = grid),                   
                  cols = colnames(for_map)[5:14], alpha=1, size = 0.5) + #pie_scale = 0.5
  geom_scatterpie_legend(for_map$radius, x=-140, y=-80, n = 3) +
  scale_fill_manual(values = tax_colours, name = "Taxonomy", labels = names(tax_colours)) +
  scale_y_continuous(name = "Latitude", breaks = seq(from = -90.0, to = 90.0, by = 30)) +
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(from = -180.0, to = 180.0, by = 40)) +
  xlab("Longitude") + ylab("Latitude") +
  coord_equal() +
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", colour = "white",fill = "white"),
        text = element_text(family = "serif", size = 10, colour = "black"),
        legend.background = element_blank(), 
        legend.box.background = element_blank())

pdf("FIGURES/RAW/Figure2.pdf", height = 12, width = 18)
 p
dev.off()

### Figure 3: Environmental distribution ----
# Formating table
aux_chord <- NULL
for (i in unique(samples$subcategory)) {
  auxfor <- abundance %>%
    select("amplicon", samples$sample[samples$subcategory %in% i]) %>%
    mutate(abundance = rowSums(.[,2:ncol(.)])) %>%
    left_join(select(taxo, amplicon, tax), by = "amplicon")  %>%
    group_by(tax) %>%
    summarise(subenv = i, ASVs = sum(abundance>0))
  if(is.null(aux_chord)){
    aux_chord <- auxfor
  } else {aux_chord <- add_case(aux_chord, auxfor)}
}

aux_chordf <- pivot_wider(data = aux_chord, names_from = tax, values_from = ASVs)

# Relative percentaje
for (i in 2:ncol(aux_chordf)) {
  aux_chordf[,i] <- aux_chordf[,i]/sum(aux_chordf[,i])*100
}

# Plot
m <- as.data.frame(aux_chordf[,-c(1,10)])
rownames(m) <- aux_chordf$subenv
m <- m[c(14,1,4,9,7,6,12,10,2,8,13,11,5,3), ]
m <- m[, c(1,4,5,6,7,2,3,8,9)]
# Legends
#
chord_colors <- c(tax_colours[-10],sub_colours)

col_fun <- function(x){sub_colours}
col_fun2 <- function(x){tax_colours}
lgd_category <- Legend(at = rownames(m), title_position = "topleft",
                       title = "Sub-environment", legend_gp = gpar(fill = col_fun(at)),
                       ncol = 1)
lgd_tax <- Legend(at = colnames(m),  title_position = "topleft",
                  title = "Taxonomy", legend_gp = gpar(fill = col_fun2(at)), ncol = 2)



circos.par(clock.wise = T)
chordDiagram(as.matrix(m), order = c(rev(colnames(m)), rev(rownames(m))),
             grid.col = chord_colors,
             grid.border = "black",
             annotationTrack =  c("grid"),
             annotationTrackHeight = mm_h(3))

draw(lgd_category, x = unit(0.09, "npc"), y = unit(0.94, "npc"),
     just = c("left", "top"))

draw(lgd_tax, x =  unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), 
     just = c("right", "bottom"))

circos.clear()
dev.copy(pdf,'FIGURES/RAW/Figure3.pdf', width=12, height=8)
dev.off()

# small one for Perkinsea unclassified
m <- tibble(aux_chordf[,10])
rownames(m) <- aux_chordf$subenv
mf <- m[c(14,1,4,9,7,6,12,10,2,8,13,11,5,3), ]
rownames(mf)<- aux_chordf$subenv[c(14,1,4,9,7,6,12,10,2,8,13,11,5,3)]
#
chord_colors <- c(tax_colours[10],sub_colours)
circos.par(clock.wise = T)
chordDiagram(as.matrix(mf), order = c(rev(colnames(mf)), rev(rownames(mf))),
             grid.col = chord_colors,
             grid.border = "black",
             annotationTrack =  c("grid"),
             annotationTrackHeight = mm_h(3))
circos.clear()
dev.copy(pdf,'FIGURES/RAW/Figure3bis.pdf', width=12, height=8)
dev.off()

### Figure 4: Novelty analysis ----
lwr <- read_csv("data/PERKINSEA_TABLES/lwr_list.csv")
similarity <- read_tsv("data/PERKINSEA_TABLES/Perkinsea_ASVs.blastn", col_names = F)

aux <- select(similarity, X1, X3) %>%
  left_join(select(lwr, Pquery, `LWR 1`), by = c("X1" = "Pquery")) %>%
  mutate(amplicon =  str_extract(X1, "\\w+")) %>%
  left_join(select(taxo, amplicon, LWR, tax), by = c("amplicon" = "amplicon"))
names(aux) <- c("amplicon", "similarity", "LWR", "id", "LWR_classification", "taxfact")

p <- ggplot(aux, aes(x = similarity, y = LWR))  +
  geom_point(aes(fill = taxfact), pch = 21, size = 2)+
  scale_fill_manual(values = tax_colours) +
  geom_hline(yintercept = mean(aux$LWR), linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = mean(aux$similarity), linetype = "dashed", color = "gray30") +
  facet_wrap(.~taxfact, nrow = 4) +
  xlab(label = "% of similarity") + ylab(label = "Likelihood Weight Ratio (LWR)") + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(color = "black", colour = "white",fill = "white"),
        legend.position = "none",
        text = element_text(family = "serif", size = 12, colour = "black"), 
        title = element_text(family = "serif", size = 14, colour = "black"))

pdf("FIGURES/RAW/Figure4.pdf", height = 9, width = 9)
p
dev.off()

### Figure S3: LWR vs % of ASVs ----

# read the LWR of each ASVs obtained by GAPPA 
data <- read_csv("data/PERKINSEA_TABLE/lwr_histogram.csv") #
data$Start <- format(data$Start, digits=2, nsmall=2)
data$Start <- paste0("> ", data$Start)
data_long <- gather(data, Group, LWR, `Percentage 1`:`Percentage 3`, factor_key=TRUE)
data_long$Group <- gsub('Percentage', 'LWR', data_long$Group)
data_long$LWR <- data_long$LWR*100
# Bar Plot
p <- ggplot(data_long, aes(x=Start, y=LWR, fill=Group)) +
  geom_bar(alpha=0.8, stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("#f60000", "#f6db00", "#4080fb")) +
  xlab("LWR") +
  ylab("% ASVs") +
  theme_bw() + 
  tm +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust = 1.2))

pdf("FIGURES/RAW/FigureS3.pdf", height = 5, width = 7)
p
dev.off()

### Figure S4: Principal Coordinate Analysis (PCoA) of weighted UniFrac distance between sample ----
weighted <- readRDS("data/PERKINSEA_TABLE/PICANTE/weighted.rds")
aux_w <- broom::tidy(weighted)

aux_wf <- filter(aux_w, item1 %in% sampF$sample & item2 %in% sampF$sample)
nams <- with(aux_wf, unique(c(as.character(item1), as.character(item2))))
dis_w <- with(aux_wf, structure(distance,
                                Size = length(nams),
                                Labels = nams,
                                Diag = FALSE,
                                Upper = FALSE,
                                method = "user",
                                class = "dist"))
pcres <- pcoa(dis_w)
pctvar <- round(100*pcres$values$Relative_eig, 1)


plotdat = data.frame(Axis1 = pcres$vectors[,1],
                     Axis2 = pcres$vectors[,2],
                     sample = names(pcres$vectors[,1]))

plotdat <- plotdat %>% left_join(select(samples, sample, environment))

p <- ggplot(data = plotdat, aes(x = Axis1, y = Axis2, color = environment)) +
  geom_point() +
  xlab(paste("Axis 1 (", pctvar[1], "%)", sep = "")) +
  ylab(paste("Axis 2 (", pctvar[2], "%)", sep = "")) +
  scale_color_manual(values = env_colours) +
  theme_bw()+
  tm

pdf("FIGURES/RAW/FigureS4.pdf", height = 5, width = 6)
p
dev.off()

### Figure S5: Diversity analysis of the Perkinsea samples ----

# PD and NRI
phy <- read.tree("data/PERKINSEA_TABLE//PICANTE/RaxML.raxml.bestTree") # Phylogeny V4
table <- read_tsv("data/PERKINSEA_TABLES/Perkinsea_filtered_rarrefy_10k.table") # Rarefied table
metadata <- read_tsv("data/PERKINSEA_TABLES/Perkinsea_metadata.table") # metadata

table$amplicon_n <- NA
for (i in 1:length(table$amplicon)) {
  table$amplicon_n[i] <- taxonomy$amplicon[grep(T, str_detect(taxonomy$amplicon, table$amplicon[i]))]
}
table$amplicon_n <- str_replace(table$amplicon_n, ";", "_")
table$amplicon_n <- str_replace(table$amplicon_n, "=", "_")
table$amplicon <- table$amplicon_n

colnames(table)[(ncol(table)-1):ncol(table)]

tablef <- table[,-c((ncol(table)-1):ncol(table))]

comm <- as.matrix(t(tablef))
rownames(comm) <- colnames(table[,1:(ncol(table)-2)])
colnames(comm) <- table$amplicon

pruneedphy <- prune.sample(comm, phy)
pd.result <- pd(comm, phy, include.root = TRUE)
pd.result$sample <- rownames(pd.result)
pd.result %>%
  left_join(select(samples, sample, environment)) %>%
  group_by(environment) %>%
  summarise(max = max(PD), mean = mean(PD))


# ses mpd = -NRI: mpd.obs.z Standardized effect size of mpd vs. null communities (= (mpd.obs - mpd.rand.mean) / mpd.rand.sd, equivalent to -NRI)
aux_cop <- cophenetic(phy)
sesmntd <- ses.mpd(comm, aux_cop, null.model = "taxa.labels")
sesmntd$sample <- rownames(sesmntd)
sesmntd$NRI <- -1*sesmntd$mpd.obs.z
sesmntdF <- select(sesmntd,sample, NRI) %>%
  left_join(select(samples, sample, environment, subcategory)) %>%
  filter(environment == "Marine")

aux_plot <- div %>%
  left_join(select(pd.result, sample, PD)) %>%
  left_join(select(samples, sample, subcategory, environment), 
            by = c("sample" = "sample")) %>%
  left_join(select(sesmntd, sample, mpd.obs.z ))

aux_plot$NRI <- -1 * aux_plot$mpd.obs.z

summary_div <- aux_plot %>%
  group_by(environment) %>%
  summarise(mean =  mean(shannon), sd = sd(shannon), 
            mean_rich = mean(richness), max_rich = max(richness))

aux_plot %>% group_by(subcategory) %>% summarise(mena = mean(PD))


aux <- select(aux_plot, shannon, richness, PD, NRI, subcategory, environment) %>%
  pivot_longer(cols = 1:4, names_to = "index", values_to = "value")

aux$order <- factor(aux$index, levels = c("shannon","richness", "PD", "NRI"))
labs.index = aux$order
p <- aux %>%
  ggplot(aes(x = factor(environment, levels = c("Marine", "Land water", "Soil")), 
             y = value)) +
  geom_boxplot(fill = "gray90")+
  geom_jitter(aes(fill = subcategory), colour = "black", size = 2,
              pch = 21, alpha = 0.7, width = 0.3) +
  scale_fill_manual(values = sub_colours) +
  labs(fill = "Subenvironment category")+
  xlab(NULL) + ylab("Index") +
  facet_wrap(. ~ order, scales = "free_y", labeller = labeller(index = labs.index)) +
  theme_bw()+
  tm

pdf("FIGURES/RAW/FigureS5.pdf", height = 5, width = 6)
p
dev.off()


### Figure S6: Rarefaction curves: Sampling effort ----
aux <- tibble(Marine = rowSums(abundance[, colnames(abundance) %in% filter(samples, environment == "Marine")$sample]),
              `Land water` = rowSums(abundance[, colnames(abundance) %in% filter(samples, environment == "Land water")$sample]),
              Soil = rowSums(abundance[, colnames(abundance) %in% filter(samples, environment == "Soil")$sample]))

rarecurve(t(aux))
dev.copy(pdf,'FIGURES/RAW/FigureS6.pdf', width=6, height=3.5)
dev.off()

### Figure S7: Phylogenetic placement of potential novel ASVs and environmental distribution ----
# NOTE: This is a composition of a heat map and a tree with the placement collapsed in groups. Here I only show the code for the heat map.

novel <- read_tsv("data/PERKINSEA_TABLE/NOVELTY/taxonomy_per_query.tsv")

ncol <- max(str_count(novel$taxopath, ";")) + 1
colnm <- str_c("TAX", 1:ncol)

novel <- novel %>%
  separate(col = taxopath, sep = ";", into = colnm, remove = F) %>%
  mutate(tax = ifelse(is.na(TAX5) | TAX5 == "Perkinsea_X", "unclassified Perkinsea", TAX5))


list <- read_tsv("data/PERKINSEA_TABLE/NOVELTY/clusters.txt")
list <- list[!is.na(str_match(list$amplicon,";size")),]
list$amplicon <- str_extract(list$amplicon, "\\w+")

aux_nov <- list %>%
  left_join(abundance, by = c("amplicon" = "amplicon")) %>%
  select(amplicon, cluster, colnames(abundance)[2:ncol(abundance)]) %>%
  pivot_longer(cols = colnames(.)[3:ncol(.)], names_to = "sample", values_to = "reads") %>%
  filter(reads > 0) %>% 
  left_join(select(samples, sample, environment, subcategory)) %>%
  group_by(amplicon, cluster, environment, subcategory) %>%
  summarise(value = n()) %>% 
  group_by(cluster, environment, subcategory) %>%
  summarise(value = n()) %>%
  left_join(group_by(list, cluster) %>% summarise(max = n())) %>%
  mutate(relative = value/max*100)

aux_nov2 <- aux_nov %>%
  group_by(cluster, subcategory) %>% 
  summarise(relativeF = sum(relative))  %>%
  pivot_wider(names_from = subcategory, values_from = relativeF)

x <- as.matrix(aux_nov2[,c(8,2,9,7,6,3,4,15,5,14,10,12,13,11)]) #c(12,5,7,14,15,13,6,9,3,4,5,2,8,10,11)
rownames(x) <- str_c("Group ", aux_nov2$cluster)
x <- x[c(30,28,27,26,25,24,23,22,21,20,19,18,17,16,14,13,12,11,10,9,8,7,6,5,4,3,36,35,34,33,32,31,29,15,2,1),]

pdf("FIGURES/RAW/FigureS7.pdf", width = 5, height = 12)
  pheatmap(x, cluster_cols = F, cluster_rows = F)
dev.off()

### FigureS8: Xcellidae rRNA/rDNA ratio in Malaspina dataset ----
#Read tables
metadata <- read_tsv("data/MALASPINA_TABLES/MPN_sf_vp_metadata.txt")
ASVs <- read_tsv("data/MALASPINA_TABLES/taxonomy_per_query.tsv")
count <- read_tsv("data/MALASPINA_TABLES/MPN_sf_vp_picoASV_counts_v3.txt")

# number of reads in metadata
metadata <- metadata %>%
  left_join(tibble(Sample = colnames(count[,3:ncol(count)]), 
                   reads = colSums(count[,3:ncol(count)])))

# Filter Xcellidae ASVs
aux <- as_tibble(t(count[,3:ncol(count)]))
colnames(aux) <- count$ASVid
auxF <- add_column(tibble(Sample = colnames(count[,3:ncol(count)])), aux)
auxF <- auxF[,c(1, grep(T, colnames(auxF) %in% ASVs$name))]

# Table with all the data, Metadata + ASVs count
table <- metadata %>%
  left_join(auxF) %>%
  filter(!is.na(reads))

# Contribution of the ASVs per sample
for (i in 1:nrow(table)) {
   for (j in 28:ncol(table)) {
    table[i,j] <- table[i,j] / table[i,27]  * 100
  }
}

# group Samples from same station and depth and calculate ASV RNA/DNA ratio
aux <- table %>%
  pivot_longer(28:ncol(table), names_to = "ASV", values_to = "Contribution") %>%
  filter(Dataset != "Surface") %>%
  select(ASV, Extract, Station, Depth, Ocean, lat, long, Contribution) %>%
  pivot_wider(names_from = Extract, values_from = Contribution) %>%
  filter(DNA != 0 | RNA != 0) %>%
  filter(!is.na(DNA)) %>%
  mutate(rna_dna = ifelse(DNA == 0, "RNA", ifelse(RNA == 0, "DNA", ifelse(RNA>DNA, "RNA>DNA", "DNA>RNA"))))



# Plot
aux$rna_dna <- as.factor(aux$rna_dna)


aux$DepthF <- as.factor(aux$Depth)

p <-ggplot(filter(aux, between(Depth, 200, 1400))) +
  geom_point(aes(x = ASV, y = DepthF, pch = rna_dna), 
             color = "black")+
  scale_shape_manual(name = "", labels = c("DNA soley",  
                                "Hypoactive (RNA<DNA)", 
                                "RNA solely", 
                                "Hyperactive (RNA>DNA)"), 
                     values = c(23, 24, 4, 21))+
  scale_fill_manual(name = "", labels = c("", ""), values = c("#0077b6"))+
  xlab("Station")+
  scale_y_discrete(name = "Depth [m]", limits = rev)+
  scale_size_continuous(name = "RNA/DNA ratio", breaks = seq(from = 0, to = 10, by = 2.5))+
  facet_wrap(~ Ocean, strip.position = "top", scales = "free_x")+
  theme_bw()+
  tm

pdf("FIGURES/RAW/FigureS8.pdf", width = 10, height = 6)
p
dev.off()
