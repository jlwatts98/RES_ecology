##### Data Analysis #####


### This script is the analysis associated with the INSERT JOURNAL TITLE
### DOI: INSERT DOI.

# source header to load packages and source functions
source("header.R")

# load in data produced by 01_process_data.R script
all_complete_traits = readr::read_csv("objects/all_complete_traits.csv")

comm_gard = readr::read_csv("objects/cleaned.csv")

traits_bysp = readr::read_csv("objects/traits_bysp.csv") %>%
  na.omit()

##### Figure 1: Common Garden RES Data #####

### figure 1A: PCA of common garden data
# omit missing values
comm_gard = comm_gard %>%
  na.omit()
comm_gard_pca = prcomp(comm_gard[, c(5:8, 13)], center = T, scale = T)

# Extract PC axes
PCAvalues <- data.frame(comm_gard_pca$x)

# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(comm_gard_pca$rotation), comm_gard_pca$rotation)

# Calculate the angles and the label offset
PCAloadings$Angle = ((180/pi) * atan(PCAloadings$PC2/PCAloadings$PC1))
PCAloadings$Offset <- c(2.1, -.5, 1.4, 1.75, -.75)

# Plot
fig1a = ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, color = "pink", alpha = 1) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (1.5 * PC1), yend = (1.5 * PC2)),
               arrow = arrow(length = unit(1/2, "picas")), color = "black", size = 1) +
  geom_text(data = PCAloadings, aes(label = Variables, x = (PC1), y = (PC2)), 
            color = "black", size = 5, angle = PCAloadings$Angle, hjust = 
              PCAloadings$Offset) +
  theme_bw() + coord_equal() +
  xlab("PC 1 (39.5% expl.)") +
  ylab("PC 2 (25.0%)") +
  theme(axis.title = element_text(size = 15))
fig1a


## figure 1B: D and %M
fig1b = ggplot(comm_gard, mapping = aes(x = D, y = `%M`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "D (mm)", y = "%M (%)") +
  theme(axis.title = element_text(size = 15))
fig1b

## figure 1C: RTD and %N
fig1c = ggplot(comm_gard, mapping = aes(x = RTD, y = `%N`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = expression(paste("RTD (g/", cm^3, ")")), y = "%N (%)") +
  theme(axis.title = element_text(size = 15))
fig1c


# linear models
d_m_lm = lm(D ~ `%M`, data = comm_gard)
summary(d_m_lm)
rtd_n_lm = lm(RTD ~ `%N`, data = comm_gard)
summary(rtd_n_lm)

fig1 = ggarrange(fig1a, fig1b, fig1c,
                 nrow = 1,
                 labels = c("A", "B", "C"))
fig1

# save the plot
ggsave(filename = "objects/fig1.jpg",
       plot = fig1, 
       width = 7, 
       height = 2.5, 
       units = "in", 
       dpi = 600)

##### Supplementary Figure 1: d_15N compared to PCs #####

# bind pc values to dataframe
comm_gard_PC = cbind(comm_gard, PCAvalues)

# delta 15N vs PC1
comm_gard_d_15N_pc1 = ggplot(comm_gard_PC, mapping = aes(x = PC1, y = d_15N)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "PC 1 (Collaboration Gradient)", y = expression(paste(delta^15, "N"))) +
  theme(axis.title = element_text(size = 15))
comm_gard_d_15N_pc1

# delta 15N vs. PC2
comm_gard_d_15N_pc2 = ggplot(comm_gard_PC, mapping = aes(x = PC2, y = d_15N)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "PC 2 (Conservation Gradient)", y = "") +
  theme(axis.title = element_text(size = 15))
comm_gard_d_15N_pc2

supp_fig1 = ggarrange(comm_gard_d_15N_pc1, comm_gard_d_15N_pc2, labels = c("A", "B"))
supp_fig1

# linear models
d_15n_pc1_lm = lm(PC1 ~ d_15N, data = comm_gard_PC)
summary(d_15n_pc1_lm)

d_15n_pc2_lm = lm(PC2 ~ d_15N, data = comm_gard_PC)
summary(d_15n_pc2_lm)

# save the plot
ggsave(filename = "objects/supp_fig1.jpg",
       plot = supp_fig1, 
       width = 7, 
       height = 4, 
       units = "in", 
       dpi = 600)

##### Supplementary Table 1: Common Garden Counts #####

# count number of samples per species of the common gardn
supp_table1 = comm_gard |>
  na.omit() |>
  group_by(species) |>
  tally()

# write to a csv file
write.csv(supp_table1, "objects/supp_table1.csv")


##### Figure 2: GrooT + Common Garden Root Economic Spectrum #####

### panel figure: Root Economic Spectrum

## figure 2A: PCA
all_complete_pca = prcomp(all_complete_traits[, c(4:8)], center = T, scale = T)

# Extract PC axes
PCAvalues <- data.frame(all_complete_pca$x)

# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(all_complete_pca$rotation), all_complete_pca$rotation)

# Calculate the angles and the label offset
PCAloadings$Angle = ((180/pi) * atan(PCAloadings$PC2/PCAloadings$PC1))
PCAloadings$Offset <- c(2.4, 1.5, 1.4, -.5, 1.5)

# plot
fig2 = ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, color = "pink", alpha = 1) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (1.5 * PC1), yend = (1.5 * PC2)),
               arrow = arrow(length = unit(1/2, "picas")), color = "black", size = 1) +
  geom_text(data = PCAloadings, aes(label = Variables, x = (PC1), y = (PC2)), 
            color = "black", size = 5, angle = PCAloadings$Angle, hjust = 
              PCAloadings$Offset) +
  theme_bw() + coord_equal() +
  xlab("PC 1 (41.5% expl.)") +
  ylab("PC 2 (22.6% expl.)") +
  theme(axis.title = element_text(size = 15))
fig2

## figure 2B: Density Plot

#fig2b = ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
  #stat_density_2d(aes(fill = -after_stat(density)), geom = "raster", 
  #               contour = FALSE,
  #               h = c(2.5, 2.5)) +
#  stat_density2d(h = c(2.5, 2.5), bins = 4, color = "black") +
#  scale_fill_gradientn(colours = c("darkred", "orange", 
#                                   "yellow", "white")) +
#  theme_pubr() + coord_equal() +
#  xlab("Collaboration Gradient (PC 1)") +
#  ylab("Conservation Gradient (PC 2)") +
#  theme(legend.position = "none",
 #       axis.ticks.length = unit(0, "cm"),
 #       axis.text = element_text(size = 0),
 #       axis.title = element_text(size = 12))
#fig2b

#fig2 = ggarrange(fig2a, fig2b, labels = c("A", "B"))

#fig2

# save the plot
ggsave(filename = "objects/fig2.jpg",
       plot = fig2, 
       width = 4, 
       height = 4, 
       units = "in", 
       dpi = 600)

##### Table 1: PCA loadings #####

# round the loading values to 3 digits
table1 = as.data.frame(all_complete_pca$rotation) |>
  mutate_if(is.numeric, round, digits=3) |>
  as.matrix()


# Write this table to a comma separated .txt file:    
write.table(table1, 
            file = "objects/table1.txt", 
            sep = ",", quote = FALSE, row.names = T)


##### Supplementary Figure 2 Correlation Matrix #####

# libraries
library(corrplot)
library(psych)

# correlation matrix
jpeg("objects/supp_fig2.jpeg", quality = 100, width = 7.5,
     height = 5, units = "in", res = 600)
pairs.panels(all_complete_traits[,c(8, 4, 7, 6, 5)], 
             method = "pearson", # correlation method
             hist.col = "blue",
             density = TRUE,  # show density plots
             ellipses = F, # show correlation ellipses
             scale = F,
             cex = 3.5,
             cex.cor=1,
             stars = T,
             gap = 0,
             smoother = T,
             cex.axis = 2
)

dev.off()

##### Supplementary Figure 3: Whittaker Biome System Classification #####

#library(devtools)
#library(plotbiomes)

#species_list = all_complete_traits$species

# fill in your gbif.org credentials 
#user <- "jwatts" # your gbif.org username 
#pwd <- "Jw017138@" # your gbif.org password
#email <- "jwatts@colgate.edu" # your email 

#library(dplyr)
#library(purrr)
#library(readr)  
#library(magrittr) # for %T>% pipe
#library(rgbif) # for occ_download
#library(taxize) # for get_gbifid_

#file_url <- "https://data-blog.gbif.org/post/2019-07-11-downloading-long-species-lists-on-gbif_files/global_tree_search_trees_1_3.csv"

# match the names 
#gbif_taxon_keys <- 
#  all_complete_traits %>% 
#  pull("species") %>% # use fewer names if you want to just test 
#  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
#  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
#  bind_rows() %T>% # combine all data.frames into one
#  readr::write_tsv(file = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
#  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
#  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
#  pull(usagekey) # get the gbif taxonkeys


# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!

#occ_download(
#  pred_in("taxonKey", gbif_taxon_keys),
#  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN')),
#  pred("hasCoordinate", TRUE),
#  pred("hasGeospatialIssue", FALSE),
#  pred_gte("year", 1900),
#  format = "SIMPLE_CSV",
#  user=user,pwd=pwd,email=email
#)

# check status of download
#occ_download_wait('0120153-220831081235567')
# retrieve download
#whittaker_data <- occ_download_get('0120153-220831081235567') %>%
#  occ_download_import()

# citation
# GBIF Occurrence Download https://doi.org/10.15468/dl.zvbvh3 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2022-10-25


# clean occurrence records
#library(CoordinateCleaner)
#library(rnaturalearthdata)

#flag problems
#dat <- data.frame(whittaker_data)
#flags <- clean_coordinates(x = dat, 
#                           lon = "decimalLongitude", 
#                           lat = "decimalLatitude",
#                           #countries = "countryCode",
#                           species = "species",
#                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
#                                     "zeros")) # most test are on by default
#Exclude problematic records
#whittaker_clean <- dat[flags$.summary,] |>
#  dplyr::select(species, decimalLatitude, decimalLongitude) |>
#  dplyr::rename(lat = decimalLatitude,
#                lon = decimalLongitude)

# remove na from lat and lon, make numeric
#whittaker_clean = whittaker_clean[!is.na(whittaker_clean$lat),]
#whittaker_clean = whittaker_clean[!is.na(whittaker_clean$lon),]
#whittaker_clean$lat = as.numeric(whittaker_clean$lat)
#whittaker_clean$lon = as.numeric(whittaker_clean$lon)

#whittaker_points = whittaker_clean |>
#  dplyr::select(lon, lat)
# get mean annual temperature and total precipitation rasters

# bioclim data
#bioclim = raster::getData("worldclim", var="bio", res = 2.5)

#temp_precip = raster::stack(bioclim[[1]], bioclim[[12]])

# extract value to point
#tp_vals = raster::extract(temp_precip, whittaker_points)

# combine the values and occurrences
#whittaker_tp = cbind(whittaker_clean, tp_vals)


# summarise by species
#whittaker_tp = 
#  whittaker_tp |>
#  dplyr::group_by(species) |>
#  dplyr::summarise(bio1 = mean(bio1, na.rm = T),
#                   bio12 = mean(bio12, na.rm = T))

# make whittaker plot

#whittaker_plot = whittaker_base_plot() +
#  geom_point(data = whittaker_tp, aes(x = bio1/10, y = bio12/10)) +
#  theme_bw() +
#  labs(x = "Mean Annual Temperature (\u00B0C)", y = "Total Annual Precipitation (cm)")

# call plot
#whittaker_plot

# save
#ggsave(filename = "objects/whittaker_biomes.jpg", 
#       plot = whittaker_plot, 
#       height = 5, width = 7, 
#       units = "in", dpi = 600)

##### Figure 3: Phylogenetic Signal across groups #####

# load in phylogeny
comm_gard_phylo = readr::read_rds("objects/comm_gard_phylo.rds")
RES_phylo = readr::read_rds("objects/RES_phylo.rds")

all_complete_traits$all = rep("all", nrow(all_complete_traits))
all_complete_traits$genus = stringr::word(all_complete_traits$species, 1)

all_complete_traits = all_complete_traits[, c(3,2, 15,4:14, 1)]
str(all_complete_traits)
all_complete_traits$genus = as.factor(all_complete_traits$genus)
all_complete_traits$family = as.factor(all_complete_traits$family)
all_complete_traits$group1 = as.factor(all_complete_traits$group1)
all_complete_traits = as.data.frame(all_complete_traits)
str(all_complete_traits)

# run moran's I correlogram across genus, family, and group
group_correlogram = correlogram.formula(PC1 + PC2 ~ group1/family/genus,
                                        data = all_complete_traits)

# assemble data frame for plotting
corr = rbind(group_correlogram$PC1, 
             group_correlogram$PC2)

corr$trait = c(rep("PC1", 3),
               rep("PC2", 3))
corr$labels = ifelse(corr$labels == "group1", "group", corr$labels)
corr$labels1 = c(1,2,3,1,2,3)

# plot
fig3 = corr |> 
  ggplot(aes(x = labels1, y = obs, group = trait)) +
  geom_point(stat = "identity", size = 4) + 
  geom_line(aes(linetype = trait), size = 1) +
  theme_bw() +
  scale_x_continuous(breaks=c(1,2,3),
                     labels=c("Genus", "Family", "Group")) +
  scale_linetype_discrete(name = "Gradient",
                          labels = c("Collaboration",
                                     "Conservation")) +
  theme(axis.title=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 15),
        legend.position= "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x=element_text(angle=0, vjust=0, hjust=0.5, size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        strip.text.x = element_text(size = 12)) +
  labs(
    x = "Taxonomic or Functional Group",
    y = "Moran's I Statistic"
  )
fig3

# save
ggsave("objects/fig3.jpg", plot = fig3,
       dpi = 600, width = 5, height = 3, units = "in")


##### Figure 4: Phylogeny + PCA #####

### fig 4a: PCA with groups
library(ggbiplot)

# assign new scale_color_discrete function
scale_color_discrete <- function(...) {
  scale_color_manual(..., values = viridis::viridis(6))
}

# plot pca
fig4a = ggbiplot(all_complete_pca, 
                              groups = all_complete_traits$group1, 
                              ellipse = T,
                              varname.size = 5) +
  theme_bw() +
  xlab("PC 1 (41.5% expl.)") +
  ylab("PC 2 (22.6% expl.)") +
  scale_color_discrete(name = "Group",
                       labels=c("Dicot","Graminoid","Gymnosperm", "Magnoliid",
                                "Monocot", "Fern")) +
  theme(legend.position = "right",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15))
fig4a

# save
ggsave(filename = "objects/fig4a.jpg", plot = fig4a,
       width = 5, height = 4, dpi = 600,
       units = "in")



### fig 4b: PC1 + PC2 phylogeny colored by local signal

# library phylogenetic packages
library(ape)
library(phytools)
library(phylosignal)
library(phylobase)

# modify all_complete_traits dataset to match phylogeny tip names
all_complete_traits$species = gsub(" ", "_", all_complete_traits$species)
rownames(all_complete_traits) = all_complete_traits$species
all_complete_traits$species = as.factor(all_complete_traits$species)

# make a phylo4d object
p4d_pc = phylo4d(RES_phylo, all_complete_traits[, c(9:10)])

# locating signal (Moran's I test for local phylogenetic signal)
# p-values computed via iteration
set.seed(1)
PC.lipa <- lipaMoran(p4d_pc)

jpeg("objects/fig4b.jpeg", quality = 100, width = 4,
     height = 4, units = "in", res = 300)
barplot.phylo4d(p4d_pc, bar.col=(PC.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE,
                bar.lwd = 1,
                show.tip = F,
                trait.bg.col = c("grey90", "white"),
                data.xlim = c(-3, 3),
                trait.cex = 1.75,
                show.data.axis = F, edge.width = .6)

dev.off()


### fig 4c: individual traits phylogeny colored by local signal

# make sure the RES_phylo object is a binary phylogenetic tree
is.binary(RES_phylo)
RES_binary = multi2di(RES_phylo)

# make a phylo4d object
p4d = phylo4d(RES_phylo, all_complete_traits[, c(4:8)])

# locating signal
RES.lipa <- lipaMoran(p4d)

# plot
jpeg("objects/fig4c.jpg", quality = 100, width = 8.5,
     height = 5, units = "in", res = 600)
barplot.phylo4d(p4d, bar.col=(RES.lipa$p.value < 0.05) + 1, center = T , scale = T,
                bar.lwd = 1.2,
                show.tip = F,
                trait.bg.col = c("grey90", "grey90",
                                 "grey90", "white",
                                 "white"),
                trait.labels = c("%M", "D", "SRL", "RTD", "%N"),
                data.xlim = c(-2, 3),
                trait.cex = 1.75,
                show.data.axis = F,
                edge.width = 0.75)

dev.off()

##### Supplementary Table 2: Phylogenetic Signal #####

# measure phylogenetic signal of individual traits
traits_signal = phyloSignal(p4d = p4d, method = "all")
traits_signal

# measure phylogenetic signal of PCs
pc_signal = phyloSignal(p4d = p4d_pc, method = "all")
pc_signal

# make into a table
traits_stat = traits_signal$stat
pc_stat = pc_signal$stat

supp_table2 = as.data.frame(t(rbind(pc_stat, traits_stat))) |>
  mutate_if(is.numeric, round, digits=2) |>
  as.matrix()

# save
write.table(supp_table2, 
            file = "objects/supp_table2.txt", 
            sep = ",", quote = FALSE, row.names = T)

##### Supplementary Table 3: Chi Squared Test of Independence on Local Phylogenetic signal #####

# extract p-values from local signal test
p_values = as.data.frame(PC.lipa$p.value)

# change to binary based on p-value 0.05
p_values$local_signal_PC1 = ifelse(p_values$PC1 < 0.05, 
                                   "signficant",
                                   "not significant")
p_values$local_signal_PC2 = ifelse(p_values$PC2 < 0.05, 
                                   "signficant",
                                   "not significant")

# make a table
supp_table3 = table(p_values$local_signal_PC1, p_values$local_signal_PC2)

# do a chi-squared test with no correction (because four unique states)
chisq.test(p_values$local_signal_PC1, p_values$local_signal_PC2, 
           correct=F)

p_values |>
  group_by(local_signal_PC1, local_signal_PC2) |>
  tally()

# write to a .txt file
write.table(supp_table3, 
            file = "objects/supp_table3.txt", 
            sep = ",", quote = FALSE, row.names = T)
##### Figure 5: Divergence Time of Family vs. %M #####

# libraries
library(stringr)
library(V.PhyloMaker)

# node information
tips_df = as.data.frame(RES_phylo$tip.label)
names(tips_df)[1] = "species"
tips_df$species = sub("_", " ", tips_df$species)
tips_df = merge(tips_df, traits_bysp, by = "species") %>%
  select(species, genus, family)
tips_df$species = sub(" ", "_", tips_df$species)

# extract divergence time information at each family node
comm_gard_nodes = build.nodes.2(RES_phylo, tips = tips_df) %>% 
  filter(level == "F")

# merge
comm_gard_family = merge(comm_gard, traits_bysp, by = "species")

comm_gard_age = merge(comm_gard_family, comm_gard_nodes, by = "family")

# family groups
fam_group = all_complete_traits[, c(1,2)]

# summarise by species, calculate standard error
diverge_comm_gard = comm_gard_age %>%
  dplyr::group_by(species, family) %>%
  dplyr::summarise(number = n(),
                   `%M` = mean(`%M.x`),
                   sd_M = sd(`%M.x`),
                   fam_divergence = mean(rn.bl)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(se = sd_M / sqrt(number))

# merge to final dataset
diverge_comm_gard = merge(diverge_comm_gard, fam_group, by = "family")

scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = viridis::viridis(6))
}

# plot
fig5a = ggplot(data = diverge_comm_gard, mapping = aes(
  x = fam_divergence, y = `%M`
)) +
  geom_errorbar(aes(ymin=`%M`-se, ymax=`%M`+se), width=1.2, lwd = 1.2) +
  geom_point(aes(size = number, fill = group1), pch = 21, alpha = .7) +
  theme_bw() +
  scale_fill_manual(name = "Group",
                      labels=c("Dicot","Graminoid", "Magnoliid",
                               "Monocot", "Fern"),
                    values = c("#440154FF", "#414487FF", 
                               "#22A884FF", "#7AD151FF",
                               "#FDE725FF")) +
  labs(y = "Mean %M (%)", 
       x = "") +
  theme(legend.position = "none",
        axis.title = element_text(size = 15))
fig5a

### Figure 5 B: Divergence Time vs. %M for all species

# repeat for all species
# node information

tips_df = as.data.frame(RES_phylo$tip.label)
names(tips_df)[1] = "species"

tips_df = merge(tips_df, all_complete_traits, by = "species") %>%
  select(species, genus, family)

RES_nodes = build.nodes.2(RES_phylo, tips = tips_df) %>% 
  filter(level == "F")

diverge_all = merge(all_complete_traits, RES_nodes, by = "family")

fig5b = ggplot(data = diverge_all, mapping = aes(
  x = rn.bl, y = `%M`
)) +
  geom_point(aes(fill = group1), pch = 21, size = 2, alpha = .7) +
  theme_bw() +
  scale_fill_discrete(name = "Group",
                       labels=c("Dicot","Graminoid","Gymnosperm", "Magnoliid",
                                "Monocot", "Fern")) +
  labs(y = "Mean %M (%)", 
       x = "Divergence Time (MY)") +
  theme(axis.title = element_text(size = 15))
fig5b
library(patchwork)

fig5 = fig5a + fig5b + plot_layout(nrow = 2, guides = 'collect') + 
  plot_annotation(tag_levels = "A")

fig5

# save
ggsave("objects/fig5.jpg", plot = fig5,
       dpi = 600, width = 5, 
       height = 4, units = "in")
##### Figure 6: Dominant CSR #####

# libraries
library(RColorBrewer)
# load in data
CSR_RES = readr::read_csv("objects/CSR_RES.csv")

### Figure 6a: PCA by Dominant CSR Strategy

# pca
CSR_RES_pca = prcomp(CSR_RES[8:12], center = T, scale = T)

pcs = as.data.frame(CSR_RES_pca$x)

# create ellipses around each max type
fig6a = ggbiplot(CSR_RES_pca,
                      groups = CSR_RES$max,
                      ellipse = T,
                      varname.size = 5) +
  theme_bw() +
  geom_point(aes(color = CSR_RES$max), size = 1.5) +
  scale_color_manual(values = brewer.pal(n = 8, name = "Dark2"),
                     name="CSR Strategy",
                     labels=c("Competitive","CSR","Ruderal", "Stress-Tolerant")) +
  xlab("PC 1 (39.4% expl.)") +
  ylab("PC 2 (24.7% expl.)") +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.margin=margin(c(1,1,1,1)))
fig6a

### Figure 6b: Boxplots
library(forcats)
PCs = stack(CSR_RES[,c(13,14)])
PCs$CSR = rep(CSR_RES$max, 2)

# reorder factor to change way it appears on the plot
PCs$CSR = fct_relevel(PCs$CSR, "C", "S", "R", "CSR")

# perform analysis of variance and post-hoc tukey tests
aov = aov(PC1 ~ max, data = CSR_RES)
summary(aov)
tukey = TukeyHSD(aov)
tukey

aov = aov(PC2 ~ max, data = CSR_RES)
summary(aov)
tukey = TukeyHSD(aov)
tukey

# create annotation dataframe based on tukey results
x = c(.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2)
y = rep(3.5, 8)
label = c("a", "a", "b", "a", "b", "b", "ab", "ab")

# plot
fig6b = ggplot(data = PCs, aes(x = CSR, y = values,
                                fill = ind)) +
  geom_boxplot(color = "black") + theme_bw() +
  scale_fill_manual(values = c("darkgrey", "white"),
                    name="Fine Root Strategies",
                    labels=c("Collaboration Gradient","Conservation Gradient")) +
  xlab("CSR Strategy") +
  ylab("PC") +
  annotate(geom = "text", x = x, y = y,
                          label = label, size = 5) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.margin=margin(c(1,1,1,1)))
fig6b

fig6 = ggarrange(fig6a, fig6b, nrow = 2, labels = c("A", "B"))
fig6

# save
ggsave("objects/fig6.jpg", plot = fig6,
       dpi = 600, width = 7, height = 4.5, units = "in")

##### Figure 7: CSR Space vs. Collaboration #####

# libraries
library(Ternary)

# Use an inverse distance weighting to interpolate between measured points
Predict <- function (predXY) {
  Distance <- function (a, b) {
    apply(a, 2, function (pt) sqrt(colSums((pt - b) ^ 2)))
  }
  dists <- Distance(xy, predXY)
  id <- 1 / dists
  idw <- id / rowSums(id)
  
  # Return:
  colSums(response * t(idw))
}

### Figure 7a: Collaboration Gradient and CSR Space

# Configure plotting area
par(mfrow = c(1, 1), mar = rep(2, 4))

# points
C = CSR_RES[,3]
S = CSR_RES[,4]
R = CSR_RES[,5]

abc <- t(cbind(C, S, R))

# Then we'll measure the response variable at each of those points:
response <- CSR_RES$PC1

jpeg("objects/fig7a.jpeg", quality = 100, width = 4,
     height = 4, units = "in", res = 600)
# Now we must start a plot, to define the coordinate system
par(mar = rep(0.2, 4))
TernaryPlot(alab = "C \u2192", blab = "\u2192 S", clab = "R \u2190",
            lab.col = c("red", "darkgreen", "blue"),
            main = "", # Title
            point = "up", lab.cex = 2, grid.minor.lines = 0,
            grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white", 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.1)

# Convert measured points to XY
xy <- TernaryToXY(abc)

# Predict at triangle centres
tri <- TriangleCentres(resolution = 9L) 
# Adjust the resolution to suit your own dataset

# Now we interpolate between our known values to generate a colour for each
# of our tiles
predicted <- Predict(tri[1:2, ])
map <- rbind(x = tri['x', ], y = tri['y', ], z = predicted,
             down = tri['triDown', ])

# Place a semitransparent colour fill over grid lines:
ColourTernary(map)

dev.off()

### calculate moran's I #####
xy = t(xy)
xy = as.data.frame(xy)
xy$PC1 = response

dists = as.matrix(dist(cbind(xy$x, xy$y)))

dists.inv <- 1/dists
diag(dists.inv) <- 0

morans_I = Moran.I(xy$PC1, dists.inv)
morans_I

### Figure 7b: Conservation Gradient in CSR space

response <- CSR_RES$PC2

spectrum = rev(brewer.pal(9, "Oranges"))
jpeg("objects/fig7b.jpeg", quality = 100, width = 4,
     height = 4, units = "in", res = 300)
# Now we must start a plot, to define the coordinate system
par(mar = rep(0.2, 4))
TernaryPlot(alab = "C \u2192", blab = "\u2192 S", clab = "R \u2190",
            lab.col = c("red", "darkgreen", "blue"),
            main = "", # Title
            point = "up", lab.cex = 2, grid.minor.lines = 0,
            grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white", 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.1)

# Convert measured points to XY
xy <- TernaryToXY(abc)


# Predict at triangle centres
tri <- TriangleCentres(resolution = 9L) 
# Adjust the resolution to suit your own dataset

# Now we interpolate between our known values to generate a colour for each
# of our tiles
predicted <- Predict(tri[1:2, ])
map <- rbind(x = tri['x', ], y = tri['y', ], z = predicted,
             down = tri['triDown', ])

# Place a semitransparent colour fill over grid lines:
ColourTernary(map, spectrum = spectrum)

dev.off()

### calculate moran's I #####
xy = t(xy)
xy = as.data.frame(xy)
xy$PC1 = response

dists = as.matrix(dist(cbind(xy$x, xy$y)))

dists.inv <- 1/dists
diag(dists.inv) <- 0

morans_I = Moran.I(xy$PC1, dists.inv)
morans_I

##### Supplementary Figure 4: Common Garden Individual Traits and CSR Boxplots#####

# libraries
library(forcats)
library(ggplot2)
# load in CSR data
CSR = readr::read_csv("objects/CSR.csv")

# add common garden species
CSR_comm_gard = merge(traits_bysp, CSR, by = "species")

# reduce for stacking
comm_gard_CSR_reduced = CSR_comm_gard |>
  dplyr::select(`%M`, D, SRL, RTD, `%N`)

# stack
cg_CSR_stack = stack(comm_gard_CSR_reduced)

# add max column
cg_CSR_stack$max = rep(CSR_comm_gard$max, 5)

# reorder factors
cg_CSR_stack$max = fct_relevel(cg_CSR_stack$max, "C", "S", "R", "CSR")
cg_CSR_stack$ind = fct_relevel(cg_CSR_stack$ind, "%M", "D", "SRL", "RTD", "%N")

# ANOVAs and post hoc tukeys
perc_M_aov = aov(`%M` ~ max, data = CSR_comm_gard)
summary(perc_M_aov)
perc_M_tukey = TukeyHSD(perc_M_aov)
perc_M_tukey

perc_N_aov = aov(`%N` ~ max, data = CSR_comm_gard)
summary(perc_N_aov)
perc_N_tukey = TukeyHSD(perc_N_aov)
perc_N_tukey


D_aov = aov(D ~ max, data = CSR_comm_gard)
summary(D_aov)
D_tukey = TukeyHSD(D_aov)
D_tukey


SRL_aov = aov(SRL ~ max, data = CSR_comm_gard)
summary(SRL_aov)
SRL_tukey = TukeyHSD(SRL_aov)
SRL_tukey

RTD_aov = aov(RTD ~ max, data = CSR_comm_gard)
summary(RTD_aov)
RTD_tukey = TukeyHSD(RTD_aov)
RTD_tukey

# plot
supp_fig4 = ggplot(data = cg_CSR_stack, 
                   mapping = aes(x = max,
                                 y = values)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ ind, ncol = 3,
             scales = "free_y") +
  labs(x = "CSR Strategy", y = "Fine Root Trait Value") +
  theme(axis.title = element_text(size = 15))
supp_fig4

# save
ggsave(filename = "objects/supp_fig4.jpg", plot = supp_fig4,
       width = 7, height = 5, dpi = 600,
       units = "in")

##### Supplementary Figure 5: Individual Traits and CSR Boxplots #####


# add groot species
groot_indiv_traits = readr::read_csv("objects/groot_indiv_traits.csv")
groot_indiv_traits$species = gsub("_", " ", groot_indiv_traits$species)

groot_indiv_traits$myc_type = ifelse(groot_indiv_traits$mycorrhizalAssociationTypeFungalRoot == "", "AM",
                                 groot_indiv_traits$mycorrhizalAssociationTypeFungalRoot)
groot_indiv_traits = groot_indiv_traits |>
  dplyr::filter(myc_type == "AM" | myc_type == "NM-AM")

CSR_traits = merge(groot_indiv_traits, CSR, by = "species")

# rbind
CSR_traits = rbind(CSR_traits[, c(1:7, 12:15)], CSR_comm_gard[, c(1, 6, 8:10, 12, 14, 19:22)]) |>
  group_by(species, max, family) |>
  summarise_all(.funs = "mean", na.rm = T) |>
  ungroup()

CSR_traits$D[is.nan(CSR_traits$D)]<-NA
CSR_traits$`%M`[is.nan(CSR_traits$`%M`)]<-NA
CSR_traits$SRL[is.nan(CSR_traits$SRL)]<-NA
CSR_traits$RTD[is.nan(CSR_traits$RTD)]<-NA
CSR_traits$`%N`[is.nan(CSR_traits$`%N`)]<-NA

CSR_traits_reduced = CSR_traits |>
  dplyr::select(`%M`, D, SRL, RTD, `%N`)
  
# filter out extreme values for better visualization
CSR_traits_reduced$D = ifelse(CSR_traits_reduced$D > 2, NA, CSR_traits_reduced$D)
CSR_traits_reduced$SRL = ifelse(CSR_traits_reduced$SRL > 500, NA, CSR_traits_reduced$SRL)
CSR_traits_reduced$RTD =ifelse(CSR_traits_reduced$RTD > 1.5, NA, CSR_traits_reduced$RTD)
CSR_traits_reduced$`%N` = ifelse(CSR_traits_reduced$`%N` > 4, NA, CSR_traits_reduced$`%N`)

# stack
CSR_stack = stack(CSR_traits_reduced)

# add max column
CSR_stack$max = rep(CSR_traits$max, 5)

# reorder factors
CSR_stack$max = fct_relevel(CSR_stack$max, "C", "S", "R", "CSR")
CSR_stack$ind = fct_relevel(CSR_stack$ind, "%M", "D", "SRL", "RTD", "%N")
CSR_stack$ind_f = factor(CSR_stack$ind, levels=c('%M','D','SRL','RTD', '%N'))

# look at structure of dataframe
str(CSR_stack)

# ANOVAs and post hoc tukeys
perc_M_aov = aov(`%M` ~ max, data = CSR_traits)
summary(perc_M_aov)
perc_M_tukey = TukeyHSD(perc_M_aov)
perc_M_tukey

perc_N_aov = aov(`%N` ~ max, data = CSR_traits)
summary(perc_N_aov)
perc_N_tukey = TukeyHSD(perc_N_aov)
perc_N_tukey


D_aov = aov(D ~ max, data = CSR_traits)
summary(D_aov)
D_tukey = TukeyHSD(D_aov)
D_tukey


SRL_aov = aov(SRL ~ max, data = CSR_traits)
summary(SRL_aov)
SRL_tukey = TukeyHSD(SRL_aov)
SRL_tukey

RTD_aov = aov(RTD ~ max, data = CSR_traits)
summary(RTD_aov)
RTD_tukey = TukeyHSD(RTD_aov)
RTD_tukey

# plot
supp_fig5 = ggplot(data = CSR_stack, 
                   mapping = aes(x = max,
                                 y = values)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ ind_f, ncol = 3,
             scales = "free_y") +
  labs(x = "CSR Strategy", y = "Fine Root Trait Value") +
  theme(axis.title = element_text(size = 15))
supp_fig5

# save
ggsave(filename = "objects/supp_fig5.jpg", plot = supp_fig5,
       width = 7, height = 5, dpi = 600,
       units = "in")