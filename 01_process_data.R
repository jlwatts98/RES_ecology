##### Process Raw Data #####
# source header 
source("header.R")

# load in data

scans = readr::read_csv("raw-data/scanR.csv")
scans$Code = str_sub(scans$file_name, end=-5)

weights_PLMC = readr::read_csv("raw-data/PLMC.csv")
percent_N_C = readr::read_csv("raw-data/percent_N_C.csv")

# merge dataframes by species code
data = merge(scans, weights_PLMC, by = "Code")
data = merge(data, percent_N_C, by = "Code")

setdiff(data$Code, percent_N_C$Code)
intersect(data$Code, percent_N_C$Code)
# get unique species
species_list = unique(data$Species_gbif) |>
  as.data.frame()
names(species_list) = "Species"
species_list$Genus = stringr::word(species_list$Species, 1)

# add SRL (m/g)
data$length_m = data$length_mm / 1000
data$SRL = data$length_m / data$Weight_g

# add RTD (g/cm^3)
data$volume_cm3 = data$volume_mm3 / 1000
data$RTD = data$Weight_g / data$volume_cm3


# clean unnecessary columns
cleaned = data[, c("Code", "Species", "Species_gbif",
                  "branching_freq_mm-1","mean_diameter_mm", 
                  "SRL", "RTD", "mycorrhizal_colonization", "coll_location",
                  "hair_length_raw", "hair_density_raw", "hair_factor",
                  "percent_N", "d_15N", "percent_C", "d_13C")] |>
  rename(
    species = Species_gbif,
    D = mean_diameter_mm,
    BF = `branching_freq_mm-1`,
    `%M` = mycorrhizal_colonization,
    HF = hair_factor,
    `%N` = percent_N,
    `%C` = percent_C
  )

# write a csv
readr::write_csv(cleaned, "objects/cleaned.csv")

# summarise traits by species
traits_bysp = cleaned |>
  dplyr::group_by(species) |>
  dplyr::summarise(SRL = mean(SRL),
                   D = mean(mean_diameter_mm),
                   RTD = mean(RTD),
                   BF = mean(`branching_freq_mm-1`),
                   `%M` = mean(mycorrhizal_colonization),
                   HF = mean(hair_factor),
                   `%N` = mean(percent_N),
                   percent_C = mean(percent_C),
                   d_15N = mean(d_15N),
                   d_13C = mean(d_13C))

# load in species classifications
root_spec = readr::read_csv("raw-data/root_spec.csv")

traits_bysp = merge(root_spec, traits_bysp, by = "species")
traits_bysp$group = traits_bysp$class
# save
readr::write_csv(traits_bysp, "objects/traits_bysp.csv")

##### GROOT DATA #####
# read in data
groot = read.delim("raw-data/GRooTFullVersion.csv", sep = ",") |>
  tidyr::unite("genus_species", 10:11, remove = F)
colnames(groot)[1] = "grootID"

# select only traits that I have data for; calculate mean and standard devation per trait
groot_reduced_traits = groot[groot$traitName %in% c("Root_C_concentration", "Root_C_N_ratio", "Root_N_concentration",
                                                    "Root_mycorrhizal colonization", "Mean_Root_diameter",
                                                    "Root_tissue_density", "Specific_root_length",
                                                    "Root_branching_density"), ] |>
  dplyr::group_by(genus_species, traitName, familyTNRS, mycorrhizalAssociationTypeFungalRoot, group, nitrogenFixationNodDB) |>
  dplyr::summarise(mean_trait_value = mean(traitValue))

# transform long to wide
groot_pivot = tidyr::pivot_wider(data = groot_reduced_traits, names_from = traitName, values_from = mean_trait_value)

groot_complete = groot_pivot |>
  ungroup() |>
  dplyr::select(genus_species,
                familyTNRS, 
                Mean_Root_diameter, 
                Root_N_concentration, 
                Root_tissue_density, 
                Specific_root_length, 
                `Root_mycorrhizal colonization`,
                mycorrhizalAssociationTypeFungalRoot,
                group,
                nitrogenFixationNodDB)

# rename columns to match other columns
names(groot_complete)[names(groot_complete) == "Root_mycorrhizal colonization"] = "%M"
names(groot_complete)[names(groot_complete) == "Mean_Root_diameter"] = "D"
names(groot_complete)[names(groot_complete) == "LRR"] = "mean_LRR"
names(groot_complete)[names(groot_complete) == "Root_N_concentration"] = "%N"
names(groot_complete)[names(groot_complete) == "Root_tissue_density"] = "RTD"
names(groot_complete)[names(groot_complete) == "Specific_root_length"] = "SRL"
names(groot_complete)[names(groot_complete) == "familyTNRS"] = "family"
names(groot_complete)[names(groot_complete) == "genus_species"] = "species"
names(groot_complete)[names(groot_complete) == "nitrogenFixationNodDB"] = "nitrogen_fix"

# fix nitrogen variable
groot_complete$`%N` = groot_complete$`%N` / 10


# write csv
readr::write_csv(groot_complete, "objects/groot_indiv_traits.csv")

groot_complete = groot_complete |>
  na.omit()

groot_complete$myc_type = ifelse(groot_complete$mycorrhizalAssociationTypeFungalRoot == "", "AM",
                                 groot_complete$mycorrhizalAssociationTypeFungalRoot)
groot_complete = groot_complete |>
  dplyr::filter(myc_type == "AM" | myc_type == "NM-AM")


# remove outliers in root tissue density
groot_complete = groot_complete[!groot_complete$RTD > 1, ]

# remove SRL outlier
max(groot_complete$SRL)
groot_complete = groot_complete[-which.max(groot_complete$SRL),]
max(groot_complete$SRL)

# remove AM columns
groot_complete = groot_complete |>
  dplyr::select(species, family, D, `%M`, SRL, RTD, `%N`, group)

# write csv
readr::write_csv(groot_complete, "objects/groot_complete.csv")

##### Combine GROOT and Common Garden Data #####
all_complete_traits = rbind(groot_complete, traits_bysp[, c(1, 6, 8:10, 12, 14, 18)])
all_complete_traits$source = c(rep("groot", nrow(groot_complete)),
                               rep("common_garden", nrow(traits_bysp)))

# change species names to spaces
all_complete_traits$species = gsub("_", " ", all_complete_traits$species)

all_complete_traits$group1 = ifelse(all_complete_traits$group == "Magnoliopsida", "Eudicot",
                                    ifelse(all_complete_traits$group == "Polypodiopsida", "Pteridophytes",
                                           ifelse(all_complete_traits$family == "Poaceae", "Graminoid",
                                                  ifelse(all_complete_traits$family == "Cyperaceae", "Graminoid",
                                                         ifelse(all_complete_traits$family == "Juncaceae", "Graminoid",
                                                                ifelse(all_complete_traits$group == "Liliopsida", "Monocot",
                                                                       ifelse(all_complete_traits$group == "Angiosperm_Monocotyl", "Monocot",
                                                                              ifelse(all_complete_traits$family == "Magnoliaceae", "Magnoliid",
                                                                                     ifelse(all_complete_traits$family == "Lauraceae", "Magnoliid",
                                                                                            ifelse(all_complete_traits$group == "Angiosperm_Magnoliid", "Magnoliid",
                                                                                                   ifelse(all_complete_traits$group == "Angiosperm_Eudicotyl", "Eudicot",
                                                                                                          all_complete_traits$group)))))))))))

all_complete_traits$group1 = ifelse(is.na(all_complete_traits$group1), "Eudicot",
                                    ifelse(all_complete_traits$family == "Magnoliaceae", "Magnoliid",
                                           ifelse(all_complete_traits$family == "Lauraceae", "Magnoliid",
                                                  ifelse(all_complete_traits$group == "", "Eudicot",
                                                  all_complete_traits$group1))))

all_complete_traits = all_complete_traits |>
  dplyr::group_by(species, family, group1) |>
  dplyr::summarise(D = mean(D),
                   `%N` = mean(`%N`),
                   RTD = mean(RTD),
                   SRL = mean(SRL),
                   `%M` = mean(`%M`)) |>
  dplyr::ungroup() |>
  na.omit() |>
  as.data.frame()


##### PCA on complete traits #####
RES_pca = prcomp(all_complete_traits[, 4:8], center = T, scale = T)
PCs = as.data.frame(RES_pca[["x"]])
all_complete_traits = cbind(all_complete_traits, PCs)
readr::write_csv(all_complete_traits, "objects/all_complete_traits.csv")

##### CSR Data #####
# load in data
CSR_raw = read.csv("raw-data/CSR_paper_supplementary_information.csv")
names(CSR_raw)[1] = "species_author"
names(CSR_raw)[9] = "C"
names(CSR_raw)[10] = "S"
names(CSR_raw)[11] = "R"
# create species column
CSR_raw$species = stringr::word(CSR_raw$species_author, 1, 2, sep = " ")

# define which CSR category is maximum
CSR_raw$max = apply(X = CSR_raw[, c(9:11)], MARGIN=1, FUN=which.max)
CSR_raw$max = ifelse(CSR_raw$CSR == "CSR", "CSR",
                     ifelse(CSR_raw$max == 1, "C",
                            ifelse(CSR_raw$max == 2, "S", "R")))
CSR = CSR_raw |>
  dplyr::select(species, max, C, S, R)

readr::write_csv(CSR, "objects/CSR.csv")
  
##### Merge CSR with RES #####
CSR_RES = merge(CSR, all_complete_traits, by = "species")

readr::write_csv(CSR_RES, "objects/CSR_RES.csv")

