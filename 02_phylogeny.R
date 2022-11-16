##### Phylogeny #####

# source header
source("header.R")

library(devtools)
library(V.PhyloMaker)
library(ape)

### load in common garden data
traits_bysp = read.csv(file = "objects/traits_bysp.csv") %>%
  na.omit()

# make a new dataframe to input into phylo.maker
phylo = data.frame(species = traits_bysp$species, 
                   genus = traits_bysp$genus, 
                   family =traits_bysp$family)
phylo_out <- phylo.maker(phylo, scenarios=c("S1","S2","S3"))
par(mfrow = c(1, 1))
plot.phylo(phylo_out$scenario.1, cex = 1.5, main = "scenario.1")
nodelabels(round(branching.times(phylo_out$scenario.1), 1), cex = 1)
scen1 = phylo_out[[1]]
plotTree(scen1,ftype="i",fsize=0.6,lwd=1)

write_rds(scen1, "objects/comm_gard_phylo.rds")


##### Phylogeny of GrooT Species + Common Garden #####

all_complete_traits = readr::read_csv("objects/all_complete_traits.csv")

all_complete_traits$genus = stringr::word(all_complete_traits$species)

phylo_input = all_complete_traits[,c(1,14,2)]
phylo_input$family = ifelse(phylo_input$species == "Nyssa sylvatica", "Nyssaceae", phylo_input$family)

phylo_out = phylo.maker(phylo_input, scenarios=c("S1","S2","S3"))
plot.phylo(phylo_out$scenario.1, cex = .75, main = "scenario.1", type = "fan")

# save phylogeny as an object
RES_phylo = phylo_out$scenario.1

readr::write_rds(RES_phylo, "objects/RES_phylo.rds")