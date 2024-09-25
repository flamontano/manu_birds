## Code to replicate results from

## Dimensions of passerine biodiversity along an elevational gradient: a nexus for historical biogeography and contemporary ecology
## Burgio et al. 2024

## Created by Flavia Montano-Centellas

rm(list=ls()) 

## required packages

library(ape)
library(picante)
library(dplyr)


#### Loading Data ### --------------------------------------------------

manu_data = read.csv("Data/Manu_Dataset.csv", header = T)

spp_manu<-unique(manu_data$spp_phylo)

### PHYLOGENETI DIVERSITY ### --------------

# downloading data from Pulido-Santa Cruz and Weir (2016)

library(rdryad)

data = dryad_download(doi = "10.5061/dryad.2v462")#provides location of files

# concensus tree from Pulido-Weir et al. XX. for plotting

Consensus_tree<- read.tree(data$`10.5061/dryad.2v462`[1])
manu_tree = ape::drop.tip(Consensus_tree,Consensus_tree$tip.label[-match(spp_manu, Consensus_tree$tip.label)])

# 1000 trees for analyses

#extract files in folder "JETZ TREES" from data$`10.5061/dryad.2v462`[2] in a known folder, in our case Data/PulidoWeir_trees

Pulidotrees_names<-list.files('Data/PulidoWeir_trees',full.names=TRUE)
PulidoTrees = lapply(Pulidotrees_names, read.tree)


#### Phylogenetic dispersion ##### --------------------------------------------------------

library(tidyverse)
#library(dplyr)
library(plyr)
#devtools::install_github("daijiang/lirrr")

## Community data for phylo
passer_comm = as.data.frame(t(manu_data[, c(10:22)])); colnames(passer_comm)<-manu_data$spp_phylo
oscines_comm = as.data.frame(t(manu_data[manu_data$suborder == "Passeri", c(10:22)]));colnames(oscines_comm)<-manu_data[manu_data$suborder == "Passeri", 43]
suboscines_comm = as.data.frame(t(manu_data[!manu_data$suborder == "Passeri", c(10:22)]));colnames(suboscines_comm)<-manu_data[!manu_data$suborder == "Passeri", 43]


# MPD #

if(!file.exists("Output/PD_null_df.rds")){
  PD_null = list()
  
  for (i in 1:length(PulidoTrees)){
    phy = PulidoTrees[[i]]
    PD_null[[i]]<-picante::ses.mpd(passer_comm, cophenetic(phy), null.model =  "independentswap")#
  }
  saveRDS(PD_null, file = "Output/PD_passer_random.rds")# results for 1000 phylogenies
  
  PD_null_df = plyr::aaply(laply(PD_null, as.matrix), c(2, 3), mean)
  
  saveRDS(PD_null_df, file = "Output/PD_null_df.rds")
}else {
  PD_null_df = readRDS("Output/PD_null_df.rds")
}


if(!file.exists("Output/PD_oscines_df.rds")){
  PD_null = list()
  
  for (i in 1:length(PulidoTrees)){
    phy = PulidoTrees[[i]]
    PD_null[[i]]<-picante::ses.mpd(oscines_comm, cophenetic(phy), null.model =  "independentswap")#
  }
  saveRDS(PD_null, file = "Output/PD_oscines_random.rds")# results for 1000 phylogenies
  PD_null_df = plyr::aaply(laply(PD_null, as.matrix), c(2, 3), mean)
  
  saveRDS(PD_null_df, file = "Output/PD_oscines_df.rds")
}else {
  PD_oscines_df = readRDS("Output/PD_oscines_df.rds")
}


if(!file.exists("Output/PD_suboscines_df.rds")){
  PD_null = list()
  
  for (i in 1:length(PulidoTrees)){
    phy = PulidoTrees[[i]]
    PD_null[[i]]<-picante::ses.mpd(suboscines_comm, cophenetic(phy), null.model =  "independentswap")#
  }
  saveRDS(PD_null, file = "Output/PD_suboscines_random.rds")# results for 1000 phylogenies
  PD_null_df = plyr::aaply(laply(PD_null, as.matrix), c(2, 3), mean)
  
  saveRDS(PD_null_df, file = "Output/PD_suboscines_df.rds")
}else {
  PD_suboscines_df = readRDS("Output/PD_suboscines_df.rds")
}

#

## FUNCTIONAL DIVERSITY ## -----------------------------------------------------------

## trait data

passer_traits<-manu_data[, c(23:42, 44:50)]; rownames(passer_traits)<-manu_data$species
oscines_traits<-manu_data[manu_data$suborder == "Passeri", c(23:42, 44:50)]; rownames(oscines_traits)<-manu_data[manu_data$suborder == "Passeri", 5]
suboscines_traits<-manu_data[!manu_data$suborder == "Passeri", c(23:42, 44:50)]; rownames(suboscines_traits)<-manu_data[!manu_data$suborder == "Passeri", 5]

# community data for FD

passer_comm = as.data.frame(t(manu_data[, c(10:22)])); colnames(passer_comm)<-manu_data$species
oscines_comm = as.data.frame(t(manu_data[manu_data$suborder == "Passeri", c(10:22)]));colnames(oscines_comm)<-manu_data[manu_data$suborder == "Passeri", 5]
suboscines_comm = as.data.frame(t(manu_data[!manu_data$suborder == "Passeri", c(10:22)]));colnames(suboscines_comm)<-manu_data[!manu_data$suborder == "Passeri", 5]

## 1. ALL TRAITS ##

## trait distance matrices
passer_all_dist<-cluster::daisy(passer_traits,  metric="gower", 
                                weights = c(1/2, 1/2, 1/5, 1/5, 1/5, 1/5, 1/5, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/5, 1/5, 1/5, 1/5, 1/5, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7))

oscines_all_dist<-cluster::daisy(oscines_traits,  metric="gower", 
                                 weights = c(1/2, 1/2, 1/5, 1/5, 1/5, 1/5, 1/5, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/5, 1/5, 1/5, 1/5, 1/5, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7))

suboscines_all_dist<-cluster::daisy(suboscines_traits,  metric="gower", 
                                    weights = c(1/2, 1/2, 1/5, 1/5, 1/5, 1/5, 1/5, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/5, 1/5, 1/5, 1/5, 1/5, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7))


if(!file.exists("Output/FD_passer_null_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(passer_all_dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(passer_all_dist, as.matrix(comm))
  }
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_null_df.rds")
}else {
  FD_passer_null_df = readRDS("Output/FD_passer_null_df.rds")
}


# oscines 

if(!file.exists("Output/FD_oscines_null_df.rds")){
  # Observed FD
  FD_oscines <- FD::fdisp(oscines_all_dist, as.matrix(oscines_comm))
  # null models
  FD_oscines_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_oscines_null[[i]]<-FD::fdisp(oscines_all_dist, as.matrix(comm))
  }
  FD_oscines_null_df = plyr::ldply(FD_oscines_null, function(x) x[[1]])
  FD_oscines_null_df = as.data.frame(t(FD_oscines_null_df))
  FD_oscines_null_df$fdis_obs = FD_oscines$FDis
  # null mean and sd
  FD_oscines_null_df$null_mean = select(FD_oscines_null_df, starts_with("V")) %>% rowMeans() 
  FD_oscines_null_df$null_sd = apply(select(FD_oscines_null_df, starts_with("V")), 1, sd)
  FD_oscines_null_df$ses = (FD_oscines_null_df$fdis_obs - FD_oscines_null_df$null_mean)/FD_oscines_null_df$null_sd
  saveRDS(FD_oscines_null_df, file = "Output/FD_oscines_null_df.rds")
}else {
  FD_oscines_null_df = readRDS("Output/FD_oscines_null_df.rds")
}

## suboscines ##
if(!file.exists("Output/FD_suboscines_null_df.rds")){
  # Observed FD
  FD_suboscines <- FD::fdisp(suboscines_all_dist, as.matrix(suboscines_comm))
  # null models
  FD_suboscines_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_suboscines_null[[i]]<-FD::fdisp(suboscines_all_dist, as.matrix(comm))
  }
  FD_suboscines_null_df = plyr::ldply(FD_suboscines_null, function(x) x[[1]])
  FD_suboscines_null_df = as.data.frame(t(FD_suboscines_null_df))
  FD_suboscines_null_df$fdis_obs = FD_suboscines$FDis
  # null mean and sd
  FD_suboscines_null_df$null_mean = select(FD_suboscines_null_df, starts_with("V")) %>% rowMeans() 
  FD_suboscines_null_df$null_sd = apply(select(FD_suboscines_null_df, starts_with("V")), 1, sd)
  FD_suboscines_null_df$ses = (FD_suboscines_null_df$fdis_obs - FD_suboscines_null_df$null_mean)/FD_suboscines_null_df$null_sd
  saveRDS(FD_suboscines_null_df, file = "Output/FD_suboscines_null_df.rds")
}else {
  FD_suboscines_null_df = readRDS("Output/FD_suboscines_null_df.rds")
}


## 2.  SIZE TRAITS ##

## trait distance matrices
dist<-cluster::daisy(as.data.frame(scale(passer_traits[,-c(3:27)])),  metric="euclidean")


if(!file.exists("Output/FD_passer_size_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_passer_size_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_size_df.rds")
}else {
  FD_passer_size_df = readRDS("Output/FD_passer_size_df.rds")
}

# oscines 
dist<-cluster::daisy(as.data.frame(scale(oscines_traits[,-c(3:27)])),  metric="euclidean")

if(!file.exists("Output/FD_oscines_size_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(oscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_oscines_size_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_oscines_size_df.rds")
}else {
  FD_oscines_size_df = readRDS("Output/FD_oscines_size_df.rds")
}



## suboscines ##
dist<-cluster::daisy(as.data.frame(scale(suboscines_traits[,-c(3:27)])),  metric="euclidean")

if(!file.exists("Output/FD_suboscines_size_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(suboscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_suboscines_size_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_suboscines_size_df.rds")
}else {
  FD_suboscines_size_df = readRDS("Output/FD_suboscines_size_df.rds")
}


## 3. DIET TRAITS

dist<-cluster::daisy(passer_traits[,-c(1,2, 8:27)],  metric="gower")


if(!file.exists("Output/FD_passer_diet_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_passer_diet_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_diet_df.rds")
}else {
  FD_passer_diet_df = readRDS("Output/FD_passer_diet_df.rds")
}

# oscines 
dist<-cluster::daisy(oscines_traits[,-c(1,2, 8:27)],  metric="gower")

if(!file.exists("Output/FD_oscines_diet_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(oscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_oscines_diet_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_oscines_diet_df.rds")
}else {
  FD_oscines_diet_df = readRDS("Output/FD_oscines_diet_df.rds")
}



## suboscines ##
dist<-cluster::daisy(suboscines_traits[,-c(1,2, 8:27)],  metric="gower")

if(!file.exists("Output/FD_suboscines_diet_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(suboscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_suboscines_diet_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_suboscines_diet_df.rds")
}else {
  FD_suboscines_diet_df = readRDS("Output/FD_suboscines_diet_df.rds")
}

## 4. FORAGING STRATEGY

dist<-cluster::daisy(passer_traits[,-c(1:7, 16:27)],  metric="gower")


if(!file.exists("Output/FD_passer_foraging_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_passer_foraging_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_foraging_df.rds")
}else {
  FD_passer_foraging_df = readRDS("Output/FD_passer_foraging_df.rds")
}



# oscines 
dist<-cluster::daisy(oscines_traits[,-c(1:7, 16:27)],  metric="gower")

if(!file.exists("Output/FD_oscines_foraging_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(oscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_oscines_foraging_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_oscines_foraging_df.rds")
}else {
  FD_oscines_foraging_df = readRDS("Output/FD_oscines_foraging_df.rds")
}



## suboscines ##
dist<-cluster::daisy(suboscines_traits[,-c(1:7, 16:27)],  metric="gower")

if(!file.exists("Output/FD_suboscines_foraging_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(suboscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_suboscines_foraging_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_suboscines_foraging_df.rds")
}else {
  FD_suboscines_foraging_df = readRDS("Output/FD_suboscines_foraging_df.rds")
}



## 5. FORAGING surface

dist<-cluster::daisy(passer_traits[,-c(1:15, 21:27)],  metric="gower")


if(!file.exists("Output/FD_passer_surface_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_passer_surface_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_surface_df.rds")
}else {
  FD_passer_surface_df = readRDS("Output/FD_passer_surface_df.rds")
}



# oscines 
dist<-cluster::daisy(oscines_traits[,-c(1:15, 21:27)],  metric="gower")

if(!file.exists("Output/FD_oscines_surface_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(oscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_oscines_surface_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_oscines_surface_df.rds")
}else {
  FD_oscines_surface_df = readRDS("Output/FD_oscines_surface_df.rds")
}



## suboscines ##
dist<-cluster::daisy(suboscines_traits[,-c(1:15, 21:27)],  metric="gower")

if(!file.exists("Output/FD_suboscines_surface_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(suboscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_suboscines_surface_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_suboscines_surface_df.rds")
}else {
  FD_suboscines_surface_df = readRDS("Output/FD_suboscines_surface_df.rds")
}



## 6. Morphology

dist<-cluster::daisy(as.data.frame(scale(passer_traits[,-c(1:20)])),  metric="euclidean")


if(!file.exists("Output/FD_passer_morpho_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(passer_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(passer_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_passer_morpho_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_passer_morpho_df.rds")
}else {
  FD_passer_morpho_df = readRDS("Output/FD_passer_morpho_df.rds")
}




# oscines 
dist<-cluster::daisy(as.data.frame(scale(oscines_traits[,-c(1:20)])),  metric="euclidean")

if(!file.exists("Output/FD_oscines_morpho_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(oscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(oscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_oscines_morpho_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_oscines_morpho_df.rds")
}else {
  FD_oscines_morpho_df = readRDS("Output/FD_oscines_morpho_df.rds")
}



## suboscines ##
dist<-cluster::daisy(as.data.frame(scale(suboscines_traits[,-c(1:20)])),  metric="euclidean")

if(!file.exists("Output/FD_suboscines_morpho_df.rds")){
  # Observed FD
  FD_passer <- FD::fdisp(dist, as.matrix(suboscines_comm))
  # null models
  FD_passer_null = list()
  for (i in 1:1000){
    set.seed(i)
    comm = picante::randomizeMatrix(suboscines_comm, null.model = c("independentswap"),iterations = 10000)
    FD_passer_null[[i]]<-FD::fdisp(dist, as.matrix(comm))
  }
  saveRDS(FD_passer_null, file = "Output/FD_subscines_morpho_null.rds")
  FD_passer_null_df = plyr::ldply(FD_passer_null, function(x) x[[1]])
  FD_passer_null_df = as.data.frame(t(FD_passer_null_df))
  FD_passer_null_df$fdis_obs = FD_passer$FDis
  # null mean and sd
  FD_passer_null_df$null_mean = select(FD_passer_null_df, starts_with("V")) %>% rowMeans() 
  FD_passer_null_df$null_sd = apply(select(FD_passer_null_df, starts_with("V")), 1, sd)
  FD_passer_null_df$ses = (FD_passer_null_df$fdis_obs - FD_passer_null_df$null_mean)/FD_passer_null_df$null_sd
  saveRDS(FD_passer_null_df, file = "Output/FD_suboscines_morpho_df.rds")
}else {
  FD_suboscines_morpho_df = readRDS("Output/FD_suboscines_morpho_df.rds")
}

## One dataset ##

FD<-FD_passer_null_df[,c(1001, 1002)]
FD$size<-FD_passer_size_df[,1001]
FD$diet<-FD_passer_diet_df[,1001]
FD$foraging<-FD_passer_foraging_df[,1001]
FD$surface<-FD_passer_surface_df[,1001]
FD$morpho<-FD_passer_morpho_df[,1001]

FD$os_all<-FD_oscines_null_df[,c(1001)]
FD$os_size<-FD_oscines_size_df[,1001]
FD$os_diet<-FD_oscines_diet_df[,1001]
FD$os_foraging<-FD_oscines_foraging_df[,1001]
FD$os_surface<-FD_oscines_surface_df[,1001]
FD$os_morpho<-FD_oscines_morpho_df[,1001]

FD$subos_all<-FD_suboscines_null_df[,c(1001)]
FD$subos_size<-FD_suboscines_size_df[,1001]
FD$subos_diet<-FD_suboscines_diet_df[,1001]
FD$subos_foraging<-FD_suboscines_foraging_df[,1001]
FD$subos_surface<-FD_suboscines_surface_df[,1001]
FD$subos_morpho<-FD_suboscines_morpho_df[,1001]
FD<-FD[,-2]


## predicted by random
FD$pred_all<-FD_passer_null_df[,c(1002)]## 
FD$pred_size<-FD_passer_size_df[,1002]
FD$pred_diet<-FD_passer_diet_df[,1002]
FD$pred_foraging<-FD_passer_foraging_df[,1002]
FD$pred_surface<-FD_passer_surface_df[,1002]
FD$pred_morpho<-FD_passer_morpho_df[,1002]

FD$pred_os_all<-FD_oscines_null_df[,c(1002)]
FD$pred_os_size<-FD_oscines_size_df[,1002]
FD$pred_os_diet<-FD_oscines_diet_df[,1002]
FD$pred_os_foraging<-FD_oscines_foraging_df[,1002]
FD$pred_os_surface<-FD_oscines_surface_df[,1002]
FD$pred_os_morpho<-FD_oscines_morpho_df[,1002]

FD$pred_subos_all<-FD_suboscines_null_df[,c(1002)]
FD$pred_subos_size<-FD_suboscines_size_df[,1002]
FD$pred_subos_diet<-FD_suboscines_diet_df[,1002]
FD$pred_subos_foraging<-FD_suboscines_foraging_df[,1002]
FD$pred_subos_surface<-FD_suboscines_surface_df[,1002]
FD$pred_subos_morpho<-FD_suboscines_morpho_df[,1002]

## FDPD
FDPD = FD

FDPD$PDpasser<-as.data.frame(PD_null_df)$mpd.obs #observed values, average of 1000 trees
FDPD$PDoscines<-as.data.frame(PD_oscines_df)$mpd.obs #observed values, average of 1000 trees
FDPD$PDsuboscines<-as.data.frame(PD_suboscines_df)$mpd.obs #observed values, average of 1000 trees

FDPD$PDpasser_pred<-as.data.frame(PD_null_df)$mpd.rand.mean #predicted values by 1000 randomizations for each tree, average of 1000 trees
FDPD$PDoscines_pred<-as.data.frame(PD_oscines_df)$mpd.rand.mean #predicted values by 1000 randomizations for each tree, average of 1000 trees
FDPD$PDsuboscines_pred<-as.data.frame(PD_suboscines_df)$mpd.rand.mean #predicted values by 1000 randomizations for each tree, average of 1000 trees


FDPD$n_passer = as.data.frame(PD_null_df)$ntaxa
FDPD$n_oscines = as.data.frame(PD_oscines_df)$ntaxa
FDPD$n_suboscines = as.data.frame(PD_suboscines_df)$ntaxa

#saveRDS(FDPD, "FDPD.rds")



### Polynomial regressions #### -----------------------------------------------------

## Raw polynomial and orthogonal polynomial regressions ##

elev = c(500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500)

elev = scale(elev)

### I. REGRESSIONS WITH OBSERVED DATA ###

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(FDPD)));names(coef_raw) = c("linear", "squared")
coef_raw_SE = as.data.frame(matrix(ncol  = 2, nrow = length(FDPD)));names(coef_raw_SE) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(FDPD)));names(coef_ortho) = c("linear", "squared")
coef_ortho_SE = as.data.frame(matrix(ncol  = 2, nrow = length(FDPD)));names(coef_ortho) =  c("linear", "squared")

for(i in 1:length(FDPD)){
  coef_raw[i,] = summary(lm(FDPD[,i] ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(FDPD[,i] ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
  coef_raw_SE[i,] = summary(lm(FDPD[,i] ~poly(elev, 2, raw = T)))$coefficients[c(2:3),2]
  coef_ortho_SE[i,] = summary(lm(FDPD[,i] ~poly(elev, 2, raw = F)))$coefficients[c(2:3),2]
}

coef_raw_df = coef_raw
coef_ortho_df = coef_ortho

coef_raw_SE_df = coef_raw_SE
coef_ortho_SE_df = coef_raw_SE

rownames(coef_raw_df)<-colnames(FDPD)
rownames(coef_ortho_df)<-colnames(FDPD)

rownames(coef_raw_SE_df)<-colnames(FDPD)
rownames(coef_ortho_SE_df)<-colnames(FDPD)

names(FDPD)


## REGRESSIONS OF RANDOMIZED ASSEMBLAGES ###

## 1. PHYLOGENETIC DIVERSITY

## Passeriformes 

data4reg<-readRDS("Output/PD_passer_random.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_PDpasser = coef_raw
coef_ortho_PDpasser = coef_ortho

## Oscines 

data4reg<-readRDS("Output/PD_oscines_random.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_PDoscines = coef_raw
coef_ortho_PDoscines = coef_ortho

# Suboscines

data4reg<-readRDS("Output/PD_suboscines_random.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$mpd.rand.mean ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_PDsuboscines = coef_raw
coef_ortho_PDsuboscines = coef_ortho

## 2. FUNCTIONAL DIVERSITY --ALL TRAITS

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser = coef_raw
coef_ortho_FDpasser = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines = coef_raw
coef_ortho_FDoscines = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_suboscines_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines = coef_raw
coef_ortho_FDsuboscines = coef_ortho


## 3. FUNCTIONAL DIVERSITY --size 

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_size_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser_size = coef_raw
coef_ortho_FDpasser_size = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_size_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines_size = coef_raw
coef_ortho_FDoscines_size = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_suboscines_size_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines_size = coef_raw
coef_ortho_FDsuboscines_size = coef_ortho


## 3. FUNCTIONAL DIVERSITY -- diet

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_diet_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser_diet = coef_raw
coef_ortho_FDpasser_diet = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_diet_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines_diet = coef_raw
coef_ortho_FDoscines_diet = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_suboscines_diet_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines_diet = coef_raw
coef_ortho_FDsuboscines_diet = coef_ortho

## 3. FUNCTIONAL DIVERSITY -- foraging strategy

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_foraging_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser_foraging = coef_raw
coef_ortho_FDpasser_foraging = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_foraging_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines_foraging = coef_raw
coef_ortho_FDoscines_foraging = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_suboscines_foraging_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines_foraging = coef_raw
coef_ortho_FDsuboscines_foraging = coef_ortho


## 3. FUNCTIONAL DIVERSITY -- foraging surface 

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_surface_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser_surface = coef_raw
coef_ortho_FDpasser_surface = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_surface_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines_surface = coef_raw
coef_ortho_FDoscines_surface = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_suboscines_surface_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines_surface = coef_raw
coef_ortho_FDsuboscines_surface = coef_ortho

## 3. FUNCTIONAL DIVERSITY -- morphology

## Passeriformes 

data4reg<-readRDS("Output/FD_passer_morpho_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDpasser_morpho = coef_raw
coef_ortho_FDpasser_morpho = coef_ortho


## Oscines

data4reg<-readRDS("Output/FD_oscines_morpho_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDoscines_morpho = coef_raw
coef_ortho_FDoscines_morpho = coef_ortho

## Suboscines

data4reg<-readRDS("Output/FD_subscines_morpho_null.rds")# 1000 dataframes, with mean of 1000 randomization each

coef_raw = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_raw) = c("linear", "squared")
coef_ortho = as.data.frame(matrix(ncol  = 2, nrow = length(data4reg)));names(coef_ortho) = c("linear", "squared")

for(i in 1:length(data4reg)){
  df = data4reg[[i]]
  coef_raw[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = T)))$coefficients[c(2:3),1]
  coef_ortho[i,] = summary(lm(df$FDis ~poly(elev, 2, raw = F)))$coefficients[c(2:3),1]
}

coef_raw_FDsuboscines_morpho = coef_raw
coef_ortho_FDsuboscines_morpho = coef_ortho


##dataframes for comparison ###

coef_df<-coef_raw_df[c(1:18, 37:39),]

coef_df$linear_0.025 <-NA
coef_df$linear_0.975<-NA
coef_df$squared_0.025 <-NA
coef_df$squared_0.975<-NA


coef_df[1,3:4]<-quantile(coef_raw_FDpasser$linear,c(0.025, 0.975))
coef_df[2,3:4]<-quantile(coef_raw_FDpasser_size$linear,c(0.025, 0.975))
coef_df[3,3:4]<-quantile(coef_raw_FDpasser_diet$linear,c(0.025, 0.975))
coef_df[4,3:4]<-quantile(coef_raw_FDpasser_foraging$linear,c(0.025, 0.975))
coef_df[5,3:4]<-quantile(coef_raw_FDpasser_surface$linear,c(0.025, 0.975))
coef_df[6,3:4]<-quantile(coef_raw_FDpasser_morpho$linear,c(0.025, 0.975))

coef_df[7,3:4]<-quantile(coef_raw_FDoscines$linear,c(0.025, 0.975))
coef_df[8,3:4]<-quantile(coef_raw_FDoscines_size$linear,c(0.025, 0.975))
coef_df[9,3:4]<-quantile(coef_raw_FDoscines_diet$linear,c(0.025, 0.975))
coef_df[10,3:4]<-quantile(coef_raw_FDoscines_foraging$linear,c(0.025, 0.975))
coef_df[11,3:4]<-quantile(coef_raw_FDoscines_surface$linear,c(0.025, 0.975))
coef_df[12,3:4]<-quantile(coef_raw_FDoscines_morpho$linear,c(0.025, 0.975))

coef_df[13,3:4]<-quantile(coef_raw_FDsuboscines$linear,c(0.025, 0.975))
coef_df[14,3:4]<-quantile(coef_raw_FDsuboscines_size$linear,c(0.025, 0.975))
coef_df[15,3:4]<-quantile(coef_raw_FDsuboscines_diet$linear,c(0.025, 0.975))
coef_df[16,3:4]<-quantile(coef_raw_FDsuboscines_foraging$linear,c(0.025, 0.975))
coef_df[17,3:4]<-quantile(coef_raw_FDsuboscines_surface$linear,c(0.025, 0.975))
coef_df[18,3:4]<-quantile(coef_raw_FDsuboscines_morpho$linear,c(0.025, 0.975))


coef_df[19,3:4]<-quantile(coef_raw_PDpasser$linear,c(0.025, 0.975))
coef_df[20,3:4]<-quantile(coef_raw_PDoscines$linear,c(0.025, 0.975))
coef_df[21,3:4]<-quantile(coef_raw_PDsuboscines$linear,c(0.025, 0.975))


coef_df[1,5:6]<-quantile(coef_raw_FDpasser$squared,c(0.025, 0.975))
coef_df[2,5:6]<-quantile(coef_raw_FDpasser_size$squared,c(0.025, 0.975))
coef_df[3,5:6]<-quantile(coef_raw_FDpasser_diet$squared,c(0.025, 0.975))
coef_df[4,5:6]<-quantile(coef_raw_FDpasser_foraging$squared,c(0.025, 0.975))
coef_df[5,5:6]<-quantile(coef_raw_FDpasser_surface$squared,c(0.025, 0.975))
coef_df[6,5:6]<-quantile(coef_raw_FDpasser_morpho$squared,c(0.025, 0.975))

coef_df[7,5:6]<-quantile(coef_raw_FDoscines$squared,c(0.025, 0.975))
coef_df[8,5:6]<-quantile(coef_raw_FDoscines_size$squared,c(0.025, 0.975))
coef_df[9,5:6]<-quantile(coef_raw_FDoscines_diet$squared,c(0.025, 0.975))
coef_df[10,5:6]<-quantile(coef_raw_FDoscines_foraging$squared,c(0.025, 0.975))
coef_df[11,5:6]<-quantile(coef_raw_FDoscines_surface$squared,c(0.025, 0.975))
coef_df[12,5:6]<-quantile(coef_raw_FDoscines_morpho$squared,c(0.025, 0.975))

coef_df[13,5:6]<-quantile(coef_raw_FDsuboscines$squared,c(0.025, 0.975))
coef_df[14,5:6]<-quantile(coef_raw_FDsuboscines_size$squared,c(0.025, 0.975))
coef_df[15,5:6]<-quantile(coef_raw_FDsuboscines_diet$squared,c(0.025, 0.975))
coef_df[16,5:6]<-quantile(coef_raw_FDsuboscines_foraging$squared,c(0.025, 0.975))
coef_df[17,5:6]<-quantile(coef_raw_FDsuboscines_surface$squared,c(0.025, 0.975))
coef_df[18,5:6]<-quantile(coef_raw_FDsuboscines_morpho$squared,c(0.025, 0.975))


coef_df[19,5:6]<-quantile(coef_raw_PDpasser$squared,c(0.025, 0.975))
coef_df[20,5:6]<-quantile(coef_raw_PDoscines$squared,c(0.025, 0.975))
coef_df[21,5:6]<-quantile(coef_raw_PDsuboscines$squared,c(0.025, 0.975))

## within CI or not

coef_df$linear_sig<-NA
coef_df$squared_sig<-NA

coef_df$linear_sig<-as.numeric(!dplyr::between(coef_df$linear, coef_df$linear_0.025, coef_df$linear_0.975))# 0 not significant, 1 significant
coef_df$squared_sig<-as.numeric(!dplyr::between(coef_df$squared, coef_df$squared_0.025, coef_df$squared_0.975))# 0 not significant, 1 significant

coef_df_raw<-coef_df

#write.csv(coef_df_raw, "coef_df_raw.csv")

### orthogonal

coef_df<-coef_ortho_df[c(1:18, 37:39),]

coef_df$linear_0.025 <-NA
coef_df$linear_0.975<-NA
coef_df$squared_0.025 <-NA
coef_df$squared_0.975<-NA


coef_df[1,3:4]<-quantile(coef_ortho_FDpasser$linear,c(0.025, 0.975))
coef_df[2,3:4]<-quantile(coef_ortho_FDpasser_size$linear,c(0.025, 0.975))
coef_df[3,3:4]<-quantile(coef_ortho_FDpasser_diet$linear,c(0.025, 0.975))
coef_df[4,3:4]<-quantile(coef_ortho_FDpasser_foraging$linear,c(0.025, 0.975))
coef_df[5,3:4]<-quantile(coef_ortho_FDpasser_surface$linear,c(0.025, 0.975))
coef_df[6,3:4]<-quantile(coef_ortho_FDpasser_morpho$linear,c(0.025, 0.975))

coef_df[7,3:4]<-quantile(coef_ortho_FDoscines$linear,c(0.025, 0.975))
coef_df[8,3:4]<-quantile(coef_ortho_FDoscines_size$linear,c(0.025, 0.975))
coef_df[9,3:4]<-quantile(coef_ortho_FDoscines_diet$linear,c(0.025, 0.975))
coef_df[10,3:4]<-quantile(coef_ortho_FDoscines_foraging$linear,c(0.025, 0.975))
coef_df[11,3:4]<-quantile(coef_ortho_FDoscines_surface$linear,c(0.025, 0.975))
coef_df[12,3:4]<-quantile(coef_ortho_FDoscines_morpho$linear,c(0.025, 0.975))

coef_df[13,3:4]<-quantile(coef_ortho_FDsuboscines$linear,c(0.025, 0.975))
coef_df[14,3:4]<-quantile(coef_ortho_FDsuboscines_size$linear,c(0.025, 0.975))
coef_df[15,3:4]<-quantile(coef_ortho_FDsuboscines_diet$linear,c(0.025, 0.975))
coef_df[16,3:4]<-quantile(coef_ortho_FDsuboscines_foraging$linear,c(0.025, 0.975))
coef_df[17,3:4]<-quantile(coef_ortho_FDsuboscines_surface$linear,c(0.025, 0.975))
coef_df[18,3:4]<-quantile(coef_ortho_FDsuboscines_morpho$linear,c(0.025, 0.975))


coef_df[19,3:4]<-quantile(coef_ortho_PDpasser$linear,c(0.025, 0.975))
coef_df[20,3:4]<-quantile(coef_ortho_PDoscines$linear,c(0.025, 0.975))
coef_df[21,3:4]<-quantile(coef_ortho_PDsuboscines$linear,c(0.025, 0.975))


coef_df[1,5:6]<-quantile(coef_ortho_FDpasser$squared,c(0.025, 0.975))
coef_df[2,5:6]<-quantile(coef_ortho_FDpasser_size$squared,c(0.025, 0.975))
coef_df[3,5:6]<-quantile(coef_ortho_FDpasser_diet$squared,c(0.025, 0.975))
coef_df[4,5:6]<-quantile(coef_ortho_FDpasser_foraging$squared,c(0.025, 0.975))
coef_df[5,5:6]<-quantile(coef_ortho_FDpasser_surface$squared,c(0.025, 0.975))
coef_df[6,5:6]<-quantile(coef_ortho_FDpasser_morpho$squared,c(0.025, 0.975))

coef_df[7,5:6]<-quantile(coef_ortho_FDoscines$squared,c(0.025, 0.975))
coef_df[8,5:6]<-quantile(coef_ortho_FDoscines_size$squared,c(0.025, 0.975))
coef_df[9,5:6]<-quantile(coef_ortho_FDoscines_diet$squared,c(0.025, 0.975))
coef_df[10,5:6]<-quantile(coef_ortho_FDoscines_foraging$squared,c(0.025, 0.975))
coef_df[11,5:6]<-quantile(coef_ortho_FDoscines_surface$squared,c(0.025, 0.975))
coef_df[12,5:6]<-quantile(coef_ortho_FDoscines_morpho$squared,c(0.025, 0.975))

coef_df[13,5:6]<-quantile(coef_ortho_FDsuboscines$squared,c(0.025, 0.975))
coef_df[14,5:6]<-quantile(coef_ortho_FDsuboscines_size$squared,c(0.025, 0.975))
coef_df[15,5:6]<-quantile(coef_ortho_FDsuboscines_diet$squared,c(0.025, 0.975))
coef_df[16,5:6]<-quantile(coef_ortho_FDsuboscines_foraging$squared,c(0.025, 0.975))
coef_df[17,5:6]<-quantile(coef_ortho_FDsuboscines_surface$squared,c(0.025, 0.975))
coef_df[18,5:6]<-quantile(coef_ortho_FDsuboscines_morpho$squared,c(0.025, 0.975))


coef_df[19,5:6]<-quantile(coef_ortho_PDpasser$squared,c(0.025, 0.975))
coef_df[20,5:6]<-quantile(coef_ortho_PDoscines$squared,c(0.025, 0.975))
coef_df[21,5:6]<-quantile(coef_ortho_PDsuboscines$squared,c(0.025, 0.975))

## within CI or not

coef_df$linear_sig<-NA
coef_df$squared_sig<-NA

coef_df$linear_sig<-as.numeric(!dplyr::between(coef_df$linear, coef_df$linear_0.025, coef_df$linear_0.975))# 0 not significant, 1 significant
coef_df$squared_sig<-as.numeric(!dplyr::between(coef_df$squared, coef_df$squared_0.025, coef_df$squared_0.975))# 0 not significant, 1 significant

coef_df_ortho<-coef_df

#write.csv(coef_df_ortho, "coef_df_ortho.csv")

