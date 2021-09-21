# Process data for the Great Lakes

# Wildlife data
data_dir <- "D:/NFWF_PhaseIII/Great_Lakes/Data - wildlife"
state_list <- list(
  'Illinois', 'Indiana', 'Michigan', 'Minnesota', 'New_York', 'Ohio',
  'Pennsylvania', 'Wisconsin')

# threatened and endangered spp for each state
esl_dir <- paste(data_dir, 'FWS_Endangered_Spp_Database', sep='/')
df_list <- list()
for (state in state_list) {
  esl_path <- paste(
    esl_dir,
    paste(state, '_species-listings-by-state-report.csv', sep=''), sep='/')
  state_df <- read.csv(esl_path)
  state_df$state <- state
  df_list[[state]] <- state_df
}
esl_df <- do.call(rbind, df_list)
colnames(esl_df) <- c(
  'scientificName', 'common_name', 'whereListed', 'region', 'status', 'group',
  'state')
esl_df <- esl_df[
  (esl_df$group != "Flowering Plants") &
    (esl_df$group != "Ferns and Allies") &
  !duplicated(esl_df$scientificName), 
  c('scientificName', 'common_name', 'whereListed', 'region', 'status',
    'group')]

# species list from State Wildlife Action Plans for each state
swap_dir <- paste(data_dir, 'SWAP', sep='/')
df_list <- list()
for (state in state_list) {
  swap_path <- paste(swap_dir, paste('SWAP_', state, '.csv', sep=''), sep='/')
  state_df <- read.csv(swap_path)
  state_df$state <- state
  df_list[[state]] <- state_df
}
combined_df <- do.call(rbind, df_list)
colnames(combined_df) <- c(
  'scientificName', 'common_name', 'in_2005', 'in_2015', 'taxon', 'match',
  'state')

# in which states does each spp appear?
spp_list <- unique(combined_df$scientificName)
spp_by_state_list <- list()
for (sp in spp_list) {
  sp_subs <- combined_df[(
    combined_df$scientificName == sp) & 
      !duplicated(combined_df$scientificName),
    c('scientificName', 'common_name', 'taxon')]
  sp_subs$states <- paste(
    combined_df[combined_df$scientificName == sp, 'state'], collapse=', ')
  spp_by_state_list[[sp]] <- sp_subs
}
spp_by_state <- do.call(rbind, spp_by_state_list)

# find threatened/endangered species not appearing in SWAP
esl_df$in_SWAP <- 0
for (r in 1:NROW(esl_df)) {
  if (esl_df[r, 'scientificName'] %in% spp_by_state$scientificName) {
    esl_df[r, 'in_SWAP'] <- 1
  }
}
missing_from_swap <- esl_df[esl_df$in_SWAP == 0, ]
missing_from_swap <- missing_from_swap[, c(
  "scientificName", "common_name", "group", "status")]
missing_from_swap$states <- 'federal'
colnames(missing_from_swap)[3] <- 'taxon'

# combine with spp by state
spp_by_state$status <- 'SWAP'
spp_by_state <- spp_by_state[, colnames(missing_from_swap)]
full_spp_list <- rbind(spp_by_state, missing_from_swap)

# aquatic species only
aquatic_spp <- full_spp_list[(full_spp_list$taxon == 'Fish') |
                               (full_spp_list$taxon == 'Fishes') |
                               (full_spp_list$taxon == 'Clams') |
                               (full_spp_list$taxon == 'Mollusks'), ]
write.csv(aquatic_spp, paste(data_dir, "aquatic_spp_list.csv", sep='/'),
          row.names=FALSE)
aquatic_spp <- read.csv(paste(data_dir, "aquatic_spp_list.csv", sep='/'))

# iucn habitat data
iucn_habitat_df <- read.csv(
  paste(data_dir, 'IUCN_habitat_threat', 'habitats.csv', sep='/'))
iucn_suitable <- iucn_habitat_df[
  iucn_habitat_df$suitability != 'Marginal',
  c('scientificName', 'code', 'name')]
colnames(iucn_suitable) <- c('scientificName', 'habitatCode', 'habitatName')

# match iucn species with SWAP and T&E species
all_spp_habitats <- merge(
  full_spp_list, iucn_suitable, by='scientificName', all.x=TRUE)
write.csv(
  all_spp_habitats,
  paste(data_dir, 'IUCN_all_spp_habitat_list.csv', sep='/'), row.names=FALSE)
aquatic_habitats <- merge(
  aquatic_spp, iucn_suitable, by='scientificName', all.x=TRUE)
length(unique(aquatic_habitats$habitatCode))  # 39 unique habitat types for aquatic spp

# match IUCN habitat classification with NWI classification
aquatic_habitat_list <- aquatic_habitats[
  (!duplicated(aquatic_habitats[,  c("habitatCode", "habitatName")])) &
    (!is.na(aquatic_habitats$habitatCode)),
  c("habitatCode", "habitatName")]
write.csv(
  aquatic_habitat_list,
  paste(data_dir, 'IUCN_aquatic_habitat_list.csv', sep='/'), row.names=FALSE)
nwi_codes_df <- read.csv("D:/Datasets/FWS/NWI/NWI-Code-Definitions-Oct-19/NWI_Code_Definitions/NWI_Code_Definitions.csv")

# match IUCN habitat classifications with C-CAP classification
iucn_classification <- iucn_suitable[
  !duplicated(iucn_suitable$habitatCode), c('habitatCode', 'habitatName')]
iucn_classification$habitatLargeClass <- sapply(
  strsplit(as.character(iucn_classification$habitatName),'-', fixed=TRUE),
  "[", 1)
write.csv(iucn_classification,
          "D:/Current/Great_Lakes/Data - wildlife/IUCN_habitat_threat/IUCN_habitat_classes.csv",
          row.names=FALSE)

# threats data from IUCN
iucn_threat_df <- read.csv(
  paste(data_dir, 'IUCN_habitat_threat', 'threats.csv', sep='/'))
iucn_threat_df <- iucn_threat_df[, c(
  c("scientificName", "name", "stressCode", "stressName"))]
colnames(iucn_threat_df) <- c(
  'scientificName', 'threatName', 'stressCode', 'stressName')

# match threats classifications with C-CAP classification
swap_threats <- merge(
  spp_by_state, iucn_threat_df, by='scientificName')
threat_classification <- swap_threats[
  !duplicated(swap_threats$threatName),
              c('threatName', 'stressCode', 'stressName')]
write.csv(threat_classification,
          "D:/Current/Great_Lakes/Data - wildlife/IUCN_habitat_threat/IUCN_threat_classes.csv",
          row.names=FALSE)

# collect all C-CAP landcover codes by taxonomic group
swap_match_df <- read.csv("D:/Current/Great_Lakes/Data - wildlife/IUCN_habitat_threat/IUCN_habitat_classes.csv")
swap_habitats <- read.csv('D:/Current/Great_Lakes/Data - wildlife/SWAP/swap_iucn_habitats.csv')
swap_ccap_habitats <- merge(swap_habitats, swap_match_df, all.x=TRUE)
swap_ccap_habitats <- swap_ccap_habitats[, c(
  "scientificName", 'common_name', 'taxon', "C.CAP_value",
  "C.CAP_description")]
swap_ccap_habitats$source <- 'SWAP'
esl_match_df <- read.csv("D:/Current/Great_Lakes/Data - wildlife/FWS_Endangered_Spp_Database/missing_from_IUCN_habitats.csv")
esl_match_df <- esl_match_df[, c(
  "scientificName", 'common_name', 'group', "C.CAP_value",
  "C.CAP_description")]
colnames(esl_match_df)[3] <- 'taxon'
esl_match_df[(esl_match_df$taxon == 'Clams'), 'taxon'] <- 'Mollusks'
esl_match_df[(esl_match_df$taxon == 'Snails'), 'taxon'] <- 'Mollusks'
esl_match_df[(esl_match_df$taxon == 'Fishes'), 'taxon'] <- 'Fish'
esl_match_df[(esl_match_df$taxon == 'Crustaceans'), 'taxon'] <- 'Other Invertebrates'
esl_match_df[(esl_match_df$taxon == 'Insects'), 'taxon'] <- 'Other Invertebrates'
esl_match_df$source <- 'ESL'
habitat_spp_df <- rbind(swap_ccap_habitats, esl_match_df)
habitat_taxon_summary <- unique(
  habitat_spp_df[c("taxon", "C.CAP_value", "C.CAP_description")])
write.csv(
  habitat_taxon_summary,
  paste(data_dir, 'C-CAP_landcover_code_by_taxon.csv', sep='/'),
  row.names=FALSE)