# Process data for the Great Lakes

# Wildlife data
data_dir <- "D:/Current/Great_Lakes/Data - wildlife"
state_list <- list(
  'Illinois', 'Indiana', 'Michigan', 'Minnesota', 'New_York', 'Ohio',
  'Pennsylvania', 'Wisconsin')

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

# leftover remnants
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

spp_by_state <- combined_df[
  !duplicated(combined_df$scientificName),
  c('scientificName', 'common_name', 'taxon')]

# iucn habitat data
iucn_habitat_df <- read.csv(
  paste(data_dir, 'IUCN_habitat_threat', 'habitats.csv', sep='/'))
iucn_suitable <- iucn_habitat_df[
  iucn_habitat_df$suitability != 'Marginal',
  c('scientificName', 'code', 'name')]
colnames(iucn_suitable) <- c('scientificName', 'habitatCode', 'habitatName')

# match iucn species with SWAP species
swap_habitats <- merge(spp_by_state, iucn_suitable, by='scientificName')
swap_habitats$habitatPrefix <- sapply(
  strsplit(as.character(swap_habitats$habitatCode),'.', fixed=TRUE), "[", 1)
write.csv(swap_habitats,
          'D:/Current/Great_Lakes/Data - wildlife/SWAP/swap_iucn_habitats.csv',
          row.names=FALSE)

# find threatened/endangered species not appearing in IUCN habitats data
esl_df$in_IUCN_habitats <- 0
esl_df$in_SWAP <- 0
for (r in 1:NROW(esl_df)) {
  if (esl_df[r, 'scientificName'] %in% iucn_suitable$scientificName) {
    esl_df[r, 'in_IUCN_habitats'] <- 1
  }
  if (esl_df[r, 'scientificName'] %in% spp_by_state$scientificName) {
    esl_df[r, 'in_SWAP'] <- 1
  }
}
missing_from_iucn <- esl_df[esl_df$in_IUCN_habitats == 0, ]
write.csv(missing_from_iucn,
          "D:/Current/Great_Lakes/Data - wildlife/FWS_Endangered_Spp_Database/missing_from_IUCN_habitats.csv",
          row.names=FALSE)

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


