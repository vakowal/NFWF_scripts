# compare species included in terrestrial index with spp in 2018 revised
# Guam Wildlife Action Plan (GWAP)

# species appearing in the 2018 revised GWAP
spp_2018 <- read.csv("C:/Users/Ginger/Documents/Guam/Data - Wildlife Index/GWAP/2019_report_spp.csv")

# species appearing in the list downloaded from USGS SWAP Species Conservation Analysis Tool
spp_download <- read.csv("C:/Users/Ginger/Documents/Guam/Data - Wildlife Index/GWAP/species_greatest_conservation_need_Guam_download_2.3.21.csv")
# leave out plants
spp_download <- spp_download[spp_download$Taxonomic.Group != 'Plants', ]
# Kim's species table
spp_kim <- read.csv("C:/Users/Ginger/Documents/Guam/Data - Wildlife Index/GWAP/GWAP_2015_vs_2019.csv")

# species habitat list from IUCN
habitat_iucn <- read.csv("C:/Users/Ginger/Documents/Guam/Data - Wildlife Index/GWAP/IUCN_download/habitats.csv")
taxon_iucn <- read.csv("C:/Users/Ginger/Documents/Guam/Data - Wildlife Index/GWAP/IUCN_download/taxonomy.csv")

# THROWAWAY: make sure the spp in habitat_iucn and taxon_iucn correspond
spp_taxon <- unique(taxon_iucn$scientificName)
spp_habitat <- unique(habitat_iucn$scientificName)
diff1 <- setdiff(spp_taxon, spp_habitat)
diff2 <- setdiff(spp_habitat, spp_taxon)

diffk <- setdiff(scinames_kim, spp_taxon)

# sets: look for spp occurring in one source but not another
# setdiff(x, y)  # in x, not in y
scinames_kim <- spp_kim$Scientific.Name
scinames_download <- spp_download$Scientific.Name.Reported.in.State.SWAP
scinames_2018 <- spp_2018$Scientific.Name
scinames_iucn <- unique(spp_iucn$scientificName)

spp_download$not_in_kim <- 1
for(r in c(1:NROW(spp_download))) {
  if (spp_download[r, 'Scientific.Name.Reported.in.State.SWAP'] %in%
      spp_kim$Scientific.Name) {
    spp_download[r, 'not_in_kim'] <- 0
  }
}
spp_download$not_in_iucn <- 1
for(r in c(1:NROW(spp_download))) {
  if (spp_download[r, 'Scientific.Name.Reported.in.State.SWAP'] %in%
      scinames_iucn) {
    spp_download[r, 'not_in_iucn'] <- 0
  }
}
spp_iucn$not_in_download <- 1
for(r in c(1:NROW(spp_download))) {
  if (spp_iucn[r, 'scientificName'] %in%
      scinames_download) {
    spp_iucn[r, 'not_in_download'] <- 0
  }
}

spp_kim$not_in_download <- 1
for(r in c(1:NROW(spp_kim))) {
  if (spp_kim[r, 'Scientific.Name']  %in%
      spp_download$Scientific.Name.Reported.in.State.SWAP) {
    spp_kim[r, 'not_in_download'] <- 0
  }
}
in_download_not_kim <- setdiff(scinames_download, scinames_kim)
in_kim_not_download <- setdiff(scinames_kim, scinames_download)

in_2018_not_download <- setdiff(scinames_2018, scinames_download)
in_download_not_2018 <- setdiff(scinames_download, scinames_2018)
