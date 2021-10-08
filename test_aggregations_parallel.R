library(tidyverse)
library(multidplyr)
library(data.table)
library(dtplyr)

#tables_folder <- 'D:/Data/Impacts_Data/tables'
tables_folder <- '/home/rstudio/data/impacts_data/tables'


###############################################################################
# Join various tables together
sls <- data.table(readRDS(file.path(tables_folder, 'sls.rds')))
countries <- data.table(readRDS(file.path(tables_folder, 'countries.rds')))
divisions <- data.table(readRDS(file.path(tables_folder, 'divisions.rds')))
setkey(sls, id)
setkey(countries, id)
setkey(divisions, id)

slsbysite <- data.table(readRDS(file.path(tables_folder, 'slsbysite.rds')))
countrybysite <- data.table(readRDS(file.path(tables_folder, 'countrybysite.rds')))
divisionbysite <- data.table(readRDS(file.path(tables_folder, 'divisionbysite.rds')))
setkey(slsbysite, sls_id)
setkey(countrybysite, country_id)
setkey(divisionbysite, division_id)

sls <- sls[slsbysite, ]
divisions<- divisions[divisionbysite, ]
countries <- countries[countrybysite, ]
sls[, id:=NULL]
divisions[, id:=NULL]
countries[, id:=NULL]
setnames(sls, "name", "sls_name")
setnames(countries, "name", "country_name")
setnames(divisions, "name", "division_name")

countries <- countries[, .(site_id, country_name)]

sites <- data.table(readRDS(file.path(tables_folder, 'sites.rds')))
# Remove shape column since it is large and not needed
sites[, shape:=NULL]

setkey(sites, id)
setkey(sls, site_id)
setkey(divisions, site_id)
setkey(countries, site_id)

sites <- sites[sls, nomatch=0]
sites <- sites[countries, nomatch=0]
sites <- sites[divisions, nomatch=0]

# Filter out blue nature alliance
bna_name <- 'C4O - Blue Nature Alliance'
table(sites$division_name == bna_name)
sites <- sites[sites$division_name != bna_name, ]

pixels <- data.table(readRDS(file.path(tables_folder, 'pixels_with_seq.rds')))
pixelsbysite <- data.table(readRDS(file.path(tables_folder, 'pixelsbysite.rds')))


# Merge pixels and pixelsbysite
setkey(pixels, id)
setkey(pixelsbysite, pixel_id)
pixels <- pixels[pixelsbysite, nomatch=0]

# Add site columns
setkey(pixels, site_id)
setkey(sites, id)
pixels <- pixels[sites, nomatch=0, allow.cartesian=TRUE]

###############################################################################
###  Partition across CPUs

cluster  <- new_cluster(12)
# Partition by site
pixels %>%
  as_tibble() %>%
  group_by(id, site_id, biome, new_or_continued_1) %>%
  partition(cluster) -> partied_data


###############################################################################
###  Run stats by new/continued and by biome


nrow(pixels)
partied_data %>%
    group_by(id) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(biome, new_or_continued_1) %>%
    mutate(high_irr_c = c_tstor_ic > 25,
           any_irr_c = c_tstor_ic > .01,
           under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))) %>%
    summarise(
        total_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
        restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
        woody_carbon_stored=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
        soil_carbon_stored=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
        irr_carbon_stored=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
        irr_carbon_ha_any=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
        irr_carbon_ha_high=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
        carbon_seq=sum(area_ha * c_potl_seq, na.rm=TRUE),
        population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    collect() -> global_stats_broken_down

global_stats_broken_down %>%
    group_by(biome, new_or_continued_1) %>%
    summarise(
        total_area_ha=sum(total_area_ha),
        restoration_area=sum(restoration_area),
        woody_carbon_stored=sum(woody_carbon_stored),
        soil_carbon_stored=sum(soil_carbon_stored),
        irr_carbon_stored=sum(irr_carbon_stored),
        irr_carbon_ha_any=sum(irr_carbon_ha_any),
        irr_carbon_ha_high=sum(irr_carbon_ha_high),
        carbon_seq=sum(carbon_seq),
        population=sum(population)) %>%
    write_csv('global_stats_broken_down.csv')


###############################################################################
###  Run stats by site

# Run stats by site
partied_data %>%
    group_by(site_id, id) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(site_id) %>%
    mutate(high_irr_c = c_tstor_ic > 25,
           any_irr_c = c_tstor_ic > .01,
           under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))) %>%
    summarise(
        total_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
        restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
        woody_carbon_stored=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
        soil_carbon_stored=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
        irr_carbon_stored=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
        irr_carbon_ha_any=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
        irr_carbon_ha_high=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
        carbon_seq=sum(area_ha * c_potl_seq, na.rm=TRUE),
        population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    collect() -> site_stats
site_stats %>%
    group_by(site_id) %>%
    summarise(
        total_area_ha=sum(total_area_ha),
        restoration_area=sum(restoration_area),
        woody_carbon_stored=sum(woody_carbon_stored),
        soil_carbon_stored=sum(soil_carbon_stored),
        irr_carbon_stored=sum(irr_carbon_stored),
        irr_carbon_ha_any=sum(irr_carbon_ha_any),
        irr_carbon_ha_high=sum(irr_carbon_ha_high),
        carbon_seq=sum(carbon_seq),
        population=sum(population)) %>%
    write_csv('site_stats.csv')


# Species
species <- data.table(readRDS(file.path(tables_folder, 'species.rds')))
speciesbysite <- data.table(readRDS(file.path(tables_folder, 'speciesbysite.rds')))

setkey(species, id)
setkey(speciesbysite, species_id)
species <- species[speciesbysite, nomatch=0]

# Add site columns
setkey(species, site_id)
setkey(sites, id)
species <- species[sites, nomatch=0, allow.cartesian=TRUE]

# Run in parallel across CPUs with multidplyr
species %>%
    distinct(id, .keep_all=TRUE) %>%
    summarise(n_species=n()) %>%
    as_tibble() %>%
    write_csv('species.csv')
