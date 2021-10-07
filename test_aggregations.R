library(tidyverse)
library(multidplyr)
library(data.table)

tables_folder <- 'D:/Data/Impacts_Data/tables'


# hist(pixels$c_tstor_soil[pixels$c_tstor_soil > 0] / 100 * pixels$area_ha[pixels$c_tstor_soil > 0])
# hist(pixels$c_tstor_woody[pixels$c_tstor_woody > 0] / 100 *  / pixels$area_ha[pixels$c_tstor_woody > 0])
#
# hist(pixels$c_tstor_soil[pixels$c_tstor_soil > 0] / 100 / pixels$area_ha[pixels$c_tstor_woody > 0])
# hist(pixels$c_tstor_woody[pixels$c_tstor_woody > 0] / 100)
#
#
#
# pixels <- readRDS(file.path(tables_folder, 'pixels_with_seq.rds'))
# pixelsbysite <- readRDS(file.path(file.path(tables_folder, 'pixelsbysite.rds')))
#
# nrow(pixelsbysite)
# nrow(pixels)
#
#
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

bna_name <- 'C4O - Blue Nature Alliance'
table(sites$division_name == bna_name)

sites <- sites[sites$division_name != bna_name, ]


pixels <- data.table(readRDS(file.path(tables_folder, 'pixels_with_seq.rds'))
pixelsbysite <- data.table(readRDS(file.path(tables_folder, 'pixelsbysite.rds')))

setkey(pixels, id)
setkey(pixelsbysite, pixel_id)
pixels <- pixels[pixelsbysite, nomatch=0]
saveRDS(pixels, 'tables/pixels_with_site_codes.rds')

setkey(pixels, site_id)
setkey(sites, id)
pixels <- pixels[sites, nomatch=0, allow.cartesian=TRUE]
saveRDS(pixels, 'tables/pixels_with_sites.rds')

#distinct(id, .keep_all=TRUE) %>%

setwd('/home/rstudio/data/impacts_data')
pixels <- readRDS(tables_folder, 'pixels_with_sites.rds')
pixels <- tbl_df(pixels)

# cluster <- new_cluster(14)
# cluster_library(cluster, "dplyr")
# pixels_par <- pixels %>% group_by(site_id, sls_name) %>% partition(cluster)
    
pixels %>%
    group_by(sls_name, id) %>%
    summarise(
        which_max = which.max(coverage_fraction),
        coverage_fraction=coverage_fraction[which_max],
        area_ha=area_ha[which_max],
        woody_carbon_stored=c_tstor_woody[which_max],
        soil_carbon_stored=c_tstor_soil[which_max],
        irr_carbon_stored=c_tstor_ic[which_max],
        population=population[which_max]
    ) %>%
    group_by(sls_name) %>%
    summarise(
        area_ha=sum(area_ha * coverage_fraction),
        woody_carbon_stored=sum(woody_carbon_stored * coverage_fraction),
        soil_carbon_stored=sum(soil_carbon_stored * coverage_fraction),
        irr_carbon_stored=sum(irr_carbon_stored * coverage_fraction),
        population=sum(population * coverage_fraction)
    ) %>%
    write_csv('stats_by_sls.csv')

pixels %>%
    group_by(site_id, id) %>%
    summarise(
        which_max = which.max(coverage_fraction),
        coverage_fraction=coverage_fraction[which_max],
        area_ha=area_ha[which_max],
        woody_carbon_stored=c_tstor_woody[which_max],
        soil_carbon_stored=c_tstor_soil[which_max],
        irr_carbon_stored=c_tstor_ic[which_max],
        population=population[which_max]
    ) %>%
    group_by(site_id) %>%
    summarise(
        area_ha=sum(area_ha * coverage_fraction),
        woody_carbon_stored=sum(woody_carbon_stored * coverage_fraction),
        soil_carbon_stored=sum(soil_carbon_stored * coverage_fraction),
        irr_carbon_stored=sum(irr_carbon_stored * coverage_fraction),
        population=sum(population * coverage_fraction)
    ) %>%
    write_csv('stats_by_site.csv')


sites_sp_save <- readRDS(tables_folder, 'sites_sp.rds')
sls <- readRDS(file.path(tables_folder, 'sls.rds'))
slsbysite <- readRDS(file.path(tables_folder, 'slsbysite.rds'))
countries <- readRDS(file.path(tables_folder, 'countries.rds'))
interventions <- readRDS(file.path(tables_folder, 'interventions.rds'))
interventionbysite <- readRDS(file.path(tables_folder, 'interventionbysite.rds'))
startags <- readRDS(file.path(tables_folder, 'startags.rds'))
startagbysite <- readRDS(file.path(tables_folder, 'startagbysite.rds'))
divisions <- readRDS(file.path(tables_folder, 'divisions.rds'))

species <- readRDS(file.path(tables_folder, 'species.rds'))
speciesbysite <- readRDS(file.path(tables_folder, 'speciesbysite.rds'))
sequestration_potential <- readRDS(file.path(data_folder_local, 'tables/c_potl_seq.rds'))
