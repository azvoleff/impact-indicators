library(tidyverse)
library(data.table)

sls <- data.table(readRDS('tables/sls.rds'))
countries <- data.table(readRDS('tables/countries.rds'))
divisions <- data.table(readRDS('tables/divisions.rds'))
setkey(sls, id)
setkey(countries, id)
setkey(divisions, id)

slsbysite <- data.table(readRDS('tables/slsbysite.rds'))
countrybysite <- data.table(readRDS('tables/countrybysite.rds'))
divisionbysite <- data.table(readRDS('tables/divisionbysite.rds'))
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

sites <- data.table(readRDS('tables/sites.rds'))
# Remove shape column since it is large and not needed
sites[, shape:=NULL]

setkey(sites, id)
setkey(sls, site_id)
setkey(divisions, site_id)
setkey(countries, site_id)

sites <- sites[sls, nomatch=0]
sites <- sites[countries, nomatch=0]
sites <- sites[divisions, nomatch=0]

pixels <- data.table(readRDS('tables/pixels.rds'))
pixelsbysite <- data.table(readRDS('tables/pixelsbysite.rds'))

setkey(pixels, id)
setkey(pixelsbysite, pixel_id)
pixels <- pixels[pixelsbysite, nomatch=0]
saveRDS(pixels, 'tables/pixels_with_site_codes.rds')

setkey(pixels, site_id)
setkey(sites, id)
pixels <- pixels[sites, nomatch=0, allow.cartesian=TRUE]
saveRDS(pixels, 'tables/pixels_with_sites.rds')

sitespixels %>%
    group_by(sls_id) %>%
    summarise(area_ha=sum(area_ha * coverage_fraction))

sitespixels %>%
    group_by(id) %>%
    summarise(area_ha=sum(area_ha * coverage_fraction))

sitespixels %>%
    group_by(id) %>%
    summarise(area_ha=sum(area_ha * coverage_fraction))

names(sitespixels)


sites_sp_save <- readRDS('tables/sites_sp.rds')
sls <- readRDS('tables/sls.rds')
slsbysite <- readRDS('tables/slsbysite.rds')
countries <- readRDS('tables/countries.rds')
interventions <- readRDS('tables/interventions.rds')
interventionbysite <- readRDS('tables/interventionbysite.rds')
startags <- readRDS('tables/startags.rds')
startagbysite <- readRDS('tables/startagbysite.rds')
divisions <- readRDS('tables/divisions.rds')

species <- readRDS('tables/species.rds')
speciesbysite <- readRDS('tables/speciesbysite.rds')
sequestration_potential <- readRDS(file.path(data_folder_local, 'tables/c_potl_seq.rds'))
