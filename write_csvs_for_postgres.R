library(tidyverse)
library(multidplyr)
library(data.table)
library(dtplyr)
library(aws.s3)
library(aws.ec2metadata) # Needed to use IAM roles on EC2

tables_folder <- 'D:/Data/Impacts_Data/tables'
csv_folder <- 'D:/Data/Impacts_Data/csv_out'
data_folder_avoided_emissions <- 'D:/Data/Impacts_Data/avoided_emissions_data'
#tables_folder <- '/home/rstudio/data/impacts_data/tables'
#csv_folder <- '/home/rstudio/data/impacts_data/csv_out'

###############################################################################
# Join various tables together

sls <- data.table(readRDS(file.path(tables_folder, 'sls.rds')))
countries <- data.table(readRDS(file.path(tables_folder, 'countries.rds')))
divisions <- data.table(readRDS(file.path(tables_folder, 'divisions.rds')))
ecosystems <- data.table(readRDS(file.path(tables_folder, 'ecosystems.rds')))
setkey(sls, id)
setkey(countries, id)
setkey(divisions, id)
setkey(ecosystems, id)

slsbysite <- data.table(readRDS(file.path(tables_folder, 'slsbysite.rds')))
countrybysite <- data.table(readRDS(file.path(tables_folder, 'countrybysite.rds')))
divisionbysite <- data.table(readRDS(file.path(tables_folder, 'divisionbysite.rds')))
setkey(slsbysite, sls_id)
setkey(countrybysite, country_id)
setkey(divisionbysite, division_id)

sls <- sls[slsbysite, ]
divisions<- divisions[divisionbysite, ]
countries <- countries[countrybysite, ]

sls[, name:=NULL]
divisions[, name:=NULL]
countries[, name:=NULL]
setnames(sls, "id", "sls_id")
setnames(countries, "id", "country_id")
setnames(divisions, "id", "division_id")

# sls[, id:=NULL]
# divisions[, id:=NULL]
# countries[, id:=NULL]
# setnames(sls, "name", "sls_name")
# setnames(countries, "name", "country_name")
# setnames(divisions, "name", "division_name")
# countries <- countries[, .(site_id, country_name)]

countries <- countries[, .(site_id, country_id)]

sites <- data.table(readRDS(file.path(tables_folder, 'sites.rds')))
# Remove shape column since it is large and not needed
sites[, shape:=NULL]
sites[, reporting_year:=NULL]

setkey(sites, id)
setkey(sls, site_id)
setkey(divisions, site_id)
setkey(countries, site_id)

sites <- sites[sls, nomatch=0]
sites <- sites[countries, nomatch=0]
sites <- sites[divisions, nomatch=0]
setkey(sites, id)

###############################################################################
###  Avoided emissions and forest change trends

pixels_ae_full_investment_period <- data.table(readRDS(file.path(tables_folder, 'pixels_ae_full_investment_period.rds')))
pixels_ae_by_year <- data.table(readRDS(file.path(tables_folder, 'pixels_ae_by_year.rds')))
pixels_ae_bysite <- data.table(readRDS(file.path(tables_folder, 'pixels_ae_bysite.rds')))
pixels_ae_sampledfraction <- data.table(readRDS(file.path(tables_folder, 'pixels_ae_sampledfraction.rds')))

pixels_ae_sampledfraction %>%
    distinct(site_id, sampled_fraction) %>%
    filter(!(grepl('BNA', site_id))) %>%
    as_tibble() -> pixels_ae_sampledfraction_filtered

inner_join(
    pixels_ae_full_investment_period, pixels_ae_bysite, by=c('id'='cell_id')
) %>%
    inner_join(pixels_ae_sampledfraction_filtered) %>%
    group_by(site_id) %>%
    summarise(
        forest_loss_avoided_ha=sum(forest_loss_avoided_ha),
        emissions_avoided_mgco2e = sum(emissions_avoided_mgco2e),
        area_sampled_ha = sum(area_sampled_ha)
    ) %>%
    mutate(reporting_year = 'FY2022') %>%
    inner_join(pixels_ae_sampledfraction) %>%
    mutate(
        forest_loss_avoided_ha = forest_loss_avoided_ha * 1 / sampled_fraction[1],
        emissions_avoided_mgco2e = emissions_avoided_mgco2e * 1 / sampled_fraction,
        area_ha = area_sampled_ha * 1 / sampled_fraction
    ) %>%
    select(-sampled_fraction) %>%
    as_tibble() -> avoided_emissions_full_period
avoided_emissions_full_period %>%
    select(-area_ha) %>%
    write_csv(file.path(csv_folder, 'avoided_emissions_full_period.csv'))

inner_join(
    pixels_ae_by_year, pixels_ae_bysite, by=c('id'='cell_id')
) %>%
    inner_join(pixels_ae_sampledfraction_filtered) %>%
    group_by(site_id, year) %>%
    summarise(
        forest_loss_avoided_ha=sum(forest_loss_avoided_ha) * 1 / sampled_fraction[1],
        emissions_avoided_mgco2e = sum(emissions_avoided_mgco2e) * 1 / sampled_fraction[1],
        area_ha = sum(area_sampled_ha) * 1 / sampled_fraction[1]
    ) %>%
    mutate(reporting_year = 'FY2022') %>%
    as_tibble() -> avoided_emissions_by_year
avoided_emissions_by_year %>%
    select(-area_ha) %>%
    write_csv(file.path(csv_folder, 'avoided_emissions_by_year.csv'))

avoided_emissions_by_year %>%
    group_by(year) %>%
    summarise(forest_loss_avoided_ha = sum(forest_loss_avoided_ha),
              emissions_avoided_mgco2e = sum(emissions_avoided_mgco2e),
              area_ha = sum(area_ha)) %>%
    pivot_longer(forest_loss_avoided_ha:area_ha) %>%
    ggplot() +
    geom_line(aes(year, value)) +
    facet_wrap(name~., scales='free_y')


# Forest change data (annual change in forest cover)
treatments_and_controls <- data.table(readRDS(file.path(data_folder_avoided_emissions, 'treatments_and_controls.RDS')))
treatment_key <- data.table(readRDS(file.path(data_folder_avoided_emissions, 'treatment_cell_key.RDS')))

treatment_key[, c("reporting_year", "region"):=NULL]
treatments_and_controls <- treatments_and_controls[, .SD, .SDcols = patterns("(cell)|(fc_[0-9]{2})")]

setkey(treatment_key, cell)
setkey(treatments_and_controls, cell)

fc <- treatment_key[treatments_and_controls, nomatch=0]

fc# Convert pixel areas into pixels areas that are actually within site, using
# coverage_fraction. Then drop that col as no longer needed
fc[, area_ha:=area_ha*coverage_fraction]
fc[, coverage_fraction:=NULL]

# Multiply each fc col by the area in hectares and coverage fraction, to
# convert the fc cols from fractions to forest areas. Neex to convert fc cols
# from percentages (/ 100)
indx <- grep('fc_', colnames(fc))
for(j in indx){
    set(fc, i=NULL, j=j, value=(fc[[j]]/ 100 ) * fc[['area_ha']])
}

# Sum across all rows to get forest cover in hectares for each year
fc[, cell:=NULL]
fc <- fc[, lapply(.SD, sum, na.rm=TRUE), by=id_numeric]

site_id_numeric_key <- data.table(read_csv(file.path(data_folder_avoided_emissions, 'site_id_key.csv')))
setkey(fc, id_numeric)
setkey(site_id_numeric_key, id_numeric)
fc <- site_id_numeric_key[fc, nomatch=0]
fc[, id_numeric:=NULL]

write_csv(fc, file.path(csv_folder, 'forest_change.csv'))
put_object(
  file.path(csv_folder, 'forest_change.csv'),
  object='impacts_data/csv_out/forest_change.csv',
  bucket='landdegradation'
)

###############################################################################
###  Main pixels data

pixels <- data.table(readRDS(file.path(tables_folder, 'pixels_with_seq.rds')))
pixelsbysite <- data.table(readRDS(file.path(tables_folder, 'pixelsbysite.rds')))

setkey(pixels, id)
setkey(pixelsbysite, pixel_id)
pixels <- pixels[pixelsbysite, nomatch=0]

# Add site columns
setkey(pixels, site_id)
setkey(sites, id)
pixels_join <- pixels[sites, nomatch=0, allow.cartesian=TRUE]

saveRDS(pixels_join, file.path(tables_folder, 'pixels_joined.rds'))


pixels_join <- readRDS(file.path(tables_folder, 'pixels_joined.rds'))

cluster  <- new_cluster(24)
# Partition by site
pixels_join %>%
	as_tibble() %>%
	group_by(id, site_id) %>%
	partition(cluster) -> partied_data


###############################################################################
### Country summaries

partied_data %>%
	group_by(id, country_id, reporting_year, new_or_continued_1, ecosystem) %>%
	slice(which.max(coverage_fraction)) %>%
	group_by(country_id, reporting_year, new_or_continued_1, ecosystem) %>%
	mutate(
		   high_irr_c = c_tstor_ic > 25,
		   any_irr_c = c_tstor_ic > .01,
		   under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
		   ) %>%
	summarise(
			  calc_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
			  restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
			  calc_tstor_woody=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_soil=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
			  calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
			  calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
			  calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
			  calc_population=sum(population * coverage_fraction, na.rm=TRUE)
			  ) %>%
  mutate(
    calc_tstor_blue=sum(calc_tstor_ic[ecosystem %in% c(5, 6, 7)])  # Calculate as sum of ic in mangroves, salt marsh, and seagrasses
  ) %>%
	collect() -> country_data
country_data %>%
	group_by(country_id, reporting_year, new_or_continued_1) %>%
	summarise(
			  ttl_area_ha=sum(calc_area_ha),
			  ttl_restoration_area=sum(restoration_area),
			  ttl_tstor_woody=sum(calc_tstor_woody),
			  ttl_tstor_soil=sum(calc_tstor_soil),
			  ttl_tstor_blue=sum(calc_tstor_blue),
			  ttl_tstor_ic=sum(calc_tstor_ic),
			  ttl_ha_ic=sum(calc_ha_ic),
			  ttl_ha_high_ic=sum(calc_ha_high_ic),
			  ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
			  ttl_carbon_sequestered=sum(calc_carbon_sequestered),
			  ttl_population=sum(calc_population)
			  ) %>%
	mutate(
	     ttl_tstor=sum(ttl_tstor_woody, ttl_tstor_soil)
		   ) %>%
	write_csv(file.path(csv_folder, 'country_data.csv'))
put_object(
		   file.path(csv_folder, 'country_data.csv'),
		   object='impacts_data/csv_out/country_data.csv',
		   bucket='landdegradation'
)

###############################################################################
### Scape summaries
partied_data %>%
	group_by(id, sls_id, reporting_year, new_or_continued_1, ecosystem) %>%
	slice(which.max(coverage_fraction)) %>%
	group_by(sls_id, reporting_year, new_or_continued_1, ecosystem) %>%
	mutate(
		   high_irr_c = c_tstor_ic > 25,
		   any_irr_c = c_tstor_ic > .01,
		   under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
		   ) %>%
	summarise(
			  calc_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
			  restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
			  calc_tstor_woody=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_soil=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
			  calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
			  calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
			  calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
			  calc_population=sum(population * coverage_fraction, na.rm=TRUE)
			  ) %>%
  mutate(
    calc_tstor_blue=sum(calc_tstor_ic[ecosystem %in% c(5, 6, 7)])  # Calculate as sum of ic in mangroves, salt marsh, and seagrasses
  ) %>%
	collect() -> sls_data
sls_data %>%
	group_by(sls_id, reporting_year, new_or_continued_1) %>%
	summarise(
			  ttl_area_ha=sum(calc_area_ha),
			  ttl_restoration_area=sum(restoration_area),
			  ttl_tstor_woody=sum(calc_tstor_woody),
			  ttl_tstor_soil=sum(calc_tstor_soil),
			  ttl_tstor_blue=sum(calc_tstor_blue),
			  ttl_tstor_ic=sum(calc_tstor_ic),
			  ttl_ha_ic=sum(calc_ha_ic),
			  ttl_ha_high_ic=sum(calc_ha_high_ic),
			  ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
			  ttl_carbon_sequestered=sum(calc_carbon_sequestered),
			  ttl_population=sum(calc_population)
			  ) %>%
	mutate(
	     ttl_tstor=sum(ttl_tstor_woody, ttl_tstor_soil)
		   ) %>%
	write_csv(file.path(csv_folder, 'sls_data.csv'))
put_object(
		   file.path(csv_folder, 'sls_data.csv'),
		   object='impacts_data/csv_out/sls_data.csv',
		   bucket='landdegradation'
)

###############################################################################
### Division summaries
partied_data %>%
	group_by(id, division_id, reporting_year, new_or_continued_1, ecosystem) %>%
	slice(which.max(coverage_fraction)) %>%
	group_by(division_id, reporting_year, new_or_continued_1, ecosystem) %>%
	mutate(
		   high_irr_c = c_tstor_ic > 25,
		   any_irr_c = c_tstor_ic > .01,
		   under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
		   ) %>%
	summarise(
			  calc_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
			  restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
			  calc_tstor_woody=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_soil=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
			  calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
			  calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
			  calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
			  calc_population=sum(population * coverage_fraction, na.rm=TRUE)
			  ) %>%
  mutate(
    calc_tstor_blue=sum(calc_tstor_ic[ecosystem %in% c(5, 6, 7)])  # Calculate as sum of ic in mangroves, salt marsh, and seagrasses
  ) %>%
	collect() -> division_data
division_data %>%
	group_by(division_id, reporting_year, new_or_continued_1) %>%
	summarise(
			  ttl_area_ha=sum(calc_area_ha),
			  ttl_restoration_area=sum(restoration_area),
			  ttl_tstor_woody=sum(calc_tstor_woody),
			  ttl_tstor_soil=sum(calc_tstor_soil),
			  ttl_tstor_blue=sum(calc_tstor_blue),
			  ttl_tstor_ic=sum(calc_tstor_ic),
			  ttl_ha_ic=sum(calc_ha_ic),
			  ttl_ha_high_ic=sum(calc_ha_high_ic),
			  ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
			  ttl_carbon_sequestered=sum(calc_carbon_sequestered),
			  ttl_population=sum(calc_population)
			  ) %>%
	mutate(
	     ttl_tstor=sum(ttl_tstor_woody, ttl_tstor_soil)
		   ) %>%
	write_csv(file.path(csv_folder, 'division_data.csv'))
put_object(
		   file.path(csv_folder, 'division_data.csv'),
		   object='impacts_data/csv_out/division_data.csv',
		   bucket='landdegradation'
)

###############################################################################
### Sites summaries
partied_data %>%
	group_by(id, site_id, reporting_year, new_or_continued_1, ecosystem) %>%
	slice(which.max(coverage_fraction)) %>%
	group_by(site_id, reporting_year, new_or_continued_1, ecosystem) %>%
	mutate(
		   high_irr_c = c_tstor_ic > 25,
		   any_irr_c = c_tstor_ic > .01,
		   under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
		   ) %>%
	summarise(
			  calc_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
			  restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
			  calc_tstor_woody=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_soil=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
			  calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
			  calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
			  calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
			  calc_population=sum(population * coverage_fraction, na.rm=TRUE)
			  ) %>%
  mutate(
    calc_tstor_blue=sum(calc_tstor_ic[ecosystem %in% c(5, 6, 7)])  # Calculate as sum of ic in mangroves, salt marsh, and seagrasses
  ) %>%
	collect() -> sites_data
sites_data %>%
	group_by(country_id, reporting_year, new_or_continued_1) %>%
	summarise(
			  ttl_area_ha=sum(calc_area_ha),
			  ttl_restoration_area=sum(restoration_area),
			  ttl_tstor_woody=sum(calc_tstor_woody),
			  ttl_tstor_soil=sum(calc_tstor_soil),
			  ttl_tstor_blue=sum(calc_tstor_blue),
			  ttl_tstor_ic=sum(calc_tstor_ic),
			  ttl_ha_ic=sum(calc_ha_ic),
			  ttl_ha_high_ic=sum(calc_ha_high_ic),
			  ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
			  ttl_carbon_sequestered=sum(calc_carbon_sequestered),
			  ttl_population=sum(calc_population)
			  ) %>%
	mutate(
		   ttl_tstor=sum(ttl_tstor_woody, ttl_tstor_soil)
		   ) %>%
	write_csv(file.path(csv_folder, 'site_data.csv'))
put_object(
		   file.path(csv_folder, 'site_data.csv'),
		   object='impacts_data/csv_out/site_data.csv',
		   bucket='landdegradation'
)

###############################################################################
### Ecosystems summaries
partied_data %>%
	group_by(id, country_id, ecosystem, reporting_year, new_or_continued_1) %>%
	slice(which.max(coverage_fraction)) %>%
	group_by(country_id, ecosystem, reporting_year, new_or_continued_1) %>%
	mutate(
		   high_irr_c = c_tstor_ic > 25,
		   any_irr_c = c_tstor_ic > .01,
		   under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
		   ) %>%
	summarise(
			  calc_area_ha=sum(area_ha * coverage_fraction, na.rm=TRUE),
			  restoration_area=sum(area_ha[under_restoration] * coverage_fraction[under_restoration], na.rm=TRUE),
			  calc_tstor_woody=sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_soil=sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
			  calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
			  calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
			  calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
			  calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
			  calc_population=sum(population * coverage_fraction, na.rm=TRUE)
			  ) %>%
	collect() -> ecosystem_data
ecosystem_data %>%
	group_by(country_id, ecosystem, reporting_year, new_or_continued_1) %>%
	summarise(
			  ttl_area_ha=sum(calc_area_ha),
			  ttl_restoration_area=sum(restoration_area),
			  ttl_tstor_woody=sum(calc_tstor_woody),
			  ttl_tstor_soil=sum(calc_tstor_soil),
			  ttl_tstor_ic=sum(calc_tstor_ic),
			  ttl_ha_ic=sum(calc_ha_ic),
			  ttl_ha_high_ic=sum(calc_ha_high_ic),
			  ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
			  ttl_carbon_sequestered=sum(calc_carbon_sequestered),
			  ttl_population=sum(calc_population)
			  ) %>%
	mutate(
	     ttl_tstor=sum(ttl_tstor_woody, ttl_tstor_soil)
		   ) %>%
	write_csv(file.path(csv_folder, 'ecosystem_data.csv'))
put_object(
		   file.path(csv_folder, 'ecosystem_data.csv'),
		   object='impacts_data/csv_out/ecosystem_data.csv',
		   bucket='landdegradation'
)


#pixels_sliced <- pixels_join[
#  pixels_join[,
#              .I[pt == max(coverage_fraction)],
#              by=list(id, country_id, reporting_year, new_or_continued_1)
#  ]$V1
#]

#pixels_sliced <- pixels_join[, .SD[which.max(coverage_fraction)], by=list(id, country_id, reporting_year, new_or_continued_1)]

#pixels_sliced[, high_irr_c := c_tstor_ic > 25]
#pixels_sliced[, any_irr_c := c_tstor_ic > .01]
#pixels_sliced[, under_restoration := (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type)]
