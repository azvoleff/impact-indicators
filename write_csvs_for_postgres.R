library(tidyverse)
library(multidplyr)
library(data.table)
library(dtplyr)
library(aws.s3)
library(aws.ec2metadata) # Needed to use IAM roles on EC2

#tables_folder <- 'D:/Data/Impacts_Data/tables'
tables_folder <- '/home/rstudio/data/impacts_data/tables'
csv_folder <- '/home/rstudio/data/impacts_data/csv_out'

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



###############################################################################
###  Partition across CPUs

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


###############################################################################
### Avoided emissions summaries


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
