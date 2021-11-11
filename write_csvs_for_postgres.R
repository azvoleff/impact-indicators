library(tidyverse)
library(multidplyr)
library(data.table)
library(dtplyr)
library(aws.s3)
library(aws.ec2metadata) # Needed to use IAM roles on EC2

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

saveRDS(pixels_join, 'pixels_join.rds')


###############################################################################
###  Partition across CPUs

cluster  <- new_cluster(24)
# Partition by site
pixels_join %>%
  as_tibble() %>%
  group_by(id, site_id) %>%
  partition(cluster) -> partied_data
  

###############################################################################
### Country summaries

partied_data %>%
    group_by(id, country_id, reporting_year, new_or_continued_1) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(country_id, reporting_year, new_or_continued_1) %>%
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
        carbon_sequestered=sum(area_ha * c_potl_seq, na.rm=TRUE),
		calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
        calc_population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    mutate(
		calc_tstor=calc_tstor_woody + calc_tstor_soil
    ) %>%
    collect() -> country_data
country_data %>%
    group_by(country_id, reporting_year, new_or_continued_1) %>%
    summarise(
        ttl_area_ha=sum(calc_area_ha),
        ttl_restoration_area=sum(restoration_area),
        ttl_tstor_woody=sum(calc_tstor_woody),
        ttl_tstor_soil=sum(calc_tstor_soil),
		ttl_tstor=sum(calc_tstor),
        ttl_tstor_ic=sum(calc_tstor_ic),
        ttl_ha_ic=sum(calc_ha_ic),
        ttl_ha_high_ic=sum(calc_ha_high_ic),
		ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
        ttl_carbon_seq=sum(carbon_sequestered),
		ttl_carbon_sequestered=sum(calc_carbon_sequestered),
        ttl_population=sum(calc_population)
    ) %>% write_csv('country_data.csv')
put_object('country_data.csv', object='impacts_data/csv_out/country_data.csv', bucket='landdegradation')

###############################################################################
### Scape summaries
partied_data %>%
    group_by(id, sls_id, reporting_year, new_or_continued_1) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(sls_id, reporting_year, new_or_continued_1) %>%
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
		calc_tstor=sum(sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE), sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE), na.rm=TRUE),
        calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
        calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
        calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
		calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
        carbon_sequestered=sum(area_ha * c_potl_seq, na.rm=TRUE),
		calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
        calc_population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    collect() -> sls_data
sls_data %>%
    group_by(sls_id, reporting_year, new_or_continued_1) %>%
    summarise(
        ttl_area_ha=sum(calc_area_ha),
        ttl_restoration_area=sum(restoration_area),
        ttl_tstor_woody=sum(calc_tstor_woody),
        ttl_tstor_soil=sum(calc_tstor_soil),
		ttl_tstor=sum(calc_tstor),
        ttl_tstor_ic=sum(calc_tstor_ic),
        ttl_ha_ic=sum(calc_ha_ic),
        ttl_ha_high_ic=sum(calc_ha_high_ic),
		ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
        ttl_carbon_seq=sum(carbon_sequestered),
		ttl_carbon_sequestered=sum(calc_carbon_sequestered),
        ttl_population=sum(calc_population)
    ) %>% write_csv('sls_data.csv')
put_object('sls_data.csv', object='impacts_data/csv_out/sls_data.csv', bucket='landdegradation')

###############################################################################
### Division summaries
partied_data %>%
    group_by(id, division_id, reporting_year, new_or_continued_1) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(division_id, reporting_year, new_or_continued_1) %>%
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
		calc_tstor=sum(sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE), sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE), na.rm=TRUE),
        calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
        calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
        calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
		calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
        carbon_sequestered=sum(area_ha * c_potl_seq, na.rm=TRUE),
		calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
        calc_population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    collect() -> division_data
division_data %>%
    group_by(division_id, reporting_year, new_or_continued_1) %>%
    summarise(
        ttl_area_ha=sum(calc_area_ha),
        ttl_restoration_area=sum(restoration_area),
        ttl_tstor_woody=sum(calc_tstor_woody),
        ttl_tstor_soil=sum(calc_tstor_soil),
		ttl_tstor=sum(calc_tstor),
        ttl_tstor_ic=sum(calc_tstor_ic),
        ttl_ha_ic=sum(calc_ha_ic),
        ttl_ha_high_ic=sum(calc_ha_high_ic),
		ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
        ttl_carbon_seq=sum(carbon_sequestered),
		ttl_carbon_sequestered=sum(calc_carbon_sequestered),
        ttl_population=sum(calc_population)
    ) %>% write_csv('division_data.csv')
put_object('division_data.csv', object='impacts_data/csv_out/division_data.csv', bucket='landdegradation')

###############################################################################
### Sites summaries
partied_data %>%
    group_by(id, site_id, reporting_year, new_or_continued_1) %>%
    slice(which.max(coverage_fraction)) %>%
    group_by(site_id, reporting_year, new_or_continued_1) %>%
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
		calc_tstor=sum(sum(c_tstor_woody / 100 * area_ha * coverage_fraction, na.rm=TRUE), sum(c_tstor_soil / 100 * area_ha * coverage_fraction, na.rm=TRUE), na.rm=TRUE),
        calc_tstor_ic=sum(c_tstor_ic * area_ha * coverage_fraction, na.rm=TRUE),
        calc_ha_ic=sum(area_ha[any_irr_c] * coverage_fraction[any_irr_c], na.rm=TRUE),
        calc_ha_high_ic=sum(area_ha[high_irr_c] * coverage_fraction[high_irr_c], na.rm=TRUE),
		calc_tonneshectare_ic=sum((c_tstor_ic * area_ha * coverage_fraction) / (area_ha * coverage_fraction), na.rm=TRUE),
        carbon_sequestered=sum(area_ha * c_potl_seq, na.rm=TRUE),
		calc_carbon_sequestered=sum(area_ha * c_potl_seq * coverage_fraction, na.rm=TRUE),
        calc_population=sum(population * coverage_fraction, na.rm=TRUE)
    ) %>%
    collect() -> site_data
site_data %>%
    group_by(site_id, reporting_year, new_or_continued_1) %>%
    summarise(
        ttl_area_ha=sum(calc_area_ha),
        ttl_restoration_area=sum(restoration_area),
        ttl_tstor_woody=sum(calc_tstor_woody),
        ttl_tstor_soil=sum(calc_tstor_soil),
		ttl_tstor=sum(calc_tstor),
        ttl_tstor_ic=sum(calc_tstor_ic),
        ttl_ha_ic=sum(calc_ha_ic),
        ttl_ha_high_ic=sum(calc_ha_high_ic),
		ttl_tonneshectare_ic=sum(calc_tonneshectare_ic),
        ttl_carbon_seq=sum(carbon_sequestered),
		ttl_carbon_sequestered=sum(calc_carbon_sequestered),
        ttl_population=sum(calc_population)
    ) %>% write_csv('site_data.csv')
put_object('site_data.csv', object='impacts_data/csv_out/site_data.csv', bucket='landdegradation')





















pixels_sliced <- pixels_join[
  pixels_join[,
              .I[pt == max(coverage_fraction)],
              by=list(id, country_id, reporting_year, new_or_continued_1)
  ]$V1
]

pixels_sliced <- pixels_join[, .SD[which.max(coverage_fraction)], by=list(id, country_id, reporting_year, new_or_continued_1)]

pixels_sliced[, high_irr_c := c_tstor_ic > 25]
pixels_sliced[, any_irr_c := c_tstor_ic > .01]
pixels_sliced[, under_restoration := (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type)]




group_by(country_id, reporting_year, new_or_continued_1) %>%
  mutate(
    high_irr_c = c_tstor_ic > 25,
    any_irr_c = c_tstor_ic > .01,
    under_restoration = (restoration_type != ' ' & restoration_type != 'Not Applicable' & !is.na(restoration_type))
  ) %>%
  
  
  group_slice_max = group[group[, .I[which.max(pt)], by=Subject]$V1]




ecosystem_data %>%
    group_by(country_id, reporting_year, new_or_continued_1, ecosystem) %>%
    summarise(
        ttl_ha_ic=sum(calc_ha_ic)
    ) %>%
write_csv('country_ecosystem_carbon.csv')



########################################
# We also need the output of a pair of crosstab queries to generate additional data. If it's simpler, we could store the pixel results in one table and the crosstab outputs in two more tables. I can join them just as easily in a SQL query inside the database.
# This might be one way to script the crosstab query in R, see what you think.

# SQL

(
	SELECT
		ct.country_id,
		ct.new_or_continued_1,
		ct.reporting_year,
		ct."1" AS ttl_tstor_pfrst,
		ct."2" AS ttl_tstor_ofrst,
		ct."3" AS ttl_tstor_grass,
		ct."4" AS ttl_tstor_wet,
		ct."5" AS ttl_tstor_mang,
		ct."6" AS ttl_tstor_sltmrsh,
		ct."7" AS ttl_tstor_sgrass,
		ct."8" AS ttl_tstor_peat
	FROM crosstab(
		'SELECT country_id, new_or_continued_1, reporting_year, ecosystem AS ecosystemid, SUM(calc_tstor) AS ttl_tstor
		FROM m_countrypixels
		WHERE (ecosystem is not null AND ecosystem > 0)
		AND new_or_continued_1 = ''Continued this year''
		GROUP BY country_id, new_or_continued_1, reporting_year, ecosystem
		ORDER BY 1, 2'::text, 'VALUES (''1''::double precision),(''2''::double precision),(''3''::double precision),(''4''::double precision),(''5''::double precision),(''6''::double precision),(''7''::double precision),(''8''::double precision)'::text)
	AS ct(country_id text, new_or_continued_1 text, reporting_year text, "1" double precision, "2" double precision, "3" double precision, "4" double precision, "5" double precision, "6" double precision, "7" double precision, "8" double precision)
	UNION ALL
	SELECT
		ct.country_id,
		ct.new_or_continued_1,
		ct.reporting_year,
		ct."1" AS ttl_tstor_pfrst,
		ct."2" AS ttl_tstor_ofrst,
		ct."3" AS ttl_tstor_grass,
		ct."4" AS ttl_tstor_wet,
		ct."5" AS ttl_tstor_mang,
		ct."6" AS ttl_tstor_sltmrsh,
		ct."7" AS ttl_tstor_sgrass,
		ct."8" AS ttl_tstor_peat
	FROM crosstab(
		'SELECT country_id, new_or_continued_1, reporting_year, ecosystem AS ecosystemid, SUM(calc_tstor) AS ttl_tstor
		FROM m_countrypixels
		WHERE (ecosystem is not null AND ecosystem > 0)
		AND new_or_continued_1 = ''New this year''
		GROUP BY country_id, new_or_continued_1, reporting_year, ecosystem
		ORDER BY 1, 2'::text, 'VALUES (''1''::double precision),(''2''::double precision),(''3''::double precision),(''4''::double precision),(''5''::double precision),(''6''::double precision),(''7''::double precision),(''8''::double precision)'::text)
	AS ct(country_id text, new_or_continued_1 text, reporting_year text, "1" double precision, "2" double precision, "3" double precision, "4" double precision, "5" double precision, "6" double precision, "7" double precision, "8" double precision)
) AS csquery ON csquery.country_id = cp.country_id AND csquery.new_or_continued_1 = cp.new_or_continued_1 AND csquery.reporting_year = cp.reporting_year

# R script (assumes a source table named CountryData, which could be created by loading the country_data.csv created in the first R script segment above, or otherwise using it directly if it's still in memory)

partied_data %>%
    aggregate(CountryData.calc_tstor, list(CountryData.country_id, CountryData.ecosystem), sum) %>%
    collect() -> country_cstor_eco
country_cstor_eco %>%
    write_csv('country_carbonstore_ecosystem.csv')
