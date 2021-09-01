library(tidyverse)
library(sf)
library(units)
library(lubridate)

sites <- readRDS('tables/sites_sp.RDS')
# Add a field to use as a flag for sites under/over 100 hectares
sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'
sites %>% left_join(
    sites_cea %>%
        as_tibble() %>%
        select(id, area_cea)
    ) -> sites

table(is.na(sites$ci_start_date))
table(is.na(sites$ci_end_date))

sites %>%
    ggplot() +
    geom_histogram(aes(ci_start_date))

# Set all start dates that are missing to 2018 (around the most common year)
sites$ci_start_date[is.na(sites$ci_start_date)] <- mdy('1/1/2018')
sites$ci_start_year <- year(sites$ci_start_date)

table(is.na(sites$ci_end_date))
sites %>%
    ggplot() +
    geom_histogram(aes(ci_end_date))
# Set all end dates that are greater than 12/31/2020 to NA, so they are treated 
# as ongoing
sites$ci_end_date[sites$ci_end_date > mdy('12/31/2020')] <- NA
sites$ci_end_year <- year(sites$ci_end_date)
table(is.na(sites$ci_end_year))

sites$id_numeric <- 1:nrow(sites)

sites %>%
    as_tibble() %>%
    select(id, id_numeric, -shape) %>%
    write_csv('site_id_key.csv')
saveRDS(sites, 'avoided_emissions_data/sites_cleaned_for_avoided_emissions.RDS')
