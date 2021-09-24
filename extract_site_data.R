library(raster)
library(gdalUtils)
library(data.table)
library(tidyverse)
library(lubridate)
library(tidyjson)
library(ISOcodes)
library(sf)
library(exactextractr)
library(foreach)
library(Rcpp)


data_folder_onedrive <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'
#data_folder_local <- 'E:/Data/Impacts_Data'
data_folder_local <- 'D:/Data/Impacts_Data'
data_folder_deg_paper <- 'E:/Data'

data_folder_local <- '/home/rstudio/data/impacts_data'


unjoin_table <- function(starts_with_text, id_name) {
    variablebysite <- sites %>%
        as_tibble() %>%
        dplyr::select(id, starts_with(starts_with_text)) %>%
        pivot_longer(starts_with(starts_with_text)) %>%
        dplyr::select(-name) %>%
        filter(!is.na(value)) %>%
        filter(!(value == " ")) %>%
        rename(name := value) %>%
        distinct()
    id_table <- variablebysite %>%
        summarise(name = unique(name)) %>%
        mutate(id = 1:n()) %>%
        relocate(id)
    variablebysite <- variablebysite %>%
        rename(site_id=id) %>%
        left_join(rename(id_table, !!quo_name(id_name) :=id)) %>%
        dplyr::select(-name)
    return(list(variablebysite, id_table))
}

# Need to split up polygons that cross the meridian
split_on_meridian <- function(p) {
    w_hemi <- st_sfc(st_polygon(list(matrix(c(-180, -180, 0, 0, -180, 90, -90, -90, 90, 90), ncol=2))), 
                       crs=4326)
    e_hemi <- st_sfc(st_polygon(list(matrix(c(0, 0, 180, 180, 0, 90, -90, -90, 90, 90), ncol=2))), 
                       crs=4326)
    e_side <- st_intersection(p, e_hemi)
    w_side <- st_intersection(p, w_hemi)
    rbind(e_side, w_side)
}


load_as_vrt <- function(folder, pattern, band=FALSE, raster=TRUE) {
    vrt_file <- tempfile(fileext='.vrt')
    files <- list.files(folder, pattern=pattern)
    if (length(files) == 0) {
        print('No files found')
        return(FALSE)
    }
    if (band) {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file, b=band)
        r <- raster(vrt_file)
    } else {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file)
        r <- stack(vrt_file)
    }
    if (raster) {
        return(r)
    } else {
        return(vrt_file)
    }
}
###############################################################################
### Load sites

sites_2020 <- st_read(paste0(data_folder_onedrive, "/Impact_Sites/CI_Online_Global_Impact_Sites.gdb"))
sites_2020 <- st_zm(sites_2020, drop=TRUE)
sites_2020 %>%
    rename(id = CI_ID) %>%
    rename_all(tolower) -> sites_2020
sites_2020 <- st_transform(sites_2020, 'EPSG:4326')

sites_2021 <- st_read(paste0(data_folder_onedrive, "/Impact_Sites/Final_PreVetting_FY2021.gdb"))
sites_2021 <- st_zm(sites_2021) # Drop zm coords
sites_2021 %>%
    rename(id = CI_ID) %>%
    rename_all(tolower) -> sites_2021
sites_2021 <- st_transform(sites_2021, 'EPSG:4326')

sites_2020 %>%
    mutate(validity = st_is_valid(., reason=TRUE)) %>%
    filter(validity != 'Valid Geometry') %>%
    as_tibble() %>%
    select(-shape) %>%
    write_csv('site_validity_2020.csv')

sites_2021 %>%
    mutate(validity = st_is_valid(., reason=TRUE)) %>%
    filter(validity != 'Valid Geometry') %>%
    as_tibble() %>%
    select(-shape) %>%
    write_csv('site_validity_2021.csv')


# Load 2021 data into sites table
sites_2021 %>%
    mutate(validity = st_is_valid(.)) %>%
    filter(validity == TRUE) %>%
    select(-validity) -> sites

# Fix error in SLS names
# table
sites$ci_sls_1 <- as.character(sites$ci_sls_1)
sites$ci_sls_1 <- gsub("Cardamom-Tonle Sap Landscape in Cambodia", "Cardamom-Tonle Sap Landscape", sites$ci_sls_1)
sites$ci_sls_2 <- as.character(sites$ci_sls_2)
sites$ci_sls_2 <- gsub("Cardamom-Tonle Sap Landscape in Cambodia", "Cardamom-Tonle Sap Landscape", sites$ci_sls_2)

# Save the spatial data into the sites_sp dataframe, and create a sites frame 
# that has the spatial data stored as WKT in a column
sites_sp <- sites
sites <- as_tibble(sites_sp)
sites$shape <- st_as_text(sites_sp$shape)

sites_sp <- split_on_meridian(sites_sp)

###################################################
# Split out data that was joined to the main tables
sls_unjoin <- unjoin_table('ci_sls', 'sls_id')
slsbysite <- sls_unjoin[[1]]
sls <- sls_unjoin[[2]]

division_unjoin <- unjoin_table('ci_division', 'division_id')
divisionbysite <- division_unjoin[[1]]
divisions <- division_unjoin[[2]]

startag_unjoin <- unjoin_table('star_tag', 'startag_id')
startagbysite <- startag_unjoin[[1]]
startags <- startag_unjoin[[2]]

intervention_unjoin <- unjoin_table('intervention_type', 'intervention_id')
interventionbysite <- intervention_unjoin[[1]]
interventions <- intervention_unjoin[[2]]

countrybysite <- dplyr::select(sites,
                        site_id=id,
                        country_id=country_iso_num_code) %>%
    mutate(country_id = country_id)
countries <- ISO_3166_1 %>%
    filter(Numeric %in% countrybysite$country_id) %>%
    dplyr::select(-Common_name) %>%
    rename(id = Numeric) %>%
    relocate(id) %>%
    rename_all(tolower)

sites <- dplyr::select(sites,
                -starts_with('ci_sls'),
                -starts_with('ci_division'),
                -starts_with('star_tag'),
                -starts_with('intervention_type'),
                -country_iso_num_code,
                -area_ha,
                -country)

sites_sp %>%
    select(id, shape) %>%
    full_join(
        select(sites, -shape)
    ) -> sites_sp_save

saveRDS(sites, 'tables/sites.rds')
saveRDS(sites_sp_save, 'tables/sites_sp.rds')
saveRDS(sls, 'tables/sls.rds')
saveRDS(slsbysite, 'tables/slsbysite.rds')
saveRDS(countries, 'tables/countries.rds')
saveRDS(countrybysite, 'tables/countrybysite.rds')
saveRDS(interventions, 'tables/interventions.rds')
saveRDS(interventionbysite, 'tables/interventionbysite.rds')
saveRDS(startags, 'tables/startags.rds')
saveRDS(startagbysite, 'tables/startagbysite.rds')
saveRDS(divisions, 'tables/divisions.rds')
saveRDS(divisionbysite, 'tables/divisionbysite.rds')

###############################################################################
### Load indicator data

#########
# Species
sp_raw <- read_csv('species_all_sites_CR-EN-VU.csv')
sp_raw %>%
    select(-site_id) %>%
    distinct() %>%
    mutate(id = 1:n()) %>%
    relocate(id) -> species
saveRDS(species, 'tables/species.rds')
left_join(sp_raw, species) %>%
    rename(species_id=id) %>%
    select(site_id, species_id) -> speciesbysite 
saveRDS(speciesbysite, 'tables/speciesbysite.rds')

#############
# Soil carbon
c_tstor_soil_rast <- load_as_vrt(data_folder_local, 'soil_carbon_stored_2020[-.0-9]*tif$')
names(c_tstor_soil_rast) <- 'c_tstor_soil'
crs(c_tstor_soil_rast) <- crs('EPSG:4326')

##############
# Woody carbon
#
# Note: can't use the load_as_vrt function as for this case there are just way 
# too many files for the woody carbon data, so load them manually into a vrt. 

# Need to reproject and crop to match the other layers
ext <- extent(c_tstor_soil_rast)
output_extent <- c(ext[1], ext[3], ext[2], ext[4])
output_dim <- c(ncol(c_tstor_soil_rast), nrow(c_tstor_soil_rast))
c_tstor_woody_vrt_file <- tempfile(fileext='.vrt')
gdalbuildvrt(file.path(data_folder_local, 
                       'woody_carbon_stored_2020*.tif'), 
             c_tstor_woody_vrt_file,
             te=output_extent, tr=c(xres(c_tstor_soil_rast), yres(c_tstor_soil_rast)))
c_tstor_woody_rast <- brick(c_tstor_woody_vrt_file)
crs(c_tstor_woody_rast) <- crs('EPSG:4326')
names(c_tstor_woody_rast) <-c(paste0("fc", 2000:2020),
                              paste0("fl", 2001:2020),
                              paste0("cb", 2000:2020),
                              paste0("ce", 2001:2020))

######################
# Irrecoverable carbon
c_tstor_ic_vrt_file <- tempfile(fileext='.vrt')
gdalbuildvrt(file.path(data_folder_local, "irrecoverable_carbon_250m_2018*.tif"),
             c_tstor_ic_vrt_file,
             te=output_extent,
             tr=c(xres(c_tstor_soil_rast), yres(c_tstor_soil_rast)))
c_tstor_ic_rast <- raster(c_tstor_ic_vrt_file)
names(c_tstor_ic_rast) <- 'c_tstor_ic'

##########################################
# Potential sequestration from restoration
rest_rast_patterns <- c('c_potl_seq_agfor0020_250m*',
                        'c_potl_seq_agfor2060_250m*',
                        'c_potl_seq_natre0020_250m*',
                        'c_potl_seq_natre2060_250m*',
                        'c_potl_seq_mtrer0020_250m*',
                        'c_potl_seq_mtrer2060_250m*',
                        'c_potl_seq_mshrr0020_250m*',
                        'c_potl_seq_mshrr2060_250m*',
                        'c_potl_seq_pwtea0020_250m*',
                        'c_potl_seq_pwpin0020_250m*',
                        'c_potl_seq_pwoco0020_250m*',
                        'c_potl_seq_pwobr0020_250m*',
                        'c_potl_seq_pwoak0020_250m*',
                        'c_potl_seq_pweuc0020_250m*')
rest_vrts <- foreach(p=rest_rast_patterns, .combine=c) %do% {
    vrt_file <- tempfile(fileext='.vrt')
    gdalbuildvrt(file.path(data_folder_local, p), vrt_file,
                 te=output_extent, tr=c(xres(c_tstor_soil_rast), yres(c_tstor_soil_rast)))
    return(vrt_file)
}
rest_c_rast <- stack(rest_vrts)
names(rest_c_rast) <- gsub('_250m\\*', '', rest_rast_patterns)
crs(rest_c_rast) <- crs('EPSG:4326')

############
# Ecosystems
eco_vrt <- tempfile(fileext='.vrt')
eco_in_files <- list.files(
        file.path(
            data_folder_local
        ),
        pattern='ecosystems_250m[-.0-9]*tif$',
        full.names=TRUE
)
gdalbuildvrt(eco_in_files,
             eco_vrt,
             te=output_extent,
             tr=c(xres(c_tstor_soil_rast),
                  yres(c_tstor_soil_rast)))
ecosystems_key <- data.frame(
    id=c(1, 2, 3, 4, 5, 6, 7, 8),
    ecosystem_name=c("Primary forest",
                     "Secondary forest",
                     "Grassland",
                     "Wetlands",
                     "Mangroves",
                     "Salt marsh",
                     "Seagrass",
                     "Peatland")
)
eco_rast <- raster(eco_vrt)
names(eco_rast) <- 'ecosystem'
    
##########################################
# Population
pop_vrt <- tempfile(fileext='.vrt')
pop_in_files <- list.files(
    file.path(
        data_folder_deg_paper,
        'Degradation_Paper',
        'GEE_Rasters'),
    pattern='pop_count_2000_05_10_15_20[-.0-9]*tif$',
    full.names=TRUE
)
gdalbuildvrt(pop_in_files,
             pop_vrt,
             te=output_extent,
             tr=c(xres(c_tstor_soil_rast),
                  yres(c_tstor_soil_rast)))
population <- stack(pop_vrt)
names(population) <- c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015', 'pop_2020')
#population <- crop(population, regions)
NAvalue(population) <- -32768
population <- population$pop_2020
names(population) <- 'population'

###############################################################################
### Extract values
rasts <- stack(c_tstor_soil_rast,
               c_tstor_woody_rast$cb2020,
               c_tstor_ic_rast,
               eco_rast,
               population,
               rest_c_rast)
names(rasts)[names(rasts) == 'cb2020'] <- 'c_tstor_woody'

exact_extract(rasts,
              sites_sp,
              include_cell=TRUE,
              include_cols=c('id', 'reporting_year')) %>%
    rbindlist() %>%
    rename(site_id=id,
           pixel_id=cell) %>%
    relocate(pixel_id, site_id, reporting_year, coverage_fraction) -> pixels
print(nrow(pixels))
pixels %>%
    select(pixel_id, site_id, coverage_fraction) -> pixelsbysite
pixels %>%
    select(-site_id, -coverage_fraction) %>%
    rename(id=pixel_id) %>%
    distinct(id, reporting_year, .keep_all=TRUE) -> pixels
print(nrow(pixels))

# Get pixel areas in hectares, and join to pixels table
sourceCpp('area_ha.cpp')
exact_extract(
    rasts[[1]],
    sites_sp,
    fun=area_ha,
    summarize_df=TRUE,
    include_xy=TRUE,
    include_cell=TRUE,
    xres=xres(rasts),
    yres=yres(rasts),
    use_cov_frac=FALSE) %>%
    rename(id=cell) %>%
    distinct(id, .keep_all=TRUE) %>%
    right_join(pixels) -> pixels
print(nrow(pixels))

# rescale the population data to account for change in resolution
population_original <- raster(pop_in_files[1])
scaling <- (xres(population) * yres(population)) /
           (xres(population_original) * yres(population_original))
pixels$population <- pixels$population * scaling

# Topcode population
pixels$population[pixels$population > 5000] <- 5000
pixels$population[is.na(pixels$population)] <- 0

saveRDS(pixels, 'tables/pixels.rds')
saveRDS(pixelsbysite, 'tables/pixelsbysite.rds')

sites <- readRDS(file.path(data_folder_local, 'tables/sites.rds'))
pixels <- readRDS(file.path(data_folder_local, 'tables/pixels.rds'))
pixelsbysite <- readRDS(file.path(data_folder_local, 'tables/pixelsbysite.rds'))

# Recode sequestration potential from pixel data
lengthyr = 1 # length of intervention in years
left_join(
    pixels,
    pixelsbysite,
    by=c('id'='pixel_id')
) %>%
    left_join(
        (
         sites %>%
            select(id, restoration_type)
        )
        , by=c('site_id'='id')
    ) %>%
    mutate(
        c_potl_seq = recode(
            restoration_type,
            'Agroforestry' = c_potl_seq_agfor0020 * lengthyr,
            'Enrichment Planting/Assisted Natural Regeneration' = c_potl_seq_natre0020 * lengthyr,
            'Natural Regeneration' = c_potl_seq_natre0020 * lengthyr,
            'Mangrove Tree Restoration' = c_potl_seq_mtrer0020 * lengthyr,
            'Mangrove Shrub Restoration' = c_potl_seq_mshrr0020 * lengthyr,
            'Silvopasture' = c_potl_seq_agfor0020 * .15 * lengthyr,
            'Seed dispersal' = c_potl_seq_natre0020 * .60 * lengthyr,
            'Rangeland Restoration - Planned Grazing' = 3.67 * lengthyr,
            .default=0
        ),
        c_potl_seq_type_verification = recode(
            restoration_type,
            'Agroforestry' = 1,
            'Enrichment Planting/Assisted Natural Regeneration' = 2,
            'Natural Regeneration' = 3,
            'Mangrove Tree Restoration' = 4,
            'Mangrove Shrub Restoration' = 5,
            'Silvopasture' = 6,
            'Seed dispersal' = 7, 
            'Rangeland Restoration - Planned Grazing' = 8,
            .default=0
        )
    ) -> sequestration_potential
saveRDS(sequestration_potential, file.path(data_folder_local, 'tables/c_potl_seq.rds'))
pixels %>%
    select(
        -starts_with('c_potl_seq')
    ) %>%
    left_join(
        sequestration_potential %>%
        select(id, c_potl_seq)
    ) -> pixels_with_seq
# Carry over maximum value of potential sequestraton within pixels
pixels_with_seq %>% 
    group_by(id) %>%
    mutate(c_potl_seq=max(c_potl_seq)) %>%
    distinct(id, .keep_all=TRUE) %>%
    ungroup() -> pixels_seq_unique
saveRDS(
    select(
        pixels_seq_unique,
        id,
        area_ha,
        c_potl_seq
    ),
    file.path(
        data_folder_local,
        'tables/pixels_with_seq.rds'
    )
)

nrow(pixels_seq_unique)
table(pixels_seq_unique$c_potl_seq)
