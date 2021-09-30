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


data_folder_input <- '/home/rstudio/data/impacts_data/inputs_for_pixels'
sites_data_folder <- '/home/rstudio/data/impacts_data/sites'
tables_folder <- '/home/rstudio/data/impacts_data/tables'
code_folder <- '/home/rstudio/data/code/impact-indicators'

# data_folder_input <- 'D:/Data/Impacts_Data'
# sites_data_folder <- 'D:/Data/Impacts_Data/sites'
# tables_folder <- 'D:/Data/Impacts_Data/tables'
# code_folder <- 'D:/Code/LandDegradation/impact_indicators/extract-indicators'

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


quarter_split <- function(polys) {
    foreach(n=1:nrow(polys), .combine=rbind) %do% {
        p <- polys[n, ]
        xmn <- extent(p)[1]
        xmx <- extent(p)[2]
        ymn <- extent(p)[3]
        ymx <- extent(p)[4]
        mid_x <- (xmx - xmn) / 2 + xmn
        mid_y <- (ymx - ymn) / 2 + ymn
        nw_bounds <- st_sfc(st_polygon(list(matrix(c(xmn, xmn, mid_x, mid_x, xmn,
                                                     ymx, mid_y, mid_y, ymx, ymx), ncol=2))),
                           crs=4326)
        ne_bounds <- st_sfc(st_polygon(list(matrix(c(mid_x, mid_x, xmx, xmx, mid_x,
                                                     ymx, mid_y, mid_y, ymx, ymx), ncol=2))),
                           crs=4326)
        se_bounds <- st_sfc(st_polygon(list(matrix(c(mid_x, mid_x, xmx, xmx, mid_x,
                                                     mid_y, ymn, ymn, mid_y, mid_y), ncol=2))),
                           crs=4326)
        sw_bounds <- st_sfc(st_polygon(list(matrix(c(xmn, xmn, mid_x, mid_x, xmn,
                                                     mid_y, ymn, ymn, mid_y, mid_y), ncol=2))),
                           crs=4326)
        sf_use_s2(FALSE)
        rbind(st_intersection(p, nw_bounds),
              st_intersection(p, ne_bounds),
              st_intersection(p, se_bounds),
              st_intersection(p, sw_bounds))
    }
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

sites <- st_read(file.path(sites_data_folder, "Final_FY21.gdb"))
sites <- st_zm(sites) # Drop zm coords
sites %>%
    rename(id = CI_ID) %>%
    rename_all(tolower) -> sites

sites <- st_transform(sites, 'EPSG:4326')

sites %>%
    mutate(validity = st_is_valid(., reason=TRUE)) %>%
    filter(validity != 'Valid Geometry') %>%
    as_tibble() %>%
    dplyr::select(-shape) %>%
    write_csv(file.path(sites_data_folder, 'site_validity.csv'))

sites %>%
    mutate(validity = st_is_valid(.)) %>%
    filter(validity == TRUE) %>%
    dplyr::select(-validity) -> sites

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

sites_sp_large <- filter(sites_sp, area_ha > 1e7)

sites_sp_large <- suppressWarnings(suppressMessages(quarter_split(sites_sp_large)))
sites_sp_large <- suppressWarnings(suppressMessages(quarter_split(sites_sp_large)))

sites_sp <- filter(sites_sp, area_ha <= 1e7) %>%
    bind_rows(sites_sp_large)

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

sites <- dplyr::select(
    sites,
    -starts_with('ci_sls'),
    -starts_with('ci_division'),
    -starts_with('star_tag'),
    -starts_with('intervention_type'),
    -country_iso_num_code,
    -area_ha,
    -country
    -join_count,
    -target_fid,
    -join_fid,
    -number,
    -numbertext,
    -area_km2,
    -name,
    -iso3,
    -name_vmap,
    -gaborone,
    -notes
)

sites_sp %>%
    dplyr::select(id, shape) %>%
    full_join(
        dplyr::select(sites, -shape)
    ) -> sites_sp_save

saveRDS(sites, file.path(tables_folder, 'sites_cleaned_cols.rds'))
saveRDS(sites_sp_save, file.path(tables_folder, 'sites_sp.rds'))
saveRDS(sls, file.path(tables_folder, 'sls.rds'))
saveRDS(slsbysite, file.path(tables_folder, 'slsbysite.rds'))
saveRDS(countries, file.path(tables_folder, 'countries.rds'))
saveRDS(countrybysite, file.path(tables_folder, 'countrybysite.rds'))
saveRDS(interventions, file.path(tables_folder, 'interventions.rds'))
saveRDS(interventionbysite, file.path(tables_folder, 'interventionbysite.rds'))
saveRDS(startags, file.path(tables_folder, 'startags.rds'))
saveRDS(startagbysite, file.path(tables_folder, 'startagbysite.rds'))
saveRDS(divisions, file.path(tables_folder, 'divisions.rds'))
saveRDS(divisionbysite, file.path(tables_folder, 'divisionbysite.rds'))
saveRDS(sites_sp_save, file.path(tables_folder, 'sites_sp.rds'))
# Save sites for ingestion to GEE (filtering out all columns except for id)
st_write(
    dplyr::select(sites_sp, id),
    file.path(sites_data_folder, "Final_PostVetting_FY2021_for_gee.shp"),
    delete_dsn=TRUE
)

###############################################################################
### Load indicator data

#############
# Carbon (woody, soil, irrecoverable)
soil_c_rast <- load_as_vrt(data_folder_input, 'soc_bc_250m_scaled100_2020[-.0-9]*tif$')
names(soil_c_rast) <- 'c_tstor_soil'
biomass_c_rast <- load_as_vrt(data_folder_input, 'biomass_carbon_gfw_bc_250m_scaled100_2020[-.0-9]*tif$')
names(biomass_c_rast) <- 'c_tstor_woody'
c_tstor_ic_rast <- load_as_vrt(data_folder_input, "irrecoverable_carbon_250m_2018[-.0-9]*tif$")
names(c_tstor_ic_rast) <- 'c_tstor_ic'
c_rasts <- stack(soil_c_rast, biomass_c_rast, c_tstor_ic_rast)
crs(c_rasts) <- crs('EPSG:4326')

# # Get values to use for reprojecting and croping other layers to match
ext <- extent(c_rasts)
output_extent <- c(ext[1], ext[3], ext[2], ext[4])
output_dim <- c(ncol(c_rasts), nrow(c_rasts))

##########################################
# Potential sequestration from restoration
rest_rast_patterns <- c('c_potl_seq_agfor0020_250m[-.0-9]*',
                        'c_potl_seq_agfor2060_250m[-.0-9]*',
                        'c_potl_seq_natre0020_250m[-.0-9]*',
                        'c_potl_seq_natre2060_250m[-.0-9]*',
                        'c_potl_seq_mtrer0020_250m[-.0-9]*',
                        'c_potl_seq_mtrer2060_250m[-.0-9]*',
                        'c_potl_seq_mshrr0020_250m[-.0-9]*',
                        'c_potl_seq_mshrr2060_250m[-.0-9]*',
                        'c_potl_seq_pwtea0020_250m[-.0-9]*',
                        'c_potl_seq_pwpin0020_250m[-.0-9]*',
                        'c_potl_seq_pwoco0020_250m[-.0-9]*',
                        'c_potl_seq_pwobr0020_250m[-.0-9]*',
                        'c_potl_seq_pwoak0020_250m[-.0-9]*',
                        'c_potl_seq_pweuc0020_250m[-.0-9]*')
rest_c_rast <- foreach(p=rest_rast_patterns, .combine=raster::stack) %do% {
    load_as_vrt(data_folder_input, p)
}
names(rest_c_rast) <- gsub('_250m\\[-\\.0-9\\]\\*', '', rest_rast_patterns)
crs(rest_c_rast) <- crs('EPSG:4326')

############
# Ecosystems
eco_rast <- load_as_vrt(data_folder_input, 'ecosystems_250m[-.0-9]*tif$')
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
saveRDS(ecosystems_key, file.path(tables_folder, 'ecosystems.rds'))
names(eco_rast) <- 'ecosystem'
    
##########################################
# Population
pop_vrt <- tempfile(fileext='.vrt')
gdalbuildvrt(file.path(data_folder_input, 'ppp_2020_1km_Aggregated.tif'),
             pop_vrt,
             te=output_extent,
             tr=c(xres(c_rasts),
                  yres(c_rasts)))
population <- raster(pop_vrt)
names(population) <- c('population')

###############################################################################
### Extract values
rasts <- stack(c_rasts,
               eco_rast,
               population,
               rest_c_rast)

exact_extract(rasts,
              sites_sp,
              include_cell=TRUE,
              include_cols=c('id', 'reporting_year')) %>%
    rbindlist() %>%
    rename(site_id=id,
           pixel_id=cell) %>%
    relocate(pixel_id, site_id, reporting_year, coverage_fraction) -> pixels_no_seq
print(nrow(pixels_no_seq))


# Take account of fact that some pixels within some sites might have been split 
# across the polygons created by the quarter_split lines above, and that these 
# pixels therefore need to be joined together within polygons so that there 
# aren't duplicate pixel IDs within the same polygon
#pixels_no_seq %>%
#    dplyr::select(pixel_id, site_id, coverage_fraction) %>%
#    group_by(site_id, pixel_id) %>%
#    summarise(coverage_fraction = sum(coverage_fraction)) -> pixelsbysite

# TODO: fix the above problem if there is a way to process faster. For now
# ignore the above - recognizing it is an insignificant issue overall as it only
# affects a small number of pixels and only in very large sites
pixels_no_seq %>%
    dplyr::select(pixel_id, site_id, coverage_fraction) %>%
    distinct(pixel_id, site_id, .keep_all=TRUE) -> pixelsbysite
print(nrow(pixelsbysite))    

pixels_no_seq %>%
    dplyr::select(-site_id, -coverage_fraction) %>%
    rename(id=pixel_id) %>%
    distinct(id, reporting_year, .keep_all=TRUE) -> pixels_no_seq
print(nrow(pixels_no_seq))

# Get pixel areas in hectares, and join to pixels_no_seq table
sourceCpp(file.path(code_folder, 'area_ha.cpp'))
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
    right_join(pixels_no_seq) -> pixels_no_seq
print(nrow(pixels_no_seq))

# rescale the population data to account for change in resolution
population_original <- raster(
    file.path(data_folder_input, 'ppp_2020_1km_Aggregated.tif')
)
scaling <- (xres(population) * yres(population)) /
           (xres(population_original) * yres(population_original))
pixels_no_seq$population <- pixels_no_seq$population * scaling

saveRDS(pixels_no_seq, file.path(tables_folder, 'pixels_no_seq.rds'))
saveRDS(pixelsbysite, file.path(tables_folder, 'pixelsbysite.rds'))


#########
# Sequestration_potential

sites <- readRDS(file.path(tables_folder, 'sites_cleaned_cols.rds'))
pixels_no_seq <- readRDS(file.path(tables_folder, 'pixels_no_seq.rds'))
pixelsbysite <- readRDS(file.path(tables_folder, 'pixelsbysite.rds'))

# Recode sequestration potential from pixel data
lengthyr = 1 # length of intervention in years
left_join(
    pixels_no_seq,
    pixelsbysite,
    by=c('id'='pixel_id')
) %>%
    left_join(
        (
         sites %>%
             dplyr::select(id, restoration_type)
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
saveRDS(sequestration_potential, file.path(tables_folder, 'c_potl_seq.rds'))

pixels_no_seq %>%
    dplyr::select(
        -starts_with('c_potl_seq')
    ) %>%
    left_join(
        sequestration_potential %>%
        select(id, c_potl_seq)
    ) -> pixels_with_seq_preclean

# Carry over maximum value of potential sequestration within pixels (could 
# happen if multiple sites code different restoration approaches for the 
# same pixel
#pixels_with_seq_preclean %>% 
#    group_by(id) %>%
#    mutate(c_potl_seq=max(c_potl_seq, na.rm=TRUE)) %>%
#    distinct(id, .keep_all=TRUE) %>%
#    ungroup() -> pixels_seq_unique

# Below doesn't seem to work right, otherwise would be a faster way to do the
# above carry over
#pixels_with_seq_preclean_dt <- data.table(pixels_with_seq_preclean)
#pixels_seq_unique <- pixels_with_seq_preclean_dt[
#    pixels_with_seq_preclean_dt[, .I[which.max(c_potl_seq)], by=id]$V1
#]

# Since above two attempts at removing dupes fail, use a simple distinct on
# id for now
pixels_with_seq_preclean %>%
    distinct(id, .keep_all=TRUE) -> pixels_seq_unique
saveRDS(
    pixels_seq_unique,
    file.path(
        tables_folder,
        'pixels_with_seq.rds'
    )
)

nrow(pixels_seq_unique)
table(pixels_seq_unique$c_potl_seq)

#########
# Species
sp_raw <- read_csv(file.path(data_folder_input, 'species', 'species_all_sites_CR-EN-VU.csv'))
sp_raw %>%
    dplyr::select(-site_id) %>%
    distinct() %>%
    mutate(id = 1:n()) %>%
    relocate(id) -> species
saveRDS(species, file.path(tables_folder, 'species.rds'))
left_join(sp_raw, species) %>%
    rename(species_id=id) %>%
    dplyr::select(site_id, species_id) -> speciesbysite 
saveRDS(speciesbysite, file.path(tables_folder, 'speciesbysite.rds'))
