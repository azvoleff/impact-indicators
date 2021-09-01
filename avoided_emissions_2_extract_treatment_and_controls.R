library(raster)
library(fasterize)
library(sf)
library(tidyverse)
library(units)
library(foreach)
library(mapview)
library(exactextractr)
library(raster)
library(tictoc)
library(gdalUtils)

options("optmatch_max_problem_size"=Inf)

data_folder <- 'D:/Data/Impacts_Data'
data_folder_avoided_emissions <- 'avoided_emissions_data'

gdal_crop <- function(r, s) {
    if (filename(r) == '') {
        stop('input raster must have a valid filename')
    }
    e <- extent(s)
    band_names <- names(r)
    out_file <- tempfile(fileext='.tif')
    sf::gdal_utils('warp', filename(r), out_file,
                   options=c('-te', e@xmin, e@ymin, e@xmax, e@ymax, '-multi', 
                             '-co', 'COMPRESS=LZW'))
    out_r <- brick(out_file)
    names(out_r) <- names(r)
    return(out_r)
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}

###############################################################################
### Final data setup

f <- treatment ~ lc_2015_agriculture + precip + temp + elev + slope + 
    dist_cities + dist_roads + crop_suitability + pop_2015 + pop_growth + 
    total_biomass
saveRDS(f, file.path(data_folder_avoided_emissions, 'formula.RDS'))

covariates <- brick(file.path(data_folder_avoided_emissions, 'covariates_covariates.tif'))
names(covariates) <- read_csv(file.path(data_folder_avoided_emissions, 'covariates_covariates.csv'))$names
lc_2015 <- brick(file.path(data_folder_avoided_emissions, 'covariates_lc_2015.tif'))
names(lc_2015) <- read_csv(file.path(data_folder_avoided_emissions, 'covariates_lc_2015.csv'))$names

##############
# Woody carbon
#
# Note: can't use the load_as_vrt function as for this case there are just way 
# too many files for the woody carbon data, so load them manually into a vrt. 

# Need to reproject and crop to match the other layers
ext <- extent(covariates)
output_extent <- c(ext[1], ext[3], ext[2], ext[4])
output_dim <- c(ncol(covariates), nrow(covariates))
c_tstor_woody_vrt_file <- tempfile(fileext='.vrt')
gdalbuildvrt(
    file.path(
        data_folder, 
        'woody_carbon_stored_2020*.tif'
    ), 
    c_tstor_woody_vrt_file,
    te=output_extent,
    tr=c(xres(covariates), yres(covariates))
)
c_tstor_woody_rast <- brick(c_tstor_woody_vrt_file)
crs(c_tstor_woody_rast) <- crs('EPSG:4326')
names(c_tstor_woody_rast) <-c(
    paste0("fc_", 2000:2020),
    paste0("fl_", 2001:2020),
    paste0("cb_", 2000:2020),
    paste0("ce_", 2001:2020)
)
fc <- c_tstor_woody_rast[[which(grepl('fc_', names(c_tstor_woody_rast)))]]

# gdal_translate(
#     c_tstor_woody_vrt_file,
#     file.path(data_folder_avoided_emissions, 'fc_2000_2020.tif'),
#     co="COMPRESS=LZW"
# )

d <- stack(covariates, lc_2015, fc)
# Ensure only layers in the formula are included (so extra data isn't being 
# passed around)
d <- d[[c(get_names(f),
          'region',
          'ecoregion',
          'pa',
           paste0('fc_', 2000:2020)
         )]]
write_csv(data.frame(names=names(d)), file='all_covariates_names.csv')

###############################################################################
###  Load sites and covariates

sites <- readRDS('avoided_emissions_data/sites_cleaned_for_avoided_emissions.RDS')
dim(sites)

# Drop sites with no overlap with GADM (since they'd throw errors later during 
# the extraction) - these are marine sites
sites <- filter(sites, !(id %in% c('242002', '242114')))
dim(sites)

# Filter to only sites over 100 ha
sites <- sites[!sites$area_cea < as_units(100, 'hectares'), ]
dim(sites)

regions <- readRDS(file.path(data_folder_avoided_emissions, 'regions.RDS'))

regions_rast <- fasterize(regions, raster(d[[1]]), field='level1_ID')
names(regions_rast) <- 'region'

d <- stack(d, regions_rast)

###############################################################################
###  Load sites and covariates

# Run extractions of treatment points individually to make catching any polygon 
# errors easier
exact_extract(
    d$region,
    sites,
    include_cell=TRUE, 
    include_cols=c('id_numeric', 'reporting_year'),
    force_df=TRUE
) %>%
    bind_rows() %>%
    rename(region=value) %>%
    filter(!is.na(region)) -> treatment_key
saveRDS(treatment_key, 'avoided_emissions_data/treatment_cell_key.RDS')

regions %>%
    filter(level1_ID %in% unique(treatment_key$region)) -> regions_filtered

quarter <- function(polys) {
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
regions_filtered_q <- foreach(n=1:nrow(regions_filtered), .combine=rbind) %do% {
    suppressWarnings(suppressMessages(quarter(regions_filtered[n, ])))
}

# Fix Fiji, which has a geometrycollection - convert it to a set of polygons 
# (it is level1_ID 908, row 441)
geom_coll <- regions_filtered_q %>% filter(st_is(., 'GEOMETRYCOLLECTION'))
regions_filtered_q <- regions_filtered_q %>% filter(!st_is(., 'GEOMETRYCOLLECTION'))
geom_coll <- st_cast(geom_coll)
regions_filtered_q <- bind_rows(regions_filtered_q, geom_coll)

geom_coll[[1,]]

# Run extraction of control and treatment data by region to make the problem 
# tractable in-memory
n <- 1
out <- foreach(
    this_region_ID=unique(regions_filtered_q$level1_ID)[44], 
    .packages=c('exactextractr', 'sf')
) %do% {
    print(paste0('Processing region ', n, ' of ', length(unique(treatment_key$region)), ' (id ',  this_region_ID, ')...'))
    this_file <- paste0(
        'avoided_emissions_data/extracted_covariates/treatments_and_controls_',
        this_region_ID,
        '.RDS'
    )
    if (file.exists(this_file)) {
        print(paste0('Skipping ', this_region_ID, '. Already processed.'))
    } else {
        this_region <- filter(regions_filtered_q, level1_ID == this_region_ID)
        covariate_values <- exact_extract(
            d,
            this_region,
            include_cell=TRUE
        )
        saveRDS(covariate_values, this_file)
    }
    n <- n + 1
}
