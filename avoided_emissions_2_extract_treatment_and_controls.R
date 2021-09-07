library(fasterize)
library(raster)
library(sf)
library(tidyverse)
library(units)
library(foreach)
library(mapview)
library(exactextractr)
library(raster)
library(tictoc)
library(gdalUtils)
library(Rcpp)

options(rasterTmpDir='/data/tmp/')

options("optmatch_max_problem_size"=Inf)

#data_folder <- 'D:/Data/Impacts_Data'
#data_folder_avoided_emissions <- 'avoided_emissions_data'

data_folder <- '/home/rstudio/data'
data_folder_avoided_emissions <- '/home/rstudio/data/impacts_data'
code_folder <- '/data/code/impact-indicators'

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
# Forest cover
# Below vrt was created with command:
# /usr/bin/gdalbuildvrt -te -179.92492803 -60.14164261 179.84992806 82.79163355 -tr 0.00833333 0.00833333 -overwrite /data/Hansen_250m/Hansen_GFC-2020-v1.8_250m_coverbyyear.vrt /home/rstudio/data/Hansen_250m/*.tif
fc_rast <- brick(file.path(
    data_folder,
    'Hansen_250m',
    'Hansen_GFC-2020-v1.8_250m_coverbyyear.vrt'
))
crs(fc_rast) <- crs('EPSG:4326')
names(fc_rast) <- paste0("fc_", 2000:2020)

# gdal_translate(
#     fc_vrt_file,
#     file.path(data_folder_avoided_emissions, 'fc_2000_2020.tif'),
#     co="COMPRESS=LZW"
# )

d <- stack(covariates, lc_2015, fc_rast)
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

sites <- readRDS(file.path(data_folder_avoided_emissions, 'sites_cleaned_for_avoided_emissions.RDS'))
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
sourceCpp(file.path(code_folder, 'area_ha.cpp'))
exact_extract(
    d[[1]],
    sites,
    fun=area_ha,
    summarize_df=TRUE,
    include_cell=TRUE,
    include_xy=TRUE,
    xres=xres(d),
    yres=yres(d),
    use_cov_frac=FALSE
) %>%
    bind_rows() %>%
    distinct() -> areas
treatment_key %>%
    left_join(areas) -> treatment_key
saveRDS(treatment_key, file.path(data_folder_avoided_emissions, 'treatment_cell_key.RDS'))

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

cl <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl)
regions_filtered_q <- foreach(
    n=1:nrow(regions_filtered),
    .combine=rbind,
    .packages=c('foreach', 'sf', 'raster')
) %dopar% {
    suppressWarnings(suppressMessages(quarter(regions_filtered[n, ])))
}

# Remove Fiji, which has a geometrycollection that will fail in extraction. Need
# to find a fix (it is level1_ID 908, row 441)
regions_filtered_q <- regions_filtered_q %>% filter(!st_is(., 'GEOMETRYCOLLECTION'))
#geom_coll <- regions_filtered_q %>% filter(st_is(., 'GEOMETRYCOLLECTION'))
#geom_coll <- st_cast(geom_coll)
#regions_filtered_q <- bind_rows(regions_filtered_q, geom_coll)

exact_extract(
    d,
    regions_filtered_q,
    include_cell=TRUE
) %>%
    bind_rows() %>%
    filter(coverage_fraction >= .99) %>%
    select(-coverage_fraction) %>%
    distinct() -> covariate_values
exact_extract(
    d[[1]],
    regions_filtered_q,
    fun=area_ha,
    summarize_df=TRUE,
    include_cell=TRUE,
    include_xy=TRUE,
    xres=xres(d),
    yres=yres(d),
    use_cov_frac=FALSE
) %>%
    bind_rows() %>%
    distinct() -> areas
dim(covariate_values)
covariate_values %>%
    left_join(areas) -> covariate_values
saveRDS(
    covariate_values,
    file.path(
        data_folder_avoided_emissions,
        'treatments_and_controls.RDS'
    )
)