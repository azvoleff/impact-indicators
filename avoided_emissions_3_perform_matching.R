library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyverse)
library(foreach)
library(optmatch)
library(lubridate)
library(biglm)
library(tictoc)
library(doParallel)

data_folder <- '/home/rstudio/data/impacts_data'

options("optmatch_max_problem_size"=Inf)

MAX_TREATMENT <- 1000
CONTROL_MULTIPLIER <- 50
    
# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
foreach_rbind <- function(d1, d2) {
    if (is.null(d1) & is.null(d2)) {
        return(NULL)
    } else if (!is.null(d1) & is.null(d2)) {
        return(d1)
    } else if (is.null(d1) & !is.null(d2)) {
        return(d2)
    } else  {
        return(bind_rows(d1, d2))
    }
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}


get_matches <- function(d, dists) {
    # If the controls are too far from the treatments (due to a caliper) then 
    # the matching may fail. Can test for this by seeing if subdim runs 
    # successfully
    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                             error=function(e) return(FALSE))
    if (subdim_works) {
        #m <- pairmatch(dists, data=d, remove.unmatchables = TRUE)
        m <- fullmatch(dists, min.controls=1, max.controls=1, data=d)
        d <- d[matched(m), ]
    } else {
        d <- data.frame()
    }
    return(d)
}

match_ae <- function(d, f) {
    m <- foreach(this_group=unique(d$group), .combine=foreach_rbind) %do% {
        this_d <- filter(d, group == this_group)
        # Calculate propensity scores with a GLM, or else use Mahalanobis 
        # distance if there aren't enough points to run a glm
        if (sum(this_d$treatment) < 1) {
            return(NULL)
        } else if (sum(this_d$treatment) < 15) {
            dists <- match_on(f, data=this_d)
        } else {
            
            model <- glm(f, data=this_d, family=binomial())
            dists <- match_on(model, data=this_d)
        }
        return(get_matches(this_d, dists))
    }
    # Need to handle the possibility that there were no matches for this 
    # treatment, meaning d will be an empty data.frame
    if (nrow(m) == 0) {
        return(NULL)
    } else {
        return(m)
    }
}

###############################################################################
###  Load sites and covariates
treatment_key <- readRDS(file.path(data_folder, 'treatment_cell_key.RDS'))
readRDS(file.path(data_folder, 'sites_cleaned_for_avoided_emissions.RDS')) %>%
    as_tibble() -> sites

# Filter to include only values of group that appear in the treatment pixels, 
# and to not include values that appear only in the treatment pixels
filter_groups <- function(vals) {
    vals$group <- interaction(vals$region, vals$ecoregion, vals$pa)
    vals <- filter(vals, group %in% unique(filter(vals, treatment)$group))
    treatment_groups <- unique(filter(vals, treatment)$group)
    control_groups <- unique(filter(vals, !treatment)$group)
    vals <- filter(vals, group %in% treatment_groups[treatment_groups %in% control_groups])
    # Filter out values of group that appear ONLY in the treatment pixels
    vals$group <- droplevels(vals$group)
    return(vals)
}

this_id <- unique(treatment_key$id_numeric)[2]

###############################################################################
###  Run matching
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
set.seed(31)
ae <- foreach(
    this_id=unique(treatment_key$id_numeric)[1:2],
    .combine=c,
    .inorder=FALSE,
    .packages=c('tidyverse', 'optmatch', 'sf', 'foreach')
) %dopar% {
    
    ###############
    # Load datasets
    
    site <- filter(sites, id_numeric == this_id)
    match_path <- file.path(data_folder, 'matches', paste0('m_', this_id, '.RDS'))
    if (file.exists(match_path)) {
        print(paste0('Skipping ', this_id, '. Already processed.'))
        return(NULL)
    } else {
        print(paste0('Processing ', this_id, '.'))
    }

    treatment_cell_IDs <- filter(treatment_key,
                                 id_numeric == this_id,
                                 !is.na(region))
    n_treatment_cells_total <- nrow(treatment_cell_IDs)
    if (n_treatment_cells_total == 0) {
        print(paste0('Skipping ', this_id, '. No treatment cells.'))
        return(NULL)
    }
    vals <- foreach(this_region = unique(treatment_cell_IDs$region),
                    .combine=rbind) %do% {
        readRDS(
            file.path(
                data_folder,
                'extracted_covariates',
                paste0(
                    'treatments_and_controls_',
                    this_region,
                    '.RDS'
                    )
                )
            ) %>%
            bind_rows() %>%
            filter(region == this_region)
    }

    vals %>% full_join(
            treatment_cell_IDs %>%
                select(cell) %>%
                mutate(treatment=TRUE)
            , by='cell') -> vals
    vals$treatment <- as.logical(vals$treatment)
    vals$treatment[is.na(vals$treatment)] <- FALSE

    # Remove areas falling within another CI site from the control sample 
    # (but DON'T remove those areas falling within this site)
    vals %>%
        filter(cell %in% filter(treatment_key, id_numeric == this_id)$cell) -> treatment_vals
    vals %>%
        filter(!(cell %in% treatment_key$cell)) -> potential_control_vals
    bind_rows(treatment_vals, potential_control_vals) -> vals

    ################
    # Setup grouping
    
    # Eliminate any pixels with NAs in group variables (happens occasionally 
    # where polygons overlap some ocean, leading to undefined ecoregion, for 
    # example)
    n_filtered <- nrow(vals)
    vals <- filter(vals,
                   !is.na(region),
                   !is.na(ecoregion),
                   !is.na(pa))
    n_filtered <- n_filtered - nrow(vals)
    if (n_filtered > 0) {
        print(paste0(this_id, ': Filtered ', n_filtered, ' rows due to missing data in grouping variables.'))
    }
    # Eliminate any groups that are only in the control pixels, or only in the 
    # treatment pixels
    vals <- filter_groups(vals)
    sample_sizes <- vals %>%
        count(treatment, group)
    # Sample the treatment cells if there are more than MAX_TREATMENT pixels, 
    # and the control cells if there are more than CONTROL_MULTIPLIER * 
    # MAX_TREATMENT pixels
    bind_rows(
            filter(vals, treatment)  %>%
                group_by(group) %>%
                sample_n(min(MAX_TREATMENT, n())),
            filter(vals, !treatment)  %>%
                group_by(this_group=group) %>%
                sample_n(min(CONTROL_MULTIPLIER * filter(sample_sizes,
                                                         treatment == TRUE,
                                                         group == this_group[1])$n,
                             n()))
            ) %>% 
        ungroup() %>%
        select(-this_group) -> vals
    # Refilter in case any groups were lost due to the sampling
    vals <- filter_groups(vals)

    # Project all items to cylindrical equal area
    # d_crop <- projectRaster(d_crop, crs=CRS('+proj=cea'), method='ngb')
    
    ################
    # Add defor data
    
    # For sites that were established in or after 2005, match on the five years 
    # of deforestation data preceding the year of establishment. For sites 
    # estab prior to 2005, don't match on defor rate
    estab_year <- site$ci_start_year
    f <- readRDS(file.path(data_folder, 'formula.RDS'))
    if (estab_year >= 2005) {
        init <- vals[, grepl(paste0('fc_20', substr(estab_year - 5, 3, 4)), names(vals))]
        final <- vals[,grepl(paste0('fc_20', substr(estab_year, 3, 4)), names(vals))]
        defor_pre_intervention <- ((final - init) / init) * 100
        # Correct for division by zero in places that had no forest cover in 
        # year 0
        defor_pre_intervention[init == 0] <- 0
        # Remove pixels that had no forest cover in year zero
        vals <- filter(vals, defor_pre_intervention != 0, !(is.na(defor_pre_intervention)))
        # Refilter in case any groups were lost due to the filtering for
        # pixels with zero forest cover
        vals <- filter_groups(vals)
        names(defor_pre_intervention) <- 'defor_pre_intervention'
        vals <- cbind(vals, defor_pre_intervention)
        f <- update(f, ~ . + defor_pre_intervention)
    }
    #vals <- vals %>% select(-starts_with('fc_'), -starts_with('fcc_'))
    
    sample_sizes <- vals %>%
        count(treatment, group)
    print(paste0(this_id, ': ', paste(filter(sample_sizes, treatment)$n, collapse=', '), ' treatment pixels'))
    print(paste0(this_id, ': ', paste(filter(sample_sizes, !treatment)$n, collapse=', '), ' control pixels'))

    ##############
    # Run matching
    
    if (nrow(filter(vals, treatment)) == 0) {
        print(paste0(this_id, ': No treatment values remaining after filtering'))
        return(NULL)
    } else {
        m <- match_ae(vals, f)
        print(paste0(this_id, ': Formatting output'))
        if (is.null(m)) {
            print(paste0(this_id, ': no matches'))
        } else {
            m$id <- this_id
            m <- m %>% dplyr::select(id, everything())
            print(paste0(this_id, ': saving output'))
            m$sampled_fraction <- sum(vals$treatment) / n_treatment_cells_total
            saveRDS(m, match_path)
        }
    }
    return(1)
}

###############################################################################
# Load all output and resave in one file
m <- foreach(f=list.files(file.path(data_folder, 'matches'), pattern ='^m_[0-9]*.RDS$'), 
              .combine=foreach_rbind) %dopar% {
    readRDS(file.path(data_folder, 'matches', f))
}
saveRDS(m, file.path(data_folder, 'matches', 'm_ALL.RDS'))


###############################################################################
# Summarize results by site

readRDS('sites.RDS') %>%
    select(id,
           CI_Start_Date_clean, 
           CI_End_Date_clean) %>%
    mutate(CI_Start_Year=year(CI_Start_Date_clean),
           CI_End_Year=ifelse(is.na(year(CI_End_Date_clean)), 2099, year(CI_End_Date_clean))) %>%
    select(-CI_Start_Date_clean, -CI_End_Date_clean) -> sites

get_chunk <- function(d, n, n_chunks=10) {
    start_ind <- round(seq(1, length(d), length.out=n_chunks + 1))[1:(n_chunks)]
    end_ind <- c(start_ind[2:n_chunks] - 1, length(d))
    return(d[start_ind[n]:end_ind[n]])
}

# Process in chunks to save memory
n_chunks <- 20
data_files <- list.files('Output', pattern ='^m_[0-9]*_[0-9]{4}.RDS$')
m_site <- foreach (i=1:n_chunks, .combine=bind_rows) %do% {
#m_site <- foreach (i=10:11, .combine=bind_rows) %do% {
    print(paste0('Progress: ', ((i-1)/n_chunks)*100, '%'))
    tic()
    this_m <- foreach(f=get_chunk(data_files, i, n_chunks),
                  .combine=bind_rows, .inorder=FALSE) %do% {

        readRDS(paste0('Output/', f)) %>%
            select(cell,
                   id,
                   treatment,
                   sampled_fraction,
                   total_biomass,
                   starts_with('fc_')) %>%
            left_join(sites, by=c('id')) %>%
            gather(year, forest_at_year_end, starts_with('fc_')) %>%
            mutate(year=2000 + as.numeric(str_replace(year, 'fc_', ''))) %>%
            group_by(id, cell, treatment) %>%
            filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
            arrange(cell, year) %>%
            mutate(forest_loss_during_year=c(NA, diff(forest_at_year_end)),
                   forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
                   biomass_at_year_end = total_biomass * forest_frac_remaining,
                   #  to convert biomass to carbon * .5
                   C_change=c(NA, diff(biomass_at_year_end)) * .5,
                   #  to convert change in C to CO2e * 3.67
                   Emissions_MgCO2e=C_change * -3.67) %>%
            filter(between(year, CI_Start_Year[1], CI_End_Year[1])) %>% # drop year prior to project start as no longer needed
            as_tibble() -> x
    }
    this_m %>%
        group_by(id, year, treatment) %>%
        summarise(CI_Start_Year=CI_Start_Year[1],
                  CI_End_Year=CI_End_Year[1],
                  # correct totals for areas where only a partial sample was used 
                  # by taking into account the fraction sampled
                  forest_loss_ha=sum(forest_loss_during_year, na.rm=TRUE) * (1 / sampled_fraction[1]),
                  Emissions_MgCO2e=sum(Emissions_MgCO2e, na.rm=TRUE) * (1 / sampled_fraction[1]),
                  n_pixels=n()) %>%
        as_tibble() -> this_m_site
    return(this_m_site)
}

saveRDS(m_site, file='output_raw_by_site.RDS')
write_csv(m_site, 'output_raw_by_site.csv')
