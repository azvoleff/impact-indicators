library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyverse)
library(foreach)
library(optmatch)
library(lubridate)
library(biglm)
library(tictoc)

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
        if (sum(this_d$treatment) > 30) {
            model <- glm(f, data=this_d, family=binomial())
            dists <- match_on(model, data=this_d)
        } else {
            dists <- match_on(f, data=this_d)
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

treatment_key <- readRDS('Output/treatment_cell_key.RDS')
readRDS('sites.RDS') %>%
    select(-geometry) %>%
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

###############################################################################
###  Run matching

set.seed(31)
ae <- foreach(this_year=unique(treatment_key$Data_Year),
              .combine=foreach_rbind, .inorder=FALSE) %do% {
    foreach(this_CI_ID=unique(treatment_key$CI_ID),
            .combine=foreach_rbind, .inorder=FALSE) %do% {
        tic()
        ###############
        # Load datasets
        
        site <- filter(sites,
                       CI_ID == this_CI_ID,
                       Data_Year == this_year)
        if (file.exists(paste0('Output/m_', this_CI_ID, '_', this_year, '.RDS'))) {
            print(paste0('Skipping ', this_CI_ID, ' for year ', this_year, '. Already processed.'))
            return(NULL)
        } else {
            print(paste0('Processing ', this_CI_ID, ' for year ', this_year, '.'))
        }

        treatment_cell_IDs <- filter(treatment_key,
                                     CI_ID == this_CI_ID,
                                     !is.na(region),
                                     Data_Year == this_year)
        n_treatment_cells_total <- nrow(treatment_cell_IDs)
        if (n_treatment_cells_total == 0) {
            print(paste0('Skipping ', this_CI_ID, ' for year ', this_year, '. No treatment cells.'))
            return(NULL)
        }
        vals <- foreach(this_region = unique(treatment_cell_IDs$region),
                        .combine=rbind) %do% {
            v <- readRDS(paste0('Output/treatments_and_controls_', this_region, 
                                '.RDS'))
            filter(v, region == this_region)
        }

        vals %>% full_join(
                treatment_cell_IDs %>%
                    select(cell, Data_Year) %>%
                    mutate(treatment=TRUE)
                , by='cell') -> vals
        vals$treatment <- as.logical(vals$treatment)
        vals$treatment[is.na(vals$treatment)] <- FALSE
        vals$Data_Year <- this_year

        # Remove areas falling within another CI site from the control sample 
        # (but DON'T remove those areas falling within this site)
        filter(vals,
               !(cell %in% filter(treatment_key,
                                  !(cell %in% filter(treatment_key,
                                                     CI_ID == this_CI_ID)$cell))$cell)) -> vals

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
            print(paste0(this_CI_ID, ': Filtered ', n_filtered, ' rows due to missing data in grouping variables.'))
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
        estab_year <- year(site$CI_Start_Date_clean)
        f <- readRDS('Output/formula.RDS')
        if (estab_year >= 2005) {
            init <- vals[, grepl(paste0('fc_', substr(estab_year - 5, 3, 4)), names(vals))]
            final <- vals[,grepl(paste0('fc_', substr(estab_year, 3, 4)), names(vals))]
            defor_pre_intervention <- ((final - init) / init) * 100
            # Correct for division by zero in places that had no forest cover in 
            # year 0
            defor_pre_intervention[init == 0] <- 0
            names(defor_pre_intervention) <- 'defor_pre_intervention'
            vals <- cbind(vals, defor_pre_intervention)
            f <- update(f, ~ . + defor_pre_intervention)
        }
        #vals <- vals %>% select(-starts_with('fc_'), -starts_with('fcc_'))
        
        sample_sizes <- vals %>%
            count(treatment, group)
        print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, treatment)$n, collapse=', '), ' treatment pixels'))
        print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, !treatment)$n, collapse=', '), ' control pixels'))

        ##############
        # Run matching
        
        if (nrow(filter(vals, treatment)) == 0) {
            print(paste0(this_CI_ID, ': No treatment values remaining after filtering'))
            return(NULL)
        } else {
            m <- match_ae(vals, f)
            print(paste0(this_CI_ID, ': Formatting output'))
            if (is.null(m)) {
                print(paste0(this_CI_ID, ': no matches'))
            } else {
                m$CI_ID <- this_CI_ID
                m$Data_Year <- this_year
                m <- m %>% dplyr::select(CI_ID, everything())
                print(paste0(this_CI_ID, ': saving output'))
                m$sampled_fraction <- sum(vals$treatment) / n_treatment_cells_total
                saveRDS(m, paste0('Output/m_', this_CI_ID, '_', this_year, '.RDS'))
            }
        }
        toc()
        return(m)
    }
}

###############################################################################
# Load all output and resave in one file
m <- foreach(f=list.files('Output', pattern ='^m_[0-9]*_[0-9]{4}.RDS$'), 
              .combine=foreach_rbind) %do% {
    readRDS(paste0('Output/', f))
}
saveRDS(m, 'Output/m_ALL.RDS')


###############################################################################
# Summarize results by site

readRDS('sites.RDS') %>%
    select(CI_ID,
           Data_Year,
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
                   CI_ID,
                   Data_Year, 
                   treatment,
                   sampled_fraction,
                   total_biomass,
                   starts_with('fc_')) %>%
            left_join(sites, by=c('CI_ID', 'Data_Year')) %>%
            gather(year, forest_at_year_end, starts_with('fc_')) %>%
            mutate(year=2000 + as.numeric(str_replace(year, 'fc_', ''))) %>%
            group_by(CI_ID, Data_Year, cell, treatment) %>%
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
        group_by(CI_ID, Data_Year, year, treatment) %>%
        summarise(CI_Start_Year=CI_Start_Year[1],
                  CI_End_Year=CI_End_Year[1],
                  # correct totals for areas where only a partial sample was used 
                  # by taking into account the fraction sampled
                  forest_loss_ha=sum(forest_loss_during_year, na.rm=TRUE) * (1 / sampled_fraction[1]),
                  Emissions_MgCO2e=sum(Emissions_MgCO2e, na.rm=TRUE) * (1 / sampled_fraction[1]),
                  n_pixels=n()) %>%
        as_tibble() -> this_m_site

    toc()
    gc()
    return(this_m_site)
}

saveRDS(m_site, file='output_raw_by_site.RDS')
write_csv(m_site, 'output_raw_by_site.csv')
