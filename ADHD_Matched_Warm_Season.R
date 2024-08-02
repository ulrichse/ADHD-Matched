library(dplyr)
library(arrow)
library(lubridate)
library(data.table)
library(haven)
library(sjPlot)
library(MASS)
library(grid)
library(forestploter)
library(dlnm)
library(lme4)



setwd( "C:/Users/ulrichse/OneDrive - Appalachian State University/Documents/RStudio/ADHD Matched")

heat <- read.csv("Data/heat.csv") %>%
  mutate(Date=as.Date(date),
         Year=year(Date),
         heatwave=ifelse(Heatwave_days > 0, 1, 0))%>%
  filter(Year >= 2008 & Year <= 2021)%>%
  dplyr::select(Date, region_code, heatwave)

data <- read.csv("Data/sample.csv") 
data[is.na(data)] <- 0

data <- data %>%
  mutate(Date=as.Date(date, format = "%m/%d/%Y"),
         Year=year(Date),
         ADHD_Black = ifelse(Black == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Black = ifelse(Black == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Black = ifelse(Black == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Black = ifelse(Black == 1 & overlap == 1, 1, 0),
         ADHD_White = ifelse(White == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_White = ifelse(White == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_White = ifelse(White == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_White = ifelse(White == 1 & overlap == 1, 1, 0),
         ADHD_POC = ifelse(POC == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_POC = ifelse(POC == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_POC = ifelse(POC == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_POC = ifelse(POC == 1 & overlap == 1, 1, 0),
         ADHD_Male = ifelse(Male == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Male = ifelse(Male == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Male = ifelse(Male == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Male = ifelse(Male == 1 & overlap == 1, 1, 0),
         ADHD_Female = ifelse(Female == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Female = ifelse(Female == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Female = ifelse(Female == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Female = ifelse(Female == 1 & overlap == 1, 1, 0),
         ADHD_Age5to11 = ifelse(Age5to11 == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Age5to11 = ifelse(Age5to11 == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Age5to11 = ifelse(Age5to11 == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Age5to11 = ifelse(Age5to11 == 1 & overlap == 1, 1, 0),
         ADHD_Age12to18 = ifelse(Age12to18 == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Age12to18 = ifelse(Age12to18 == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Age12to18 = ifelse(Age12to18 == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Age12to18 = ifelse(Age12to18 == 1 & overlap == 1, 1, 0),
         ADHD_Age18plus = ifelse(Age18plus == 1 & ADHD_Prisec == 1, 1, 0),
         MDD_Age18plus = ifelse(Age18plus == 1 & MDD_Prisec == 1, 1, 0),
         Selfharm_Age18plus = ifelse(Age18plus == 1 & Selfharm_Prisec == 1, 1, 0),
         Overlap_Age18plus = ifelse(Age18plus == 1 & overlap == 1, 1, 0)
  ) 

datagrp <- data %>%
  group_by(Date, region_code) %>%
  summarize(across(starts_with(c("ADHD", "MDD", "Selfharm", "Overlap")), sum, na.rm = TRUE))

data <- heat %>%
  left_join(datagrp, by=c("Date", "region_code"))%>%
  mutate(zip=as.numeric(region_code),
         date=Date,
         year=year(date),
         doy=yday(date),
         dow=wday(date))

data <- data %>%
  filter(month(Date) >= 5 & month(Date) <= 9)

data[is.na(data)] <- 0

# Function to exclude ZCTAs with no heatwave days during the study period from matching
zctas_without_hw <- function(data){
  
  no_hw <- data %>%
    group_by(zip)%>%
    summarize(hw=sum(heatwave))%>% # Create list of ZCTAs without heatwaves
    filter(hw==0)%>%
    dplyr::select(zip)
  
  data <- data %>% # Remove ZCTAs without heatwaves from the data
    dplyr::filter(!zip %in% no_hw)
  
}

# Create matched dataframe with lags 0 to 7
matched_df_lag7 <- function(data, n_controls = 3, lag_range = 0:7, control_doy_range = -3:3) {
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  for (i in 1:length(zip_list)) {
    df <- subset(dat, zip == zip_list[i])
    
    # Exclude the 3 days within any other heatwave
    df$time <- 1:nrow(df)
    cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
    df$cand_control <- TRUE
    df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE
    
    case_dates <- subset(df, heatwave == 1)
    control_dates <- subset(df, heatwave == 0)
    
    for (j in 1:nrow(case_dates)) {
      # Choose lags (lagged 0 to 3)
      lag_dates <- case_dates[j, ]$date + lag_range
      lag_case <- subset(df, date %in% lag_dates)
      
      # Choose 10 comparable unexposed days for each heatwave-exposed day
      control_range <- case_dates[j, ]$doy + control_doy_range
      control_subset <- subset(control_dates,
                               control_dates$year != case_dates[j, ]$year &
                                 doy %in% control_range &
                                 cand_control)
      controls <- dplyr::sample_n(control_subset, n_controls)
      
      # Choose lagged days for selected unexposed days
      la_con <- c(1:7)
      for (p in 1:length(la_con)) {
        lag_control_dates <- controls$date + la_con[p]
        lag_control_each <- subset(df, date %in% lag_control_dates)
        
        if (p == 1) {
          lag_control <- lag_control_each
        } else {
          lag_control <- rbind(lag_control, lag_control_each)
        }
      }
      j_stratum <- rbind(lag_case, controls, lag_control)
      stratum <- paste("stratum", j, sep = ".")
      j_stratum$stratum <- stratum
      status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
      j_stratum$status <- status
      lag <- c(rep(0:7, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:7), length.out = nrow(lag_control)))
      j_stratum$lag <- lag
      
      if (j == 1) {
        new_df <- j_stratum
      } else {
        new_df <- rbind(new_df, j_stratum)
      }
    }
    if (i == 1) {
      matched_df <- new_df
    } else {
      matched_df <- rbind(matched_df, new_df)
    }
  }
  
  return(matched_df)
  
  gc()
  
}

# Create matched dataframe with lags 0 to 5
matched_df_lag5 <- function(data, n_controls = 3, lag_range = 0:5, control_doy_range = -3:3) {
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  for (i in 1:length(zip_list)) {
    df <- subset(dat, zip == zip_list[i])
    
    # Exclude the 3 days within any other heatwave
    df$time <- 1:nrow(df)
    cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
    df$cand_control <- TRUE
    df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE
    
    case_dates <- subset(df, heatwave == 1)
    control_dates <- subset(df, heatwave == 0)
    
    for (j in 1:nrow(case_dates)) {
      # Choose lags (lagged 0 to 3)
      lag_dates <- case_dates[j, ]$date + lag_range
      lag_case <- subset(df, date %in% lag_dates)
      
      # Choose 10 comparable unexposed days for each heatwave-exposed day
      control_range <- case_dates[j, ]$doy + control_doy_range
      control_subset <- subset(control_dates,
                               control_dates$year != case_dates[j, ]$year &
                                 doy %in% control_range &
                                 cand_control)
      controls <- dplyr::sample_n(control_subset, n_controls)
      
      # Choose lagged days for selected unexposed days
      la_con <- c(1:5)
      for (p in 1:length(la_con)) {
        lag_control_dates <- controls$date + la_con[p]
        lag_control_each <- subset(df, date %in% lag_control_dates)
        
        if (p == 1) {
          lag_control <- lag_control_each
        } else {
          lag_control <- rbind(lag_control, lag_control_each)
        }
      }
      j_stratum <- rbind(lag_case, controls, lag_control)
      stratum <- paste("stratum", j, sep = ".")
      j_stratum$stratum <- stratum
      status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
      j_stratum$status <- status
      lag <- c(rep(0:5, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:5), length.out = nrow(lag_control)))
      j_stratum$lag <- lag
      
      if (j == 1) {
        new_df <- j_stratum
      } else {
        new_df <- rbind(new_df, j_stratum)
      }
    }
    if (i == 1) {
      matched_df <- new_df
    } else {
      matched_df <- rbind(matched_df, new_df)
    }
  }
  
  return(matched_df)
  
  gc()
  
}

# Create matched dataframe with lags 0 to 3
matched_df_lag3 <- function(data, n_controls = 3, lag_range = 0:3, control_doy_range = -3:3) {
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  for (i in 1:length(zip_list)) {
    df <- subset(dat, zip == zip_list[i])
    
    # Exclude the 3 days within any other heatwave
    df$time <- 1:nrow(df)
    cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
    df$cand_control <- TRUE
    df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE
    
    case_dates <- subset(df, heatwave == 1)
    control_dates <- subset(df, heatwave == 0)
    
    for (j in 1:nrow(case_dates)) {
      # Choose lags (lagged 0 to 3)
      lag_dates <- case_dates[j, ]$date + lag_range
      lag_case <- subset(df, date %in% lag_dates)
      
      # Choose 10 comparable unexposed days for each heatwave-exposed day
      control_range <- case_dates[j, ]$doy + control_doy_range
      control_subset <- subset(control_dates,
                               control_dates$year != case_dates[j, ]$year &
                                 doy %in% control_range &
                                 cand_control)
      controls <- dplyr::sample_n(control_subset, n_controls)
      
      # Choose lagged days for selected unexposed days
      la_con <- c(1:3)
      for (p in 1:length(la_con)) {
        lag_control_dates <- controls$date + la_con[p]
        lag_control_each <- subset(df, date %in% lag_control_dates)
        
        if (p == 1) {
          lag_control <- lag_control_each
        } else {
          lag_control <- rbind(lag_control, lag_control_each)
        }
      }
      j_stratum <- rbind(lag_case, controls, lag_control)
      stratum <- paste("stratum", j, sep = ".")
      j_stratum$stratum <- stratum
      status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
      j_stratum$status <- status
      lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
      j_stratum$lag <- lag
      
      if (j == 1) {
        new_df <- j_stratum
      } else {
        new_df <- rbind(new_df, j_stratum)
      }
    }
    if (i == 1) {
      matched_df <- new_df
    } else {
      matched_df <- rbind(matched_df, new_df)
    }
  }
  
  return(matched_df)
  
  gc()
  
}

# Create crossbasis for lags -2 to 7
matched_cb_lagn2_lag7 <- function(data, n_controls = 3, lag_range = -2:7, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(-2, 7),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to generate individual lags from -2 to 7
get_daily_lags <- function(matched_df, matched_cb){
  
  datasets <- list()
  combined_data_list <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit_am5 <- gnm::gnm(formula,
                        eliminate = factor(zip), family = quasipoisson(link = "log"),
                        data = matched_df)
    pred_am5 <- dlnm::crosspred(matched_cb, fit_am5, at = 1)
    
    # Print matRRfit, matRRlow, and matRRhigh
    print(pred_am5$matRRfit)
    print(pred_am5$matRRlow)
    print(pred_am5$matRRhigh)
    
    datasets[[paste0(outcome, "_pred_am5")]] <- pred_am5
    
    combined_data_name <- paste0(outcome, "_combined_data")
    
    # Create an empty data frame with dynamic column names
    max_columns <- max(sapply(datasets, function(x) ncol(x$matRRfit), simplify = TRUE),
                       sapply(datasets, function(x) ncol(x$matRRlow), simplify = TRUE),
                       sapply(datasets, function(x) ncol(x$matRRhigh), simplify = TRUE))
    
    combined_data <- data.frame(matrix(NA, nrow = 0, ncol = max_columns + 2))
    colnames(combined_data) <- c(names(datasets[[1]]$matRRfit), "Type", "Pred_Model")
    
    # Populate the combined data frame
    for (name in names(datasets)) {
      matRRfit <- datasets[[name]]$matRRfit
      matRRlow <- datasets[[name]]$matRRlow
      matRRhigh <- datasets[[name]]$matRRhigh
      
      combined_data <- rbind(combined_data, cbind(matRRfit, Type = "matRRfit", Pred_Model = name))
      combined_data <- rbind(combined_data, cbind(matRRlow, Type = "matRRlow", Pred_Model = name))
      combined_data <- rbind(combined_data, cbind(matRRhigh, Type = "matRRhigh", Pred_Model = name))
    }
    
    # Assign the combined data frame to the list with dynamic name
    combined_data_list[[combined_data_name]] <- combined_data
  }
  
  return(combined_data)
  
  gc()
  
}

# Function to create crossbasis for lags 0 to 3
matched_cb_lag3 <- function(data, n_controls = 3, lag_range = 0:3, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0, 3),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to create crossbasis for lags 0 to 7
matched_cb_lag7 <- function(data, n_controls = 3, lag_range = 0:7, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0, 7),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to create crossbasis for lags 0 to 5
matched_cb_lag5 <- function(data, n_controls = 3, lag_range = 0:5, control_doy_range = -3:3){
  
  data <- zctas_without_hw(data)
  dat <- data
  zip_list <- unique(dat$zip)
  setorder(dat, zip, date)
  
  # Use "dlnm" package to generate the distributed lag function for "heatwave"
  for (i in 1:length(zip_list)) {
    orig_dat <- subset(dat, zip == zip_list[i])
    match_dat <- subset(matched_df, zip == zip_list[i])
    
    orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0, 5),
                                argvar = list(fun = "lin"),
                                arglag = list(fun = "integer"))
    obs_n <- nrow(orig_dat)
    orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
    orig_cb_matr$date <- orig_dat$date
    matched_date <- match_dat %>% dplyr::select(date)
    matched_cb_matr <- orig_cb_matr %>%
      dplyr::right_join(matched_date, by = "date") %>%
      dplyr::select(-date) %>% as.matrix()
    
    if (i == 1) {
      matched_cb_matrix <- matched_cb_matr
    } else {
      matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
    }
    # Add attributes
    matched_dim <- dim(matched_cb_matrix)
    attr <- attributes(orig_cb)
    attr$dim <- matched_dim
    matched_cb <- matched_cb_matrix
    attributes(matched_cb) <- attr
  }
  
  return(matched_cb)
  
  gc()
  
}

# Function to generate cumulative RR from lag0 to lag7
cumulative_lag0_lag7 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 8
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6) + exp(x7) + exp(x8)) / 8, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

# Function to generate cumulative RR from lag0 to lag5
cumulative_lag0_lag5 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 6
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6)) / 6, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

# Cumulative RR for lag0 to lag3
cumulative_lag0_lag3 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 4
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}


# Function to generate cumulative RR from lag0 to lag7
cumulative_lagn2_lag7 <- function(matched_df, matched_cb){
  
  outcome_data <- list()
  
  for (outcome in outcomes) {
    formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
    
    fit <- gnm::gnm(formula,
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df)
    
    pred <- dlnm::crosspred(matched_cb, fit, at = 1)
    
    over_rr <- sum(pred$matRRfit) / 10
    
    library(msm)
    estvar <- pred$vcov
    estmean <- c(pred$coefficients)
    
    over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6) + exp(x7) + exp(x8) + exp(x9) + exp(x10)) / 10, 
                                   estmean, estvar)
    over_rr_low <- over_rr / exp(1.96 * over_rr_se)
    over_rr_high <- over_rr * exp(1.96 * over_rr_se)
    
    # Create a data frame for the current outcome
    outcome_df <- data.frame(outcome = outcome,
                             over_rr = over_rr,
                             over_rr_se = over_rr_se,
                             over_rr_low = over_rr_low,
                             over_rr_high = over_rr_high,
                             estvar = estvar,
                             estmean = estmean)
    
    outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
  }
  
  data_frames <- unname(outcome_data)
  merged_df <- do.call(rbind, data_frames)
  return(merged_df)
  
  gc()
  
}

outcomes <- c("ADHD_Prisec",
              "ADHD_Black", 
              "ADHD_White", 
              "ADHD_POC", 
              "ADHD_Male", 
              "ADHD_Female", 
              "ADHD_Age5to11", 
              "ADHD_Age12to18", 
              "ADHD_Age18plus",
              "MDD_Prisec",
              "MDD_Black", 
              "MDD_White", 
              "MDD_POC", 
              "MDD_Male", 
              "MDD_Female", 
              "MDD_Age5to11", 
              "MDD_Age12to18", 
              "MDD_Age18plus",
              "Selfharm_Prisec",
              "Selfharm_Black", 
              "Selfharm_White", 
              "Selfharm_POC", 
              "Selfharm_Male", 
              "Selfharm_Female", 
              "Selfharm_Age5to11", 
              "Selfharm_Age12to18", 
              "Selfharm_Age18plus",
              "overlap",
              "Overlap_Black", 
              "Overlap_White", 
              "Overlap_POC", 
              "Overlap_Male", 
              "Overlap_Female", 
              "Overlap_Age5to11", 
              "Overlap_Age12to18", 
              "Overlap_Age18plus")

#Run the match design for the warm season for reach of the four outcomes separately: 
# ADHD_Prisec, MDD_Prisec, Selfharm_Prisec, Overlap
# for each lag0 to lag7
# for cumulative lag 0-3, 0-5, and 0-7

matched_df <- read_parquet("matched_df_ADHD_warm_lag7.parquet")
matched_cb <- matched_cb_lag7(data)

combined_data <- get_daily_lags(matched_df, matched_cb)
write.csv(combined_data, "Results/daily_lags_warm.csv")

merged_df <- cumulative_lag0_lag7(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
write.csv(merged_df, "Results/Warm_Cumulative_Lag7.csv")
gc()

matched_cb <- matched_cb_lag5(data)
matched_df <- matched_df_lag5(data)
merged_df <- cumulative_lag0_lag5(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
write.csv(merged_df, "Results/Warm_Cumulative_Lag5.csv")
gc()

matched_df <- matched_df_lag3(data) 
matched_cb <- matched_cb_lag3(data)
merged_df <- cumulative_lag0_lag3(matched_df, matched_cb)
merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
write.csv(merged_df, "Results/Warm_Cumulative_Lag3.csv")
gc()

