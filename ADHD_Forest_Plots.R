
library(tidyr)
library(forestploter)
library(dplyr)
library(grid)
library(tidyverse)

setwd( "C:/Users/ulric/OneDrive - Appalachian State University/Documents/RStudio/ADHD Matched")

forest_data <- read.csv("Results/Annual_Cumulative_Lag3.csv")

dt <- forest_data %>%
  dplyr::select(-X, -estmean)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, Outcome, Category)

dt <- pivot_wider(dt, names_from = Outcome, values_from = c(over_rr_se, over_rr, over_rr_low, over_rr_high))

# Add two blank columns for CI
dt$`ADHD` <- paste(rep(" ", 20), collapse = " ")
dt$`MDD` <- paste(rep(" ", 20), collapse = " ")
dt$`Selfharm` <- paste(rep(" ", 20), collapse = " ")
dt$`Overlap` <- paste(rep(" ", 20), collapse = " ")

# Generate point estimation and 95% CI. Paste two CIs together and separate by line break.
dt$`ADHD RR (95% CI)` <- ifelse(is.na(dt$over_rr_se_ADHD), "",
                                    sprintf("%.2f (%.2f to %.2f)",
                                            dt$over_rr_ADHD, dt$over_rr_low_ADHD, dt$over_rr_high_ADHD))
dt$`MDD RR (95% CI)` <- ifelse(is.na(dt$over_rr_se_MDD), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$over_rr_MDD, dt$over_rr_low_MDD, dt$over_rr_high_MDD))
dt$`Selfharm RR (95% CI)` <- ifelse(is.na(dt$over_rr_se_Selfharm), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$over_rr_Selfharm, dt$over_rr_low_Selfharm, dt$over_rr_high_Selfharm))
dt$`Overlap RR (95% CI)` <- ifelse(is.na(dt$over_rr_se_Overlap), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$over_rr_Overlap, dt$over_rr_low_Overlap, dt$over_rr_high_Overlap))

est = list(dt$over_rr_ADHD,
           dt$over_rr_MDD,
           dt$over_rr_Selfharm,
           dt$over_rr_Overlap)
lower = list(dt$over_rr_low_ADHD,
             dt$over_rr_low_MDD,
             dt$over_rr_low_Selfharm,
             dt$over_rr_low_Overlap) 
upper = list(dt$over_rr_high_ADHD,
             dt$over_rr_high_MDD,
             dt$over_rr_high_Selfharm,
             dt$over_rr_high_Overlap)

dt <- dt %>%
  dplyr::select(Category, 
                `ADHD RR (95% CI)`, 
                'ADHD',
                'MDD RR (95% CI)',
                'MDD',
                'Selfharm RR (95% CI)',
                'Selfharm',
                'Overlap RR (95% CI)',
                'Overlap')

tm <- forest_theme(base_size = 10)
p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(3, 5, 7, 9),
            ticks_at = c(1),
            xlim = c(0.95, 1.05),
            ref_line = 1,
            title = "Annual Lag0-Lag3",
            theme = tm)
plot(p)
