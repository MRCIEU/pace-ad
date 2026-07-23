library(tidyverse)
library(dplyr)
library(purrr)
library(scales)

models <- with(expand.grid(num = 1:3, letter = c("a","b","c")),
               paste0("m", num, letter))


files<-list.files(pattern="samplesizesrdsf")

read_one <- function(f) {
r<-read.table(f,sep=" ",he=F)
names(r)<-c("cohort","N","N_cases","N_controls")
w<-which(r$cohort%in%c("CHS","nIOW","UKIDS"))
r<-r[-w,]
r$N <- round(as.numeric(r$N),0)
r$N_cases <- round(as.numeric(r$N_cases),0)
r$N_controls <- round(as.numeric(r$N_controls),0)
r$prevalence <- r$N_cases / r$N
fname<-gsub(".txt","",f)
fname<-gsub("samplesizesrdsf","",fname)
r$model <- fname
r
}

all_df <- do.call(rbind, lapply(files, read_one))

# Split into a list of data.frames by model label in `cohort`
samplesizes <- split(all_df, all_df$model)


totals_by_model <- map_dfr(samplesizes, ~ {
  .x %>%
    summarise(
      total_N         = sum(N, na.rm = TRUE),
      total_cases     = sum(N_cases, na.rm = TRUE),
      total_controls  = sum(N_controls, na.rm = TRUE),
      prevalence_overall = total_cases / total_N,
      imbalance       = sum(N_cases + N_controls - N, na.rm = TRUE),
      n_cohorts       = n()
    ) %>%
    mutate(model = unique(.x$model))
})

samplesizes_with_totals <- map(samplesizes, function(df) {
  # Compute totals for this model
  tot <- df %>%
    summarise(
      cohort      = "Total",
      N           = sum(N, na.rm = TRUE),
      N_cases     = sum(N_cases, na.rm = TRUE),
      N_controls  = sum(N_controls, na.rm = TRUE),
      model       = unique(model),
      .groups     = "drop"
    ) %>%
    mutate(prevalence = ifelse(N > 0, N_cases / N, NA_real_)) %>%
    select(names(df))  # ensure same column order as original

  bind_rows(df, tot)
})