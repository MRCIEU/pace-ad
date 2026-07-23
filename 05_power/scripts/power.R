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

#iow<-which(samplesizes$m1a[,"cohort"]=="IOW") 
#samplesizes$m1a[iow,"N_cases"]<-179
#samplesizes$m1a[iow,"N_controls"]<-437
#samplesizes$m1a[iow,"prevalence"]<-samplesizes$m1a[iow,"N_cases"]/samplesizes$m1a[iow,"N"]

#iow<-which(samplesizes$m2a[,"cohort"]=="IOW") 
#samplesizes$m2a[iow,"N_cases"]<-99
#samplesizes$m2a[iow,"N_controls"]<-569
#samplesizes$m2a[iow,"prevalence"]<-samplesizes$m2a[iow,"N_cases"]/samplesizes$m2a[iow,"N"]

chs_c<-which(samplesizes$m3c[,"cohort"]=="CHS") 
samplesizes$m3c[chs_c,"N"]<-197
samplesizes$m3c[chs_c,"N_cases"]<-35
samplesizes$m3c[chs_c,"N_controls"]<-162
samplesizes$m3c[chs_c,"prevalence"]<-samplesizes$m3c[chs_c,"N_cases"]/samplesizes$m3c[chs_c,"N"]

#iow<-which(samplesizes$m3a[,"cohort"]=="IOW") 
#samplesizes$m3a[iow,"N_cases"]<-25
#samplesizes$m3a[iow,"N_cases"]<-25
#samplesizes$m3a[iow,"N_controls"]<-437
#samplesizes$m3a[iow,"prevalence"]<-samplesizes$m3a[iow,"N_cases"]/samplesizes$m3a[iow,"N"]


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


###OLD VERSION
#samp_size_old <- load("./03_qc/data/samplesizes.RData")

#inma<-which(samplesizes$m1a[,"cohort"]=="INMA") 
#samplesizes$m1a[inma,"N_cases"]<-146
#samplesizes$m1a[inma,"N_controls"]<-162
#samplesizes$m1a[inma,"N"]<-sum(samplesizes$m1a[inma,"N_cases"],samplesizes$m1a[inma,"N_controls"])
#samplesizes$m1a[inma,"prevalence"]<-100*(samplesizes$m1a[inma,"N_cases"]/samplesizes$m1a[inma,"N"])

#nestw<-which(samplesizes$m1a[,"cohort"]=="NEST_W")
#samplesizes$m1a<-samplesizes$m1a[-(nestw),]

#total<-which(samplesizes$m1a[,"cohort"]=="Total")
#samplesizes$m1a[total,"N"]<-sum(samplesizes$m1a[-total,"N"])
#samplesizes$m1a[total,"N_cases"]<-sum(samplesizes$m1a[-total,"N_cases"])
#samplesizes$m1a[total,"N_controls"]<-sum(samplesizes$m1a[-total,"N_controls"])
#samplesizes$m1a[total,"prevalence"]<-100*(samplesizes$m1a[total,"N_cases"]/samplesizes$m1a[total,"N"])

#samplesizes$m1b<-samplesizes$m1a
#samplesizes$m1c<-samplesizes$m1a
######################################################################################################

#inma<-which(samplesizes$m2a[,"cohort"]=="INMA") 
#samplesizes$m2a[inma,"N_cases"]<-83
#samplesizes$m2a[inma,"N_controls"]<-220
#samplesizes$m2a[inma,"N"]<-sum(samplesizes$m2a[inma,"N_cases"],samplesizes$m2a[inma,"N_controls"])
#samplesizes$m2a[inma,"prevalence"]<-100*(samplesizes$m2a[inma,"N_cases"]/samplesizes$m2a[inma,"N"])

#total<-which(samplesizes$m2a[,"cohort"]=="Total")
#samplesizes$m2a[total,"N"]<-sum(samplesizes$m2a[-total,"N"])
#samplesizes$m2a[total,"N_cases"]<-sum(samplesizes$m2a[-total,"N_cases"])
#samplesizes$m2a[total,"N_controls"]<-sum(samplesizes$m2a[-total,"N_controls"])
#samplesizes$m2a[total,"prevalence"]<-100*(samplesizes$m2a[total,"N_cases"]/samplesizes$m2a[total,"N"])

#samplesizes$m2b<-samplesizes$m2a
#samplesizes$m2c<-samplesizes$m2a

########################################################################################################
#chs<-which(samplesizes$m3c[,"cohort"]=="CHS") 
#samplesizes$m3c[chs,"N_cases"]<-35
#samplesizes$m3c[chs,"N_controls"]<-162
#samplesizes$m3c[chs,"N"]<-sum(samplesizes$m3c[chs,"N_cases"],samplesizes$m3c[chs,"N_controls"])
#samplesizes$m3c[chs,"prevalence"]<-100*(samplesizes$m3c[chs,"N_cases"]/samplesizes$m3c[chs,"N"])

#total<-which(samplesizes$m3c[,"cohort"]=="Total")
#samplesizes$m3c[total,"N"]<-sum(samplesizes$m3c[-total,"N"])
#samplesizes$m3c[total,"N_cases"]<-sum(samplesizes$m3c[-total,"N_cases"])
#samplesizes$m3c[total,"N_controls"]<-sum(samplesizes$m3c[-total,"N_controls"])
#samplesizes$m3c[total,"prevalence"]<-100*(samplesizes$m3c[total,"N_cases"]/samplesizes$m3c[total,"N"])

#samplesizes$m2b<-samplesizes$m2a
#samplesizes$m2c<-samplesizes$m2a

###############################################

#m3a has different numbers for IOW in m3b and m3c.
#new_load function?
#samp_size<-samplesizes
samp_size<-samplesizes_with_totals
models <- names(samp_size)
top_hit_files <- list.files("./05_power/data", full.names=T)

#top_hit_files <- top_hit_files[-grep("info", top_hit_files)]
#godmc_methchar <- load("mean_allcpgs.Robj")
load("/user/home/epzjlm/scratch/repo/godmc_phase2_analysis/07_enrichments/mean_allcpgs.Robj")
godmc_methchar<-df.all
godmc_methchar <- godmc_methchar %>%
	dplyr::select(cpg, sdcpg)


# ------------------------------------------------
# Extract data from meta-analysis
# ------------------------------------------------

#r<-read.table("childhood_model-A.tsv",he=T)
#head(r)
#         CpG       Beta        SE         P     Direction Het_ISq  Het_P
#1 cg03030757  -5.773720 1.2063859 1.702e-06 --+----------     0.0 0.7355
#2 cg03738669  -3.168814 0.6737079 2.557e-06 --+----++-+--    30.9 0.1365
#3 cg00888758  -3.332836 0.7128703 2.936e-06 ----------+--     0.0 0.9077
#4 cg09768093 -18.895025 4.1168986 4.440e-06 --------+--+-     0.0 0.8255
#5 cg03303325   4.331360 0.9831669 1.055e-05 ++++++++++--+     0.0 0.8356
#6 cg09730123  -3.722404 0.8563897 1.382e-05 --------+-++-     0.0 0.4521

# get per model total sample sizes
ma_samp_size <- lapply(samp_size, function(x) {
	y <- x[x$cohort == "Total",]
	df <- tibble(n = y$N, 
				 n_case = y$N_cases, 
				 n_control = y$N_controls)
	return(df)
})

## Get top hits
get_model <- function(res_file)
{
    #stringr::str_extract(res_file, "m[1-3][a-c]")
     
}

convert_units <- function(OR, sd)
{
	beta <- log(OR)
	return(exp(beta * sd))
}

x=top_hit_files[1]
top_hits_res <- lapply(top_hit_files, function(x) {
	read_tsv(x) %>%
		#dplyr::select(CpG = MarkerName, beta = Effect, SE = StdErr, P = Pvalue) %>%
		dplyr::select(CpG = CpG, beta = Beta, SE = SE, P = P) %>%
		arrange(P) %>%
		head(n = 30)
        
})
#names(top_hits_res) <- get_model(top_hit_files)
n.out<-gsub("./05_power/data/","",top_hit_files)
n.out<-gsub(".tsv","",n.out)
n.out<-gsub("model-","",n.out)
n.out<-gsub("_"," ",n.out)
n.out<-gsub("childhood","m1",n.out)
n.out<-gsub("early-onset","m2",n.out)
n.out<-gsub("persistent","m3",n.out)
n.out<-gsub(" A","a",n.out)
n.out<-gsub(" B","b",n.out)
n.out<-gsub(" C","c",n.out)

names(top_hits_res) <- n.out

### setup for report
model_ad_nams <- tibble(ad = c(rep("Childhood AD", 3), 
							   rep("Early-onset AD", 3), 
							   rep("Persistent AD", 3)), 
						model = c("m1a", "m1b", "m1c",
								  "m2a", "m2b", "m2c", 
								  "m3a", "m3b", "m3c"))
study_data_report <- map_dfr(models, function(mod) {
	if (!grepl("c", mod)) return(NULL)
	ad_nam <- model_ad_nams[model_ad_nams$model == mod, "ad", drop = T]
	ss_df <- ma_samp_size[[mod]] %>%
		mutate(trait = ad_nam, K = scales::number(n_case / n,accuracy=0.001))
	thr_mods <- grep(mod, names(top_hits_res))
	thr <- bind_rows(top_hits_res[thr_mods])
	## ADD METH SD TO THE RESULTS
	thr <- left_join(thr, godmc_methchar, by = c("CpG" = "cpg"))
	OR <- ifelse(sign(thr$beta) == 1, exp(thr$beta), 1/exp(thr$beta))
	OR_sd <- convert_units(OR, thr$sdcpg)
	OR_sd <- OR_sd[!is.na(OR_sd)]
	out <- ss_df %>%
		mutate(`OR per SD median (range)` = paste0(round(median(OR_sd), 2), " (", round(min(OR_sd), 2), ", ", round(max(OR_sd), 2), ")"), 
			   `Lowest P` = scales::label_scientific(digits = 2)(min(thr$P))) %>%
		dplyr::select(-n_case, -n_control) %>%
		dplyr::select(trait, everything())
	return(out)	    
})

write.table(study_data_report, file = "./05_power/report/report_data/study_data.tsv", 
		 	col.names = T, row.names = F, quote = F, sep = "\t")

#study_data_report
# A tibble: 3 x 5
#  trait              n K     `OR per SD median (range)` `Lowest P`
#  <chr>          <dbl> <chr> <chr>                      <chr>     
#1 Childhood AD    5110 0.325 1.05 (1.02, 1.14)          4.3e-06   
#2 Early-onset AD  5847 0.263 1.06 (1.01, 1.14)          2.5e-07   
#3 Persistent AD   3189 0.125 1.09 (1.05, 1.14)          1.6e-06   


# ------------------------------------------------
# Convert ORs to difference in meth
# ------------------------------------------------

# convery log OR to OR

hit_res <- lapply(top_hits_res, function(tab) {
	tab %>%
		mutate(OR = exp(beta)) %>%
		dplyr::select(CpG, OR)
})

# function to convert ORs to having only a +ve effect
convert_or <- function(df) ifelse(df$OR < 1, 1/df$OR, df$OR)

min_or <- min(unlist(lapply(hit_res, convert_or)))
max_or <- max(unlist(lapply(hit_res, convert_or)))


# ------------------------------------------------
# run simulations just using logistic regression as normal
# ------------------------------------------------

params <- expand.grid(meth_mean = seq(0.3, 0.7, 0.1),
					  meth_diff = seq(0.01, 0.05, 0.01),  
					  meth_sd = seq(0.05, 0.2, 0.05), 
					  n = seq(signif(min(study_data_report$n), 1), signif(max(study_data_report$n), 1), 500), 
					  case_prop = as.numeric(c(min(study_data_report$K), max(study_data_report$K))), 
					  sim = 1:1000) %>% as_tibble()

## This params below is to check results are consistent with this paper: 
## https://link.springer.com/article/10.1186/s12864-019-5761-7

# params <- expand.grid(meth_mean = seq(0.3, 0.7, 0.1), 
# 					  meth_diff = 0.02,  
# 					  meth_sd = 0.05, 
# 					  n = 1000, 
# 					  case_prop = 0.5, 
# 					  sim = 1:100) %>% as_tibble()

# params <- params[params$sim == 1, ]
x=1

prog_rows <- seq(0, nrow(params), nrow(params)/10)

msg_progress <- function(prog_rows, row) 
{
	if (row %in% prog_rows) {
		row_prog <- which(prog_rows %in% row)
		percent_done <- (row_prog - 1) * 10
		message("Progress = ", percent_done, "%")
	}
}

#' Calculate power to detect an EWAS association 
#'
#' @param meth_mean mean methylation of CpG site
#' @param meth_diff methylation difference between cases and controls
#' @param meth_sd SD of methylation of CpG site
#' @param n sample size
#' @param case_prop proportion of the samples that are cases
#'
#' @return tibble of input parameters, OR and P from `glm()` 
calc_power <- function(meth_mean, meth_diff, meth_sd, n, case_prop)
{
	## Generate cases and controls and methylation levels for each
	n_case <- n * case_prop
	mean_m_case <- meth_mean + meth_diff/2
	mean_m_control <- meth_mean - meth_diff/2
	tab <- tibble(case = c(rep(1, n_case), rep(0, n - n_case)))
	m_case <- rnorm(n_case, mean=mean_m_case, sd=meth_sd)
	m_control <- rnorm(n - n_case, mean=mean_m_control, sd=meth_sd)
	tab$meth <- c(m_case, m_control)
	
	## Run the logistic regression and extract OR and P
	fit <- glm(case ~ meth, data = tab, family = "binomial")
	or <- exp(coef(fit)[2])
	or_sd <- convert_units(or, meth_sd)
	p <- summary(fit)$coefficients[2,4]
	out <- tibble(meth_mean, meth_diff, meth_sd, n, case_prop, or, or_sd, p)
	return(out)
}

start_time <- proc.time()
out_res <- lapply(1:nrow(params), function(x) {
	msg_progress(prog_rows, row = x)
	set.seed(x)
	# print(x)
	df <- params[x, ]
	res <- with(df, calc_power(meth_mean, meth_diff, meth_sd, n, case_prop))
	return(res)
})
time_taken <- proc.time() - start_time
time_taken
## WILL TAKE AROUND 5H TO COMPLETE!!!


fin_res <- bind_rows(out_res)

write.table(fin_res, file = "./05_power/results/power_simulation_data.tsv", 
			col.names = T, row.names = F, quote = F, sep = "\t")


### do downstream analysis
fin_res <- read_tsv("./05_power/results/power_simulation_data.tsv")
fin_res$or_sd <- convert_units(fin_res$or, fin_res$meth_sd)

sig_level <- 1e-7

## differences in mean meth??

power_levels_mm <- fin_res %>%
	group_by(meth_mean, meth_diff, meth_sd, n, case_prop) %>%
	summarise(power = sum(p < sig_level) / max(fin_res$sim), 
			  avg_or_persd = mean(or_sd)) %>%
	mutate(K = comma(case_prop))

meth_mean_power_box <- ggplot(power_levels_mm, aes(x = as.factor(meth_mean), y = power)) + 
	geom_boxplot() + 
	theme_bw()

meth_mean_or_box <- ggplot(power_levels_mm, aes(x = as.factor(meth_mean), y = avg_or)) + 
	geom_boxplot() + 
	theme_bw()

outp <- cowplot::plot_grid(meth_mean_or_box, meth_mean_power_box,
				  labels = c("OR change", "Power change"),
             	  ncol = 1, nrow = 2)

ggsave("./05_power/results/mean_meth_boxplots.pdf", plot = outp)

ggsave("./05_power/results/mean_meth_power_boxplot.pdf", plot = meth_mean_power_box)

## THIS IS ASSUMING NO DIFFERENCES IN POWER ACROSS MEAN METH!!

power_levels <- fin_res %>%
	group_by(meth_diff, meth_sd, n, case_prop) %>%
	summarise(power = sum(p < sig_level) / (max(fin_res$sim) * length(unique(fin_res$meth_mean))), 
			  avg_or_persd = median(or_sd), min_or_persd = min(or_sd), max_or_persd = max(or_sd)) %>%
	mutate(K = comma(case_prop))

power_levels

## ORs look a bit much... but let's see what they look like when converting to OR/SD increase in DNAm
convert_units <- function(OR, sd)
{
	beta <- log(OR)
	return(exp(beta * sd))
}

conv_or <- power_levels %>%
	transmute(power = power, 
			  avg_or = convert_units(avg_or, meth_sd), 
			  min_or = convert_units(min_or, meth_sd), 
			  max_or = convert_units(max_or, meth_sd))

summary(conv_or) # normal looking ORs

## Checking out how OR changes over different categories
power_levels %>% 
	dplyr::filter(meth_diff == 0.01, meth_sd == 0.05, case_prop == 0.12)

power_levels %>% 
	dplyr::filter(meth_diff == 0.01, meth_sd == 0.05, n == 3000)

power_levels %>% 
	dplyr::filter(meth_diff == 0.01, n == 3000, case_prop == 0.12)

# power.t.test(n=300, delta=0.05, sd=0.1,  sig.level=0.05/485000)
# n 
# cases
# meth sd
# OR (or meth_diff)

write.table(power_levels, file = "./05_power/report/report_data/power-levels.tsv", 
			row.names = F, col.names = T, sep = "\t")

power_plots <- ggplot(power_levels, aes(x = n, y = power, colour = K, group = K)) + 
	geom_point() + 
	geom_line() + 
	# labs(legend = "K") + 
	facet_grid(meth_diff ~ meth_sd, 
			   labeller = label_both) + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90))

ggsave("./05_power/results/pwr_plots.pdf", plot = power_plots)

median_or_range <- comma(min(power_levels$avg_or)) 

or_plots <- ggplot(power_levels, aes(x = n, y = avg_or, colour = K, group = K)) + 
	geom_point() + 
	geom_line() + 
	facet_grid(meth_diff ~ meth_sd, 
			   labeller = label_both) + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90))

ggsave("./05_power/results/or_plots.pdf", plot = or_plots)

power.t.test(n=1250, delta=0.01, sd=0.05, sig.level = 1e-7)


