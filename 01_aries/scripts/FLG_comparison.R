# FLG comparison
# ------------------------------------------------------------------------
library(haven)
library(tidyverse)
library(dplyr)
library(janitor)

## data
phen_file <- "./data/childhood-ad-doctor-only.tsv"
new_phen_file <- "./data/childhood-ad-doctor-only-NEW.tsv"
phen_dat <- read_tsv(phen_file)
new_phen_dat <- read_tsv(new_phen_file)

## ad = all are controls if not cases
## ad2 = controls = when "no" answered to eczema questions, cases = when "yes and doctor diagnosed"

ecz_81m <- replace(new_phen_dat$kq035, new_phen_dat$kq035 == 3, 0)
ecz_81m <- ifelse(!ecz_81m %in% c(1, 0), NA, ecz_81m)

ecz_91m <- replace(new_phen_dat$kr042, new_phen_dat$kr042 == 3, 0)
ecz_91m <- ifelse(!ecz_91m %in% c(1, 0), NA, ecz_91m)

ecz_103m <- replace(new_phen_dat$ks1042, new_phen_dat$ks1042 == 3, 0)
ecz_103m <- ifelse(!ecz_103m %in% c(1, 0), NA, ecz_103m)

ecz_128m <- replace(new_phen_dat$kv1060, new_phen_dat$kv1060 == 3, 0)
ecz_128m <- ifelse(!ecz_128m %in% c(1, 0), NA, ecz_128m)

ecz_diag <- replace(new_phen_dat$kv1070, new_phen_dat$kv1070 %in% c(1, 4), 0)
ecz_diag <- replace(ecz_diag, ecz_diag %in% c(2, 3), 1)
ecz_diag <- ifelse(!ecz_diag %in% c(1, 0), NA, ecz_diag)

ecz_all <- rowSums(cbind(ecz_81m, ecz_91m, ecz_103m, ecz_128m, ecz_diag), na.rm=TRUE)
ecz_all_na <- rowSums(cbind(is.na(ecz_81m), is.na(ecz_91m), is.na(ecz_103m), is.na(ecz_128m), is.na(ecz_diag)), na.rm=TRUE)
ecz_all[ecz_all_na == 5] <- NA
new_phen_dat$childhood_ad2 <- replace(ecz_all, ecz_all > 0, 1)


old_phen_dat_all <- read_dta("./data/dd_PHENO_CLEAN.dta")
flg_vars <- read_dta("./data/children_FLG_variables.dta")
flg_vars <- flg_vars[, c("aln", "qlet", "FLG_comb")]

comb_dat <- old_phen_dat_all %>%
	dplyr::select(aln, qlet, childhood_AD, childhood_AD_dd) %>%
	left_join(new_phen_dat[, c("aln", "qlet", "childhood_ad")]) %>%
	left_join(flg_vars)

ad_vars <- grep("childhood", colnames(comb_dat), value = T)

#' Extract res from glm() function
#' 
#' @param glm_obj object obtained from running the glm() function
#' @return table containing summary stats
summ_glm_res <- function(glm_obj)
{
	summ_res <- summary(glm_obj)
	out <- tibble(Beta = summ_res$coef[2, 1], 
				  SE = summ_res$coef[2, 2], 
				  P = summ_res$coef[2, 4],
				  OR=exp(summ_res$coef[2, 1]),
				 ci_low =exp(confint(glm_obj)[2,1]),
				 ci_high =exp(confint(glm_obj)[2,2]))

	return(out)
}


assoc_res <- lapply(ad_vars, function(ad) {
	cat(ad,"\n")
	test_dat <- comb_dat %>%
		dplyr::select(aln, qlet, one_of(ad), FLG_comb) %>%
		na.omit()
	form <- as.formula(paste0(ad, " ~ ", "FLG_comb"))
	tb<-tabyl(test_dat,FLG_comb,.data[[ad]])
	glm_res <- glm(form, data = test_dat, family = "binomial")
	summ_stats <- summ_glm_res(glm_res)
	out <- summ_stats %>%
		mutate(N = nrow(test_dat), N_cases = sum(test_dat[[ad]]), N_controls = N - N_cases)
	return(list(summary=out, table=tb))
})
names(assoc_res) <- ad_vars

## Removing own childhood AD definition
assoc_res <- assoc_res[names(assoc_res) != "childhood_ad"]

out_res <- bind_rows(assoc_res, .id = "definition") %>%
	mutate(definition = ifelse(definition == "childhood_AD", "doc diagnosis or rash", "doc diagnosis only"))

write.table(out_res, file = "flg-eczema-case-assoc_revision.tsv", sep="\t", quote=F, row.names=F, col.names = T)

comb_dat2 <- old_phen_dat %>%
	dplyr::select(aln, qlet, childhood_AD, childhood_AD_dd, childhood_AD_dd2) %>%
	left_join(new_phen_dat[, c("aln", "qlet", "childhood_ad")]) %>%
	left_join(flg_vars)

summary(flg_vars)

