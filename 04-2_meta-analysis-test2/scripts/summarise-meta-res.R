# ------------------------------------------------------------
# Summarise meta-analysis results
# ------------------------------------------------------------

## Aim: To take meta-analysis results and summarise them in a succint table

## Date: 2022-03-23

# srun --job-name "InteractiveJob" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=5GB --pty bash

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # QQs etc.
library(cowplot) # organising plots nicely
library(RColorBrewer) # Colour in the plots
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
meta_files <- args[1]
old_file <- args[2]
plots_outfile <- args[3]
summary_outfile <- args[4]
comp_p_outfile <- args[5]

# meta_files <- "results/metal-res/m1c-hr.txt"
# old_file <- "../04_meta-analysis/results/metal-res/m1c.txt"
# plots_outfile <- "results/meta-summary/man-qqs-m1c-hr.png"
# summary_outfile <- "results/meta-summary/comparison-summary.RData"
# comp_p_outfile <- "results/meta-summary/hr-main-meta-comp-plot.png"

meta_files <- unlist(str_split(meta_files, " "))

# ------------------------------------------------------------
# Summarise the results to fit the paper
# ------------------------------------------------------------

## general data functions
get_model <- function(res_file)
{
    gsub(".txt", "", basename(res_file))
    # stringr::str_extract(res_file, "m[1-3][a-c]-?.*")
}

read_meta_file <- function(res_file)
{
	read_tsv(res_file) %>%
		dplyr::select(name = MarkerName, beta = Effect, SE = StdErr, P = Pvalue, Isq = HetISq, het_p = HetPVal)
}

## Get inflation stats
get_lambda <- function(res_file) {
	res <- read_meta_file(res_file) %>%
        dplyr::select(name, beta, SE, P)
    lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    # get top hit as well
    out <- res %>%
        arrange(P) %>%
        head(n = 1) %>%
        mutate(lambda = lamb)
	return(out)
}

lambda_list <- lapply(meta_files, get_lambda)
names(lambda_list) <- get_model(meta_files)
lambda_tib <- bind_rows(lambda_list, .id = "model")
# lambda_tib <- tibble(model = names(lambda_list), lambda = unlist(lambda_list))

## Get heterogeneity stats for top hits
get_het_stats <- function(res_file)
{
	res <- read_meta_file(res_file) %>%
		arrange(P) %>%
		head(n = 30)
	out <- tibble(model = get_model(res_file), name = res$name, Isq = res$Isq, het_p = res$het_p)
	return(out)
}

het_stats <- map_dfr(meta_files, get_het_stats)
summary(het_stats)

# ## Write it out
# write.table(lambda_tib, file = lambda_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
# write.table(het_stats, file = het_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## QQ plots
make_qq <- function(res_file)
{
	res <- read_meta_file(res_file)
	lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", comma(lamb))
	ewaff_qq <- ewaff.qq.plot(res$P, lambda.method = "none") + 
		theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
		labs(title = get_model(res_file)) + 
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
		theme(text = element_text(size = 8), legend.position = "none", plot.title = element_blank())
}

plot_qqs <- function(qqlist)
{
    m_qqs <- lapply(qqlist, function(x) {x + theme(legend.position = "none")})
    plots <- cowplot::plot_grid(plotlist = m_qqs, nrow=3)
    return(plots)
}

qq_plots <- lapply(meta_files, make_qq)
names(qq_plots) <- sapply(meta_files, get_model)
# ggsave(qq_outfile, plot_qqs(qq_plots))

## Manhattan plots
annotation <- meffil::meffil.get.features("450k")
annotation <- annotation %>% 
    mutate(chr = gsub("chr", "", chromosome)) %>%
    mutate(chr = gsub("X", "23", chr)) %>% 
    mutate(chr = as.numeric(gsub("Y", "24", chr)))

make_man <- function(res_file, cpg_annotations)
{
    res <- read_meta_file(res_file) %>%
        left_join(cpg_annotations)
    # to highlight
    cpg_h <- res[res$P < 1e-7, ]$name
    gg_man <- gg.manhattan(df = res, 
                           hlight = cpg_h, 
                           title = NULL, 
                           SNP = "name", 
                           CHR = "chr", 
                           BP = "position", 
                           P = "P", 
                           sig = 1e-7, 
                           sugg = 1e-5, 
                           lab = TRUE, 
                           colour = TRUE)
    gg_man <- gg_man + 
        theme(plot.title = element_blank(), text = element_text(size = 10), axis.text.x = element_text(angle = 90, size = 8))
    return(gg_man)
}

plot_mans <- function(pheno_mod, manlist) 
{
    m_man <- manlist[grep(pheno_mod, names(manlist))]
    m_man <- lapply(m_man, function(x) {x + theme(axis.title.x = element_blank(), title = element_blank())})
    plots <- cowplot::plot_grid(plotlist = m_man, labels = names(m_man), nrow = 3)
    return(plots)
}

mans <- lapply(meta_files, make_man, annotation)
names(mans) <- sapply(meta_files, get_model)
models <- sapply(meta_files, get_model)

# ggsave(m1_man_outfile, plot_mans("m1", mans))
# ggsave(m2_man_outfile, plot_mans("m2", mans))
# ggsave(m3_man_outfile, plot_mans("m3", mans))



all_plots <- cowplot::plot_grid(qq_plots[[1]], mans[[1]], 
                                nrow = 2, ncol = 1, rel_heights = c(1, 2), 
                                labels = c("Hanifin-Rajka"))

ggsave(plots_outfile, plot = all_plots)

## Correlation plot
## Get effect estimates
meta_files <- c(meta_files, old_file)
cpg_lists <- lapply(meta_files, function(x) {
    res <- data.table::fread(x)
    return(res$MarkerName)
})

common_cpgs <- intersect(cpg_lists[[1]], cpg_lists[[2]])

betas <- lapply(meta_files, function(x) {
    res <- data.table::fread(x) %>%
        dplyr::filter(MarkerName %in% common_cpgs) %>%
        arrange(MarkerName) %>%
        dplyr::select(CpG = MarkerName, Beta = Effect)
    return(res)
})
names(betas) <- sapply(meta_files, get_model)

models <- names(betas)
## correlation
beta_cors <- cor(betas[[models[1]]]$Beta, betas[[models[2]]]$Beta)

old_top_hits <- read_meta_file(old_file) %>%
    arrange(P) %>%
    head(n = 30)
## top hits correlation
old_betas <- lapply(betas, function(bee) {
    bee %>%
        dplyr::filter(CpG %in% old_top_hits$name)
})

tophit_beta_cors <- cor(old_betas[[models[1]]]$Beta, old_betas[[models[2]]]$Beta)

## top hits replication
hr_res <- read_meta_file(meta_files[grep("hr.txt", meta_files)]) %>%
    dplyr::filter(name %in% old_top_hits$name) %>%
    mutate(replicated = P < 0.05/30)

summary(hr_res)

# ------------------------------------------------------------
# Write out results summary
# ------------------------------------------------------------

summ_out <- list(correlation = list(all = beta_cors, top_dd = dd_beta_cors, top_sr = sr_beta_cors), 
                 replication = list(top_dd = sr_dd, top_sr = dd_sr))

save(summ_out, file = summary_outfile)

## What we want in report

# Short intro on what analyses were doing

# Table 1: Cohorts in each meta-analysis

# Doctor-diagnosed/"good" definition
# | Cohort | Cases | Controls | N | Definition
#  ...
# Self-report/"poor" definition
# | Cohort | Cases | Controls | N | Definition 

# Figure 1: QQs and Manhattans

# Table 2: Replication and correlation

# | Metric                  | Value |
# | Correlation of all      |       |
# | Correlation of top DD   |       |
# | Correlation of top SR   |       |
# | Replication of DD (0.05)|       |
# | Replication of SR (0.05)|       |




### FOREST PLOT

## Josine's instructions
old_top_hits_p4 <- read_meta_file(old_file) %>%
    dplyr::filter(P < 1e-4)

hr_res <- read_meta_file(meta_files[grep("hr.txt", meta_files)]) %>%
    dplyr::filter(name %in% old_top_hits_p4$name) %>%
    mutate(replicated = P < 0.05/30)

## need columns: # name Group Estimate 2.5 % 97.5 % 
forest_res <- bind_rows(hr_res, old_top_hits_p4) %>%
    mutate(`2.5 %` = beta - (1.96 * SE), 
           `97.5 %` = beta + (1.96 * SE), 
           studies = ifelse(is.na(replicated), "all", "HR")) %>%
    dplyr::select(name, studies, Estimate = beta, `2.5 %`, `97.5 %`) %>%
    mutate(studies = as.factor(studies)) %>%
    arrange(studies, name)

forest_res$facet_var <- facet_var_gen(forest_res, 2, "studies")

outplot1 <- forest_plot(forest_res, 
                       col_num = 2, 
                       group = "studies", 
                       y_axis = "name", 
                       null_at = 0)
ggsave(comp_p_outfile, plot = outplot1)

hr_res2 <- read_meta_file(meta_files[grep("hr.txt", meta_files)]) %>%
    dplyr::filter(P < 1e-4) %>%
    mutate(studies = "HR")

old_hits <- read_meta_file(old_file) %>%
    dplyr::filter(name %in% hr_res2$name) %>%
    mutate(studies = "all")

## need columns: # name Group Estimate 2.5 % 97.5 % 
forest_res <- bind_rows(hr_res2, old_hits) %>%
    mutate(`2.5 %` = beta - (1.96 * SE), 
           `97.5 %` = beta + (1.96 * SE)) %>%
    dplyr::select(name, studies, Estimate = beta, `2.5 %`, `97.5 %`) %>%
    mutate(studies = as.factor(studies)) %>%
    arrange(studies, name)

forest_res$facet_var <- facet_var_gen(forest_res, 2, "studies")

outplot2 <- forest_plot(forest_res, 
                        col_num = 2, 
                        group = "studies", 
                        y_axis = "name", 
                        null_at = 0)

ggsave("results/meta-summary/hr-hits-forest.png", plot = outplot2)


517 (93, 424)
153 (69, 84)
605 (177, 428)

605+153+517
93 + 69 + 177

## Checking for replication from IOW previous study
supp_tab1_file <- "/pace_ad/09_replication/data/pmid-26199674_supp-table-1.xlsx"
supp_tab1 <- readxl::read_xlsx(supp_tab1_file)
all_cpgs <- supp_tab1$CpGs
sig_cpgs <- supp_tab1[supp_tab1[["P-Value"]] < 0.05, "CpGs", drop=T]
hr_res3 <- read_meta_file(meta_files[grep("hr.txt", meta_files)]) %>%
    dplyr::filter(name %in% all_cpgs) %>%
    mutate(sig = ifelse(name %in% sig_cpgs, TRUE, FALSE)) %>%
    arrange(P)


hr_res3$fdr <- p.adjust(hr_res3$P, method = "fdr")


## Checking for FLG associations 
genes_file <- "/pace_ad/08_candidate-genes/data/candidate-genes-cpgs.tsv"
gene_dat <- read_tsv(genes_file)

## just take first ranked genes as that's ~27k sites
gene_dat <- gene_dat %>%
    dplyr::filter(rank == "first" | gene == "FLG")

flg_dat <- gene_dat %>%
    dplyr::filter(gene %in% "FLG")
hr_res2[hr_res2$name %in% flg_dat$name,]

