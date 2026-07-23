
# ------------------------------------------------------------
# Summarise meta-analysis results
# ------------------------------------------------------------

## Aim: To take meta-analysis results and summarise them in a succint table

## Date: 2022-03-23

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # QQs etc.
library(cowplot) # organising plots nicely
library(RColorBrewer) # Colour in the plots
#library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
meta_files <- args[1]
plots_outfile <- args[2]
summary_outfile <- args[3]

plots_outfile <- "results/meta-summary/man-qqs_m1.tiff"
summary_outfile <- "results/meta-summary/comparison-summary_m1.RData"
qq_outfile <- "/results/meta-summary/meta-qqs_m1.png"

meta_files <- "results/metal-res/m1a.txt results/metal-res/m1b.txt results/metal-res/m1c.txt"
meta_files <- unlist(str_split(meta_files, " "))

m1_man_outfile <- "results/meta-summary/meta-m1a-manhattan.png"
m2_man_outfile <- "results/meta-summary/meta-m1b-manhattan.png"
m3_man_outfile <- "results/meta-summary/meta-m1c-manhattan.png"

n<-gsub("results/metal-res/","",meta_files)
n<-gsub("m1a","Unadjusted",n)
n<-gsub("m1b","Maternal factors adj",n)
n<-gsub("m1c","Main",n)
n<-gsub(".txt","",n)
nlab<-n




# ------------------------------------------------------------
# Summarise the results to fit the paper
# ------------------------------------------------------------

## general data functions
get_model <- function(res_file)
{
    stringr::str_extract(res_file, "m[1-3][a-c]-..")
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
    
    cat(res_file,"\n")
    n<-gsub("results/metal-res/","",res_file)
    n<-gsub("m1a","Unadjusted",n)
    n<-gsub("m1b","Maternal factors adj",n)
    n<-gsub("m1c","Main",n)
    n<-gsub(".txt","",n)

    lamb <- median(qchisq(res$P, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
    lamb2 <- paste("lambda == ", sprintf("%.2f", lamb))
    ewaff_qq <- ewaff.qq.plot(res$P, lambda.method = "none", 
                              xlab = bquote(-log[10]("expected P")), 
                              ylab = bquote(-log[10]("observed P"))) + 
        theme_bw() + 
        annotate("text", x = -Inf, y = Inf, label = lamb2, hjust = 0, vjust = 1, parse = TRUE) + 
        labs(title = n) +
        #theme(plot.title = element_markdown()) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
        theme(text = element_text(size = 8)) + 
        theme(legend.position = "none")
}

plot_qqs <- function(qqlist)
{
    m_qqs <- lapply(qqlist, function(x) {x})
    plots <- cowplot::plot_grid(plotlist = m_qqs, nrow=3)
    return(plots)
}

qq_plots <- lapply(meta_files, make_qq)
names(qq_plots) <- sapply(meta_files, get_model)
ggsave(qq_outfile, plot_qqs(qq_plots))

## Manhattan plots
annotation <- meffil::meffil.get.features("450k")
annotation <- annotation %>% 
    mutate(chr = gsub("chr", "", chromosome)) %>%
    mutate(chr = gsub("X", "23", chr)) %>% 
    mutate(chr = as.numeric(gsub("Y", "24", chr)))

scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
    x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
    y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
    z.cell <- x.cell * resolution + y.cell
    frequency.table <- table(z.cell)
    frequency <- rep(0,max(z.cell, na.rm=T))
    frequency[as.integer(names(frequency.table))] <- frequency.table
    f.cell <- frequency[z.cell]
    
    big.cells <- length(which(frequency > max.per.cell))
    sort(c(which(f.cell <= max.per.cell),
           sample(which(f.cell > max.per.cell),
                  size=big.cells * max.per.cell, replace=F)),
         decreasing=F)
}

gg.manhattan <- function(df, hlight, col = "default",
                         title = "Manhattan", SNP = "SNP", CHR = "CHR", BP = "BP", P = "P",
                         sig = 5e-8, sugg = 1e-5, lab = FALSE, colour = TRUE, 
                         remove_chr_labs = NULL){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    dplyr::group_by(!! as.name(CHR)) %>% 
    dplyr::summarise(chr_len = as.numeric(max(!! as.name(BP), na.rm = TRUE))) %>% 
    
    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    dplyr::left_join(df, ., by=setNames(CHR, CHR)) %>%
    
    # Add a cumulative position of each SNP
    dplyr::arrange({{ CHR }}, {{ BP }}) %>%
    dplyr::mutate( BPcum = !! as.name(BP) + tot) %>%
    
    # Add highlight and annotation information
    dplyr::mutate( is_highlight := ifelse(!! as.name(SNP) %in% hlight, "yes", "no")) %>%
    dplyr::mutate( is_annotate := ifelse(!! as.name(SNP) %in% hlight, "yes", "no"))

  # change CHR to a factor
  df.tmp[[CHR]] <- as.factor(df.tmp[[CHR]])

  # thin the scatter plot
  selection.idx <- scatter.thinning(x = df.tmp$BPcum, y = -log10(df.tmp[[P]]), resolution=100, max.per.cell=100)
  df.select <- df.tmp[selection.idx, ]
  
  df.select$stat <- -log10(df.select[[P]])

  # sort the colour out
  if (col == "default") {
    col <- RColorBrewer::brewer.pal(9, "Greys")[c(5,9)] 
  }
  
  # for the colour later on
  chr_n <- length(unique(df.select[[CHR]]))

  # get chromosome center positions for x-axis
  axisdf <- df.select %>% 
    dplyr::group_by_at(CHR) %>% 
    dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    if (!is.null(remove_chr_labs)) {
        axisdf[[CHR]] <- ifelse(axisdf[[CHR]] %in% remove_chr_labs, "", axisdf[[CHR]])
    }

  p <- ggplot2::ggplot(df.select, aes(x = BPcum, y = stat)) +
    # Show all points
    ggplot2::geom_point(aes_string(color = CHR), alpha=0.8, size=1) +
    ggplot2::scale_color_manual(values = rep_len(col, length.out = chr_n)) +

    # custom axes:
    ggplot2::scale_x_continuous( label = axisdf[[CHR]], breaks= axisdf$center ) +
    ggplot2::scale_y_continuous(expand = c(0, 1)) + # expand=c(0,1)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggplot2::ggtitle(paste0(title)) +
    ggplot2::labs(x = "Chromosome", y = "-log10(P)") +
    
    # add genome-wide sig and sugg lines
    ggplot2::geom_hline(yintercept = -log10(sig)) +
    ggplot2::geom_hline(yintercept = -log10(sugg), linetype="dashed")# +
    
    if (colour) {
        p <- p + 
            ggplot2::geom_point(data=subset(df.select, is_highlight=="yes"), color="orange", size=2)
    } else {
        p <- p +
            ggplot2::geom_point(data=subset(df.select, is_highlight=="yes"), shape=2, size=2)
    }
    # Add highlighted points
    # geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
    p <- p +    
    # Custom the theme:
    ggplot2::theme_bw(base_size = 22) +
    ggplot2::theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
   if (lab && length(hlight) > 0) p <- p + ggrepel::geom_label_repel(data=df.select[df.select$is_annotate=="yes",],
                                      ggplot2::aes_string(label=SNP, alpha=0.7), size=5, force=1.3)
   return(p)
}

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
    m_man <- lapply(m_man, function(x) {x + theme(title = element_blank())})
    plots <- cowplot::plot_grid(plotlist = m_man, labels = names(m_man), nrow = 3)
    return(plots)
}

mans <- lapply(meta_files, make_man, annotation)
names(mans) <- nlab

#models <- sapply(meta_files, get_model)
models<-nlab
m1a_model <- which(models%in%c("Unadjusted"))
m1b_model <- grep("Maternal", models, value = T)
m1c_model <- grep("Main", models, value = T)
names(qq_plots)<-nlab
theme_set(theme_bw(base_size = 12))   # or theme_minimal(), theme_cowplot()
all_plots <- cowplot::plot_grid(qq_plots[[m1a_model]], mans[[m1a_model]], 
                                qq_plots[[m1b_model]], mans[[m1b_model]],
                                qq_plots[[m1c_model]], mans[[m1c_model]],
                                nrow = 3, ncol = 2, rel_widths = c(1, 2), 
                                
                                label_x = 0,
                                label_y = 1,
                                hjust = 0, vjust = 1
                                )

ggsave(plots_outfile, plot = all_plots,dpi=300,device="tiff",bg="white")









