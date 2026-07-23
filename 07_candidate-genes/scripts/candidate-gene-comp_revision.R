#module load languages/R/4.2.3

library(tidyverse) # tidy code and data
#library(meffil) # to get DNAm positions
library(readxl) # read in the candidate gene data
library(biomaRt) # ensembl gene IDs and positions
library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
library(org.Hs.eg.db)
library(AnnotationDbi)
# genes_file <- "data/candidate-genes-cpgs.tsv"
# meta_file <- "../04_meta-analysis/results/metal-res/m1c.txt"
# summary_out <- "results/m1c-comp-summary.RData"
# qqplot_out <- "results/plots/m1c-qq-comparison.png"

raw_datfile <- "../data/candidate-genes-raw.xlsx"
raw_dat <- read_xlsx(raw_datfile, n_max = 26)

genes_file<-"../data/GWAS_genes.txt"
meta_file <- "../../04_meta-analysis/results/childhood_model-C.tsv"
## data
raw_dat <- read_tsv(genes_file,na ="N/A")
meta_res <- read_tsv(meta_file)

colnames(raw_dat)
colnames(raw_dat) <- gsub(" |-", "_", colnames(raw_dat))

## separate out gene names and prioritisation scores
clean_dat <- raw_dat %>%
	separate(Top_Ranked_Gene, "first", "\\(") %>%
	separate(Second_Ranked_Gene, "second", "\\(") %>%
	separate(Third_Ranked_Gene, "third", "\\(") %>%
	dplyr::select(-Locus, -Variant, -`Chr:position`) %>% 
	pivot_longer(everything(), names_to = "rank", values_to = "gene") %>%
	mutate(gene = substr(gene ,1, nchar(gene) - 1))

#head(clean_dat)
# A tibble: 6 x 2
#  rank   gene  
#  <chr>  <chr> 
#1 first  HRNR  
#2 second RPTN  
#3 third  CRNN  
#4 first  IL-6R 
#5 second UBE2Q1
#6 third  ADAR 

clean_dat$gene<-gsub(" /.*", "", clean_dat$gene)

gene_dat <- clean_dat %>%
	dplyr::filter(rank == "first" | gene == "FLG")

cpg.pval<-meta_res$P
names(cpg.pval)<-meta_res$CpG

genelist1<-list(unique(gene_dat$gene))
names(genelist1) <- c("GS1")

entrez <- mapIds(
  org.Hs.eg.db,
  keys = genelist1$GS1,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"   # or "CharacterList"
)

entrez

length(genelist1$GS1)
#[1] 70
length(entrez)
#[1] 70
which(is.na(entrez))
#KIAA1109 
#      12 
w<-which(names(entrez)%in%"KIAA1109")
entrez[w]<-"84162"

sig.cpg<-meta_res[meta_res$P<0.05,"CpG"]

res <- gsameth(
  sig.cpg = sig.cpg$CpG,        # your significant CpGs
  all.cpg = meta_res$CpG,        # background CpGs
  collection = entrez,     # ← your custom gene sets
  array.type = "450K",      # or "EPIC"
  plot.bias = TRUE, 
  prior.prob = TRUE
)

topGSA(res)

res <- gsameth(
  sig.cpg = sig.cpg$CpG,        # your significant CpGs
  all.cpg = meta_res$CpG,        # background CpGs
  collection = entrez,     # ← your custom gene sets
  array.type = "450K",      # or "EPIC"
  plot.bias = TRUE, 
  prior.prob = TRUE,
  sig.genes = TRUE
)
topGSA(res)



#res
#            N DE      P.DE       FDR
#collection 70 35 0.5080048 0.5080048

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sigCpGs<-meta_res[meta_res$P<0.05,]
gst <- gometh(sig.cpg=sigCpGs, all.cpg=meta_res$P, collection="GO",
              plot.bias=TRUE, array.type=T)
topGSA(gst, n=10)









res1 = methylglm(cpg.pval = cpg.pval,array.type="450K",group="all",
                 maxsize = 500, GS.idtype="SYMBOL",GS.list=genelist1)
head(res1, 15)

data("cpg2gene450k")


res1 = methylgometh(cpg.pval = cpg.pval,array.type="450K",group="all",sig.cut="0.001")

FullAnnot = prepareAnnot(CpG2Gene)
res3<-methylRRA(cpg.pval = cpg.pval, FullAnnot=FullAnnot, method = "ORA", 
                    GS.list = genelist1, GS.idtype = "SYMBOL", 
                    minsize = 100, maxsize = 300)

# ------------------------------------------------------------
# Extract CpGs surrounding the genes
# ------------------------------------------------------------
 
mef_probe_info <- meffil.get.features()
mef_probe_info <- mef_probe_info %>%
	dplyr::filter(grepl("^cg[0-9]*", name)) %>%
	mutate(chromosome = gsub("^chr", "", chromosome)) %>%
	dplyr::select(name, chromosome, position) %>%
	as_tibble

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attri <- listAttributes(mart)
chr <- c(1:22, "x", "y")
all_genes <- lapply(chr, function(chr_num) {
	print(chr_num)
	genes <- getBM(
	    attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
	    #filters=c("chromosome_name", "start", "end"),
	    filters="chromosome_name",
	    #values=list(chromosome_name=chr_num, start="1",end=as.character(3e8)),
    	values=chr_num,
    	mart=mart
    )
    genes$chromosome_name <- as.character(genes$chromosome_name)
	return(genes)   
})
library(tidyverse) # couldn't load this before as it crashed things! 
all_genes <- bind_rows(all_genes) %>%
	distinct()

clean_dat %>%
	left_join(all_genes, by = c("gene" = "hgnc_symbol")) %>%
	dplyr::select(rank, gene, chromosome = chromosome_name, start_position, end_position) %>%
	dplyr::filter(is.na(chromosome))

## make some manual edits to clean dat to get it in line with the 
new_KIAA0391 <- "PRORP" # new hgnc name for one gene
clean_dat$gene <- ifelse(clean_dat$gene == "KIAA0391", new_KIAA0391, clean_dat$gene)

# need to remove the dash for IL-XXXX 
clean_dat_no_il <- clean_dat[!grepl("IL-", clean_dat$gene),]
clean_dat_il <- clean_dat[grepl("IL-", clean_dat$gene),]
clean_dat_il$gene <- gsub("-", "", clean_dat_il$gene)
clean_dat <- bind_rows(clean_dat_no_il, clean_dat_il)

## Add in FLG to clean_dat as not included
flg_dat <- tibble(rank = NA, gene = "FLG")
clean_dat <- bind_rows(clean_dat, flg_dat)

## Combine datasets to get positions for genes
gene_dat <- clean_dat %>%
	left_join(all_genes, by = c("gene" = "hgnc_symbol")) %>%
	dplyr::select(rank, gene, chromosome = chromosome_name, start_position, end_position)

## Write out cleaned gene data
write.table(gene_dat, file = clean_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## Think about how to actually extract the data
gene_dat <- gene_dat %>%
	mutate(new_start = start_position - 1e6, new_end = end_position + 1e6)

x=1
mapped_cpgs <- lapply(1:nrow(gene_dat), function(x) {
	gd <- gene_dat[x, ]
	cpgs_within_gene <- between(mef_probe_info$position, gd$new_start, gd$new_end)
	if (sum(cpgs_within_gene) == 0) return(NULL)
	cpg_dat <- mef_probe_info %>%
		dplyr::filter(cpgs_within_gene, chromosome == gd$chromosome) %>%
		mutate(gene = gd$gene, rank = gd$rank) %>%
		distinct()
	return(cpg_dat)
})
mapped_cpgs <- bind_rows(mapped_cpgs)

table(mapped_cpgs$rank)
 # first second  third
 # 27245  27448  26844

## There are 27K near the first ranked gene so I think keeping to that would be fine

## Write out the cpgs, genes and positions
write.table(mapped_cpgs, file = cpg_genes_output, col.names = T, row.names = F, quote = F, sep = "\t")



