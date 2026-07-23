library(tidyverse)
library(dplyr)
library(purrr)
library(scales)
library(ggforestplot)

top_hit_files <- list.files("./05_power/data", full.names=T)
top_hit_files_mc<-top_hit_files[grep("model-C",top_hit_files)]

n1<-gsub("./05_power/data/","",top_hit_files)
n1<-gsub(".tsv","",n1)
n1<-gsub("model-A","unadj",n1)
n1<-gsub("model-B","maternalfactorsadj",n1)
n1<-gsub("model-C","main",n1)
n2<-n1[grep("main",n1)]

top_hits_list <- Map(function(file, lab) {
  read_tsv(file, show_col_types = FALSE) %>%
    dplyr::select(CpG = CpG, beta = Beta, SE = SE, P = P,Isquared=Het_ISq,Het_P=Het_P) %>%
    arrange(P) %>%
    slice_head(n = 30) %>%
    #mutate(AD = lab, .before = 1)},top_hit_files_mc,n1) 
    mutate(AD = lab, .before = 1)},top_hit_files,n1)
length(top_hits_list)

top_hits_list_mc <- Map(function(file, lab) {
  read_tsv(file, show_col_types = FALSE) %>%
    dplyr::select(CpG = CpG, beta = Beta, SE = SE, P = P,Isquared=Het_ISq,Het_P=Het_P) %>%
    arrange(P) %>%
    slice_head(n = 30) %>%
    mutate(AD = lab, .before = 1)},top_hit_files_mc,n2) 
    #mutate(AD = lab, .before = 1)},top_hit_files,n1)
length(top_hits_list_mc)

# Bind
top_hits_res <- dplyr::bind_rows(top_hits_list)
nrow(top_hits_res)
#270

top_hits_res_mc <- dplyr::bind_rows(top_hits_list_mc)
nrow(top_hits_res_mc)
#90

top_hits_res_split<-top_hits_res %>%
  separate(AD, into = c("ADsubtype", "model"), sep = "_")

top_hits_res_split_mc<-top_hits_res_mc %>%
  separate(AD, into = c("ADsubtype", "model"), sep = "_")

###

top_hits_list <- Map(function(file, lab) {
  read_tsv(file, show_col_types = FALSE) %>%
    dplyr::select(CpG = CpG, beta = Beta, SE = SE, P = P,Isquared=Het_ISq,Het_P=Het_P) %>%
    arrange(P) %>%
    filter(CpG%in%unique(top_hits_res$CpG)) %>%
    mutate(AD = lab, .before = 1)},top_hit_files,n1) 

top_hits_list_mc <- Map(function(file, lab) {
  read_tsv(file, show_col_types = FALSE) %>%
    dplyr::select(CpG = CpG, beta = Beta, SE = SE, P = P,Isquared=Het_ISq,Het_P=Het_P) %>%
    arrange(P) %>%
    filter(CpG%in%unique(top_hits_res_mc$CpG)) %>%
    mutate(AD = lab, .before = 1)},top_hit_files,n1) 


top_hits_res2 <- dplyr::bind_rows(top_hits_list)
top_hits_res2<-top_hits_res2 %>%
  separate(AD, into = c("ADsubtype", "model"), sep = "_")

dim(top_hits_res2)

length(top_hits_res$CpG)
#270
length(unique(top_hits_res$CpG))
#162
length(top_hits_res2$CpG)
#1458
length(unique(top_hits_res2$CpG))
#162

top_hits_res2_mc <- dplyr::bind_rows(top_hits_list_mc)
top_hits_res2_mc<-top_hits_res2_mc %>%
  separate(AD, into = c("ADsubtype", "model"), sep = "_")

dim(top_hits_res2_mc)
#801
length(top_hits_res_mc$CpG)
#90
length(unique(top_hits_res_mc$CpG))
#89
length(top_hits_res2_mc$CpG)
#801
length(unique(top_hits_res2_mc$CpG))
#89


top_hits_res_wide <- top_hits_res2 %>%
  pivot_wider(
    id_cols = CpG,
    names_from = c(ADsubtype, model),
    values_from = c(Isquared, Het_P)
  )
dim(top_hits_res_wide)
#[1] 162  19

n<-names(top_hits_res_wide)
n<-gsub("Het_P_","Het P (",n)
n<-gsub("Isquared_","I-squared (",n)
n<-gsub(".gz",")",n)

top_hits_res_mc$AD<-gsub("_main"," AD",top_hits_res_mc$AD)
top_hits_res_mc$AD<-gsub(".gz","",top_hits_res_mc$AD)
ad_subtypes <- top_hits_res_mc %>%
group_by(CpG) %>%
summarise(
ADsubtype = paste(unique(AD), collapse = "; "),
.groups = "drop"
)

top_hits_res_wide_mc <- top_hits_res2_mc %>%
  pivot_wider(
    id_cols = CpG,
    names_from = c(ADsubtype, model),
    values_from = c(Isquared, Het_P)
) %>%
left_join(ad_subtypes, by="CpG")

dim(top_hits_res_wide_mc)
#[1] 89  19

glimpse(top_hits_res2_mc)

top_hits_res2_mc %>% count(CpG, ADsubtype, model) %>%filter(n > 1)

n_mc<-names(top_hits_res_wide_mc)
n_mc<-gsub("Het_P_","Het P (",n_mc)
n_mc<-gsub("Isquared_","I-squared (",n_mc)
n_mc<-gsub(".gz",")",n_mc)
names(top_hits_res_wide_mc)<-n_mc

write.table(top_hits_res_wide,"TableS7_longformat.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(top_hits_res_wide_mc,"TableS7_longformat_mc.txt",sep="\t",col.names=T,row.names=F,quote=F)

##Get 1e5 hits and forestplot
childhood<-top_hits_res2[which(top_hits_res2$ADsubtype=="childhood"&top_hits_res2$model=="main"&top_hits_res2$P<1e-5),"CpG"]
earlyonset<-top_hits_res2[which(top_hits_res2$ADsubtype=="early-onset"&top_hits_res2$model=="main"&top_hits_res2$P<1e-5),"CpG"]
persistent<-top_hits_res2[which(top_hits_res2$ADsubtype=="persistent"&top_hits_res2$model=="main"&top_hits_res2$P<1e-5),"CpG"]

top_hits_res3<-top_hits_res2[top_hits_res2$CpG%in%unique(c(childhood$CpG,earlyonset$CpG,persistent$CpG)),]

top_hits_res3$CpG<-factor(top_hits_res3$CpG, levels=c(childhood$CpG,earlyonset$CpG,persistent$CpG))
top_hits_res3$ADsubtype<-factor(top_hits_res3$ADsubtype, levels=c("childhood", "early-onset","persistent"))

top_hits_res3$model<-gsub("unadj","unadjusted",top_hits_res3$model)
top_hits_res3$model<-gsub("maternalfactorsadj","maternal factors",top_hits_res3$model)
#top_hits_res3$model<-gsub("main","maternal factors and cell counts",top_hits_res3$model)
#top_hits_res3$model<-factor(top_hits_res3$model, levels=c("unadjusted","maternal factors","maternal factors and cell counts"))
top_hits_res3$model<-factor(top_hits_res3$model, levels=c("unadjusted","maternal factors","main"))

top_hits_res3$group<-NA
w<-which(top_hits_res3$CpG%in%childhood$CpG)
top_hits_res3$group[w]<-"childhood AD"
w<-which(top_hits_res3$CpG%in%earlyonset$CpG)
top_hits_res3$group[w]<-"early-onset AD"
w<-which(top_hits_res3$CpG%in%persistent$CpG)
top_hits_res3$group[w]<-"persistent AD"


top_hits_res3$group <- factor(
  top_hits_res3$group,
  levels = c("childhood AD", "early-onset AD", "persistent AD")   # desired top→bottom
)

top_hits_res3$name<-paste0(top_hits_res3$ADsubtype," (",top_hits_res3$model,")")

top_hits_res3$CpG2<-paste0(top_hits_res3$CpG," (",top_hits_res3$group,")")
top_hits_res3$CpG2<-paste0(top_hits_res3$CpG," (",gsub(" AD","",top_hits_res3$group),")")
#cats<-c("childhood AD unadjusted","childhood AD maternal factors adj","childhood AD maternal factors and cell counts adj", "early-onset AD unadjusted","early-onset AD maternal factors adj","early-onset AD maternal factors and cell counts adj","persistent AD unadjusted","persistent AD maternal factors adj","persistent AD maternal factors and cell counts adj")
#cats<-c("childhood unadjusted","childhood maternal factors","childhood maternal factors and cell counts", "early-onset unadjusted","early-onset maternal factors","early-onset maternal factors and cell counts","persistent unadjusted","persistent maternal factors","persistent maternal factors and cell counts")
cats<-c("childhood (unadjusted)","childhood (maternal factors)","childhood (main)", "early-onset (unadjusted)","early-onset (maternal factors)","early-onset (main)","persistent (unadjusted)","persistent (maternal factors)","persistent (main)")

match(cats,top_hits_res3$name)
top_hits_res3$name<-factor(top_hits_res3$name, levels=cats)
head(top_hits_res3$name)

top_hits_res3[top_hits_res3$CpG%in%c("cg11772801"),"name"]

#model_rank_map<-c("childhood AD unadjusted" = 1,"childhood AD maternal factors adj" = 2 ,"childhood AD maternal factors and cell counts adj" = 3, "early-onset AD unadjusted" = 4,"early-onset AD maternal factors adj" = 5,"early-onset AD maternal factors and cell counts adj" = 6,"persistent AD unadjusted" = 7,"persistent AD maternal factors adj" = 8,"persistent AD maternal factors and cell counts adj" = 9)
#model_rank_map<-c("childhood AD"=1,"early-onset AD"=2,"persistent AD"=3)
#model_rank_map<-c("childhood unadjusted" = 1,"childhood maternal factors" = 2 ,"childhood maternal factors and cell counts" = 3, "early-onset unadjusted" = 4,"early-onset maternal factors" = 5,"early-onset maternal factors and cell counts" = 6,"persistent unadjusted" = 7,"persistent maternal factors" = 8,"persistent maternal factors and cell counts" = 9)
model_rank_map<-c("childhood (unadjusted)" = 1,"childhood (maternal factors)" = 2 ,"childhood (main)" = 3, "early-onset (unadjusted)" = 4,"early-onset (maternal factors)" = 5,"early-onset (main)" = 6,"persistent (unadjusted)" = 7,"persistent (maternal factors)" = 8,"persistent (main)" = 9)

#
top_hits_res3 <- top_hits_res3 %>%
  mutate(
    model_rank = (model_rank_map[as.character(name)])
  )

#top_hits_res3 <- top_hits_res3 %>%
#  group_by(CpG) %>%
#  mutate(name = fct_reorder(name, model_rank,.desc = TRUE)) %>%  # .desc = TRUE for reverse
#  ungroup()

top_hits_res3[top_hits_res3$CpG%in%c("cg11772801"),"model_rank"]
top_hits_res3[top_hits_res3$CpG%in%c("cg11772801"),"name"]
facet_lab<-as.character(top_hits_res3$CpG2)
names(facet_lab)<-top_hits_res3$CpG 

p<-forestplot(
  df = top_hits_res3,
  estimate = beta,
  pvalue = P,
  se=SE,
  xlab = "beta (95% CI)",
  colour = ADsubtype,
  psignif=1,
  shape=model,
  name=name
) +
facet_wrap(~ CpG, ncol = 3, scales = "free_x",labeller=labeller(CpG=as_labeller(facet_lab))) +

scale_y_discrete(limits = rev(levels(top_hits_res3$name))) +
#guides(
#    colour = guide_legend(order = 1),
#    shape  = guide_legend(order = 2)
# ) +

#theme(legend.position = c(1, -0.008),
#      legend.justification = c(1, 0),
#      legend.text          = element_text(size = 10),
#      legend.title         = element_text(size = 12)
#        )

theme(
    legend.position = c(-0.06,-0.07),
    legend.justification = c(0, 0),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  ) +
 
guides(
    colour = guide_legend(order=1),
    shape  = guide_legend(order=2)
  ) +


theme(
  plot.margin = margin(10, 10, 20, 10)   # big bottom margin
)

ggsave(p,filename="forestplot2.pdf",height=12,width=10)

