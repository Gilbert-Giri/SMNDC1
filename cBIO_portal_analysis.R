
# SMNDC1 analyses
# Part of ERK1 splice switching reprograms kinase output to drive resistance to KRAS inhibitors 
# code by Vincent L. Cannataro, Ph.D. 


# All paths are relative to the working directory

## load in packages ---- 

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(ggtext)
library(scales)

# cbioportal uses hg19 https://docs.cbioportal.org/user-guide/faq/

# options(timeout = 600)
# remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release")
library(ces.refset.hg19)



## download data  ---- 

# mRNA vs. copy number 
# this analysis: https://www.cbioportal.org/results/plots?cancer_study_list=laml_tcga_pan_can_atlas_2018%2Cacc_tcga_pan_can_atlas_2018%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cchol_tcga_pan_can_atlas_2018%2Ccoadread_tcga_pan_can_atlas_2018%2Cdlbc_tcga_pan_can_atlas_2018%2Cesca_tcga_pan_can_atlas_2018%2Cgbm_tcga_pan_can_atlas_2018%2Chnsc_tcga_pan_can_atlas_2018%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Clihc_tcga_pan_can_atlas_2018%2Cluad_tcga_pan_can_atlas_2018%2Clusc_tcga_pan_can_atlas_2018%2Cmeso_tcga_pan_can_atlas_2018%2Cov_tcga_pan_can_atlas_2018%2Cpaad_tcga_pan_can_atlas_2018%2Cpcpg_tcga_pan_can_atlas_2018%2Cprad_tcga_pan_can_atlas_2018%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_tcga_pan_can_atlas_2018%2Cstad_tcga_pan_can_atlas_2018%2Ctgct_tcga_pan_can_atlas_2018%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Cucs_tcga_pan_can_atlas_2018%2Cucec_tcga_pan_can_atlas_2018%2Cuvm_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants%2Cgistic%2Crna_seq_v2_mrna_median_Zscores&case_set_id=all&gene_list=SMNDC1&geneset_list=%20&tab_index=tab_visualize&Action=Submit&plots_horz_selection=%7B%7D&plots_vert_selection=%7B%22selectedDataSourceOption%22%3A%22rna_seq_v2_mrna_median_all_sample_Zscores%22%2C%22logScale%22%3A%22true%22%7D&plots_coloring_selection=%7B%22selectedOption%22%3A%22undefined_%7B%5C%22clinicalAttributeId%5C%22%3A%5C%22CANCER_TYPE%5C%22%2C%5C%22patientAttribute%5C%22%3Afalse%2C%5C%22studyId%5C%22%3A%5C%22acc_tcga_pan_can_atlas_2018%5C%22%7D%22%7D
# and click download data in top right
# store in input_data/cbioportal/smndc1_express_rel_all_samples_v_copynum_TCGA_all_cancers.txt

# segment data
# this analysis: https://www.cbioportal.org/results/cnSegments?cancer_study_list=laml_tcga_pan_can_atlas_2018%2Cacc_tcga_pan_can_atlas_2018%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cchol_tcga_pan_can_atlas_2018%2Ccoadread_tcga_pan_can_atlas_2018%2Cdlbc_tcga_pan_can_atlas_2018%2Cesca_tcga_pan_can_atlas_2018%2Cgbm_tcga_pan_can_atlas_2018%2Chnsc_tcga_pan_can_atlas_2018%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Clihc_tcga_pan_can_atlas_2018%2Cluad_tcga_pan_can_atlas_2018%2Clusc_tcga_pan_can_atlas_2018%2Cmeso_tcga_pan_can_atlas_2018%2Cov_tcga_pan_can_atlas_2018%2Cpaad_tcga_pan_can_atlas_2018%2Cpcpg_tcga_pan_can_atlas_2018%2Cprad_tcga_pan_can_atlas_2018%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_tcga_pan_can_atlas_2018%2Cstad_tcga_pan_can_atlas_2018%2Ctgct_tcga_pan_can_atlas_2018%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Cucs_tcga_pan_can_atlas_2018%2Cucec_tcga_pan_can_atlas_2018%2Cuvm_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants%2Cgistic%2Crna_seq_v2_mrna_median_Zscores&case_set_id=all&gene_list=SMNDC1&geneset_list=%20&tab_index=tab_visualize&Action=Submit&plots_horz_selection=%7B%7D&plots_vert_selection=%7B%22selectedDataSourceOption%22%3A%22rna_seq_v2_mrna_median_all_sample_Zscores%22%2C%22logScale%22%3A%22true%22%7D&plots_coloring_selection=%7B%22selectedOption%22%3A%22undefined_%7B%5C%22clinicalAttributeId%5C%22%3A%5C%22CANCER_TYPE%5C%22%2C%5C%22patientAttribute%5C%22%3Afalse%2C%5C%22studyId%5C%22%3A%5C%22acc_tcga_pan_can_atlas_2018%5C%22%7D%22%7D&comparison_subtab=survival&comparison_selectedEnrichmentEventTypes=%5B%22HOMDEL%22%2C%22AMP%22%5D
# click download button in top right "download copy number segment file for the selected samples"
# store in input_data/cbioportal/cna/combined_study_segments.seg

# also create the "figures", "tables", and "output_data" directories to store the outputs into.

## load in data ----- 

mRNA_v_copy <- data.table::fread(input = "input_data/cbioportal/smndc1_express_rel_all_samples_v_copynum_TCGA_all_cancers.txt")
all_seg_data <- data.table::fread("input_data/cbioportal/cna/combined_study_segments.seg")



## cancer type v. expression v. copy number plot ---- 

z_threshold <- 1.2


# take the mRNA copy data, copy gains are gain + amplification, high expression >= z_threshold
mRNA_v_copy <- mRNA_v_copy |>  
  mutate(copy_up = case_when(
    `SMNDC1: Putative copy-number alterations from GISTIC` %in% c("Gain","Amplification") ~ "copy gain",
    TRUE ~ "no copy gain")) |> 
  mutate(expression_up = 
           case_when(`SMNDC1: mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM)` >= z_threshold ~ "high expression", 
                     TRUE ~ "no high expression") )


# make a column with all combinations of copy number and expression levels 
mRNA_v_copy_exp <- mRNA_v_copy |> 
  mutate(copy_by_expression = paste(copy_up, expression_up)) |> 
  group_by(`Cancer Type`) |> 
  count(copy_by_expression) |> 
  mutate(copy_by_expression = case_when(
    copy_by_expression == "copy gain high expression" ~ "Amplification and mRNA high",
    copy_by_expression == "copy gain no high expression" ~ "Amplification",
    copy_by_expression == "no copy gain high expression" ~ "mRNA high",
    copy_by_expression == "no copy gain no high expression" ~ "No alterations")) |> 
  mutate(copy_by_expression_fct = 
           factor(copy_by_expression, 
                  levels = c("Amplification and mRNA high","Amplification","mRNA high","No alterations")))


# calculate the proportion of each cancer type in each of the above buckets
mRNA_v_copy_exp <- mRNA_v_copy_exp |> 
  mutate(prop_type = n/sum(n))

# find the descending order of "no alterations" among cancer types
mRNA_v_copy_exp_order <- mRNA_v_copy_exp |> 
  filter(copy_by_expression == "No alterations") |> 
  arrange(desc(prop_type)) |> 
  pull(`Cancer Type`)


mRNA_v_copy_exp <- mRNA_v_copy_exp |> 
  mutate(cancer_type_order = factor(`Cancer Type`, levels = mRNA_v_copy_exp_order))


# find the total number of samples that have alterations
variant_prop_df <- mRNA_v_copy_exp |> 
  group_by(`Cancer Type`) |> 
  mutate(total_samples = sum(n)) |> 
  filter(copy_by_expression == "No alterations") |> 
  mutate(variant_prop = 1 - (n/total_samples))





# plot mRNA and copy number vs. cancer type, all cancer types
mRNA_v_copy_exp |> 
  ggplot(aes(y=cancer_type_order, x=prop_type,fill=copy_by_expression_fct)) + 
  geom_col() + 
  geom_richtext(data = variant_prop_df, 
            aes(y=`Cancer Type`,x=0, label = 
                  paste0(round(variant_prop*100,1),"%", " *N* = ",total_samples),
                fill = NA, label.color = NA
                ),
            hjust=0) + 
  scale_fill_viridis_d() + 
  labs(fill = "Copy and expression status",
       y=NULL,
       x="Proportion") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse = TRUE))

ggsave("figures/supp_fig_all_cancers_mRNA_amp.pdf")
ggsave("figures/supp_fig_all_cancers_mRNA_amp.png")


# same as above, but without the labels
mRNA_v_copy_exp |> 
  ggplot(aes(y=cancer_type_order, x=prop_type,fill=copy_by_expression_fct)) + 
  geom_col() + 
  # geom_richtext(data = variant_prop_df, 
  #               aes(y=`Cancer Type`,x=0, label = 
  #                     paste0(round(variant_prop*100,1),"%", " *N* = ",total_samples),
  #                   fill = NA, label.color = NA
  #               ),
  #               hjust=0) + 
  scale_fill_viridis_d() + 
  labs(fill = "Copy and expression status",
       y=NULL,
       x="Proportion") + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse = TRUE))

ggsave("figures/supp_fig_all_cancers_mRNA_amp_nonums.pdf")
ggsave("figures/supp_fig_all_cancers_mRNA_amp_nonums.png")


cancers_to_keep <- c("Pancreatic Cancer", 
                     "Ovarian Epithelial Tumor",
                     "Hepatobiliary Cancer", 
                     "Non-Small Cell Lung Cancer",
                     "Esophagogastric Cancer")
                     
cancers_to_keep <- cancers_to_keep[c(2,5,4,3,1)]



# same as above, but certain cancer types
fig1a <- mRNA_v_copy_exp |> 
  filter(`Cancer Type` %in% cancers_to_keep) |> 
  ggplot(aes(y=cancer_type_order, x=prop_type,fill=copy_by_expression_fct)) + 
  geom_col() + 
  geom_richtext(data = variant_prop_df |> filter(`Cancer Type` %in% cancers_to_keep), 
            aes(y=`Cancer Type`,x=0, label =  
                  paste0(round(variant_prop*100,1),"%", " *N* = ",total_samples),
                fill = NA, label.color = NA),
            hjust=0) + 
  scale_fill_viridis_d(breaks = levels(mRNA_v_copy_exp$copy_by_expression_fct)[c(1,3,2,4)]) +
  labs(fill = "SMNDC1 expression and copy number status",
       y=NULL,
       x="Proportion") + 
  theme_bw() + 
  theme(legend.position = "bottom",legend.title.position = "top") + 
  guides(fill = guide_legend(reverse = TRUE,nrow = 2))

ggsave("figures/fig_1a_mRNA_amp.pdf",height = 7,width = 7,plot = fig1a)
ggsave("figures/fig_1a_mRNA_amp.png",height = 7,width = 7,plot = fig1a)



# same as above, but without the labels. 
fig1a_nonums <- mRNA_v_copy_exp |> 
  filter(`Cancer Type` %in% cancers_to_keep) |> 
  ggplot(aes(y=cancer_type_order, x=prop_type,fill=copy_by_expression_fct)) + 
  geom_col() + 
  # geom_richtext(data = variant_prop_df |> filter(`Cancer Type` %in% cancers_to_keep), 
  #               aes(y=`Cancer Type`,x=0, label =  
  #                     paste0(round(variant_prop*100,1),"%", " *N* = ",total_samples),
  #                   fill = NA, label.color = NA),
  #               hjust=0) + 
  scale_fill_viridis_d(breaks = levels(mRNA_v_copy_exp$copy_by_expression_fct)[c(1,3,2,4)]) +
  labs(fill = "SMNDC1 expression and copy number status",
       y=NULL,
       x="Proportion") + 
  theme_bw() + 
  theme(legend.position = "bottom",legend.title.position = "top") + 
  guides(fill = guide_legend(reverse = TRUE,nrow = 2))

ggsave("figures/fig_1a_mRNA_amp_nonums.pdf",height = 7,width = 7,plot = fig1a_nonums)
ggsave("figures/fig_1a_mRNA_amp_nonums.png",height = 7,width = 7,plot = fig1a_nonums)





## segment analysis ------ 

# find the regions for SMNDC1
smndc1_ranges <- ces.refset.hg19$gr_genes[which(ces.refset.hg19$gr_genes$names == "SMNDC1"),]


# all seg data as a GRanges
all_seg_data_GR <- GenomicRanges::makeGRangesFromDataFrame(df = all_seg_data, 
                                                           keep.extra.columns = T, 
                                                           seqnames.field  = "chrom",
                                                           start.field = "loc.start",
                                                           end.field = "loc.end")


# overlapping segments of the gene
overlapping_segs_w_smndc1 <- subsetByOverlaps(all_seg_data_GR, smndc1_ranges, 
                                              type = "any") # warning OK because only looking at CHROM 10






overlaps_df <- GenomicRanges::as.data.frame(overlapping_segs_w_smndc1)
overlaps_df <- overlaps_df |> 
  mutate(length = end-start) |> # length of each segment
  arrange(desc(length)) |>  # desc order
  mutate(index_val = 1:n()) # make col that is 1 --> nrow()



# find just the samples with a gain in SMNDC1
mRNA_v_copy |> 
  # filter(`Cancer Type` == "Pancreatic Cancer") |> 
  filter(`SMNDC1: Putative copy-number alterations from GISTIC` %in% c("Gain","Amplification")) |> 
  pull(`Sample Id`) -> sample_w_gain


# filter the samples with an overlap in the gene by just those with a gain in the gene
overlaps_gain <- overlaps_df |> 
  filter(ID %in% sample_w_gain)


# smndc1 range
smndc1_df <- as.data.frame(smndc1_ranges)


# chrom 10 range
all_seg_data |> 
  filter(chrom == 10) |> 
  pull(loc.end) |> 
  max() -> 
  chrom_10_max


smndc1_start <- min(smndc1_df$start)
smndc1_end <- max(smndc1_df$end)




overlaps_gain <- left_join(overlaps_gain, mRNA_v_copy, by = c("ID" = "Sample Id"))


# glimpse(overlaps_gain)

write.csv(x = overlaps_gain, file = "tables/SMNDC1_CNAsegments_expression.csv")


overlaps_gain |> 
  group_by(`Cancer Type`) |> 
  mutate(index_val = 1:n()) |> 
  ggplot() + 
  geom_segment(aes(x = start, xend=end,y=index_val,yend=index_val,
                   # color = `SMNDC1: Putative copy-number alterations from GISTIC`
                   )) + 
  facet_wrap(~`Cancer Type`,scales = "free_y") + 
  geom_vline(xintercept = smndc1_start,color="red") + 
  geom_vline(xintercept = smndc1_end, color = "red") + 
  labs(y="Segments", x="Chromosome 10 position") + 
  scale_x_continuous(labels = scales::comma,limits = c(0,chrom_10_max)) + 
  scale_y_continuous(
    labels = NULL
  ) +
  scale_color_manual(values = c("#1b9e77","#d95f02"))  +
  theme_classic() + 
  theme(axis.ticks.y = element_blank()) + 
  labs(title = "CNA segments classified as gain or amplification overlapping with SMNDC1",
       subtitle = "SMNDC1 location is highlighted in red")


ggsave(filename = "figures/supp_fig_all_cancers_SMNDC1_gain_segments.pdf", height = 10, width = 17)
ggsave(filename = "figures/supp_fig_all_cancers_SMNDC1_gain_segments.png", height = 10, width = 17)






fig1b <- overlaps_gain |> 
  filter(`Cancer Type` %in% cancers_to_keep) |> 
  mutate(`Cancer Type` = factor(`Cancer Type`, levels = cancers_to_keep)) |> 
  group_by(`Cancer Type`) |> 
  mutate(index_val = 1:n()) |> 
  ggplot() + 
  geom_segment(aes(x = start, xend=end,y=index_val,yend=index_val,
                   # color = `SMNDC1: Putative copy-number alterations from GISTIC`
  )) + 
  facet_wrap(~`Cancer Type`,scales = "free_y",ncol = 1) + 
  geom_vline(xintercept = smndc1_start,color="red") + 
  geom_vline(xintercept = smndc1_end, color = "red") + 
  labs(y="Segments", x="Chromosome 10 position") + 
  scale_x_continuous(labels = scales::comma,limits = c(0,chrom_10_max)) + 
  # scale_y_continuous(
  #   # labels = NULL
  # ) +
  scale_color_manual(values = c("#1b9e77","#d95f02"))  +
  theme_classic() 
# + 
  # theme(axis.ticks.y = element_blank()) + 
  # labs(title = "CNA segments classified as gain or amplification overlapping with SMNDC1",
       # subtitle = "SMNDC1 location is highlighted in red")


ggsave(filename = "figures/fig_1b_SMNDC1_gain_segments.pdf", height = 7, 
       width = 7,plot = fig1b)

ggsave(filename = "figures/fig_1b_SMNDC1_gain_segments.png", height = 7, 
       width = 7,plot = fig1b)









overlaps_gain <- overlaps_gain |> 
  mutate(segment_begins_before_gene = start < smndc1_start) |> 
  mutate(segment_ends_after_gene = end > smndc1_end) |> 
  mutate(segment_within_gene =  !segment_begins_before_gene & !segment_ends_after_gene)

overlaps_gain$segment_begins_before_gene |> table()
overlaps_gain$segment_ends_after_gene |> table()

overlaps_gain$segment_within_gene |> table()


overlaps_gain |> 
  filter(!segment_begins_before_gene | !segment_ends_after_gene)

# only 2/767 have segments not completely overlapping SMNDC1, both in Endometrial Cancer



# April 2026 work ------ 
# for the SMNDC1 gain segments (attached). (Fig 1b)
# 
# 1) Can you please let me know what are the genes that are commonly amplified with SMNDC1? 
# 
# 2) Also, how big is the minimal segment that it is commonly amplified?



# (1) above, what genes are commonly amplified with SMNDC1 ----- 


overlaps_gain_cancersubset <- overlaps_gain |> 
  filter(`Cancer Type` %in% cancers_to_keep)


overlaps_gain_cancersubset

# find the overlap of this with and ces.refset.hg19$gr_genes

overlaps_gain_cancersubset_GR <- GenomicRanges::makeGRangesFromDataFrame(df = overlaps_gain_cancersubset, 
                                                           keep.extra.columns = T)


gr_genes_10 <- ces.refset.hg19$gr_genes[GenomicRanges::seqnames(ces.refset.hg19$gr_genes) == 10]


all_overlaps <- findOverlaps(gr_genes_10, overlaps_gain_cancersubset_GR, 
                                              type = "within") # warning OK because only looking at CHROM 10

all_overlaps |> 
  as.data.frame() |> 
  dplyr::count(queryHits, subjectHits)

# query is the exomes... gene names... need to collapse by names? 


names(gr_genes_10) <- gr_genes_10$names


gr_genes_10_split <- split(gr_genes_10,names(gr_genes_10))

gene_starts <- min(start(gr_genes_10_split))

gene_ends <- max(end(gr_genes_10_split))

all(names(gene_starts) == names(gene_ends)) # same genes same order 

gr_genes_10_short <- GRanges(seqnames = 10, ranges = IRanges(start = gene_starts, end = gene_ends),names = names(gene_starts))


all_overlaps <- findOverlaps(gr_genes_10_short, overlaps_gain_cancersubset_GR, 
                                              type = "within") # warning OK because only looki

# query = genes, subject = segment 

gene_names_index <- gr_genes_10_short |> 
  as.data.frame() |> 
  mutate(gene_index = 1:n())

# for each segment, need to find all the genes.
# then need to count the genes by number of times found

all_overlaps_named <- all_overlaps |> 
  as.data.frame() |> 
  left_join(gene_names_index, by=c("queryHits" = "gene_index"))


# counting gene names 
gene_shared_segments <- all_overlaps_named |> 
  dplyr::count(names,name = "number_of_segments_shared") |> 
  mutate(proportion_of_segments_shared = number_of_segments_shared/length(overlaps_gain_cancersubset_GR)) |> 
  arrange(desc(proportion_of_segments_shared))

gene_shared_segments

write_csv(x = gene_shared_segments, file = "output_data/gene_shared_segments.csv")

# (2) Also, how big is the minimal segment that it is commonly amplified? ------ 

median_non_whole_chrom_length <- overlaps_gain_cancersubset |> 
  filter(length < max(length)) |> 
  summarize(median(length))


ggplot(data = overlaps_gain_cancersubset) + 
  geom_histogram(aes(x=length)) + 
  # geom_vline(xintercept = as.numeric(median_non_whole_chrom_length)) + 
  labs(x="Segment length", y= "Count", title = "All gain segments overlapping with SMNDC1")

ggsave(filename = "figures/gain_segment_lengths.pdf")


overlaps_gain_noWGD_density <- overlaps_gain_cancersubset |> 
  filter(length < max(length)) |> 
  pull(length) |> 
  density()

peak_x <- overlaps_gain_noWGD_density$x[which.max(overlaps_gain_noWGD_density$y)]
peak_y <- overlaps_gain_noWGD_density$y[which.max(overlaps_gain_noWGD_density$y)]





ggplot(data = overlaps_gain_cancersubset |> 
  filter(length < max(length))) + 
  geom_density(aes(x=length)) + 
  geom_vline(xintercept = peak_x) + 
  labs(x="Segment length", y= "Density", 
    title = "Sub-whole-chromosomal Gain segments overlapping with SMNDC1", 
    subtitle = paste("Vertical line represents peak of density plot:",round(peak_x)))

ggsave(filename = "figures/sub_whole_chrom_segment_gain_density.pdf")


overlaps_gain_cancersubset |> 
  filter(length < max(length)) |> 
  summarize(median(length))
  

#  which genes are in this smaller segment bucket with SMNDC1? 

overlaps_gain_cancersubset_sm_segments <- overlaps_gain_cancersubset |> 
  filter(length < peak_x)


overlaps_gain_cancersubset_GR_sm_segments <- GenomicRanges::makeGRangesFromDataFrame(df = overlaps_gain_cancersubset_sm_segments, 
                                                           keep.extra.columns = T)


all_overlaps_sm_segments <- findOverlaps(gr_genes_10_short, overlaps_gain_cancersubset_GR_sm_segments, 
                                              type = "within") # warning OK 

# query = genes, subject = segment 

# for each segment, need to find all the genes.
# then need to count the genes by number of times found

all_overlaps_named_sm_segments <- all_overlaps_sm_segments |> 
  as.data.frame() |> 
  left_join(gene_names_index, by=c("queryHits" = "gene_index"))


# counting gene names 
gene_shared_segments_sm_segments <- all_overlaps_named_sm_segments |> 
  dplyr::count(names,name = "number_of_segments_shared") |> 
  mutate(proportion_of_segments_shared = number_of_segments_shared/length(overlaps_gain_cancersubset_GR_sm_segments)) |> 
  arrange(desc(proportion_of_segments_shared))

gene_shared_segments_sm_segments

write_csv(x = gene_shared_segments_sm_segments, file = "output_data/gene_shared_segments_smaller_segments.csv")
