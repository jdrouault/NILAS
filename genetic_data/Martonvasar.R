evtools::install_github("hadley/readxl") # development version


library(readxl)
library(xlsx)
library(tidyverse)
library(rstatix)
library(ggpmisc)
library(agricolae)
library(multcompView)
library(emmeans)
library(pwr)
library(patchwork)

cbPalette <- c("grey", "#0072B2", "#56B4E9", "#E69F00", "#C77CFF" ,"#F0E442", "#009E73", "#56B4E9",  "#0072B2", "#D55E00", "#CC79A7", "darkblue") #red


### date today 
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)

###### ======= Create paths 
root <- ("~/Dropbox-SE Dropbox/allshare/") #C:/Users/drouault
root_alternative <-"~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/" #C:/Users/drouault

NILAS_genotypes_path <- paste0(root, "stage_eloise/NILAS_haploimputer/Data/dp6q30dr0.40maf0.01/outputs/")
haplot_imput_data <- paste0(root, "stage_eloise/NILAS_haploimputer/Data/")
save_path <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Tables/genotypic_data/consistency_replicates/")
phenotype_path <- paste0(root_alternative, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_day_dta_04Avr11/") ### need to change the path to access either the blues relatedd to gdd_dta_em or to gdd_dta 
fig_consistency <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/constistency_replicates/")
fig_genes_introgression <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/effect_genes_introgressionNB/")
fig_phenotype <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/phenotypic_data/")
file_path <- paste0(root_alternative, "/NILAS/main_folder/drouault/paper/Files/")
genetic_data_save_path <- paste0(root_alternative, "/NILAS/main_folder/drouault/data/genotypic_data/NILs/")
phenospace_path <-paste0(root_alternative, "PhotoGrid/paper/data/phenotypic_space/")
mixed_model_results <- paste0(root_alternative,"NILAS/main_folder/drouault/paper/Files/")
phenotypic_data <- paste0(root_alternative,"NILAS/main_folder/drouault/data/cleaned_data/NILs_diversity_panel/")
supp_file <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Files/")

######## ========== Only phenotypic analysis 

#### ========== Load data 

## Data with linear model and fixed effect

lines_to_select <- read.csv(paste0(genetic_data_save_path, "20240702Selected_NILs_after_filtration.csv"), header = T) ### the 490 selected lines #20240619Selected_NILs_after_filtration.csv

adjustement_for_nb_genes_effect <- lines_to_select %>%
  dplyr::select(entryid, sum_gene_pr, sum_gene_fl, Nb_ft_pr_gwa) %>%
  dplyr::rename(entryid2 = entryid)

lines_PRNFT <- read.csv(paste0(phenotype_path, "lines_slope_dta_em_LD_SD_lmregression_parents_20240723.csv"), header = T) %>%
  dplyr::mutate(qtl = as.character(qtl)) %>%
  dplyr::rename(indi.cof = indi.cof, intercept = intercept) %>%
  dplyr::filter(entryid2 %in% lines_to_select$entryid) %>%
  dplyr::full_join(adjustement_for_nb_genes_effect) %>%
  dplyr::mutate(indi.cof_corr_nb_genes = indi.cof - (0.0131 * (sum_gene_pr + Nb_ft_pr_gwa)), 
                intercept_corr_nb_genes = intercept - (0.0729 * (sum_gene_fl + Nb_ft_pr_gwa)))


parents_PRNFT <- read.csv(paste0(phenotype_path, "parents_slope_dta_em_LD_SD_lmregression_20240723.csv"), header = T) %>%
  dplyr::mutate(qtl = as.character(qtl)) %>%
  dplyr::rename(indi.cof = indi.cof_corrected, intercept = intercept_corrected) %>%
  dplyr::mutate(genotype_to_color = ifelse(is.na(dp_name), rp_name , dp_name )) %>%
  dplyr::mutate(indi.cof_corr_nb_genes = indi.cof, 
                intercept_corr_nb_genes = intercept) # consider it as dummie variable to make the plot 

## Data with mixed model and NIL random effect

prnftpt_mixed_model <- read_excel(paste0(mixed_model_results, "File_S4_new_20241008.xlsx")) %>%
  as.data.frame()

slope_int_mixed_model_haplotype <- prnftpt_mixed_model %>%
  dplyr::filter(genetic_type== "haplotype") %>%
  dplyr::mutate(across(intercept:slope, ~as.numeric(as.character(.)))) %>%
  dplyr::select(entryid, intercept, slope, rp, dp, rp_name, dp_name) %>%
  dplyr::rename(entryid2 = entryid)

## data with single introgression 

  genetic_info <- read.csv(paste0(supp_file, "File_S5.csv"), header = T, sep = ";", dec = ",")


#### Create set of 20 genotypes for Martonvasar

slope_int_mixed_model_haplotype_20 <- slope_int_mixed_model_haplotype %>%
  dplyr::filter(entryid2 %in% genetic_info$entryid ) %>%
  dplyr::arrange(desc(slope)) %>%
  dplyr::top_n(20, slope)
  
final_file <- genetic_info %>%
  dplyr::filter(entryid %in% slope_int_mixed_model_haplotype_20$entryid2 ) %>%
  dplyr::rename(entryid2= entryid) %>%
  dplyr::full_join(slope_int_mixed_model_haplotype_20) %>%
  dplyr::select(entryid2, intercept, slope, rp, dp, rp_name, dp_name, ZmPR, 
                  DP_BG:RP_FG, Nb_ZmPR1:Nb_ZmPR4, sum_introgression, 
                zmpr1_cct:zmpr4_cacta, biol_rep, tech_rep, 
                sum_gene_pr:sum_gene_fl) %>%
  dplyr::arrange(desc(slope))
  


write.csv(final_file, paste0(supp_file, "Martonvasar_subset_NILAS_", date_today, ".csv"), row.names = F)
   
