devtools::install_github("hadley/readxl") # development version


library(ggplot2)
library(readxl)
library(xlsx)
library(tidyverse)
library(ggpubr)
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

slope_int_mixed_model_NIL <- prnftpt_mixed_model %>%
  dplyr::filter(genetic_type== "NIL") %>%
  dplyr::mutate(across(intercept:slope, ~as.numeric(as.character(.)))) %>%
  dplyr::select(entryid, intercept, slope, rp, dp, rp_name, dp_name) %>%
  dplyr::rename(entryid2 = entryid)

slope_int_mixed_model_parents <- prnftpt_mixed_model %>%
  dplyr::filter(genetic_type== "parent") %>%
  dplyr::mutate(across(intercept:slope, ~as.numeric(as.character(.)))) %>%
  dplyr::select(entryid, intercept, slope, rp, dp, rp_name, dp_name) %>%
  dplyr::rename(entryid2 = entryid) %>%
  dplyr::mutate(genotype_to_color = ifelse(dp_name =="NA", rp_name , dp_name )) 
  


## Diversity panel 

diversity_panel <- read.csv(paste0(phenospace_path, "20240705Phenospace_diversity_panel_LD_SD_lmregression.csv"), header = T, sep = ",") #20240530Phenospace_diversity_panel_LD_SD_lmregression.csv

#### ========== Mean common set genotype 

# parents_PRNFT1 <- parents_PRNFT %>%
#   dplyr::rename(genotype = genotype_to_color )

parents_PRNFT1 <- slope_int_mixed_model_parents %>%
  dplyr::rename(genotype = genotype_to_color )

common_geno <- diversity_panel %>%
  dplyr::filter(genotype %in% unique(parents_PRNFT1$genotype)) %>%
  dplyr::select(genotype, indi.cof_corrected, intercept_corrected)

diff <- parents_PRNFT1 %>%
  dplyr::select(genotype, slope, intercept) %>% #indi.cof
  dplyr::full_join(common_geno) %>%
  tidyr::pivot_longer(c=(slope:intercept_corrected), #indi.cof
                      names_to   = "parameters", 
                      values_to = ".values")

correction_plot <- diff %>%
  # dplyr::filter(parameters == "indi.cof" | parameters == "indi.cof_corrected") %>%
  dplyr::filter(parameters == "intercept" | parameters == "intercept_corrected") %>%
  dplyr::rename(population = parameters) %>%
  # dplyr::mutate(population = dplyr::recode(population, indi.cof = "NILAS",   indi.cof_corrected = "diversity_panel")) %>%
  dplyr::mutate(population = dplyr::recode(population, intercept = "NILAS",   intercept_corrected = "diversity_panel")) %>%
   ggplot(aes(x = .values, y = genotype)) +
  geom_point(size = 3, aes(color = population)) +
  # xlab("Slope of sensitivity") + #
  xlab("flo per se") + #
  ylab("genotypes") + 
  scale_color_manual(values = cbPalette) +
  # ylim(-30,20) +
  theme(plot.background = element_blank()
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
        ,aspect.ratio = 1
        ,title = element_text(size=14, face="bold")
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=18, margin = margin(t=18))
        ,axis.text.x = element_text(size=18, margin = margin(t=11), angle=90, vjust = 0.5)
        ,axis.title.y = element_text(size=18, margin = margin(r=18))
        ,axis.text.y = element_text(size=18)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  )
correction_plot

# ggsave(file = paste0(fig_phenotype,"Flo_per_se_deviation_diversity_panel_NILAS" ,date_today, ".pdf"),
       # correction_plot, width = 10, height = 12)

### Get teh general mean and the deviation from the mean for slope and BVP 

diff <- parents_PRNFT1 %>%
  dplyr::select(genotype, slope, intercept) %>% #indi.cof
  dplyr::full_join(common_geno) %>%
  dplyr::filter(!is.na(indi.cof_corrected )) %>%#is going to remove all the genotypes that are not common between NILAS and diversity panel 
  dplyr::mutate(mean_slope_NILAS = mean(slope), #indi.cof
                mean_interecpt_NILAS = mean(intercept), 
                mean_slope_divpanel = mean(indi.cof_corrected), 
                mean_intercept_divpanel = mean(intercept_corrected), 
                common_mean_slope = (mean_slope_NILAS + mean_slope_divpanel)/ 2, 
                common_mean_intercept = (mean_interecpt_NILAS + mean_intercept_divpanel)/ 2, 
                deviation_from_slope_NILAS = abs(common_mean_slope - mean_slope_NILAS), 
                deviation_from_slope_divpanel = abs(common_mean_slope - mean_slope_divpanel), 
                deviation_from_intercept_NILAS = abs(common_mean_intercept - mean_interecpt_NILAS), 
                deviation_from_intercept_divpanel = abs(common_mean_intercept - mean_intercept_divpanel)) 
  

correction_plot_slope <- diff %>%
  dplyr::select(mean_slope_NILAS:deviation_from_intercept_divpanel) %>%
  dplyr::mutate(type = "slope") %>%
  dplyr::distinct() %>%
  # tidyr::pivot_longer(col = mean_slope_NILAS:deviation_from_intercept_divpanel, 
  #                     names_to = "parameters", values_to = ".values") %>%
  ggplot(aes(x = type, y = common_mean_slope)) +
  geom_pointrange(aes(ymin = common_mean_slope-deviation_from_slope_NILAS, ymax = common_mean_slope+deviation_from_slope_divpanel)) +
    xlab("parameter") + #
  ylab("Deviation from common mean") +
  # ylim(-30,20) +
  theme(plot.background = element_blank()
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
        ,aspect.ratio = 1
        ,title = element_text(size=14, face="bold")
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=18, margin = margin(t=18))
        ,axis.text.x = element_text(size=18, margin = margin(t=11), angle=90, vjust = 0.5)
        ,axis.title.y = element_text(size=18, margin = margin(r=18))
        ,axis.text.y = element_text(size=18)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  )
correction_plot_slope

#### ========== Adjust NILAS and the diversity panel with the deviatoion from the general mean to be able to superpose them 

correction <- diff %>%
  dplyr::select(deviation_from_slope_NILAS:deviation_from_intercept_divpanel) %>%
  dplyr::distinct()

diversity_panel1 <- diversity_panel %>%
  dplyr::select(genotype, indi.cof_corrected, intercept_corrected) %>%
  dplyr::rename(indi.cof_raw = indi.cof_corrected, intercept_raw = intercept_corrected ) %>%
  dplyr::mutate(indi.cof_adjusted = indi.cof_raw - correction$deviation_from_slope_divpanel,
                intercept_adjusted = intercept_raw - correction$deviation_from_intercept_divpanel) %>%
  dplyr::mutate(exp = "Diversity panel") #Panzea_MDR

NILAS1 <- slope_int_mixed_model_NIL %>% #lines_PRNFT %>%
  dplyr::full_join(slope_int_mixed_model_parents) %>% #parents_PRNFT
  dplyr::rename(genotype = entryid2, indi.cof_raw = slope, intercept_raw = intercept ) %>%
  dplyr::select(genotype_to_color, genotype, indi.cof_raw, intercept_raw) %>%
  dplyr::mutate(indi.cof_adjusted = indi.cof_raw + correction$deviation_from_slope_divpanel,
                intercept_adjusted = intercept_raw + correction$deviation_from_intercept_divpanel) %>%
  dplyr::mutate(exp = "NILAS") #Maize ATLAS

combined_full <- NILAS1 %>%
  dplyr:::full_join(diversity_panel1) %>%
  dplyr::mutate(genotype_to_color = ifelse(genotype %in% lines_to_select$entryid | exp == "Diversity panel" | !is.na(genotype_to_color), genotype_to_color, "several_introgressions"))

write.csv(combined_full, paste0(phenotypic_data, "NILs_combined_diversity_panel", "_",date_today, ".csv"), row.names = F) 


#### ========== Superpose NILAS and the diversity panel

# diversity_panel1 <- diversity_panel %>%
#   dplyr::select(genotype, indi.cof_corrected, intercept_corrected) %>%
#   dplyr::rename(indi.cof_raw = indi.cof_corrected, intercept_raw = intercept_corrected ) %>%
#   dplyr::mutate(exp = "Panzea_MDR")
# 
# NILAS1 <- lines_PRNFT %>%
#   dplyr::full_join(parents_PRNFT) %>%
#   dplyr::rename(genotype = entryid2, indi.cof_raw = indi.cof, intercept_raw = intercept ) %>%
#   dplyr::select(genotype_to_color, genotype, indi.cof_raw, intercept_raw) %>%
#   dplyr::mutate(exp = "Maize ATLAS")
# 
# combined_full <- NILAS1 %>%
#   dplyr:::full_join(diversity_panel1)
  

#### ========== Center and scale the trait space with NILAS and the diversity panel
# 
# NILAS <- lines_PRNFT %>%
#   dplyr::full_join(parents_PRNFT)
# 
# ## Get the mean for each variable for each dataset 
# mean_PS_diversity_panel <- mean(diversity_panel$indi.cof_corrected, na.rm = T)
# mean_BVP_diversity_panel <- mean(diversity_panel$intercept_corrected, na.rm = T)
# 
# mean_PS_lines_PRNFT <- mean(NILAS$indi.cof, na.rm = T)
# mean_BVP_lines_PRNFT <- mean(NILAS$intercept, na.rm = T)
# 
# # mean_PS_parents_PRNFT <- mean(parents_PRNFT$indi.cof, na.rm = T)
# # mean_BVP_parents_PRNFT <- mean(parents_PRNFT$intercept, na.rm = T)
# 
# ## Get the sd for each variable for each dataset 
# sd_PS_diversity_panel <- sd(diversity_panel$indi.cof_corrected, na.rm = T)
# sd_BVP_diversity_panel <- sd(diversity_panel$intercept_corrected, na.rm = T)
# 
# sd_PS_lines_PRNFT <- sd(NILAS$indi.cof, na.rm = T)
# sd_BVP_lines_PRNFT <- sd(NILAS$intercept, na.rm = T)
# 
# # sd_PS_parents_PRNFT <- sd(parents_PRNFT$indi.cof, na.rm = T)
# # sd_BVP_parents_PRNFT <- sd(parents_PRNFT$intercept, na.rm = T)
# 
# 
# ## Center the Variables
# #Subtract the mean from each observation in their respective datasets
# 
# diversity_panel <- diversity_panel %>%
#   dplyr::mutate(indi.cof_centered = indi.cof_corrected - mean_PS_diversity_panel,
#                 intercept_centered = intercept_corrected - mean_BVP_diversity_panel)
# 
# NILAS <- NILAS %>%
#   dplyr::mutate(indi.cof_centered = indi.cof - mean_PS_lines_PRNFT, 
#                 intercept_centered = intercept -  mean_BVP_lines_PRNFT)
# 
# # parents_PRNFT <- parents_PRNFT %>%
# #   dplyr::mutate(indi.cof_centered = indi.cof - mean_PS_parents_PRNFT, 
# #                 intercept_centered = intercept -  mean_BVP_parents_PRNFT)
#   
# 
# ## Scale the Variables
# #Divide the centered variables by their respective standard deviation
# 
# diversity_panel <- diversity_panel %>%
#   dplyr::mutate(indi.cof_scaled = indi.cof_centered / sd_PS_diversity_panel ,
#                 intercept_scaled = intercept_centered / sd_BVP_diversity_panel)  %>%
#   dplyr::select( indi.cof_scaled, intercept_scaled, genotype) %>%
#   dplyr::mutate(exp = "Panzea_MDR")
# 
# NILAS <- NILAS %>%
#   dplyr::mutate(indi.cof_scaled = indi.cof_centered / sd_PS_lines_PRNFT, 
#                 intercept_scaled = intercept_centered / sd_BVP_lines_PRNFT ) %>%
#   dplyr::select(genotype_to_color, indi.cof_scaled, intercept_scaled, entryid2) %>%
#   dplyr::rename(genotype = entryid2) %>%
#   dplyr::mutate(exp = "Maize ATLAS")
# 
# # parents_PRNFT <- parents_PRNFT %>%
# #   dplyr::mutate(indi.cof_scaled = indi.cof_centered / sd_PS_parents_PRNFT, 
# #                 intercept_scaled = intercept_centered / sd_BVP_parents_PRNFT ) %>%
# #   dplyr::select(genotype_to_color, indi.cof_scaled, intercept_scaled, entryid2) %>%
# #   dplyr::rename(genotype = entryid2) %>%
# #   dplyr::mutate(exp = "Maize ATLAS")
# 
# combined_datasets <- NILAS %>%
#   # dplyr::full_join(parents_PRNFT) %>%
#   dplyr::full_join(diversity_panel) 

#### ========== Phenotypic Space NILAS

###### Trait space with diversity panel and NILAS centered and scaled 

list <- c( "CML10", "CML258" , "CML277", "CML341", "CML373",
           "Tzi8", "Tzi9", "LH123Ht", "2369" )

pheno_space_centered_scaled <- combined_full  %>% #combined_datasets
  dplyr::mutate(#dp = as.factor(dp), 
                #rp = as.factor(rp), 
                #dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                # genotype_to_color = as.factor(genotype_to_color),
                # genotype_to_color = factor(genotype_to_color, levels = list),
                # genotype_to_color = ifelse(genotype_to_color == "LH123Ht"| genotype_to_color== "2369", "Background parents", "Donor parents"), 
                exp1 = ifelse(genotype_to_color == "several_introgressions", "Several introgressions", exp), 
                exp1 = ifelse(is.na(exp1), exp, exp1)) %>%
  dplyr::select(-exp) %>%
  dplyr::rename(exp = exp1) %>%
  dplyr::mutate(exp = factor(exp, levels = c( "Diversity panel", "NILAS", "Several introgressions"))) %>%
  dplyr::arrange(exp) %>%
  ggplot() +
  geom_point(aes(x = intercept_adjusted #_raw, #scaled
                 , y = indi.cof_adjusted #_raw, #scaled
                 ,  color = exp), alpha = 0.8, size = 3, fill = "black") 
  
pheno_space_centered_scaled <-pheno_space_centered_scaled +
    geom_point(data = combined_full %>% filter(!is.na(.$genotype_to_color) & genotype_to_color != "several_introgressions") %>% 
      dplyr::mutate(genotype_to_color = as.factor(genotype_to_color),
                    genotype_to_color = factor(genotype_to_color, levels = list),
      genotype_to_color = ifelse(genotype_to_color == "LH123Ht"| genotype_to_color== "2369" , "Background parents", "Donors parents")), 
       aes(x = intercept_adjusted #_raw, scaled
          , y = indi.cof_adjusted #_raw, scaled
          , color = (genotype_to_color)), size = 5) +  
  # geom_xsideboxplot(data = . %>% filter(is.na(genotype_to_color)), aes( x = intercept, y = rp,  group = rp), orientation = "y", stat = "boxplot") +
  # geom_ysideboxplot( data = . %>% filter(is.na(genotype_to_color)), aes(x = rp, y = indi.cof, group = rp), orientation = "x", stat = "boxplot")+
  
  # geom_point(aes(x = dp, y = slope, group = Nb_ZmPR4_peak), fill = "black",position=position_dodge(width=0.75)) + 
  # scale_fill_manual(values = c("#5DA5DA", "#B276B2")) +
  scale_color_manual(values=cbPalette, breaks = c("Diversity panel", "NILAS", "Several introgressions","Background parents", "Donors parents")) +
  theme(plot.background = element_blank()
        ,legend.position = "right"#legend.position = c(0.5,0.8)#"top"
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 20, face="bold")#, family="Courier New Bold")
        ,legend.title = element_text(size = 20, face="bold")
        ,aspect.ratio = 1
        #,legend.title = element_text(colour="blue", size=16, face="bold")
        #,legend.text = element_text(colour="blue", size=16, face="bold")
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
        #,panel.spacing = unit(1, "lines")
        ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
        ,panel.background = element_blank()
        #,strip.background = element_blank()
        #,strip.text = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=30, margin = margin(t=18), face="bold")
        ,axis.text.x = element_text(size=30, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(r=18), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  ylab((quote(bold("Photoperiod sensitivity") ~ bold("(") ~bold(GDD[hour]) ~ bold(")")))) +
  xlab(quote(bold("Flowering time") ~ bolditalic("per se") ~ bold("(GDD)"))) +
labs(color = "")



# ### Change the plot title when saving whether it is adjusted or centered and scaled 
# ggsave(paste0(fig_phenotype, "pheno_space_centered_scaled","_", date_today,".pdf"), pheno_space_centered_scaled, width = 8, height = 6)
# 
# 

# 
# library(ggside)
# 
# pheno_space <- slope_int_mixed_model_NIL %>% #lines_PRNFT
#   dplyr::full_join(slope_int_mixed_model_parents) %>% #parents_PRNFT
#   dplyr::mutate(dp = as.factor(dp), 
#                 rp = as.factor(rp), 
#                 dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
#                 genotype_to_color = as.factor(genotype_to_color),
#                 genotype_to_color = factor(genotype_to_color, levels = list)) %>%
#   dplyr::arrange(genotype_to_color) %>%
#   ggplot() +
#   geom_point(aes(x = intercept_corr_nb_genes, y = indi.cof_corr_nb_genes), color = "black", alpha = 0.5, size = 3) + 
#   geom_point(data = . %>% filter(!is.na(genotype_to_color)), aes(x = intercept, y = indi.cof, color = (genotype_to_color)), size = 3) +
#   # geom_xsideboxplot(data = . %>% filter(is.na(genotype_to_color)), aes( x = intercept, y = rp,  group = rp), orientation = "y", stat = "boxplot") +
#   # geom_ysideboxplot( data = . %>% filter(is.na(genotype_to_color)), aes(x = rp, y = indi.cof, group = rp), orientation = "x", stat = "boxplot")+
#   
#   # geom_point(aes(x = dp, y = slope, group = Nb_ZmPR4_peak), fill = "black",position=position_dodge(width=0.75)) + 
#   # scale_fill_manual(values = c("#5DA5DA", "#B276B2")) +
#   scale_color_manual(values=cbPalette) +
#   theme(plot.background = element_blank()
#         #,legend.position = "bottom"#legend.position = c(0.5,0.8)#"top"
#         ,legend.direction="vertical"#,legend.direction="vertical"
#         ,legend.key.size = unit(1, "cm")
#         ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
#         ,legend.title = element_text(size = 12, face="bold")
#         ,aspect.ratio = 1
#         #,legend.title = element_text(colour="blue", size=16, face="bold")
#         #,legend.text = element_text(colour="blue", size=16, face="bold")
#         #,panel.grid.major = element_blank()
#         ,panel.grid.minor = element_blank()
#         ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
#         ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
#         #,panel.spacing = unit(1, "lines")
#         ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
#         ,panel.background = element_blank()
#         #,strip.background = element_blank()
#         #,strip.text = element_blank()
#         ,plot.title = element_text(face="bold")
#         ,axis.title.x = element_text(size=18, margin = margin(t=18), face="bold")
#         ,axis.text.x = element_text(size=18, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) +
#   xlab(quote(bold("Flowering time") ~ bolditalic("per se"))) +
#   ylab("Photoperiod sensitivity") +
#   labs(color = "Parents")
# 
# # geom_text(aes(x = dp, y = 25,  label = Spearman_p_sign), size = 6, check_overlap = TRUE)
# pheno_space 


#ggsave(paste0(fig_phenotype, "pheno_space_corr","_", date_today,".pdf"), pheno_space, width = 8, height = 6)


###### Boxplots 

boxplot_slope <- slope_int_mixed_model_NIL %>% #lines_PRNFT
  dplyr::full_join(slope_int_mixed_model_parents) %>% #parents_PRNFT
  dplyr::mutate(dp = as.factor(dp), 
                rp = as.factor(rp), 
                dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                genotype_to_color = as.factor(genotype_to_color),
                genotype_to_color = factor(genotype_to_color, levels = list)) %>%
  dplyr::arrange(genotype_to_color) %>%
  ggplot() + 
geom_boxplot( data = . %>% filter(is.na(genotype_to_color)), aes(x = rp, y = slope, group = rp)) + 
  theme(plot.background = element_blank()
        #,legend.position = "bottom"#legend.position = c(0.5,0.8)#"top"
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
        ,legend.title = element_text(size = 12, face="bold")
        # ,aspect.ratio = 1
        #,legend.title = element_text(colour="blue", size=16, face="bold")
        #,legend.text = element_text(colour="blue", size=16, face="bold")
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
        #,panel.spacing = unit(1, "lines")
        ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
        ,panel.background = element_blank()
        #,strip.background = element_blank()
        #,strip.text = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=30, margin = margin(t=18), face="bold")
        ,axis.text.x = element_text(size=30, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(r=18), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  # xlab("Backcross parents") +
  ylab((quote(bold("Photoperiod sensitivity") ~ bold("(") ~bold(GDD[hour]) ~ bold(")")))) 
  
boxplot_intercept <- slope_int_mixed_model_NIL %>% #lines_PRNFT
  dplyr::full_join(slope_int_mixed_model_parents) %>% #parents_PRNFT
  dplyr::mutate(dp = as.factor(dp), 
                rp = as.factor(rp), 
                dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                genotype_to_color = as.factor(genotype_to_color),
                genotype_to_color = factor(genotype_to_color, levels = list)) %>%
  dplyr::arrange(genotype_to_color) %>%
  ggplot(data = . %>% filter(is.na(genotype_to_color)), aes(x = rp, y = intercept, group = rp)) + 
  geom_boxplot( ) + 
  theme(plot.background = element_blank()
        #,legend.position = "bottom"#legend.position = c(0.5,0.8)#"top"
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
        ,legend.title = element_text(size = 12, face="bold")
        # ,aspect.ratio = 1
        #,legend.title = element_text(colour="blue", size=16, face="bold")
        #,legend.text = element_text(colour="blue", size=16, face="bold")
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
        #,panel.spacing = unit(1, "lines")
        ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
        ,panel.background = element_blank()
        #,strip.background = element_blank()
        #,strip.text = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=30, margin = margin(t=18), face="bold")
        ,axis.text.x = element_text(size=30, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(r=18), face="bold")
        ,axis.text.y = element_text(30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  # xlab("Backcross parents") +
  ylab(quote(bold("Flowering time") ~ bolditalic("per se") ~ bold("(GDD)"))) 


data_boxplot <- slope_int_mixed_model_NIL %>% #lines_PRNFT
  dplyr::full_join(slope_int_mixed_model_parents) %>% #parents_PRNFT
  dplyr::mutate(dp = as.factor(dp), 
                rp = as.factor(rp), 
                dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                genotype_to_color = as.factor(genotype_to_color),
                genotype_to_color = factor(genotype_to_color, levels = list), 
                rp = ifelse(rp == 1, "2369", "LH123Ht")) %>%
  dplyr::arrange(genotype_to_color) 

boxplot_slope <-  ggplot(data = data_boxplot %>% filter(is.na(genotype_to_color)), aes(x = rp, y = slope, group = rp)) + ylab((quote(bold("Photoperiod sensitivity") ~ bold("(") ~bold(GDD[hour]) ~ bold(")"))))     
boxplot_intercept <- ggplot(data = data_boxplot %>% filter(is.na(genotype_to_color)), aes(x = rp, y = intercept, group = rp)) + ylab(quote(bold("Flowering time") ~ bolditalic("per se") ~ bold("(GDD)"))) 



p3b <- NULL
box_p <- NULL
box_p <- (boxplot_slope | boxplot_intercept)  +
  plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 28, face = "bold"))
  ) &
  geom_boxplot( )   &
  theme_classic() &
  theme(
    plot.background = element_blank()
    #,legend.position = "bottom"#legend.position = c(0.5,0.8)#"top"
    ,legend.direction="vertical"#,legend.direction="vertical"
    ,legend.key.size = unit(1, "cm")
    ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
    ,legend.title = element_text(size = 12, face="bold")
    # ,aspect.ratio = 1
    #,legend.title = element_text(colour="blue", size=16, face="bold")
    #,legend.text = element_text(colour="blue", size=16, face="bold")
    #,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
    ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
    #,panel.spacing = unit(1, "lines")
    ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
    ,panel.background = element_blank()
    #,strip.background = element_blank()
    #,strip.text = element_blank()
    ,plot.title = element_text(face="bold")
    ,axis.title.x = element_text(size=30, margin = margin(t=18), face="bold")
    ,axis.text.x = element_text(size=30, margin = margin(t=40),  angle = 45)
    ,axis.title.y = element_text(size=30, margin = margin(r=18), face="bold")
    ,axis.text.y = element_text(size=30, margin = margin(r=18), face="bold")
    ,axis.ticks.y = element_blank()
    ,axis.ticks.x = element_blank()
    # top, right, bottom, left
    ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
    
  ) &
  labs(x = "Backcross parents")
 
  
  
p3b <- pheno_space_centered_scaled + box_p + plot_layout(ncol = 2, nrow = 1, guides = 'keep', tag_level = 'new')  + plot_annotation(tag_levels = c('A'))

pdf(paste(fig_phenotype, "pheno_space_boxplots","_", date_today,".pdf", sep=""),
    width = 10, height = 11)
box_p
dev.off()

pdf(paste(fig_phenotype, "pheno_space","_", date_today,".pdf", sep=""),
    width = 12, height = 12)
pheno_space_centered_scaled
dev.off()


ggsave(paste0(fig_phenotype, "pheno_space_boxplots","_", date_today,".pdf"), p3b, width = 30, height = 20)

# geom_point(aes(x = dp, y = slope, group = Nb_ZmPR4_peak), fill = "black",position=position_dodge(width=0.75)) + 

##### Filter the genotype that are higher in slope and higher and lower BVP than the rp

slope_LH123Ht <- slope_int_mixed_model_NIL %>%
  dplyr::filter(  rp_name == "LH123Ht" & indi.cof > 1.048985) #intercept > 799
nrow(slope_LH123Ht) 
min(slope_LH123Ht$indi.cof)
max(slope_LH123Ht$indi.cof)

slope_2369 <- slope_int_mixed_model_NIL %>%
  dplyr::filter(rp_name == "2369" & indi.cof > -8.631828) #intercept > 722
nrow(slope_2369) 
min(slope_2369$indi.cof)
max(slope_2369$indi.cof)


BVP_LH123Ht <- lines_PRNFT %>%
  dplyr::filter(  rp_name == "LH123Ht") #intercept > 799
nrow(BVP_LH123Ht) 
min(BVP_LH123Ht$intercept)
max(BVP_LH123Ht$intercept)

BVP_2369 <- lines_PRNFT %>%
  dplyr::filter(rp_name == "2369") #intercept > 722
nrow(BVP_2369) 
min(BVP_2369$intercept)
max(BVP_2369$intercept)


var_BVP <- slope_int_mixed_model_NIL %>%
  dplyr::full_join(parents_PRNFT) %>%
  dplyr::mutate(dp = as.factor(dp), 
                rp = as.factor(rp), 
                dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                genotype_to_color = as.factor(genotype_to_color),
                genotype_to_color = factor(genotype_to_color, levels = list)) %>%
  dplyr::filter(is.na(genotype_to_color)) %>% # remove the rp and dp 
  dplyr::group_by(rp) %>%
  dplyr::summarise(N = length(unique(entryid2)),
                   mcount   = sum(intercept)/N,
                   varBPV = sum((intercept - mcount)^2, (N - n())*mcount^2)/(N - 1), 
                   mean_BVP = mean(intercept, na.rm =T), 
                   sd_intercept = sd(intercept, na.rm = T), 
                   Coeff_var = sd_intercept/ abs(mean_BVP))

var_SLOPE <- slope_int_mixed_model_NIL %>%
  dplyr::full_join(parents_PRNFT) %>%
  dplyr::mutate(dp = as.factor(dp), 
                rp = as.factor(rp), 
                dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7")), 
                genotype_to_color = as.factor(genotype_to_color),
                genotype_to_color = factor(genotype_to_color, levels = list), 
                slope = as.numeric(slope)) %>%
  dplyr::filter(is.na(genotype_to_color)) %>% # remove the rp and dp 
  dplyr::group_by(rp) %>%
  dplyr::summarise(N = length(unique(entryid2)),
                   mcount   = sum(slope)/N,
                   varslope = sum((slope - mcount)^2, (N - n())*mcount^2)/(N - 1), 
                   mean_slope = mean(slope, na.rm =T), 
                   sd_slope = sd(slope, na.rm = T), 
                   Coeff_var = sd_slope/ abs(mean_slope))

# ## plot photoperiod across MET 
# 
# plot_MET_photoperiod <- blues_dta_em %>%
#   ggplot()+
#   geom_histogram(aes(x= daylength_ti), fill = "#0072B2") +
#   theme(plot.background = element_blank()
#         #,legend.position = "bottom"#
#         ,legend.position = "none"
#         ,legend.direction="vertical"#,legend.direction="vertical"
#         ,legend.key.size = unit(1, "cm")
#         ,legend.text = element_text(size = 10, face="bold")#, family="Courier New Bold")
#         ,aspect.ratio = 1
#         #,legend.title = element_text(colour="blue", size=16, face="bold")
#         #,legend.text = element_text(colour="blue", size=16, face="bold")
#         #,panel.grid.major = element_blank()
#         ,panel.grid.minor = element_blank()
#         ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
#         ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
#         #,panel.spacing = unit(1, "lines")
#         ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
#         ,panel.background = element_blank()
#         #,strip.background = element_blank()
#         #,strip.text = element_blank()
#         ,plot.title = element_text(face="bold")
#         ,axis.title.x = element_text(size=18, margin = margin(t=18), face="bold")
#         ,axis.text.x = element_text(size=18, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) + xlab("photoperiod experienced at BVP-ti")
# 


#general_fig_path <- ("~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/NILAS/main_folder/drouault/paper/Figures/")
# ggsave(paste0(general_fig_path, "photoperiod_MET","_", date_today,".pdf"), plot_MET_photoperiod)
