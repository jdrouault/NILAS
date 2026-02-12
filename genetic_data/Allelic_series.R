#### This script describe the allelic series presents at ZmPR1-4

library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggpmisc)
library(agricolae)
library(multcomp)
library(multcompView)
library(emmeans)
library(pwr)


### date today 
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)

###### ======= Create paths 
root <- ("~/Dropbox-SE Dropbox/allshare/") #C:/Users/drouault
root_alternative <-"~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/" #C:/Users/drouault
NILAS_genotypes_path <- paste0(root, "stage_eloise/NILAS_haploimputer/Data/dp6q30dr0.40maf0.01/outputs/")
haplot_imput_data <- paste0(root, "stage_eloise/NILAS_haploimputer/Data/")
save_path <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Tables/genotypic_data/consistency_replicates/")
phenotype_path <- paste0(root_alternative, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_day_dta_04Avr11/") ### change to access BLUes_dta_em or BLUEs_dta
fig_consistency <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/constistency_replicates/")
fig_genes_introgression <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/effect_genes_introgressionNB/")
fig_phenotype <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/phenotypic_data/")
file_path <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Files/")
genetic_data_save_path <- paste0(root_alternative, "NILAS/main_folder/drouault/data/genotypic_data/NILs/")
table_phenotype <- paste0(root_alternative, "/NILAS/main_folder/drouault/paper/Tables/phenotypic_data/")
mixed_model_PRNTFTP <- paste0(root_alternative, "NILAS/main_folder/drouault/data/cleaned_data/from_mixed_model/")
file_supp_path <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Files/")

############ =============== Load the data 

# load BLUEs data with dta_em 

blues_dta_em <- read.csv(paste0(phenotype_path, "haplotype_blues_20240723.csv"), header = T) %>% #haplotypes_effect_gdd_dta_em_env_20240627.csv
  dplyr::mutate(qtl = as.character(qtl))

# load genetic data selected NILs 

summary_qtl_insertion <- read.csv(paste0(genetic_data_save_path, "20240702Selected_NILs_after_filtration.csv"), header = T)  %>%
  dplyr::mutate(qtl = as.character(qtl))

# Merge genetic data selected NILs and BLUEs 

blues_dta_em_summary <- summary_qtl_insertion %>% # filter on the selected lines from biol rep 
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(blues_dta_em) %>%
  dplyr::filter(entryid2 %in% unique(summary_qtl_insertion$entryid))

# PRNFT genotype filtered 

haplotype_slope <- read.csv(paste0(mixed_model_PRNTFTP, "haplotype_slope_20240928.csv"), header = T) %>%
  dplyr::rename(indi.cof = Estimate, entryid2 = haplotype, Se_slope = Std..Error, tvalue_slope = t.value, pval_slope = Pr...t.., BH_slope = BH)

haplotype_int <- read.csv(paste0(mixed_model_PRNTFTP, "haplotype_int_20240928.csv"), header = T) %>%
  dplyr::rename(intercept = Estimate, entryid2 = haplotype, Se_int = Std..Error, tvalue_int = t.value, pval_int = Pr...t.., BH_int = BH)

# phenotypic_data <- read.csv(paste0(phenotype_path, "haplotype_slope_dta_em_LD_SD_lmregression_20240723.csv"), header = T, sep = ",") %>% # haplotypes_effect_gdd_dta_em_env_20240508.csv
#   dplyr::mutate(qtl = as.character(qtl))

# sensitivity_phenotype_data <- read.csv(paste0(phenotype_path, "haplotype_slope_BVP_dta_em_20240508.csv"), header = T, sep = ",") %>%
#   dplyr::mutate(qtl = as.character(qtl))


dp_name<- c("NA","CML10", "CML258" , "CML277", "CML341" , "CML373" , "Tzi8" , "Tzi9" )
dp <- 0:7
mini_table <- cbind(dp_name, dp)
mini_table <- as.data.frame(mini_table)
mini_table$dp <- as.numeric(mini_table$dp)

rp_name<- c("NA","2369", "LH123Ht" )
rp <- 0:2
mini_table2 <- cbind(rp_name, rp)
mini_table2 <- as.data.frame(mini_table2)
mini_table2$rp <- as.numeric(mini_table2$rp)

phenotype_genotype <- summary_qtl_insertion %>%
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(haplotype_slope) %>%
  dplyr::full_join(haplotype_int) %>%
  dplyr::filter(entryid2 %in% unique(summary_qtl_insertion$entryid)) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid)), 
                qtl = as.numeric(gsub(".*q|g.*", "", entryid))) %>%
  dplyr::full_join(mini_table) %>%
  dplyr::full_join(mini_table2) 

# phenotype_genotype <- read.csv(paste0(phenotype_path,"20240619PRNFT_genotype_filtered_NILs.csv"), header = T)

### Make file S7
file_s7 <- phenotype_genotype %>%
filter_at(vars(zmpr1_cct, zmpr1_elf4, zmpr2_zcn8, zmpr2_rap2.7, zmpr3_conz1, zmpr3_cct, zmpr4_cacta), any_vars(. == 1)) %>%
  dplyr::select(Sample:RP_FG,Nb_ZmPR, Nb_ZmPR1_peak:Nb_ZmPR4_peak,Nb_ZmPR_peak, zmpr1_cct:zmpr4_cacta)


# write.csv(file_s7, paste0(file_supp_path, "File_S7", "_",date_today, ".csv"), row.names = F) 


############################################################ Detection allelic serie ############################################################ 


### === make the required tables 

dp_name <- phenotype_genotype %>%
  dplyr::select(dp_name, dp) %>%
  dplyr::mutate( dp = as.factor(dp)) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(dp_name) & dp !=0)


rp_name <- phenotype_genotype %>%
  dplyr::select(rp_name, rp) %>%
  # dplyr::mutate( rp = as.factor(rp)) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(rp_name) & rp !=0)

colnames_list <- c("dp", "rp" , "estimate", "statistic" ,"p.value" ,"parameter", "conf.low", 
                  "conf.high")
list <- c(1:4)

### ======== Dectection of allelic serie on slope 

list_genes <- dput(colnames(phenotype_genotype[,c(26:32)]))

# Make required dataset 

data2 <- NULL

for(i in 1:length(list_genes)){
   i <- 7
 i <- list_genes[[i]]
  print(paste0("ZmPR", i))
  ZmPR <- phenotype_genotype  %>%  #blues_dta_em_summary #phenotype_genotype
    # dplyr::filter(qtl == i) %>%
    dplyr::select(entryid2,dp, rp, dp_name, rp_name, qtl,i , indi.cof) #indi.cof, intercept #diff_GDD_dta_em, intercept # env, daylength_ti #paste0("Nb_ZmPR", i, "_peak")
  
  pwc2_ZmPR <-ZmPR %>%
    dplyr::filter(if_any(matches(i), ~ . == 1)) %>% #"^Nb_ZmPR[0-9]{1}_peak$"
    dplyr::mutate(dp = as.factor(dp)) %>%
    dplyr::filter(!grepl("i12", entryid2)) 
  
t_test_results_ZmPR <- pwc2_ZmPR %>% 
  dplyr::filter(!grepl("i12", entryid2)) %>%
  tidyr::unite(ID, c("dp", "rp"), sep = "_", remove = F) %>%
  dplyr::mutate(dp = as.numeric(dp), 
                rp = as.numeric(rp)) %>%
  split(.$ID) %>%
  map(function(x) if (length(x$indi.cof) > 1) x %>% group_by(dp, rp ) %>% do(                               
    t.test(x$indi.cof, mu = 0) %>%  
    broom::tidy()) %>% ungroup() else #data <- c(unique(x$dp), unique(x$rp), rep(NA, length(colnames_list) - 2)), #select(-method, -alternative)
    tibble::as_tibble(matrix(c(unique(x$dp),unique(x$rp), rep(NA, length(colnames_list) -2)), nrow = 1, ncol = length(colnames_list)), .name_repair = ~ colnames_list)) %>%
  bind_rows()

t_test_results_ZmPR$Spearman_p_sign <- "ns"
t_test_results_ZmPR$Spearman_p_sign[is.na(t_test_results_ZmPR$p.value)] <- "NA" 
t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.05] <- "*"
t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.01] <- "**"
t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.001] <- "***"
t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.0001] <- "****"

t_test_results_ZmPR <- t_test_results_ZmPR %>%
  tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) 

### === Parametric test performed with interaction between dp*rp 

lm_mod_ZmPR <- lm(indi.cof ~ dp * rp , data = pwc2_ZmPR)
anova_ZmPR  <- aov(lm_mod_ZmPR) 
data <- tidy(aov(lm_mod_ZmPR))
data$genes <- i
data2 <- rbind.data.frame(data, data2) 

vcov(lm_mod_ZmPR)

# Perform post-hoc tests using emmeans
emm_ZmPR <- emmeans(lm_mod_ZmPR, ~ dp | rp)  # Estimate marginal means for dp within each level of rp
summary(emm_ZmPR)

summary(emm_ZmPR, infer = c(TRUE, TRUE))


cld_1 <- cld(emm_ZmPR, Letters = letters)
cld_df <- unnest(cld_1)

### Result of the parametric test with dp*rp

min <- min(cld_df$emmean, na.rm =T) - cld_df$SE[which(cld_df$lower.CL ==min(cld_df$lower.CL, na.rm =T))] -1
max <-  max(cld_df$upper.CL, na.rm =T) + cld_df$SE[which(cld_df$upper.CL ==max(cld_df$upper.CL, na.rm =T))] + 3
letter_position <-round(max(cld_df$upper.CL, na.rm =T) + cld_df$SE[which(cld_df$upper.CL ==max(cld_df$upper.CL, na.rm =T))] + 1.5, 1)


plot_peak_ZmPR <- pwc2_ZmPR %>%
  tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) %>%
  dplyr::mutate(dp = as.numeric(dp)) %>%
  dplyr::full_join(t_test_results_ZmPR) %>% #changer ici
  dplyr::mutate(dp = as.factor(dp)) %>%
  # dplyr::full_join(dp_name) %>%
  # dplyr::full_join(rp_name) %>%
  dplyr::full_join(cld_df) %>%
  dplyr::filter(!grepl("i12", entryid2)) %>%
  dplyr::mutate( dp = as.factor(dp), 
                 dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
  dplyr::arrange(dp) %>%
    dplyr::filter(!is.na(ID)) %>%
  ggplot() +
  scale_color_manual(values = c("black")) +
  # ggplot(., aes(x = dp, y = emmean, color = rp)) +
  geom_point(aes(x = dp_name, y = emmean), position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(x = dp_name, ymin = emmean - SE, ymax = emmean + SE),
                width = 0.2, position = position_dodge(width = 0.5)) +
  # geom_boxplot(aes(x = dp_name, y = intercept  )) + #fill = "#B276B2"
  # geom_point(aes(x = dp_name, y = intercept)) +
  geom_text(aes(label = .group, x = dp_name, y = letter_position ), #y = emmea,
  position = position_dodge(width = 0.5), vjust = 0, size = 10) +
 # geom_text(aes(x = dp_name, y = letter_position - 1.5, label = Spearman_p_sign), color = "black", size = 8, position = position_dodge(width = .8), size = 6) +
  facet_wrap(~rp_name) +
  # ggplot() +
  # scale_color_manual(values = c("black")) +
  # geom_boxplot(aes(x = dp_name, y = indi.cof)) + #fill = "#B276B2" #indi.cof
  # geom_point(aes(x = dp_name, y = indi.cof )) + #indi.cof
  # geom_text(aes(label = .group, x = dp_name, y = letter_position), #y = emmea, #max(indi.cof)
            # position = position_dodge(width = 0.5), vjust = 0, size = 10) +
  ylim(min,max) +
  # geom_text(aes(x = dp_name, y = letter_position - 2, label = Spearman_p_sign), color = "black", position = position_dodge(width = .8), size = 6) +
  facet_wrap(~rp_name) + 
  theme_classic() +
  labs(
    x = "donor parents",
    y = (quote(bold("Haplotype effect on PS") ~ bold("(") ~bold(GDD[hour]) ~ bold(")")))
, title = paste0("ZmPR", gsub("zmpr", "", i))) + #
  theme(plot.background = element_blank()
        #,legend.position = "bottom"#
        ,legend.position = "none"
        ,legend.direction="vertical"#,legend.direction="vertical"
        ,legend.key.size = unit(1, "cm")
        ,legend.text = element_text(size = 30, face="bold")#, family="Courier New Bold")
        ,strip.text = element_text(size = 30)
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
        ,plot.title = element_text(size = 30, face="bold")
        ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
        ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(r=30), face="bold")
        ,axis.text.y = element_text(size=30) #30 instead of 18
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) 

plot_peak_ZmPR 
ggsave(paste0(fig_genes_introgression, "allelicserie_slope_effect",i,"_", date_today,".pdf"), plot_peak_ZmPR, width = 12, height = 9)

}

# write.csv(data2, paste0(table_phenotype, "p_value_test_anova_slope", "_",date_today, ".csv"), row.names = F)

ZmPR4 <- blues_dta_em_summary %>% #sensitivity_genotype #phenotype_genotype
  dplyr::filter(qtl == 4) %>%
  dplyr::select(entryid2,dp, rp, qtl, diff_gdd, Nb_ZmPR4_peak, env, daylength_TI) #diff_GDD_dta_em, intercept

pwc2_ZmPR4 <-ZmPR4 %>%
  dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
  dplyr::filter(Nb_ZmPR4_peak == 1
                & parameters == "LD") %>%
  dplyr::mutate(dp = as.factor(dp))

summary_entryid_env <- ZmPR4 %>%
  dplyr::filter(!grepl("i12", entryid2)) %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) 


pwc2_ZmPR
negative_slope_ZmPR4 <- phenotype_genotype %>%
  dplyr::filter(entryid2 == "NILASq1g11i04" | entryid2 == "NILASq1g11i05" |
                  entryid2 == "NILASq1g11i07" )

negative_slope_ZmPR4 <- blues_dta_em_summary %>%
  dplyr::filter(entryid2 == "NILASq1g11i04" | entryid2 == "NILASq1g11i05" |
                  entryid2 == "NILASq1g11i07" )

plot_negative_slope_ZmPR4 <- negative_slope_ZmPR4 %>%
  dplyr::mutate(Parameters = ifelse(daylength_ti < 14.8 , "SD", "LD"))%>% # 14.8 mean threshold 7 parents 
  dplyr::group_by(entryid2, Parameters) %>%
  dplyr::mutate(entryid2 = as.character(entryid2), 
                mean_daylength= mean(daylength_ti, na.rm =T)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::mutate(min_daylength = min(mean_daylength), 
                new_mean_daylength = mean_daylength - min_daylength)%>%
  ggplot(aes(x = new_mean_daylength , y = diff_GDD_dta_em, color = entryid2)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)

# ggsave(paste0(fig_phenotype, "Negative_slopes_CML10_LH123Ht","_", date_today,".pdf"), plot_negative_slope_ZmPR4, width = 6, height = 5)

### ======== Detection of allelic serie on BVP 

list_genes <- dput(colnames(phenotype_genotype[,c(26:32)]))

#make the required data
data1 <- NULL 

for(i in 1:length(list_genes)){
  # i <- 1
  i <- list_genes[[i]]
  print(paste0("ZmPR", i))
  ZmPR <- phenotype_genotype %>% #sensitivity_genotype #blues_dta_em_summary 
    # dplyr::filter(qtl == i) %>%
    dplyr::select(entryid2,dp, rp, dp_name, rp_name, qtl, i, indi.cof, intercept) #diff_GDD_dta_em, intercept # env, daylength_ti #paste0("Nb_ZmPR", i, "_peak")
  
  pwc2_ZmPR <-ZmPR %>%
    dplyr::filter(across(matches(i), ~ . == 1)) %>% #"^Nb_ZmPR[0-9]{1}_peak$"
    dplyr::mutate(dp = as.factor(dp)) %>%
    dplyr::filter(!grepl("i12", entryid2))
  
  t_test_results_ZmPR <- pwc2_ZmPR %>% 
    dplyr::filter(!grepl("i12", entryid2)) %>%
    tidyr::unite(ID, c("dp", "rp"), sep = "_", remove = F) %>%
    dplyr::mutate(dp = as.numeric(dp), 
                  rp = as.numeric(rp)) %>%
    split(.$ID) %>%
    map(function(x) if (length(x$intercept) > 1) x %>% group_by(dp, rp ) %>% do(                               
      t.test(x$intercept, mu = 0) %>%  
        broom::tidy()) %>% ungroup() else #data <- c(unique(x$dp), unique(x$rp), rep(NA, length(colnames_list) - 2)), #select(-method, -alternative)
          tibble::as_tibble(matrix(c(unique(x$dp),unique(x$rp), rep(NA, length(colnames_list) -2)), nrow = 1, ncol = length(colnames_list)), .name_repair = ~ colnames_list)) %>%
    bind_rows()
  
  t_test_results_ZmPR$Spearman_p_sign <- "ns"
  t_test_results_ZmPR$Spearman_p_sign[is.na(t_test_results_ZmPR$p.value)] <- "NA" 
  t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.05] <- "*"
  t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.01] <- "**"
  t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.001] <- "***"
  t_test_results_ZmPR$Spearman_p_sign[t_test_results_ZmPR$p.value < 0.0001] <- "****"
  
  t_test_results_ZmPR <- t_test_results_ZmPR %>%
    tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) 
  
  ### === Parametric test performed with interaction between dp*rp 
  
  lm_mod_ZmPR <- lm(intercept ~ dp * rp , data = pwc2_ZmPR)
  anova_ZmPR  <- aov(lm_mod_ZmPR) 
  data <- tidy(aov(lm_mod_ZmPR))
  data$genes <- i
  data1 <- rbind.data.frame(data, data1)
 
  
  # Perform post-hoc tests using emmeans
  emm_ZmPR <- emmeans(lm_mod_ZmPR, ~ dp * rp)  # Estimate marginal means for dp within each level of rp
  summary(emm_ZmPR)
  
  # add letters
  cld_1 <- cld(emm_ZmPR, Letters = letters, level = 0.01)
  cld_df <- unnest(cld_1)
  
  ### Result of the parametric test with dp*rp
  
  min <- min(cld_df$emmean, na.rm =T) - cld_df$SE[which(cld_df$lower.CL ==min(cld_df$lower.CL, na.rm =T))] - 5
  max <-  max(cld_df$upper.CL, na.rm =T) + cld_df$SE[which(cld_df$upper.CL ==max(cld_df$upper.CL, na.rm =T))] + 5
  letter_position <-round(max(cld_df$upper.CL, na.rm =T) + cld_df$SE[which(cld_df$upper.CL ==max(cld_df$upper.CL, na.rm =T))] + 4, 1)
  
  plot_peak_ZmPR <- cld_df  %>% #pwc2_ZmPR
    tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) %>%
    dplyr::mutate(dp = as.numeric(dp)) %>%
    dplyr::full_join(t_test_results_ZmPR) %>% #changer ici
    dplyr::mutate(dp = as.factor(dp)) %>%
    dplyr::full_join(dp_name) %>%
    dplyr::full_join(rp_name) %>%
    dplyr::filter(!is.na(ID)) %>%
    # dplyr::full_join(cld_df) %>%
    dplyr::mutate(dp = as.factor(dp)) %>%
    # dplyr::filter(!grepl("i12", entryid2)) %>%
    dplyr::mutate( dp = as.factor(dp), 
                   dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
    dplyr::arrange(dp) %>%
    # dplyr::filter(!is.na(ID) & !is.na(entryid2)) %>%
    ggplot() +
    scale_color_manual(values = c("black")) +
    # ggplot(., aes(x = dp, y = emmean, color = rp)) +
    geom_point(aes(x = dp_name, y = emmean), position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(x = dp_name, ymin = emmean - SE, ymax = emmean + SE),
                    width = 0.2, position = position_dodge(width = 0.5)) +
    # geom_boxplot(aes(x = dp_name, y = intercept  )) + #fill = "#B276B2"
    # geom_point(aes(x = dp_name, y = intercept)) +
    geom_text(aes(label = .group, x = dp_name, y = letter_position ), #y = emmea,
              position = position_dodge(width = 0.5), vjust = 0, size = 10) +
    # geom_text(aes(x = dp_name, y = letter_position - 3, label = Spearman_p_sign), color = "black", size = 8, position = position_dodge(width = .8), size = 6 ) +
    facet_wrap(~rp_name) +
    ylim(min,max) +
    theme_classic() +
    theme(plot.background = element_blank()
          #,legend.position = "bottom"#
          ,legend.position = "none"
          ,legend.direction="vertical"#,legend.direction="vertical"
          ,legend.key.size = unit(1, "cm")
          ,legend.text = element_text(size = 30, face="bold")#, family="Courier New Bold")
          ,strip.text = element_text(size = 30)
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
          ,plot.title = element_text(size = 30, face="bold")
          ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
          ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
          ,axis.title.y = element_text(size=30, margin = margin(r=30), face="bold")
          ,axis.text.y = element_text(size=30) #30 instead of 18
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()
          # top, right, bottom, left
          ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
    ) + 
    labs(
      x = "donor parents",
      y = quote(bold(atop(bold("Haplotype effect on"), bold("flowering time" ~ bolditalic("per se") ~ "(GDD)"))))
, title = paste0("ZmPR", gsub("zmpr", "", i)))
  
  plot_peak_ZmPR 
  ggsave(paste0(fig_genes_introgression, "allelicserie_BVP",i,"_", date_today,".pdf"), plot_peak_ZmPR, width = 12, height = 9)
  
}

# write.csv(data1, paste0(table_phenotype, "p_value_test_anova_BVP", "_",date_today, ".csv"), row.names = F)

############ =============== Heat map

##### Heat Map ZmPR4 expression across environments 

blues_dta_em_summary <- blues_dta_em_summary %>%
  dplyr::mutate(dp = as.factor(dp)) %>%
  dplyr::full_join(dp_name) %>%
  dplyr::full_join(rp_name) %>%
  dplyr::mutate(env = as.factor(env), 
                env = factor(env, levels = c("2015_PR", "2015_FL",   "2016_FL",  "2015_cTX", "2016_cTX" ,"2015_nTX", "2016_nTX", "2015_NC", "2016_NC", 
                                             "2015_DE" , "2016_DE"  ,"2015_IA" , "2016_IA" , "2015_WI" , "2016_WI" )), 
                dp_name = as.factor(dp_name),
                rp_name = as.factor(rp_name), 
                dp_name = factor(dp_name, levels = c("CML10", "CML258", "CML277", 
                                                     "CML341", "CML373", "Tzi8", "Tzi9")))

 # i <- 4

list <- c(1:4)

for (i in 1:length(list)){
  i <- 2
  x <-  subset(blues_dta_em_summary, qtl == i)
  name <- paste0("Nb_ZmPR",i,"_peak")
  x <- x[,c("entryid2","dp", "rp", "dp_name", "rp_name" , "qtl", "diff_GDD_dta_em", name ,"diff_GDD_dta_em", "env", "daylength_ti")] #diff_GDD_dta_em, intercept
  
  x <- subset(x, x[,8] == 1)

  heatmap_CACTA <- x %>%
    # dplyr::mutate(dp = as.factor(dp)) %>%
    # dplyr::full_join(dp_name) %>%
    # dplyr::full_join(rp_name) %>%
    # dplyr::filter(paste0("Nb_ZmPR",1,"_peak") == 1) %>%
    dplyr::filter(!grepl("i12", entryid2)) %>%
    # dplyr::mutate(env = as.factor(env), 
    #               env = factor(env, levels = c("2015_PR", "2015_FL",   "2016_FL",  "2015_cTX", "2016_cTX" ,"2015_nTX", "2016_nTX", "2015_NC", "2016_NC", 
    #                                            "2015_DE" , "2016_DE"  ,"2015_IA" , "2016_IA" , "2015_WI" , "2016_WI" )), 
    #               dp_name = as.factor(dp_name),
    #               rp_name = as.factor(rp_name), 
    #               dp_name = factor(dp_name, levels = c("CML10", "CML258", "CML277", 
    #                                                    "CML341", "CML373", "Tzi8", "Tzi9")))%>%
    dplyr::group_by(dp_name,rp_name, env) %>%
    dplyr::mutate(median_diff_GDD_dta_em = median(diff_GDD_dta_em, na.rm =T) ) %>%
    ggplot(aes(env, dp_name, fill= median_diff_GDD_dta_em)) + 
    scale_fill_gradientn(colors = c( "blue", "white", "red"), aes(limits = range(median_diff_GDD_dta_em)),
                        limits=c(-100,170),breaks = seq(-100, 170, by = 20), na.value = "grey50", name = "haplotype effect") +
    geom_tile() +
    theme_classic() +
    theme(plot.background = element_blank(),
          #,legend.position = "bottom"#
          ,legend.position = "right"
          ,legend.direction="vertical"#,legend.direction="vertical"
          ,legend.key.size = unit(5, "cm")
          ,legend.key.height = unit(2, 'cm') #change legend key height
          ,legend.key.width = unit(1, 'cm') #change legend key width
          ,legend.text = element_text(size = 14, face="bold")#, family="Courier New Bold")
          ,aspect.ratio = 1
          ,legend.title = element_text(colour="black", size=16, face="bold")
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
          ,axis.title.x = element_text(size=18, margin = margin(t=18), face="bold")
          ,axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, margin = margin(t=11))
          ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
          ,axis.text.y = element_text(size=18)
          ,axis.ticks.y = element_blank()
          ,axis.ticks.x = element_blank()
          ,strip.text = element_text(size = 14, face="bold") 
          # top, right, bottom, left
          ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
    ) + geom_vline(aes(xintercept = "2016_cTX"), color = "black") +
    facet_wrap(~rp_name) + ylab("donor parents") + xlab ("environments") + labs(fill = "haplotype effect \\n (GDD)", title = name) 
  heatmap_CACTA
  
  ggsave(paste0(fig_genes_introgression, "heatmap_ZmPR",i,"_", date_today,".pdf"), heatmap_CACTA, width = 12, height = 6)
  
}

############ =============== PCA genes clustering 

library(devtools)
library(ggbiplot)
library(missMDA)
library("FactoMineR")


pivoted_PCA <- phenotype_genotype %>%
  dplyr::select(entryid2, Nb_ZmPR1_peak:Nb_ZmPR4_peak, indi.cof, intercept) %>%
  dplyr::filter(if_any(matches("^Nb_ZmPR[0-9]_peak$"), ~ . == 1)) %>% 
  tidyr::pivot_longer(col = Nb_ZmPR1_peak:Nb_ZmPR4_peak, 
                      names_to = "genes", 
                      values_to = "peak") %>%
  dplyr::filter(peak == 1) %>%
  dplyr::select(-peak) %>%
  tidyr::pivot_wider(
                     names_from = genes, 
                     values_from = c(indi.cof, intercept)) %>%
  dplyr::select(-entryid2) %>%
  tidyr::pivot_longer(col = indi.cof_Nb_ZmPR1_peak:intercept_Nb_ZmPR4_peak, 
                      names_to = c("parameter", ".value"),
                      names_sep = "_Nb_",
                      values_to = "value") 

pivoted_PCA_slope <- pivoted_PCA %>%
  dplyr::filter(parameter == "indi.cof") %>%
  dplyr::parameter(-parameter)

pivoted_PCA_intercept <- pivoted_PCA %>%
  dplyr::filter(parameter == "intercept") %>%
  dplyr::select(-parameter)
  
#nb <- estim_ncpPCA(pivoted_PCA_slope[,-1],method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data

### PCA slope 
res.comp <- imputePCA(pivoted_PCA_slope[,-1], ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
res.pca <- PCA(res.comp$completeObs, ncp =5, graph=FALSE, axes = c(1,2))
plot(res.pca, choix="var")


### PCA intercept 
res.comp <- imputePCA(pivoted_PCA_intercept, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
res.pca <- PCA(res.comp$completeObs, ncp =10, graph=FALSE, axes = c(1,2))
plot(res.pca, choix="var")


# pivoted_PCA <- pivoted_PCA[,-1]
# pivoted_PCA <-na.omit(pivoted_PCA)
  
pc <- princomp(pivoted_PCA_slope,
             center = TRUE,
             scale. = TRUE, 
             na.action = na.exclude())
attributes(pc)


g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = pivoted_PCA$genes,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)


# ############ =============== Focus on ZmCCT10 allelic effect - ZmPR4
# 
# cacta_ZmPR4 <- blues_dta_em_summary %>% #sensitivity_genotype #phenotype_genotype
#   dplyr::filter(qtl == 4) %>%
#   dplyr::select(entryid2,dp, rp, qtl, diff_GDD_dta_em, Nb_ZmPR4_peak, diff_GDD_dta_em, env, daylength_ti) #diff_GDD_dta_em, intercept
# 
# pwc2 <-  cacta_ZmPR4 %>%
#   dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1
#                 & parameters == "LD") %>%
#   dplyr::mutate(dp = as.factor(dp))
# 
# summary_entryid_env <- cacta_ZmPR4 %>%
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   as.data.frame() %>%
#   dplyr::group_by(entryid2) %>%
#   dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) 
# 

# For ZmCCT10 it is known that CML277 and CML341 carry - cacta so that expression ZmCCT10 is not impeded and should delay flowering time, 
#so we compare the effect of the alleles of CML277 and CML341 on the slopes 


# ##### ==== Test if the diff_gdd_dta differs from 0, per dp (haplotype effect on PS)
# 
# #### Verification of the assumptions to choose between a parametric or non parametric test 
# 
# # Check for the assumption of Normality
# shapiro_result <- pwc2_ZmPR1 %>% 
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::group_by(dp, rp) %>% 
#   do(                               
#     shapiro.test(.$diff_GDD_dta_em) %>%
#       broom::tidy()                        # tidy into tibble
#   ) #If the p-value from the Shapiro-Wilk test is greater than 0.05, 
# #the normality assumption is met, and we proceed with the t-test.
# 
# 
# 
# ### Result of the parametric test with dp*rp
# 
# plot_peak_ZmPR1 <- pwc2_ZmPR %>%
#   tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) %>%
#   dplyr::mutate(dp = as.factor(dp)) %>%
#   # dplyr::full_join(t_test_results_ZmPR1) %>% #changer ici 
#   # dplyr::mutate(dp = as.factor(dp)) %>%
#   dplyr::full_join(dp_name) %>%
#   dplyr::full_join(rp_name) %>%
#   dplyr::full_join(cld_df_ZmPR1) %>%
#   dplyr::mutate(dp = as.character(dp), 
#                 parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::filter(Nb_ZmPR1_peak == 1 & parameters == "LD" ) %>% 
#   # dplyr::mutate(cld = cld_df$rp$dp$Letters[dp])%>%
#   dplyr::mutate( dp = as.factor(dp), 
#                  dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
#   dplyr::arrange(dp) %>%
#   ggplot() +
#   scale_color_manual(values = c("black")) +
#   geom_boxplot(aes(x = dp_name, y = diff_GDD_dta_em  )) + #fill = "#B276B2"
#   geom_point(aes(x = dp_name, y = diff_GDD_dta_em, color = parameters )) +
#   # ggplot(cld_df, aes(x = dp, y = emmean, color = rp, group = rp)) +
#   # geom_point(position = position_dodge(width = 0.5), size = 3) +
#   # geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
#   #               width = 0.2, position = position_dodge(width = 0.5)) +
#   geom_text(aes(label = .group, x = dp_name, y = max(diff_GDD_dta_em) + 2), #y = emmea,
#             position = position_dodge(width = 0.5), vjust = 0, size = 6) +
#   # geom_text(aes(x = dp_name, y = 124, label = Spearman_p_sign), color = "black", position = position_dodge(width = .8), size = 6) +
#   facet_wrap(~rp_name) + 
#   theme_classic() +
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
#         ,axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) + 
#   labs(
#     x = "donor parents",
#     y = "Haplotype phenotypic effect")
# 
# 
# plot_peak_ZmPR1
# 
# ggsave(paste0(fig_genes_introgression, "alleles_ZmPR4_rp_groups","_", date_today,".pdf"), plot_cacta_peak_2, width = 10, height = 6)


# # Pairwise comparisons
# pairwise_comparisons_ZmPR <- pairs(emm_ZmPR)
# summary(pairwise_comparisons_ZmPR)
# 
# # Extract the estimated marginal means
# emm_df_ZmPR <- as.data.frame(emm_ZmPR)
# 
# # Create a dataframe for pairwise comparisons
# pairwise_df_ZmPR <- as.data.frame(pairwise_comparisons_ZmPR)
# cld_list <- list()
# for (a in unique(pwc2_ZmPR$rp)) {
#   data_subset <- subset(pwc2_ZmPR, rp == a)
#   lm_subset <- lm(indi.cof ~ dp, data = data_subset)
#   emm_subset <- emmeans(lm_subset, ~ dp)
#   cld_1 <- cld(emm_subset, Letters = letters)
#   cld_df <- as.data.frame(cld_1)
#   cld_df$rp <- a
#   cld_list[[a]] <- cld_df
# }
# 
# cld_df_ZmPR <- do.call(rbind, cld_list) #%>%
#   dplyr::full_join(dp_name)
# cld_df_ZmPR

# ##### ==== Test if the diff_gdd_dta differs from 0, per dp (haplotype effect on PS)
# 
# #### Verification of the assumptions to choose between a parametric or non parametric test 
# 
# # Check for the assumption of Normality
# shapiro_result <- cacta_ZmPR4 %>% 
#   dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1
#                 & parameters == "LD") %>% 
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::group_by(dp, rp) %>% 
#   do(                               
#     shapiro.test(.$diff_GDD_dta_em) %>%
#       broom::tidy()                        # tidy into tibble
#   ) #If the p-value from the Shapiro-Wilk test is greater than 0.05, 
# #the normality assumption is met, and we proceed with the t-test.
# 
# ### Power t-test 
# 
# effect_size <- (mean(pwc2$diff_GDD_dta_em[pwc2$rp == 1]) - mean(pwc2$diff_GDD_dta_em[pwc2$rp == 2])) / 
#   sd(pwc2$diff_GDD_dta_em)
# 
# t_test_power <- pwr.t.test(d = effect_size, sig.level = 0.05, power = NULL, n = length(pwc2$diff_GDD_dta_em[pwc2$rp== 1]), type = c("one.sample"),  alternative = c("two.sided"))
# print(t_test_power)
# 
# ### Wilcox test power
# 
# simulate_power_wilcox <- function(n_sim = 1000, n1 = 159, n2 = 189, mu1 = 19.922, mu2 = 6.206778, sd = 25.40136, alpha = 0.05) {
#   p_values <- replicate(n_sim, {
#     group1 <- rnorm(n1, mean = mu1, sd = sd)
#     group2 <- rnorm(n2, mean = mu2, sd = sd)
#     wilcox.test(group1, group2)$p.value
#   })
#   power <- mean(p_values < alpha)
#   return(power)
# }
# 
# # library(MKpower)
# # rxy <- function(n) rnorm(n, mean = mean(group1),sd = sd(group1))
# # sim.ssize.wilcox.test(rx = rxy, n.min =  10, n.max=  30, step.size = 1, iter = 1000, type = "one.sample")
# 
# 
# #### Parametric t-test (assumption following a normal law)
# 
# t_test_results_2 <- cacta_ZmPR4 %>% 
#   dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1
#                 & parameters == "LD") %>% 
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::group_by(dp, rp) %>% 
#   do(                               
#     t.test(.$diff_GDD_dta_em, mu = 0) %>%  # compare mean to of rp 1  
#       broom::tidy()                        # tidy into tibble
#   )
# 
# t_test_results_2$Spearman_p_sign <- "ns"
# t_test_results_2$Spearman_p_sign[t_test_results_2$p.value < 0.05] <- "*"
# t_test_results_2$Spearman_p_sign[t_test_results_2$p.value < 0.01] <- "**"
# t_test_results_2$Spearman_p_sign[t_test_results_2$p.value < 0.001] <- "***"
# t_test_results_2$Spearman_p_sign[t_test_results_2$p.value < 0.0001] <- "****"
# 
# t_test_results_2 <- t_test_results_2 %>%
#   tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) 
# 
# #### Non-parametric test Wilcoxon Rank-Sum Test (Unpaired) (assumption not following a normal law)
# 
# #An unpaired Wilcoxon rank-sum test (also known as the Mann-Whitney U test) is used when you have two independent samples.
# 
# Wilcoxon_test_results <- cacta_ZmPR4 %>% 
#   dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1
#                 & parameters == "LD") %>% 
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::group_by(dp, rp) %>% 
#   do(                               
#     wilcox.test(.$diff_GDD_dta_em, mu = 0) %>%  
#       broom::tidy()                        # tidy into tibble
#   )
# 
# 
# #### ==== Test if several group (>2) exist which means there is an allelic serie 
# 
# # test if alleles of CML277 and CML341 are significantly differents from the others alleles sources as they contains -CACTA
# 
# dp_name <- phenotype_genotype %>%
#   dplyr::select(dp_name, dp) %>%
#   dplyr::mutate( dp = as.factor(dp)) %>%
#   dplyr::distinct()
# 
# rp_name <- phenotype_genotype %>%
#   dplyr::select(rp_name, rp) %>%
#   # dplyr::mutate( rp = as.factor(rp)) %>%
#   dplyr::distinct()
# 
# ### parametric test
# 
# lm_mod <- lm(diff_GDD_dta_em ~ dp  , data = pwc2)
# anova  <- aov(lm_mod)
# 
# tukey.test <- TukeyHSD(aov(lm_mod))
# plot(tukey.test)
# cld <- multcompLetters4(anova, tukey.test)
# 
# ### Parametric test performed with interaction between dp*rp 
# lm_mod_2<- lm(diff_GDD_dta_em ~ dp * rp , data = pwc2)
# anova_2  <- aov(lm_mod_2) 
# 
# # Perform post-hoc tests using emmeans
# emm <- emmeans(lm_mod_2, ~ dp | rp)  # Estimate marginal means for dp within each level of rp
# summary(emm)
# 
# # Pairwise comparisons
# pairwise_comparisons <- pairs(emm)
# summary(pairwise_comparisons)
# 
# # Extract the estimated marginal means
# emm_df <- as.data.frame(emm)
# 
# # Create a dataframe for pairwise comparisons
# pairwise_df <- as.data.frame(pairwise_comparisons)
# cld_list <- list()
# for (i in unique(pwc2$rp)) {
#   # i <- 1
#   data_subset <- subset(pwc2, rp == i)
#   lm_subset <- lm(diff_GDD_dta_em ~ dp, data = data_subset)
#   emm_subset <- emmeans(lm_subset, ~ dp)
#   cld_1 <- cld(emm_subset, Letters = letters)
#   cld_df <- as.data.frame(cld_1)
#   cld_df$rp <- i
#   cld_list[[i]] <- cld_df
# }
# 
# cld_df <- do.call(rbind, cld_list) %>%
#   # dplyr::mutate(rp = as.factor(rp)) %>%
#   dplyr::full_join(dp_name)
# 
# 
# ### non parametric test 
# # bartlett.test(subset(cacta_ZmPR4, Nb_ZmPR4_peak == 1)$slope , (subset(cacta_ZmPR4, Nb_ZmPR4_peak == 1)$dp))
# # res.kruskal <- pwc2 %>% kruskal_test(slope ~ dp)
# # res.kruskal
# 
# ##### ==== plot the results 
# 
# ### Result of the parametric test with dp only
# 
# plot_cacta_peak <- pwc2 %>%
#   dplyr::full_join(dp_name) %>%
#   dplyr::mutate(dp = as.character(dp), 
#                 parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1 & parameters == "LD" ) %>% 
#   # dplyr::mutate(cld = cld$dp$Letters[dp])%>%
#   dplyr::mutate( dp = as.factor(dp), 
#                  dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
#   dplyr::arrange(dp) %>%
#   ggplot() +
#   theme_classic() +
#   geom_boxplot(aes(x = dp_name, y = diff_GDD_dta_em  )) + #fill = "#B276B2"
#   geom_point(aes(x = dp_name, y = diff_GDD_dta_em, color = parameters )) +
#   scale_color_manual(values = c("black")) +
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
#         ,axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) +  
#   geom_text( aes(x = dp_name, y = 115, label = cld), size = 6, check_overlap = TRUE) +
#   geom_text(data= t_test_results_2,  aes(x = dp, y =120, label = Spearman_p_sign), color = "black", position = position_dodge(width = .8), size = 6) +
#   xlab("donors parents") +
#   ylab("Haplotype effect of GDD to anthesis")
# 
# plot_cacta_peak
# 
# # ggsave(paste0(fig_genes_introgression, "alleles_ZmPR4","_", date_today,".pdf"), plot_cacta_peak, width = 10, height = 6)
# 
# ### Result of the parametric test with dp*rp
# 
# plot_cacta_peak_2 <- pwc2 %>%
#   tidyr::unite(ID, c('dp', 'rp'), sep = "_", remove = F) %>%
#   dplyr::mutate(dp = as.numeric(dp)) %>%
#   dplyr::full_join(t_test_results_2) %>%
#   dplyr::mutate(dp = as.factor(dp)) %>%
#   dplyr::full_join(dp_name) %>%
#   dplyr::full_join(rp_name) %>%
#   dplyr::full_join(cld_df) %>%
#   dplyr::mutate(dp = as.character(dp), 
#                 parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>%
#   dplyr::filter(!grepl("i12", entryid2)) %>%
#   dplyr::filter(Nb_ZmPR4_peak == 1 & parameters == "LD" ) %>% 
#   # dplyr::mutate(cld = cld_df$rp$dp$Letters[dp])%>%
#   dplyr::mutate( dp = as.factor(dp), 
#                  dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
#   dplyr::arrange(dp) %>%
#   ggplot() +
#   scale_color_manual(values = c("black")) +
#   geom_boxplot(aes(x = dp_name, y = diff_GDD_dta_em  )) + #fill = "#B276B2"
#   geom_point(aes(x = dp_name, y = diff_GDD_dta_em, color = parameters )) +
#   # ggplot(cld_df, aes(x = dp, y = emmean, color = rp, group = rp)) +
#   # geom_point(position = position_dodge(width = 0.5), size = 3) +
#   # geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
#   #               width = 0.2, position = position_dodge(width = 0.5)) +
#   geom_text(aes(label = .group, x = dp_name, y = max(diff_GDD_dta_em) + 2), #y = emmea,
#             position = position_dodge(width = 0.5), vjust = 0, size = 6) +
#   geom_text(aes(x = dp_name, y = 124, label = Spearman_p_sign), color = "black", position = position_dodge(width = .8), size = 6) +
#   facet_wrap(~rp_name) + 
#   theme_classic() +
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
#         ,axis.text.x = element_text(size=18, angle = 90, vjust = 0.5, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) + 
#   labs(
#     x = "donor parents",
#     y = "Haplotype phenotypic effect") 
# 
# 
# plot_cacta_peak_2
# 
# # ggsave(paste0(fig_genes_introgression, "alleles_ZmPR4_rp_groups","_", date_today,".pdf"), plot_cacta_peak_2, width = 10, height = 6)
# 
# ############ =============== Focus on ... allelic effect - ZmPR1
# 
# 
# for(i in 1:4){
#   
#   ############ ===============  Allelic serie associated to the slope
#   
#   paste0("ZmPR", i)
#   ZmPR <- phenotype_genotype %>% #sensitivity_genotype #blues_dta_em_summary 
#     dplyr::filter(qtl == i) %>%
#     dplyr::select(entryid2,dp, rp, dp_name, rp_name, qtl, paste0("Nb_ZmPR", i, "_peak"), indi.cof, intercept) #diff_GDD_dta_em, intercept # env, daylength_ti
#   
#   pwc2_ZmPR <-ZmPR %>%
#     # dplyr::mutate(parameters = ifelse(daylength_ti < 14.8 , "SD", "LD")) %>% # no need to do that when testing the slope or intercept because is already related to SD and LD
#     dplyr::filter(across(matches("^Nb_ZmPR[0-9]{1}_peak$"), ~ . == 1)) %>%
#     # dplyr::filter( parameters == "LD") %>% ### working with LD which corresponds to the slope 
#     dplyr::mutate(dp = as.factor(dp))
#   
#   # summary_entryid_env <- ZmPR %>%
#   #   dplyr::filter(!grepl("i12", entryid2)) %>%
#   #   as.data.frame() %>%
#   #   dplyr::group_by(entryid2) %>%
#   #   dplyr::summarise_at(vars(env), ~ sum(!(is.na(.))))
#   
#   ##### ==== Test if the diff_gdd_dta differs from 0, per dp (haplotype effect on PS)
#   
#   ### === Verification of the assumptions to choose between a parametric or non parametric test 
#   
#   # Check for the assumption of Normality
#   
#   shapiro_result <- pwc2_ZmPR %>% 
#     dplyr::filter(!grepl("i12", entryid2)) %>%
#     dplyr::group_by(dp, rp) %>% 
#     do(                               
#       shapiro.test(.$diff_GDD_dta_em) %>%
#         broom::tidy()                        # tidy into tibble
#     ) #If the p-value from the Shapiro-Wilk test is greater than 0.05, 
#   #the normality assumption is met, and we proceed with the t-test.
#   
#   shapiro_result_above_pvalue <- shapiro_result %>%
#     dplyr::filter(p.value > 0.05)
#   
#   ### === Fint a parametric test or non parametric test according to the results of the shapiro test 
#   
#   if(nrow(shapiro_result_above_pvalue) / nrow(shapiro_result) > 0.5) {
#     print("t-test")
#     
#     
#   }else{
#     print("Wilcoxon test")
#     
#     #### Non-parametric test Wilcoxon Rank-Sum Test (Unpaired) (assumption not following a normal law)
#     
#     #An unpaired Wilcoxon rank-sum test (also known as the Mann-Whitney U test) is used when you have two independent samples.
#     
#     Wilcoxon_test_results <-  pwc2_ZmPR %>% 
#       dplyr::filter(!grepl("i12", entryid2)) %>%
#       dplyr::group_by(dp, rp) %>% 
#       do(                               
#         wilcox.test(.$diff_GDD_dta_em, mu = 0) %>%  
#           broom::tidy()                        # tidy into tibble
#       )
#     
#     }
#   
# }






# 
# #### old plot 
# 
# p_values <-cacta_ZmPR4 %>%
#   group_by(dp) %>%
#   group_map(~ t.test(slope ~ Nb_ZmPR4_peak, .x)$p.value, .keep = TRUE) %>%
#   map_dfr(~ tibble(p_value = .), .id = "dp") 
# 
# cacta_ZmPR4 <- cacta_ZmPR4 %>%
#   dplyr::mutate(dp = as.character(dp)) %>%
#   left_join(p_values, by = "dp")
# 
# cacta_ZmPR4$Spearman_p_sign <- "ns"
# cacta_ZmPR4$Spearman_p_sign[cacta_ZmPR4$p_value < 0.05] <- "*"
# cacta_ZmPR4$Spearman_p_sign[cacta_ZmPR4$p_value < 0.01] <- "**"
# cacta_ZmPR4$Spearman_p_sign[cacta_ZmPR4$p_value < 0.001] <- "***"
# cacta_ZmPR4$Spearman_p_sign[cacta_ZmPR4$p_value < 0.0001] <- "****"
# 
# library(ggthemes)
# scales::show_col(few_pal()(8))
# 
# 
# plot_cacta <- cacta_ZmPR4 %>%
#   dplyr::full_join(dp_name) %>%
#   dplyr::filter(!is.na(dp)) %>% 
#   dplyr::mutate(Nb_ZmPR4_peak = as.factor(Nb_ZmPR4_peak)) %>%
#   dplyr::mutate(dp = as.factor(dp), 
#                 dp = factor(dp, levels = c("1","2" ,"3","4", "5", "6", "7"))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = dp_name, y = slope, fill = Nb_ZmPR4_peak)) + 
#   geom_point(aes(x = dp_name, y = slope, group = Nb_ZmPR4_peak), fill = "black",position=position_dodge(width=0.75)) + 
#   scale_fill_manual(values = c("#5DA5DA", "#B276B2")) +
#   scale_color_manual(values = c("#5DA5DA", "#B276B2")) +
#   theme_classic() +
#   theme(plot.background = element_blank()
#         #,legend.position = "bottom"legend.position = c(0.5,0.8)#"top"
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
#         ,axis.text.x = element_text(size=18, angle = 45, margin = margin(t=11))
#         ,axis.title.y = element_text(size=18, margin = margin(r=18), face="bold")
#         ,axis.text.y = element_text(size=18)
#         ,axis.ticks.y = element_blank()
#         ,axis.ticks.x = element_blank()
#         # top, right, bottom, left
#         ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
#   ) 
# geom_text(aes(x = dp, y = 25,  label = Spearman_p_sign), size = 6, check_overlap = TRUE)
# plot_cacta 
# 
# ggsave(paste0(fig_genes_introgression, "alleles_ZmPR4_all_ILs","_", date_today,".pdf"), plot_cacta, width = 14, height = 8)
