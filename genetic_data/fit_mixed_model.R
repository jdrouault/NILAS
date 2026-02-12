library(tidyverse)
library(nlme)
library(patchwork)

### date today 
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)

### ======= Create the path 

root <- "~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/" #### replace by your path 
# root <- "~/Dropbox-SE Dropbox/Randall Wisser/inrae_jd/" #### replace by your path 
haplotypes_path_blues <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_day_dta_04Avr11/")
genetic_data_save_path <- paste0(root_alternative, "/NILAS/main_folder/drouault/data/genotypic_data/NILs/")
fig_phenotype <- paste0(root, "NILAS/main_folder/drouault/paper/Figures/phenotypic_data/")
saving_path <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/from_mixed_model/")
file_supp_path <- paste0(root, "NILAS/main_folder/drouault/paper/Files/")

### ======= Load data

#The data used here are the DLs and Ttem-ant obtained in the script PRNT_blues_dta.Rmd

## Haplotype data 
data <- read.csv(paste0(haplotypes_path_blues, "dataset_lmer_regression_20240911.csv"), header = T) %>%
  dplyr::filter(!grepl("q2g21i12",entryid2) & !grepl("q3g31i01",entryid2) ### filter the lines that have high level of heterozygosity
                & !grepl("q2g12i04",entryid2) & 
                  !grepl("q2g52i02",entryid2))

## NILs data 
regression_NILS <- read.csv(paste0(haplotypes_path_blues, "dataset_lmer_regression_NILs_20240921.csv"), header = T) %>%
  dplyr::filter(!grepl("q2g21i12",entryid2) & !grepl("q3g31i01",entryid2) ### filter the lines that have high level of heterozygosity
                & !grepl("q2g12i04",entryid2) & 
                  !grepl("q2g52i02",entryid2))


## Parents data 
regression_parents <- read.csv(paste0(haplotypes_path_blues, "dataset_lmer_regression_parents_20240921.csv"), header = T) %>%
  dplyr::filter(!grepl("q2g21i12",entryid2) & !grepl("q3g31i01",entryid2) ### filter the lines that have high level of heterozygosity
                & !grepl("q2g12i04",entryid2) & 
                  !grepl("q2g52i02",entryid2))



### ======= Fit the model for the haplotypes  

# diff_gdd = haplotype effect anthesis measured in gdd
# new_mean_daylength = daylength average in SD and LD with intercept SD on 0

data$Parameters <- as.factor(data$Parameters)
data$Parameters <- factor(data$Parameters, levels=rev(levels(data$Parameters)))
data$DL <- as.numeric(as.factor(data$Parameters))
data$DL_num <- data$new_mean_daylength
#note, I did not allow free values for DL after noting the mean DLs were the same to the tenths 
data$DL_num[which(data$DL_num>0)] <- 3.5


temp_data <- data %>%
  dplyr::filter(rp == 1 & qtl==4) %>%
  dplyr::filter(entryid2 != "NILASq4g41i09")

# test_data <- subset(data, data$entryid2=="NILASq1g11i01")
# ggplot(test_data, aes(x = DL_num, y = diff_gdd)) +
#   geom_point() +
#   stat_smooth(method = "lm", se = FALSE)
# summary(lm(diff_gdd ~ DL_num, data=test_data))

## fixed effect model, 
fits <- lmList(diff_gdd ~ DL_num | entryid2, data=data, pool=TRUE, na.action=na.exclude)
int <- as.data.frame(summary(fits)$coefficients[,,1])
int <- cbind.data.frame(haplotype=rownames(int), int)
rownames(int) <- 1:nrow(int)
head(int)
slp <- as.data.frame(summary(fits)$coefficients[,,2])
slp <- cbind.data.frame(haplotype=rownames(slp), slp)
rownames(slp) <- 1:nrow(slp)
head(slp)

#write.csv(int, paste0(saving_path, "haplotype_int", "_",date_today, ".csv"), row.names = F) 
#write.csv(slp, paste0(saving_path, "haplotype_slope", "_",date_today, ".csv"), row.names = F) 

### Haplotype detected with statistical threshold 
int$BH <- p.adjust(int$`Pr(>|t|)`, method="BH")
count(int[which(int$BH<0.05),])
count(int[which(int$BH<0.05 & int$Estimate <0),])
slp$BH <- p.adjust(slp$`Pr(>|t|)`, method="BH")
  slp[which(slp$BH<0.05),]

  
### Haplotype detected with biological threshold 
slp$biological_threshold <- (1.96*slp$`Std. Error`)
sup_biol <- slp$haplotype[which(slp$Estimate > (slp$biological_threshold))]
inf_biol <- slp$haplotype[which(slp$Estimate < -(slp$biological_threshold))]

biological_detected_haplotype_slope <- c(sup_biol, inf_biol)

int$biological_threshold <- (1.96*int$`Std. Error`)
sup_biol_int <- int$haplotype[which(int$Estimate > (int$biological_threshold))]
inf_biol_int <- int$haplotype[which(int$Estimate < -(int$biological_threshold))]

biological_detected_haplotype_int <- c(sup_biol_int, inf_biol_int)

## ======= Fit the model for the NILs

regression_NILS$Parameters <- as.factor(regression_NILS$Parameters)
regression_NILS$Parameters <- factor(regression_NILS$Parameters, levels=rev(levels(regression_NILS$Parameters)))
regression_NILS$DL <- as.numeric(as.factor(regression_NILS$Parameters))
regression_NILS$DL_num <- regression_NILS$new_mean_daylength
#note, I did not allow free values for DL after noting the mean DLs were the same to the tenths 
# regression_lines$DL_num[which(regression_lines$DL_num>0)] <- 3.5

## fixed effect model, 
fits_NILs <- lmList(GDD_dta_em ~ DL_num | entryid2, data=regression_NILS, pool=TRUE, na.action=na.exclude)
int_NILs <- as.data.frame(summary(fits_NILs)$coefficients[,,1])
int_NILs <- cbind.data.frame(NILs=rownames(int_NILs), int_NILs)
rownames(int_NILs) <- 1:nrow(int_NILs)
head(int_NILs)
slp_NILs <- as.data.frame(summary(fits_NILs)$coefficients[,,2])
slp_NILs <- cbind.data.frame(NILs=rownames(slp_NILs), slp_NILs)
rownames(slp_NILs) <- 1:nrow(slp_NILs)
head(slp_NILs)

hist(slp_NILs$Estimate)

### Haplotype detected with statistical threshold 
int_NILs$BH <- p.adjust(int_NILs$`Pr(>|t|)`, method="BH")
slp_NILs$BH <- p.adjust(slp_NILs$`Pr(>|t|)`, method="BH")

#write.csv(int_NILs, paste0(saving_path, "NILs_int", "_",date_today, ".csv"), row.names = F) 
#write.csv(slp_NILs, paste0(saving_path, "NILs_slope", "_",date_today, ".csv"), row.names = F) 

## ======= Fit the model for the parents

regression_parents$Parameters <- as.factor(regression_parents$Parameters)
regression_parents$Parameters <- factor(regression_parents$Parameters, levels=rev(levels(regression_parents$Parameters)))
regression_parents$DL <- as.numeric(as.factor(regression_parents$Parameters))
regression_parents$DL_num <- regression_parents$new_mean_daylength
#note, I did not allow free values for DL after noting the mean DLs were the same to the tenths 
# regression_parents$DL_num[which(regression_parents$DL_num>0)] <- 3.5

## fixed effect model, 
fits_parents <- lmList(GDD_dta_em ~ DL_num | entryid2, data=regression_parents, pool=TRUE, na.action=na.exclude)
int_parents <- as.data.frame(summary(fits_parents)$coefficients[,,1])
int_parents <- cbind.data.frame(Parents=rownames(int_parents), int_parents)
rownames(int_parents) <- 1:nrow(int_parents)
head(int_parents)
slp_parents <- as.data.frame(summary(fits_parents)$coefficients[,,2])
slp_parents <- cbind.data.frame(Parents=rownames(slp_parents), slp_parents)
rownames(slp_parents) <- 1:nrow(slp_parents)
head(slp_parents)

hist(slp_parents$Estimate)

#write.csv(int_parents, paste0(saving_path, "parents_int", "_",date_today, ".csv"), row.names = F) 
#write.csv(slp_parents, paste0(saving_path, "parents_slope", "_",date_today, ".csv"), row.names = F) 


#### Make file S4

## made sub file with haplotypes 
summary_haplo <- data %>%
  dplyr::select(entryid2, Parameters, mean_daylength, new_mean_daylength) %>%
  dplyr::rename(entryid = entryid2, adjsuted_mean_daylength = new_mean_daylength) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider( names_from = Parameters,
                      values_from = c(mean_daylength, adjsuted_mean_daylength))

summary_table_LD <- data %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 <mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_LD = env) %>%
  dplyr::rename(entryid = entryid2)

summary_table_SD <- data%>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 > mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_SD = env) %>%
  dplyr::rename(entryid = entryid2)

slp <- slp %>%
  dplyr::rename(slope = Estimate, entryid = haplotype, se_slope = `Std. Error`, tvalue_slope = `t value`, pval_slope = `Pr(>|t|)`, adjusted_pval_slope = BH, biological_threshold_slope = biological_threshold ) %>%
  dplyr::mutate(genetic_type = "haplotype") %>%
  dplyr::full_join(summary_table_SD) %>%
  dplyr::full_join(summary_table_LD) %>%
  dplyr::full_join(summary_haplo)

haplotype_prnftp <- int %>%
  dplyr::rename(intercept = Estimate, entryid = haplotype, se_intercept = `Std. Error`, tvalue_intercept = `t value`, pval_intercept = `Pr(>|t|)`, adjusted_pval_intercept = BH, biological_threshold_intercept = biological_threshold) %>%
  dplyr::full_join(slp)

## made sub file with NILs  

summary_NILS <- regression_NILS %>%
  dplyr::select(entryid2, Parameters, mean_daylength, new_mean_daylength) %>%
  dplyr::rename(entryid = entryid2, adjsuted_mean_daylength = new_mean_daylength) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider( names_from = Parameters,
                      values_from = c(mean_daylength, adjsuted_mean_daylength))

summary_table_LD_NILS <- regression_NILS %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 <mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_LD = env) %>%
  dplyr::rename(entryid = entryid2)

summary_table_SD_NILS <- regression_NILS %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 > mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_SD = env) %>%
  dplyr::rename(entryid = entryid2)

slp_NILs <- slp_NILs %>%
  dplyr::rename(slope = Estimate, entryid = NILs, se_slope = `Std. Error`, tvalue_slope = `t value`, pval_slope = `Pr(>|t|)`, adjusted_pval_slope = BH) %>%
  dplyr::mutate(genetic_type = "NIL") %>%
  dplyr::full_join(summary_table_LD_NILS) %>%
  dplyr::full_join(summary_table_SD_NILS) %>%
  dplyr::full_join(summary_NILS)

NILs_prnftp <- int_NILs %>%
  dplyr::rename(intercept = Estimate, entryid = NILs, se_intercept = `Std. Error`, tvalue_intercept = `t value`, pval_intercept = `Pr(>|t|)`, adjusted_pval_intercept = BH) %>%
  dplyr::full_join(slp_NILs)

## made sub file with parents 

summary_parents <- regression_parents %>%
  dplyr::select(entryid2, Parameters, mean_daylength, new_mean_daylength) %>%
  dplyr::rename(entryid = entryid2, adjsuted_mean_daylength = new_mean_daylength) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider( names_from = Parameters,
                      values_from = c(mean_daylength, adjsuted_mean_daylength))

summary_table_LD_parents <- regression_parents %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 <mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_LD = env) %>%
  dplyr::rename(entryid = entryid2)

summary_table_SD_parents <- regression_parents %>%
  as.data.frame() %>%
  dplyr::group_by(entryid2) %>%
  dplyr::filter(13.83 > mean_daylength) %>%
  dplyr::group_by(entryid2) %>%
  dplyr::summarise_at(vars(env), ~ sum(!(is.na(.)))) %>%
  dplyr::rename(env_SD = env) %>%
  dplyr::rename(entryid = entryid2)

slp_parents <- slp_parents %>%
  dplyr::rename(slope = Estimate, entryid = Parents, se_slope = `Std. Error`, tvalue_slope = `t value`, pval_slope = `Pr(>|t|)`) %>%
  dplyr::mutate(genetic_type = "parent") %>%
  dplyr::full_join(summary_table_LD_parents) %>%
  dplyr::full_join(summary_table_SD_parents) %>%
  dplyr::full_join(summary_parents)
  
parents_prnftp <- int_parents %>%
  dplyr::rename(intercept = Estimate, entryid = Parents, se_intercept = `Std. Error`, tvalue_intercept = `t value`, pval_intercept = `Pr(>|t|)`) %>%
  dplyr::full_join(slp_parents)

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



all_prnftp <- parents_prnftp %>%
  dplyr::full_join(NILs_prnftp) %>%
  dplyr::full_join(haplotype_prnftp) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid)), 
                qtl = as.numeric(gsub(".*q|g.*", "", entryid))) %>%
  dplyr::full_join(mini_table) %>%
  dplyr::full_join(mini_table2)

  

# write.csv(all_prnftp, paste0(file_supp_path, "File_S4_new", "_",date_today, ".csv"), row.names = F) 


##### === Fished exact test 


summary_table_slope <- slp %>%
  summarise(
    above_threshold = sum(slope > biological_threshold_slope),        # Count values above the threshold
    below_threshold = sum(slope < -biological_threshold_slope)         # Count values below the threshold
  )

summary_table_int <- int %>%
  summarise(
    above_threshold = sum(Estimate > biological_threshold),        # Count values above the threshold
    below_threshold = sum(Estimate < -biological_threshold)         # Count values below the threshold
  )

fisher_matrix<- matrix(c(55, 610, 172, 493), nrow = 2)

# Perform Fisher's Exact Test
fisher_result <- fisher.test(fisher_matrix)

##### === Create histogram for detected haplotypes 

## Statistical threshold 

slope <- slp %>%
ggplot() +
 geom_histogram(aes(x=Estimate, fill = (BH < 0.05)), color = "black", alpha = 0.8)+ 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF")) +
  theme_classic() +
  theme(plot.background = element_blank(),
        #,legend.position = "bottom"#
        ,legend.position = "none"
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
        ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
        ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(30), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,strip.text = element_text(size = 30, face="bold"),
        plot.tag = element_text(size = 30), 
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  labs(tag = "(A)")+
  ylab("") +
  xlab("Haplotype effect \n on photoperiod sensitivity")

intercept <- int %>%
  ggplot() +
  geom_histogram(aes(x=Estimate, fill = (BH < 0.05)), color = "black", alpha = 0.8) + 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF")) +
  theme_classic() +
  theme(plot.background = element_blank(),
        #,legend.position = "bottom"#
        ,legend.position = "none"
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
        ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
        ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(30), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,strip.text = element_text(size = 30, face="bold"),
        plot.tag = element_text(size = 30), 
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  ylab("") +
  labs(tag = "(B)")+
  xlab(quote(bold(atop("Haplotype effect on", bold("flowering time" ~ bolditalic("per se"))))))
# xlab(quote(bold("Haplotype effect on flowering time") ~ bolditalic("per se"))) 
  

multi_plot <- slope + intercept + plot_layout(ncol = 2, nrow = 1, guides = 'keep') #& 
#plot_annotation(tag_levels = "A") 
multi_plot


## Biological threshold 

slope_bio <- slp %>%
  ggplot() +
  geom_histogram(aes(x=Estimate, fill = (Estimate > biological_threshold | Estimate < - biological_threshold)), color = "black", alpha = 0.8)+ 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF")) +
  theme_classic() +
  theme(plot.background = element_blank(),
        #,legend.position = "bottom"#
        ,legend.position = "none"
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
        ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
        ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(30), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,strip.text = element_text(size = 30, face="bold"), 
        plot.tag = element_text(size = 30), 
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")) +
  ylab("") +
  labs(tag = "(C)")+
  xlab("Haplotype effect \n on photoperiod sensitivity")


summary_table_slope <- slp %>%
  summarise(
    above_threshold = sum(slope > biological_threshold_slope),        # Count values above the threshold
    below_threshold = sum(slope < -biological_threshold_slope)         # Count values below the threshold
  )

intercept_bio <- int %>%
  ggplot() +
  geom_histogram(aes(x=Estimate, fill = (Estimate > biological_threshold | Estimate < - biological_threshold)), color = "black", alpha = 0.8) + 
  scale_fill_manual(values = c("#868686FF", "#EFC000FF")) +
  theme_classic() +
  theme(plot.background = element_blank(),
        #,legend.position = "bottom"#
        ,legend.position = "none"
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
        ,axis.title.x = element_text(size=30, margin = margin(t=30), face="bold")
        ,axis.text.x = element_text(size=30, angle = 90, vjust = 0.5, margin = margin(t=11))
        ,axis.title.y = element_text(size=30, margin = margin(30), face="bold")
        ,axis.text.y = element_text(size=30)
        ,axis.ticks.y = element_blank()
        ,axis.ticks.x = element_blank()
        ,strip.text = element_text(size = 30, face="bold"),
        plot.tag = element_text(size = 30), 
        # top, right, bottom, left
        ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
  ) +
  labs(tag = "(D)")+
  ylab("") +
  xlab(quote(bold(atop("Haplotype effect on", bold("flowering time" ~ bolditalic("per se"))))))
# xlab(quote(bold("Haplotype effect on flowering time") ~ bolditalic("per se"))) 


multi_plot <- slope + intercept + slope_bio+ intercept_bio + plot_layout(ncol = 2, nrow = 2, guides = 'keep') #& 
#plot_annotation(tag_levels = "A") 
multi_plot


ggsave(file = paste0(fig_phenotype, date_today, "histogram_haplotype_effect_slope_intercept.pdf"),
       multi_plot, width = 18, height = 18)


#### Check if the lines that were marked as significantly different have high het, several intro or candidate gene

genetic_data_save_path <- paste0(root, "/NILAS/main_folder/drouault/data/genotypic_data/NILs/")


  
  ### slope 

## statistics 

  summary_too_much_intro <- summary_too_much_intro %>%
    dplyr::mutate(lines = gsub("s[0-9]p[0-9]r[0-9]$", "", Sample))
  
  slp[which(slp$BH<0.05),]$haplotype %in% summary_too_much_intro$lines 

  candidates_genes_position_lines <- summary_qtl_insertion %>%
    dplyr::filter(zmpr1_cct == 1 | zmpr1_elf4 == 1 | zmpr2_zcn8 == 1 | 
                    zmpr2_rap2.7 == 1 | zmpr3_conz1 ==1| zmpr3_cct == 1| zmpr4_cacta ==1  )

  slp[which(slp$BH<0.05),]$haplotype %in% candidates_genes_position_lines$entryid 
  
  ##without candidates and without several intro
  candidates_svlintro_biol_int <- unique(c(final_subset$haplotype, candidate_stat_inbt$haplotype))
  
  subset_int$haplotype[!(subset_int$haplotype %in% candidates_svlintro_biol_int)]  
  

  ## Biological 
  
  #Several intro
  biol_haplotypes_too_many_intro <- biological_detected_haplotype_slope[biological_detected_haplotype_slope %in% summary_too_much_intro$lines]
  
  #Candidates 
  biol_haplotypes_candidates_genes <- biological_detected_haplotype_slope[biological_detected_haplotype_slope %in% candidates_genes_position_lines$entryid]
  
  #without candidates and without several intro
  candidates_svlintro_biol_slope <- unique(c(biol_haplotypes_candidates_genes, biol_haplotypes_too_many_intro))
  
  biological_detected_haplotype_slope[!(biological_detected_haplotype_slope %in% candidates_svlintro_biol_slope)]  
  
  ### intercept 
  
  ## Statistics
  
  subset_int <- int[int$BH < 0.05, ]
  dim(subset_int)
  
  #With several intro
  final_subset <- subset_int[subset_int$haplotype %in% summary_too_much_intro$lines, ]
  dim(final_subset)
  
  #With candidate gene
  candidate_stat_inbt <- (subset_int[subset_int$haplotype %in% candidates_genes_position_lines$entryid, ])
  
  #Without candidate gene and with a single introgression 
  
  candidates_svlintro_biol_int <- unique(c(final_subset$haplotype, candidate_stat_inbt$haplotype))
  
  subset_int$haplotype[!(subset_int$haplotype %in% candidates_svlintro_biol_int)]  

  
  ## Biological 
  
  #With several introgressions
  biol_haplotypes_too_many_intro <- biological_detected_haplotype_int[biological_detected_haplotype_int %in% summary_too_much_intro$lines]
  length(biol_haplotypes_too_many_intro)
  
  #With candidate genes
  biol_haplotypes_candidates_genes <- biological_detected_haplotype_int[biological_detected_haplotype_int %in% candidates_genes_position_lines$entryid]
  length(biol_haplotypes_candidates_genes)
  
  #Without candidate gene and with a single introgression 
  
  candidates_svlintro_biol_int <- unique(c(biol_haplotypes_too_many_intro, biol_haplotypes_candidates_genes))
  
  biological_detected_haplotype_int[!(biological_detected_haplotype_int %in% candidates_svlintro_biol_int)]  
  
  