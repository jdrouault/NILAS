library(tidyverse)
library(nlme)
library(patchwork)
library(gridExtra)
library(lme4)



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


### date today 
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)

### ======= Create the path 

root <- "~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/" #### replace by your path 
# root <- "~/Dropbox-SE Dropbox/Randall Wisser/inrae_jd/" #### replace by your path 
haplotypes_path_blues <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_day_dta_04Avr11/")
genetic_data_save_path <- paste0(root, "/NILAS/main_folder/drouault/data/genotypic_data/NILs/")
fig_phenotype <- paste0(root, "NILAS/main_folder/drouault/paper/Figures/phenotypic_data/")
source_path <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/from_mixed_model/")
genetic_data_save_path <- paste0(root, "NILAS/main_folder/drouault/data/genotypic_data/NILs/")
fig_genotype <-  paste0(root,"NILAS/main_folder/drouault/paper/Figures/genotypic_data/effect_genes_introgressionNB/")


### ======= Load data

##Haplotypes 
haplotype_slope <- read.csv(paste0(source_path, "haplotype_slope_20240921.csv"), header = T)
haplotype_int <- read.csv(paste0(source_path, "haplotype_int_20240921.csv"), header = T)

##NILs (filtered for heterozygosity)
NILs_slope <- read.csv(paste0(source_path, "NILs_slope_20240921.csv"), header = T)
NILs_int <- read.csv(paste0(source_path, "NILs_int_20240921.csv"), header = T)

##Parents (filtered for heterozygosity)
parents_slope <- read.csv(paste0(source_path, "parents_slope_20240921.csv"), header = T)
parents_int <- read.csv(paste0(source_path, "parents_int_20240921.csv"), header = T)

## Filtrer avec la liste 1 
 
group1_NILs <- read.csv(paste0(genetic_data_save_path, "20240913Selected_NILs_after_filtration.csv"), header = T)

several_intro_list <- read.csv(paste0(genetic_data_save_path, "20240913lines_double_triple_quadruple_intro.csv"), header = T)

group2_NILs <- read.csv(paste0(genetic_data_save_path, "20240912Selected_NILs_after_filtration_for_heterozigosity.csv"), header = T) 

summary_several_intro_list <- read.csv(paste0(genetic_data_save_path, "20240912Selected_NILs_after_filtration_for_heterozigosity.csv"), header = T)  %>%
  dplyr::mutate(qtl = as.character(qtl)) %>%
  dplyr::filter(entryid %in% several_intro_list$lines)

### ======= Epistasis analysis 

#### With the background 

##Filtering to keep only NILs with a single introgression
NILs_slope_signle_intro <- NILs_slope %>%
  dplyr::filter(NILs %in% unique(group1_NILs$entryid)) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", NILs)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", NILs)), 
                qtl = as.numeric(gsub(".*q|g.*", "", NILs))) 

NILs_int_single_intro <- NILs_int %>%
  dplyr::filter(NILs %in% unique(group1_NILs$entryid)) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", NILs)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", NILs)), 
                qtl = as.numeric(gsub(".*q|g.*", "", NILs))) 

haplotype_slopes_single_intro <- haplotype_slope %>%
  dplyr::filter(haplotype %in% unique(group1_NILs$entryid)) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", haplotype)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", haplotype)), 
                qtl = as.numeric(gsub(".*q|g.*", "", haplotype))) 

haplotype_int_single_intro <- haplotype_int %>%
  dplyr::filter(haplotype %in% unique(group1_NILs$entryid)) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", haplotype)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", haplotype)), 
                qtl = as.numeric(gsub(".*q|g.*", "", haplotype))) 

### Slope

distribution_plot <- haplotype_slopes_single_intro %>%
  # dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", NILs)),
  #               dp = as.numeric(gsub(".*g|[0-9]i.*", "", NILs)), 
  #               qtl = as.numeric(gsub(".*q|g.*", "", NILs))) %>%
  dplyr::rename(gb = rp) %>%
  dplyr::mutate(gb = ifelse(gb == 1, "2369", "LH123hT"),
                qtl = dplyr::recode(qtl, "1" = "ZmPR1" , "2" =  "ZmPR2", "3" =  "ZmPR3", "4" = "ZmPR4" ), 
                ) %>%
  ggplot() +
  geom_histogram(aes(x = Estimate), colour = 1, fill = "#0072B2") +
  scale_fill_manual(values = c("grey", "#3399FF"))+
  facet_grid(gb~qtl, scales= "fixed") +
  xlab("slope of sensitivity") +
  # geom_vline(data = dataInt,  mapping = aes(xintercept= Int), linetype="dashed", colour = dataInt$category, linewidth=2) +
  # scale_x_continuous(breaks=seq(), 
  #                    expand = c(0.05, 0)) + 
  labs(title = '') 

distribution_plot <- distribution_plot+ theme(plot.background = element_blank(),
                                              legend.position = "none"
                                              ,legend.direction="vertical"
                                              ,legend.key.size = unit(1, "cm")
                                              ,legend.text = element_text(size = 30)
                                              ,legend.title = element_text(size = 30)
                                              ,strip.text = element_text(size = 30)
                                              ,aspect.ratio = 1
                                              ,panel.grid.minor = element_blank()
                                              ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
                                              ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
                                              ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
                                              ,panel.background = element_blank()
                                              ,plot.title = element_text(face="bold")
                                              ,axis.title.x = element_text(size=25, margin = margin(t=25), face="bold")
                                              ,axis.text.x = element_text(size=25, margin = margin(t=25))
                                              ,axis.title.y = element_text(size=25, margin = margin(r=25), face="bold")
                                              ,axis.text.y = element_text(size=25)
                                              ,axis.ticks.y = element_blank()
                                              ,axis.ticks.x = element_blank()
                                              ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
) 

distribution_plot

# ggsave(paste0(fig_phenotype, "distribution_slope_NILs_faceted","_", date_today,".pdf"), distribution_plot, width = 25, height = 22)


### intercept 
distribution_plot_intercept <- haplotype_int %>%
  # dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", NILs)),
  #               dp = as.numeric(gsub(".*g|[0-9]i.*", "", NILs)), 
  #               qtl = as.numeric(gsub(".*q|g.*", "", NILs))) %>%
  dplyr::rename(gb = rp) %>%
  dplyr::mutate(gb = ifelse(gb == 1, "2369", "LH123hT"),
                qtl = dplyr::recode(qtl, "1" = "ZmPR1" , "2" =  "ZmPR2", "3" =  "ZmPR3", "4" = "ZmPR4" ), 
  ) %>%
  ggplot() +
  geom_histogram(aes(x = Estimate), colour = 1, fill = "#0072B2") +
  scale_fill_manual(values = c("grey", "#3399FF"))+
  facet_grid(gb~qtl, scales= "fixed") +
  xlab("slope of sensitivity") +
  # geom_vline(data = dataInt,  mapping = aes(xintercept= Int), linetype="dashed", colour = dataInt$category, linewidth=2) +
  # scale_x_continuous(breaks=seq(), 
  #                    expand = c(0.05, 0)) + 
  labs(title = '') 

distribution_plot_intercept <- distribution_plot_intercept+ theme(plot.background = element_blank(),
                                              legend.position = "none"
                                              ,legend.direction="vertical"
                                              ,legend.key.size = unit(1, "cm")
                                              ,legend.text = element_text(size = 30)
                                              ,legend.title = element_text(size = 30)
                                              ,strip.text = element_text(size = 30)
                                              ,aspect.ratio = 1
                                              ,panel.grid.minor = element_blank()
                                              ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=3)
                                              ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=3)
                                              ,panel.border = element_rect(colour="grey", size = 0.5, fill=NA)
                                              ,panel.background = element_blank()
                                              ,plot.title = element_text(face="bold")
                                              ,axis.title.x = element_text(size=25, margin = margin(t=25), face="bold")
                                              ,axis.text.x = element_text(size=25, margin = margin(t=25))
                                              ,axis.title.y = element_text(size=25, margin = margin(r=25), face="bold")
                                              ,axis.text.y = element_text(size=25)
                                              ,axis.ticks.y = element_blank()
                                              ,axis.ticks.x = element_blank()
                                              ,plot.margin=unit(c(0.55,0.5,0.8,0.5), "cm")
) 

distribution_plot_intercept

# ggsave(paste0(fig_phenotype, "distribution_intercept_NILs_faceted","_", date_today,".pdf"), distribution_plot_intercept, width = 25, height = 22)


### testing the difference between variances

## Slope 
# empty_matrix <- matrix(NA, nrow = 1, ncol = 2)
data <- NULL
f_test_pval <- do.call(rbind,lapply(split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl)), function(x) {
  # x <- split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl))[[2]]
  rp1 <- x %>% filter(rp == 1)
  rp2 <- x %>% filter(rp == 2)
  f_test <- var.test(rp1$Estimate, rp2$Estimate, alternative = "two.sided")
  f_pval <- f_test$p.value
  qtl <- unique(x$qtl)
  data <- cbind(f_pval,qtl, data)
}))

## Intercept 

data <- NULL
f_test_pval_intercept <- do.call(rbind,lapply(split(NILs_int_single_intro,as.factor(NILs_int_single_intro$qtl)), function(x) {
  # x <- split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl))[[2]]
  rp1 <- x %>% filter(rp == 1)
  rp2 <- x %>% filter(rp == 2)
  f_test <- var.test(rp1$Estimate, rp2$Estimate, alternative = "two.sided")
  f_pval <- f_test$p.value
  qtl <- unique(x$qtl)
  data <- cbind(f_pval,qtl, data)
}))

### testing the difference between means

## Slope 
anova_result <- aov(Estimate ~ rp, data = NILs_slope_signle_intro)
summary(anova_result)

NILs_slope_signle_intro 

haplotype_slopes_single_intro_pivoted <- haplotype_slopes_single_intro %>%
  tidyr::pivot_wider( names_from = qtl,
                      values_from = haplotype) %>%
  dplyr::mutate(across(`1`:`2`, ~(as.factor(.))))
  

anova_result_interaction <- aov(aov(Estimate ~ as.factor(rp) + as.factor(qtl) + as.factor(rp*qtl), data = haplotype_slopes_single_intro))
  # Estimate ~ as.factor(rp) + as.factor(`1`) + as.factor(`2`) + as.factor(`3`) + as.factor(`4`), 
  #                                 # as.factor(rp):as.factor(`1`) +
  #                                 # as.factor(rp):as.factor(`2`) +
  #                                 # as.factor(rp):as.factor(`3`) +
  #                                 # as.factor(rp):as.factor(`4`),
  #                               data = na.moit(haplotype_slopes_single_intro_pivoted)

summary(anova_result_interaction)


anova_result_interaction <- lmer(Estimate ~  rp + qtl + (1 | rp:qtl), 
                                 data = haplotype_slopes_single_intro)

summary(anova_result_interaction)

fits_parents <- lmList(Estimate ~ rp + qtl | rp:qtl, data=haplotype_slopes_single_intro, pool=TRUE, na.action=na.exclude)



t_test_pval <- do.call(rbind,lapply(split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl)), function(x) {
  # x <- split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl))[[2]]
  rp1 <- x %>% filter(rp == 1)
  rp2 <- x %>% filter(rp == 2)
  t_test_result <- t.test(rp1$Estimate, rp2$Estimate, var.equal = TRUE) 
  test_pval <- t_test_result$p.value
  qtl <- unique(x$qtl)
  data <- cbind(test_pval,qtl, data)
}))

## Intercept 

anova_result <- aov(Estimate ~ rp, data = NILs_int_single_intro)
summary(anova_result)

anova_result_interaction <- aov(Estimate ~ rp *qtl , data = NILs_int_single_intro)
summary(anova_result_interaction)

t_test_pval_int <- do.call(rbind,lapply(split(NILs_int_single_intro,as.factor(NILs_int_single_intro$qtl)), function(x) {
  # x <- split(NILs_slope_signle_intro,as.factor(NILs_slope_signle_intro$qtl))[[2]]
  rp1 <- x %>% filter(rp == 1)
  rp2 <- x %>% filter(rp == 2)
  t_test_result <- t.test(rp1$Estimate, rp2$Estimate, var.equal = TRUE) 
  test_pval <- t_test_result$p.value
  qtl <- unique(x$qtl)
  data <- cbind(test_pval,qtl, data)
}))



#### Epistasis Between several ZmPR introgressions


## Slope 

haplotype_slope_1<- haplotype_slope %>%
  dplyr::rename(lines=haplotype) %>%
  dplyr::full_join(several_intro_list) %>%
  dplyr::mutate(total = ifelse(is.na(total), 1, total))%>%
  dplyr::filter(!grepl("q2g21i12s1p3r1",Sample) & !grepl("q3g31i01s2p1r1",Sample)
                & !grepl("q2g12i04s1p3r1",Sample) & 
                  !grepl("q2g52i02s1p1r1",Sample)) 
  # dplyr::mutate(total = as.factor(total)) %>%

haplotype_slope_plot <- haplotype_slope_1 %>%
  dplyr::mutate(total = as.numeric(total)) %>%
  ggplot(aes(x = total , y =Estimate)) +
  geom_point(aes(x = total , y =Estimate, group = total)) +
  # geom_abline(slope = coef(total.lm)[["total"]], 
  #             intercept = coef(total.lm)[["(Intercept)"]]) +
  stat_poly_line(aes(x = total , y =Estimate), color = "black") +
  stat_poly_eq(use_label(c("eq", "R2")), size = 6) +
  xlab("Total number of ZmPR introgressed") +
  ylab("Photoperiod sensitivity") +
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
  ) 

haplotype_slope_1 <- haplotype_slope_1 %>%
  dplyr::mutate(total = as.numeric(total))

total.lm <- lm(Estimate ~ total, haplotype_slope_1)
summary(total.lm)


emm_ZmPR <- emmeans(total.lm, ~ total)  # Estimate marginal means for dp within each level of rp
summary(emm_ZmPR)

cld_1 <- cld(emm_ZmPR, Letters = letters)
cld_df <- unnest(cld_1)

## Intercept

haplotype_int_1<- haplotype_int %>%
  dplyr::rename(lines=haplotype) %>%
  dplyr::full_join(several_intro_list) %>%
  dplyr::mutate(total = ifelse(is.na(total), 1, total))%>%
  dplyr::filter(!grepl("q2g21i12s1p3r1",Sample) & !grepl("q3g31i01s2p1r1",Sample)
                & !grepl("q2g12i04s1p3r1",Sample) & 
                  !grepl("q2g52i02s1p1r1",Sample)) 
# dplyr::mutate(total = as.factor(total)) %>%

haplotype_int_1 <- haplotype_int_1 %>%
  dplyr::mutate(total = as.numeric(total))

total.lm <- lm(Estimate ~ total, haplotype_int_1)
summary(total.lm)

emm_ZmPR <- emmeans(total.lm, ~ total)  # Estimate marginal means for dp within each level of rp
summary(emm_ZmPR)

cld_1 <- cld(emm_ZmPR, Letters = letters)
cld_df <- unnest(cld_1)

##plot 
haplotype_int_plot <- haplotype_int_1 %>%
  dplyr::mutate(total = as.numeric(total))%>%
  ggplot(aes(x = total , y =Estimate)) +
  geom_point(aes(x = total , y =Estimate, group = total)) +
  # geom_abline(slope = coef(total.lm)[["total"]], 
  #             intercept = coef(total.lm)[["(Intercept)"]]) +
  stat_poly_line(aes(x = total , y =Estimate), color = "black") +
  stat_poly_eq(use_label(c("eq", "R2")), size = 6) +
  xlab("Total number of ZmPR introgressed") +
  ylab(quote(bold("flowering time" ~ bolditalic("per se")))) +
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
  ) 




# ggarrange(haplotype_slope_plot, haplotype_int_plot, 
#           labels = c("A", "B"),
#           ncol = 2, nrow = 1)


ggsave(paste0(fig_genotype, "fig_nbZmPR_introgressed","_", date_today,".pdf"), 
       ggarrange(haplotype_slope_plot, haplotype_int_plot, 
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1), width = 12, height = 5)


group1_NILs_plot  <- haplotype_slope_single_intro %>%
  dplyr::filter(rp == 2) %>%
  dplyr::mutate(i = str_extract(haplotype, "i\\d+"), 
                i = str_extract(i, "\\d{2}")) %>%
  ggplot(aes(x = i, y = Estimate)) +
  geom_hline(yintercept = 4.51264) +
  geom_hline(yintercept = -4.51264) + 
  geom_col() +
  facet_grid(dp ~ qtl)


### Run the NILAS.fig script first and set rp = 1, dp = 4 (5 & 4 are switched in the script), ZmPR = 1

candidate_gene_coordinate <- group2_NILs  %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid)), 
                qtl = as.numeric(gsub(".*q|g.*", "", entryid))) 


  

Q <- list()
 rp <- 1:2
 dp <- 3:4

 Palette <- c("#0072B2", "#E69F00", "white")

for (r in 1:2) {
  for (d in  3:4) {

    # r <- 1
    # d <- 3
    
    candidate_gene_coordinate_sub <- candidate_gene_coordinate %>%
      dplyr::filter(rp == 1 | rp == 2 ) %>%
      dplyr::filter( dp == 4 | dp ==  3) %>% 
      dplyr::filter(qtl == 4)%>%
      dplyr::filter(zmpr4_cacta == 1) #rp == r & qtl == 4 & dp == d & 
    
    n <- paste(r, d, sep= "_")
    
    haplotype_single_intro <- group1_NILs %>%
      dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid)),
                    dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid)), 
                    qtl = as.numeric(gsub(".*q|g.*", "", entryid)))  %>%
      dplyr::filter(rp == 1 | rp == 2 ) %>%
      dplyr::filter( dp == 4 | dp ==  3) %>% 
      dplyr::filter(qtl == 4)
      
    
Q<- haplotype_slope%>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", haplotype)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", haplotype)), 
                qtl = as.numeric(gsub(".*q|g.*", "", haplotype))) %>%
  dplyr::filter(rp == 1 | rp == 2 ) %>%
  dplyr::filter( dp == 4 | dp ==  3) %>% 
  dplyr::filter(qtl == 4) %>%
  dplyr::rename(entryid2 = haplotype) %>%
  dplyr::mutate(i = str_extract(entryid2, "i\\d+"), 
                i = str_extract(i, "\\d{2}"), 
                Type = "haplotype", 
                i = factor(i, levels = c("12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01")),
                marker = ifelse(entryid2 %in% haplotype_single_intro$entryid, "filtered", "not" ),
                filtered = ifelse(marker == "filtered", "yes", "no"),
                marker = ifelse(entryid2 %in% candidate_gene_coordinate_sub$entryid, "ok", "no"), 
                marker = factor(marker, levels = c("no", "ok")), 
                filtered = factor(filtered, levels = c("yes", "no"))) %>%
  dplyr::full_join(NILAS_newid) %>% 
  # dplyr::mutate(id_new = as.factor(id_new),
  #               id_new = fct_rev(id_new)) %>%
  # dplyr::full_join(ref_parent) %>%
  ggplot(aes(x = (-id_new+49), y = Estimate)) + # xmin = id_new, xmax =  (id_new), ymin =Estimate - Std..Error, ymax = Estimate  + Std..Error)) +
  geom_hline(yintercept = 4.51264, linetype = "dashed") +
  geom_hline(yintercept = -4.51264, linetype = "dashed") + 
  geom_hline(yintercept = 0) + 
  scale_x_discrete(limits = levels(NILAS_newid$id_new)) +
  scale_y_continuous(breaks = c(-6:20), limits = c(-6, 20), expand = c(0.01, 0))  + 
  geom_errorbar( aes(x=-id_new+49, ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.5, colour="black", alpha=0.5, size=0.8) +
  # geom_errorbarh() + 
  # Set limits
  geom_col(color = "black", aes(fill = marker, alpha = filtered)) + coord_flip() +
  # geom_segment(aes(y = unique(rep(min(Estimate),48)), x = -id_new+49+0.5,
  #                  yend = unique(rep(max(Estimate),48)), xend = -id_new+49+0.5),
  #              colour="black") +
  geom_vline(aes(xintercept = -id_new+48+0.5), colour="grey", size=0.5) +
  scale_fill_manual(values = c(Palette)) +
  scale_alpha_manual(values = c(1, 0.5)) +
  ylab("Haplotype effect \n on photoperiod sensitivity") +
  xlab("") + 
  ggtitle("B") + #paste(r,d,"ZmPR4", sep = " ")
  theme(plot.background = element_blank()
        ,legend.position = "none"
        #,legend.position = "top"
        #,legend.direction="vertical"
        #,legend.key.size = unit(0.5, "cm")
        #,legend.text = element_text(size = 10)
        #,aspect.ratio = 1
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_blank()
        ,panel.grid.major.y = element_blank() #element_line(colour="white", linetype = "dashed")    # ,panel.border = element_rect(colour="grey", linewidth = 0.5, fill=NA)
        ,panel.background = element_blank() 
        ,plot.title = element_text(face="bold", vjust=1)
        ,axis.title.x = element_text(size=14, vjust=1.5)
        ,axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1)
        ,axis.title.y = element_text(size=14, angle=45, hjust=1, vjust=1)
        ,axis.text.y = element_blank() #element_text(size=14, angle=45, hjust=1, vjust=1)
        #top, right, bottom, left
        ,plot.margin=unit(c(0.1,0.5,1,0.01), "cm") 
  )
  }}
 
 Q
 
 
grid.arrange(grobs = Q, ncol = 5) 
  
  
pdf(paste(fig_genotype, "Fine_mapping_ZmPR4" , ".pdf", sep=""),
      width=20 , height=12)
 cacta_fine_map | Q
  #grid.arrange(grobs = Q, ncol = 5,  nrow=1, widths=c(0.3,0.15,0.15,0.4))
  dev.off()
  
  
# Perform paired t-test for each value

ref_parent <- subset(parents_slope, Parents == "NILASq0g01i14")
ref_parent <- ref_parent %>%
  dplyr::rename(entryid2 = Parents) %>%
dplyr::mutate(i = str_extract(entryid2, "i\\d+"), 
              i = str_extract(i, "\\d{2}"), 
              Type = "parent") 

haplotype_slopes_single_intro_filtered <- haplotype_slopes_single_intro%>%
  dplyr::filter(rp == 1 & qtl == 2 & dp == 5) %>%
  dplyr::mutate(i = str_extract(haplotype, "i\\d+"), 
                i = str_extract(i, "\\d{2}"), 
                Type = "haplotype") %>%
  dplyr::rename(entryid2 = haplotype) %>%
  dplyr::full_join(ref_parent)

model <- lm(Estimate ~ i * Type, data = haplotype_slopes_single_intro_filtered)
summary(model)

anova_result <- anova(model)
print(anova_result)

reference <- rep(unique(subset(parents_slope, Parents == "NILASq0g01i14")$Estimate), 12)

mm <- lm(Estimate ~ i, data = haplotype_slopes_single_intro_filtered)
summary(mm)

contrast_14 <- emmeans(mm, ~ i) %>% 
  contrast(method = "eff", adjust = "none", ref = reference)

t_test <- t.test(haplotype_slopes_single_intro_filtered$Estimate, reference, paired = TRUE)
summary(t_test)
t_test$p.value
