library(tidyverse)
library(stats)
library(ggpubr)
library(grDevices)
library(scales)
library(RColorBrewer)

### root <- ("~/Dropbox-SE Dropbox/Randall Wisser/inrae_jd/")
root <- ("/Users/jdrouault/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/")

### paths 
genotypic_data_path <- paste0(root,"NILAS/main_folder/drouault/data/genotypic_data/introgression_depth/")
genes_info <-  paste0(root,"NILAS/main_folder/drouault/data/genotypic_data/genes_positions/")
haplotypes_path_blues <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_dta_em/") 
locus_type <-  paste0(root,"NILAS/main_folder/drouault/data/genotypic_data/locus_type/")
plot_path <- paste0(root, "NILAS/main_folder/drouault/plots/results/genotypic_data/") 


### download files 

###### get the files woith the number of introgression per types of qtl, for each group (qtl x rp x dp)
load_data <- function(path) { 
  files <- dir(path, pattern = '*.csv', full.names = TRUE)
  tables <- lapply(files, read.csv, sep = ";", header = T)
  do.call(rbind, tables)
}

introgression_positions <- load_data(locus_type)


# tropical_genes_introgressed <-read.csv2(file = paste0(genes_info, "2369CxCML277_tropical_genes_introgressed.csv"), 
                                                            # sep = ";" , header=TRUE)

haplotype_effect <- read.csv(paste0(haplotypes_path_blues, "haplotypes_effect_slope_20230711.csv"), sep = "," , header=TRUE)

#### get the effect of the introgressed genes on photoperiod sensitivity 

introgression_positions <- introgression_positions %>%
  tidyr::separate(Sample, into = c("entryid", 'code'), "s[0-9]") %>%
  dplyr::select(-code) 
  
genes_effect <- haplotype_effect %>%
  dplyr::filter(entryid %in%  introgression_positions$entryid) %>%
  dplyr::select(entryid, slope) %>%
  dplyr::full_join(introgression_positions) %>%
  dplyr::mutate(qtl = gsub(".*q|g.*","" ,entryid)) %>%
  # dplyr::select(-code) %>%
  # dplyr::filter(Type == "A") %>%
  dplyr::group_by(entryid, slope) %>%
  # dplyr::summarise(across(c(Nb_pr_qtl, Nb_pr_gwa, Nb_ft_gwa, Nb_ft_cg), sum)) %>%
  tidyr::pivot_longer(cols= c("Nb_pr_qtl", "Nb_pr_gwa", 'Nb_ft_gwa', "Nb_ft_cg"),
                      names_to = c(".value", "detection_technic"),
                      names_pattern = "(.*)_(.*_.*)") %>%
  dplyr::mutate(detection_technic = ifelse(grepl("pr_", detection_technic), "pr", "ft")) %>%
  dplyr::group_by(entryid, slope, qtl, detection_technic) %>%
  dplyr::summarise(Nb= sum(Nb))
  

date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)
write.csv(genes_effect, paste0(genes_info, "genes_PS_effect", "_",date_today, ".csv"), row.names = F)


#### Plot regression of number of QTL (pr and ft) introgressed against photoperiod-sensitivity

# colourCount = length(unique(genes_effect$qtl))
# getPalette = colorRampPalette(brewer.pal(12, "Set2"))
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_PS_nPS_rp1 <-  genes_effect %>%
  ggplot2::ggplot( aes(x=Nb, y = slope, fill = qtl)) + 
  geom_point( shape = 21, size = 5, stroke = 2) + #shape = detection_technic
  #scale_shape_manual(values= c4)) + #0,1,
  scale_fill_manual(values=cbPalette[c(1:5)]) +
  geom_point(data=subset(genes_effect, grepl("i12", entryid)), aes(x=Nb, y = slope), color = "red",  shape = 4, size=7, stroke = 2,  show.legend = F) +
  # scale_color_manual(values = c("red")) +
  guides(color = guide_legend(override.aes = list(colour = c(NA,NA, NA, NA), shape = c(NA,NA, NA, NA)) ) ) +
  # scale_color_manual(values = c("#5DA5DA", "#FAA43A" ,"#60BD68" ,"#B276B2")) + #"#5DA5DA", "#FAA43A" , #getPalette(colourCount))
  theme(
    legend.position="right",
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 15),
    panel.grid.major.x = element_line(color = "gray", size = .5, linetype = "dashed"),# explicitly set the horizontal lines (or they will disappear too) # remove the vertical grid lines
    panel.grid.major.y = element_line(color = "gray", size = .5, linetype = "dashed"),
    panel.grid.minor.y= element_blank(), 
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size= 18, angle=70, vjust = 0.5), 
    axis.text.y = element_text( size= 18),
    axis.title.x = element_text(size=25, margin = margin(t=11)),
    axis.title.y = element_text(size=25, margin = margin(r=11)),
    axis.ticks.length=unit(.25, "cm"),
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text = element_text(size = 20)) +
  ylab ( "photoperiod_sensitivity") +
  xlab("genes_introgressed") +
  facet_wrap(~detection_technic, scales = "free_x") +
  labs(shape = "genes_control") 
# geom_vline(aes(xintercept = (mean_slope)), color = "blue", size =1) 
plot_PS_nPS_rp1 <-plot_PS_nPS_rp1 + theme_classic(base_size = 25) +theme(axis.text.x = element_text(size= 30, angle=70, vjust = 0.5), 
                                                                          axis.text.y = element_text( size= 30),
                                                                          axis.title.x = element_text(size=30, margin = margin(t=11)),
                                                                          axis.title.y = element_text(size=30, margin = margin(r=11)),
                                                                          legend.text = element_text(size=30),
                                                                          legend.title = element_text(size=30))

date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)
ggsave(paste0(plot_path, "genes_introgress_corr_PS","_", date_today,".png"), plot_PS_nPS_rp1, width = 8, height = 7)


######## deprecated ########
### pivot longer introgression_positions

# introgression_positions1 <- introgression_positions %>%
#   # dplyr::select(-X) %>%
#   stats::setNames(gsub("CHROM", "_chr", names(.))) %>%
#   stats::setNames(gsub("r1$*", "r1_introgresssiontype", names(.))) %>%
#   stats::setNames(gsub("5_Breakpoint", "Breakpoint5", names(.))) %>%
#   stats::setNames(gsub("3_Breakpoint", "Breakpoint3", names(.))) %>%
#   tidyr::pivot_longer(cols= everything(), 
#                       names_to = c("lines",".value"),
#                       names_pattern = "(.*)_(.*)")
# values_to = c("chr", "breakpoints", "introgression"))
# tropical_genes_introgressed1 <- tropical_genes_introgressed %>%
#   tidyr::separate(Sample, into = c("entryid", 'code'), "s[0-9]") ## o,nly qtl 1 shpuld be removed when eloise send the count she made for this qtl

