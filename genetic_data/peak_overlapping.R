library(tidyverse)
library(stats)
library(ggpubr)
library(grDevices)
library(scales)
library(RColorBrewer)
library(readr)
library(readxl)
library(grid)


#### get Palette 
cbPalette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
               "#920000","#924900","#db6d00","#24ff24","#ffff6d")

#### root <- ("~/Dropbox-SE Dropbox/Randall Wisser/inrae_jd/")
root <- ("/Users/jdrouault/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/")

#### paths 
peak_overlap <-  paste0(root,"NILAS/main_folder/drouault/data/genotypic_data/peak_overlap/")
haplotypes_path_blues <- paste0(root, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_dta_em/") 
plot_path <- paste0(root, "NILAS/main_folder/drouault/plots/results/genotypic_data/") 

#### Download files 
haplotype_effect <- read.csv(paste0(haplotypes_path_blues, "haplotypes_effect_slope_20230711.csv"), sep = "," , header=TRUE)
overlapping_peak <- read.csv(paste0(peak_overlap, "20230801_overlapping_peak.csv"), sep = "," , header=TRUE)

#### Merge peak overlap and PS phenotypic data 
overlapping_peak <- overlapping_peak %>% ### technique de bourin à revoir quand j'aurai le fichier d'eloise 
  tidyr::separate(Sample, c("entryid","sample"), sep = "s",  remove = T) %>%
  dplyr::select(entryid, Peak_Overlapping) %>%
  dplyr::distinct(entryid, Peak_Overlapping)

haplotype_effect <- haplotype_effect %>%
  dplyr::full_join(overlapping_peak)

#### T-test on the distribution of the overlapping / non-overlapping haplotypes 
t <- t.test(haplotype_effect$slope[which(haplotype_effect$Peak_Overlapping == "True")], 
            haplotype_effect$slope[which(haplotype_effect$Peak_Overlapping == "False")], 
            var.equal = F)

(t$p.value)

counts = haplotype_effect %>% 
  dplyr::filter(entryid != "2369CxCML10C_ZmPR1_peak.txt") %>%
  group_by(Peak_Overlapping) %>% 
  count(Peak_Overlapping)

totalcounts = haplotype_effect %>% 
  dplyr::filter(entryid != "2369CxCML10C_ZmPR1_peak.txt") %>%
  # group_by(Peak_Overlapping) %>% 
  count(Peak_Overlapping)


box_plot <- haplotype_effect %>%
  dplyr::filter(entryid != "2369CxCML10C_ZmPR1_peak.txt") %>%
ggplot() +
  geom_boxplot(aes(
    x = Peak_Overlapping,
    y = slope, fill = Peak_Overlapping))+
  geom_text(data=counts,aes(x =Peak_Overlapping, y=5.5,label=n), col="blue", size =8) +
  # geom_text(data=totalcounts,aes(x=1.5,y=slope,label=value),col="blue")
  # geom_text(aes(x = t$p.value, y = 10, label = pval), size = 8) +
  labs(x = "overlap",
       y = "PS value") +
  scale_fill_brewer(palette = "Set1") + 
  theme(
      legend.position="right",
      panel.spacing = unit(0.5, "lines"),
      strip.text.x = element_text(size = 15),
      panel.grid.major.x = element_line(color = "gray", size = .5, linetype = "dashed"),# explicitly set the horizontal lines (or they will disappear too) # remove the vertical grid lines
      panel.grid.major.y = element_line(color = "gray", size = .5, linetype = "dashed"),
      panel.grid.minor.y= element_blank(), 
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size= 25, angle=70, vjust = 0.5), 
      axis.text.y = element_text( size= 25),
      axis.title.x = element_text(size=25, margin = margin(t=11)),
      axis.title.y = element_text(size=25, margin = margin(r=11)),
      axis.ticks.length=unit(.25, "cm"),
      strip.background = element_blank(), 
      strip.placement = "outside", 
      strip.text = element_text(size = 20)) 
  # ggtitle("Cluster Gene Expression across Phenotypes")

box_plot <- box_plot+ theme_classic(base_size = 20) +theme(axis.text.x = element_text(size= 25, angle=70, vjust = 0.5),
axis.text.y = element_text( size= 25),
axis.title.x = element_text(size=25, margin = margin(t=11)),
axis.title.y = element_text(size=25, margin = margin(r=11)),
legend.text = element_text(size=20),
legend.title = element_text(size=20))

date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)
ggsave(paste0(plot_path, "peak_overlapping","_", date_today,".png"), box_plot, width = 10, height = 8)


  ### Compare PS mean of haplotypes overlapping the QTL and introgression with the entire QTL 


cat_haplotype_effect <- haplotype_effect %>%
  dplyr::filter(Peak_Overlapping == "True" | grepl("i12", entryid)) %>%
  dplyr::rename(gb = rp) %>%
  dplyr::mutate(gb = as.character(gb), 
                qtl =as.character(qtl),
                Category = ifelse(grepl("i12", entryid), "full_QTL", "peak_overlap")) %>%
  dplyr::group_by(Category, gb, qtl) %>%
  dplyr::summarise(average_PS = mean(slope, na.rm =T)) 

histo_compare_PS <- cat_haplotype_effect %>%
  ggplot(aes(x = interaction(gb, qtl), y = average_PS, fill = Category )) +
  geom_col(position = "dodge", color="black", alpha = 0.5, width = 0.8) + 
  # annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = 3.250,
  #          label = rep(c("QTL 1", "QTL 2", "QTL 3", "QTL 4"), 1)) +
  # annotate("text", x = 1:8, y = cat_haplotype_effect$average_PS, label = rep(c("gb 1", "gb 2"), 4)) +
  # annotate("segment", x = x.scale, xend = x.scale, y = 0, yend = 4, lty = 2) + 
  ylab("average PS effect") +
  xlab("genetic background") +  theme_classic() +
  # theme(
  #       axis.title.x = element_blank(),
  #       axis.text.x = element_blank()) +
 ylim(min(cat_haplotype_effect$average_PS),4)
histo_compare_PS <- histo_compare_PS +coord_cartesian(clip ="off")
histo_compare_PS <- histo_compare_PS + theme_classic(base_size = 25) 
  
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)
ggsave(paste0(plot_path, "compare_overlap_fullintro","_", date_today,".png"), histo_compare_PS, width = 10, height = 8)  

g2 <- ggplot_gtable(ggplot_build(histo_compare_PS))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)

histo_compare_PS <- histo_compare_PS + theme_classic(base_size = 15) 

ggsave(paste0(plot_path, "average_PS_effect","_", date_today,".png"), histo_mean_PS, width = 12, height = 15)
