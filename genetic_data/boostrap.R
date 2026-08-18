# Import libaries
library(boot)

library(anyLib)
anyLib(c("data.table",
         "dplyr",
         "ggplot2",
         "gridExtra",
         "patchwork",
         "plyr",
         "reshape2",
         "stringr", 
         "forcats", 
         "ggpubr", 
         "stat", 
         "ggforce"
))

root <- "~/Documents/Drouault_etal_paper2/01_paperdev/ATLAS/paper/"#"~/Dropbox-SE Dropbox/Justine Drouault/inrae_jd/" #### replace by your path 
# root <- "~/Dropbox-SE Dropbox/Randall Wisser/inrae_jd/" #### replace by your path 
genetic_data_save_path <- paste0(root, "data/genotypic_data/NILs/")
source_path <- paste0(root, "data/cleaned_data/from_mixed_model/")
genetic_data_save_path <- paste0(root, "data/genotypic_data/NILs/")
fig_genotype <-  paste0(root,"Figures/genotypic_data/effect_genes_introgressionNB/")
dataset_dir = "dp6q30dr0.40maf0.01"
directory<- paste0("~/Documents/Drouault_etal_paper2/01_paperdev/ATLAS/paper/Stage_Eloise/Data/", dataset_dir, "/outputs/") 
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today)

ZmPR_list <- 4 
rp_list <- c(1, 2) #"2369", "LH123Ht"
dp_list <-  c(3, 4) #"CML277", "CML341"

haplotype_slope <- read.csv(paste0(source_path, "haplotype_slope_20240921.csv"), header = T)

phenotype_zmpr4 <- haplotype_slope %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", haplotype)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", haplotype)),
                qtl = as.numeric(gsub(".*q|g.*", "", haplotype))) %>%
  dplyr::rename(entryid2 = haplotype) %>%
  dplyr::filter(rp %in% rp_list) %>%
  dplyr::filter(dp %in% dp_list) %>%
  dplyr::filter(qtl %in% ZmPR_list) %>%
  #  dplyr::full_join(NILAS_newid) %>%
  #  dplyr::mutate(id_new = as.factor(id_new),
  #                   id_new = fct_rev(id_new)) %>%
  #  dplyr::arrange(desc(id_new)) %>%
  dplyr::filter(entryid2 %in% unique(geno_group$entryid)) %>% ## keep only the NILs of group A
  dplyr::mutate(
    i = str_extract(entryid2, "i\\d+"),
    i = str_extract(i, "\\d{2}"),
    i = factor(
      i,
      levels = c(
        "12",
        "11",
        "10",
        "09",
        "08",
        "07",
        "06",
        "05",
        "04",
        "03",
        "02",
        "01"
      )
    )
  )


groupA_NILs <- read.csv(paste0(genetic_data_save_path, "20240913Selected_NILs_after_filtration.csv"), header = T)

### here all this should be done before this script 
geno_group <- groupA_NILs  %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid)), 
                qtl = as.numeric(gsub(".*q|g.*", "", entryid))) %>%
  dplyr::filter(rp %in% rp_list) %>%
  dplyr::filter(dp %in% dp_list) %>% 
  dplyr::filter(qtl %in% ZmPR_list)

geno_data <- read.csv(
  paste0(
    file.path(directory),
    "/",
    "_",
    "20260728",
    "_",
    "mapping_data_NILAS",
    ".csv"
  ),
  header = T,
  sep = ","
) %>%
  dplyr::select(-X) %>% ## This step, i need to see why i have duplicates and after remove this two commandes
  dplyr::filter(Sample %in% unique(geno_group$entryid)) %>%
  dplyr::group_by(Sample) %>%
  dplyr::distinct() %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", Sample)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", Sample)),
                qtl = as.numeric(gsub(".*q|g.*", "", Sample))) 


## Create a list to define the order of the samples for the rest of the analysis
only_RP_sample <- geno_data %>%
  group_by(Sample) %>%
  filter(!any(geno == "DP")) %>%
  ungroup()## some sample only have RP we still keep in the lisf it they are not introgessed elsewhere in the genome 

list_RP_geno <-  as.list(unique(only_RP_sample$Sample))
sample_order <- geno_data %>%
  dplyr::filter(geno == "DP") %>%
  dplyr::group_by(Sample, rp) %>%
  dplyr::summarise(min_i_start = min(i.start, na.rm = TRUE)) 

## Alternative (separate per rp and then order within rp)
sample_order_rp_1 <- dput(levels(reorder(sample_order$Sample[sample_order$rp == 1], (sample_order$min_i_start[sample_order$rp == 1]))))


sample_order_rp <- dput(c(list_RP_geno, sample_order_rp_1))

## Order genotypes dataset
calculation_file <- geno_data %>%
  dplyr::filter(rp == 1 & !(Sample %in% list_RP_geno[[1]])) %>%
  dplyr::group_by(Sample, marker) %>%
  dplyr::mutate(dim_interval = i.end - i.start) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Sample = as.factor(Sample), 
                Sample = factor(Sample, levels = sample_order_rp)) #or sample_order

min_value <- abs(min(phenotype_zmpr4$Estimate, na.rm = T))


phenotype_zmpr4 <- phenotype_zmpr4 %>% 
  # dplyr::select(entryid2, Estimate, marker, id_new) %>%
  dplyr::rename(Sample  = entryid2) %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Estimate_transformed = Estimate + min_value, 
                Estimate_of_origin = ifelse(Estimate_transformed < min_value, 0, Estimate),
                Addition = ifelse(Estimate < 0, min_value - abs(Estimate) , min_value), 
                Sample = as.factor(Sample)) %>%
  tidyr::pivot_longer(c(Estimate_of_origin, Addition), 
                      names_to = "Type", 
                      values_to = "values") %>%
  # dplyr::arrange(desc(Estimate_transformed)) %>%
  dplyr::distinct()

n_iter <- 100
n_sample <- 12
ttest_pval  <- vector("list", n_iter)


### Begin function
for (i in seq_len(n_iter)) {
  # indices <- sample_order_rp_1
  indices_n <- sample(unique(calculation_file$Sample), n_sample)
  d <- calculation_file %>%
    dplyr::filter(Sample %in% indices_n) %>%
    droplevels %>%
    as.data.frame()

  get_segment <- function(difference) {
  seg <- integer(length(difference))
  seg[1] <- 1L

  for (i in 2:length(difference)) {
    seg[i] <- seg[i - 1]
    if (difference[i] > 0) {
      seg[i] <- seg[i] + 1L
    }
  }
  seg
}

## Calculate interval start and end point
dataframe_interval <- do.call(rbind, lapply(split(d, as.factor(d$Sample)), function(x){
  # print(unique(x$Sample))
  # x <- split(d,as.factor(d$Sample))[["NILASq4g31i10"]]
  x <- x[x$geno == "DP",]
  if(nrow(x) == 0){ ## don't use the sample in the calcul for interval as it does not contains dp
  }else{

    result <-  x%>%
      as.data.frame()%>%
      dplyr::arrange(marker)  %>%
      dplyr::mutate(
        diff = if_else(
          dplyr::row_number() == 1,
          0,
          i.start - dplyr::lag(i.end)),
        segment = get_segment(diff)) %>%
      dplyr::group_by(segment, rp) %>%
      dplyr::summarise(
        begin = min(i.start, na.rm = TRUE),
        end   = max(i.end, na.rm = TRUE),
        .groups = "drop"
      )
  }
}))

i.start <- dataframe_interval$begin
i.end <- dataframe_interval$end
interval_list <- c(i.start, i.end)
interval_list <- unique(interval_list)
interval_list <- sort(interval_list)

# Define in which interval the segment falls
dataset_final <- NULL
dataset_to_return <- NULL
table_for_fine_mapping <- do.call(
  rbind,
  lapply(
    split(d, d$Sample),
    function(x) {
      # x <- split(d,  as.factor(d$Sample))[["NILASq4g41i11"]]
      
      ## Run only if you want to get the interval per rp
      # rp <- unique(x$rp)
      # interval_list <- if (rp == 1) interval_list_rp1 else interval_list_rp2
      # begin_interval <- if (rp == 1) start_point_interval_ttest_rp1 else start_point_interval_ttest_rp2
      # end_interval <- if (rp == 1) end_point_interval_ttest_rp1 else end_point_interval_ttest_rp2



      dataset_final <- NULL
      # a <- 1

      for (i in 1:nrow(x)) {
        # print(i)

        # i <- 5  # remove after debugging

        dataset1 <- x[i, ]
        snp_start <- x$i.start[i]
        snp_end <- x$i.end[i]
        max_snp_end <- tail(x$i.end, 1)

        a <- 1
        interval_begin <- interval_list[[a]]
        interval_end <- interval_list[[a + 1]]
        end_max <- tail(interval_list, 1)

        # if(snp_start < begin_interval | snp_start > end_interval) {
        #   a <- NA
        #   dataset1[, "interval"] <- a
        #   dataset_final <- rbind(dataset_final, dataset1)}
        if (snp_start < interval_end & snp_end > interval_end & snp_start < end_max ) {

          while (snp_end > interval_end) {
            a <- a + 1
            interval_begin <- interval_list[[a]]
            interval_end <- interval_list[[a + 1]]
          }

          dataset1[, "interval"] <- a
          dataset_final <- rbind(dataset_final, dataset1)

        } else if (snp_start >= interval_end & snp_start < end_max) {

          while (snp_start >= interval_end) {
            a <- a + 1
            interval_begin <- interval_list[[a]]
            interval_end <- interval_list[[a + 1]]
          }

          dataset1[, "interval"] <- a
          dataset_final <- rbind(dataset_final, dataset1)

        } else if (snp_end == max_snp_end & max_snp_end < end_max) {

          a <- length(interval_list)
          dataset1[, "interval"] <- a
          dataset_final <- rbind(dataset_final, dataset1)

        } else if (snp_end > end_max) {

          a <- NA
          dataset1[, "interval"] <- a
          dataset_final <- rbind(dataset_final, dataset1)

        } else {

          dataset1[, "interval"] <- a
          dataset_final <- rbind(dataset_final, dataset1)

        }

      }

      dataset_final
    }
  )
)

## get the interval of cacta
snp_start <- 94430850
snp_end <- 94433495
a <- 1

# With both rp together
interval_begin <- interval_list[[a]]
interval_end <- interval_list[[a + 1]]
end_max <- tail(interval_list, 1)

if (snp_start < interval_end & snp_end > interval_end & snp_start < end_max ) {

  while (snp_end > interval_end) {
    # print(a)
    a <- a + 1
    interval_begin <- interval_list[[a]]
    interval_end <- interval_list[[a + 1]]
  }

  a


} else if (snp_start >= interval_end & snp_start < end_max) {

  while (snp_start >= interval_end) {
    a <- a + 1
    interval_begin <- interval_list[[a]]
    interval_end <- interval_list[[a + 1]]
  }

  a

}  else if (snp_end > end_max) {

  a <- NA
  a

} else {

  a

}


#Define the interval of bp in which no segment are NA, all the sample don't start and nd at the same bp so to avoid having interval for which the sample contains no values (NA), we define an interval in which all the sample have segment attributed to dp or rp
start_point_sample<- d %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(start_point = min(i.start, na.rm =T))

start_point_interval_ttest <- max(unique(start_point_sample$start_point))

end_point_sample<- d %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(end_point = max(i.end, na.rm =T))

end_point_interval_ttest <- min(unique(end_point_sample$end_point))

# Create plotting and t test dataset
table_for_fine_mapping <- table_for_fine_mapping %>%
  dplyr::filter(i.start > start_point_interval_ttest &
                  i.end < end_point_interval_ttest)

table_for_fine_mapping_interval_content <- table_for_fine_mapping %>%
  dplyr::group_by(Sample, interval) %>%
  dplyr::mutate(RP = sum(geno == "RP") / n(), DP = 1 - RP,
                i.start = min(i.start),
                i.end = max(i.end)) %>%
  dplyr::select(Sample, RP, DP, interval, i.start, i.end) %>%
  dplyr::distinct()

interval_with_DP <- table_for_fine_mapping_interval_content %>%
  dplyr::group_by(interval) %>%
  dplyr::summarise(sum_DP = sum(DP, na.rm = T),
                   sum_RP = sum(RP, na.rm = T)) %>%
  dplyr::filter(sum_DP == 0 | sum_DP == 1 | sum_RP == 0 | sum_RP == 1 | sum_DP <1 | sum_RP <1)

interval_for_plotting <- table_for_fine_mapping_interval_content %>%
  dplyr::filter(!(interval %in% unique(interval_with_DP$interval))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(interval) %>%
  dplyr::summarise(i.start = min(i.start, na.rm =T))

phenotype_zmpr4 <- subset(phenotype_zmpr4, Sample = indices_n)


ttest_pval[[i]] <- table_for_fine_mapping_interval_content %>%
  dplyr::full_join(phenotype_zmpr4[, c("Sample", "Estimate_transformed")] %>%
                     dplyr::distinct()) %>%
  tidyr::pivot_longer(cols = RP:DP,
                      names_to = "geno",
                      values_to = "geno_value") %>%
  dplyr::filter(geno_value != 0 &
                  !(interval %in% unique(interval_with_DP$interval))) %>%
  dplyr::group_by(interval) %>%
  dplyr::summarise(pval = t.test(Estimate_transformed ~ geno)$p.value, .groups = "drop", var.equal = FALSE) %>%
  dplyr::mutate(pval_adjusted = p.adjust(pval,method="BH"), 
                run = i,
                causal_variant_interval = a) %>%
  dplyr::full_join(interval_for_plotting)
# return(ttest_pval)

}


ttest_tot <- bind_rows(ttest_pval, .id = "column_label")

for(i in 1:3){
interval_value_plot <-  ttest_tot %>% #ttest_pval
  ggplot(aes( x= i.start, y = -log10(pval_adjusted))) + #value_inerval_tot #i.start
  geom_point() + 
  facet_wrap_paginate(~run, ncol = 6, nrow = 6, page =i) + 
  # geom_vline(data = ttest_tot %>% dplyr::group_by(run), aes(xintercept = i.start[interval == causal_variant]), color = "red") +
  # geom_vline(data = interval_value_plot %>%
  #              dplyr::group_by(rp) %>%
  #              dplyr::summarise(a_interval = i.start[interval == a]),  aes(group=rp, xintercept = i.start[interval == a]), color = "red") +
  geom_vline(aes(xintercept =  94430850), color = "red") +
  geom_vline(aes(xintercept =  94433495), color = "red") +
  geom_hline(yintercept = -log10(0.05), color = "blue") +
  ylab("-log10(p-value) of the interval") +
  xlab("Interval in \n(quantile positions, bp)") +
  theme(plot.background = element_blank()
        ,legend.position = "top"
        #,legend.direction="vertical"
        ,legend.key.size = unit(0.5, "cm")
        ,legend.text = element_text(size = 10)
        #,aspect.ratio = 1
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", linetype = "dashed")
        ,panel.grid.major.y = element_line(colour="grey")
        ,panel.border = element_rect(colour="grey", linewidth = 0.5, fill=NA)
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold", vjust=-12.375, margin=margin(b=43))
        ,axis.title.x = element_text(size=14, vjust=-4.5)
        ,axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)
        ,axis.title.y = element_text(size=14)
        ,axis.text.y =element_text(size=12)
        #top, right, bottom, left
        #,plot.margin=unit(c(0.125,0.5,0.72,0.5), "cm")
        ,plot.margin=unit(c(0.02,0.5,0.72,0.5), "cm")
  )
ggsave(plot = interval_value_plot, filename = paste0('/Users/lps-drouault/Documents/Drouault_etal_paper2/01_paperdev/ATLAS/paper/plots/sample_test_intreval_page_', i, '.pdf'), width = 14, height = 10)}


significative_run <- ttest_tot %>%
  dplyr::filter(pval_adjusted < 0.05) %>%
  dplyr::select(run) %>%
  dplyr::distinct() %>%
  dplyr::summarise(count = n())
  

