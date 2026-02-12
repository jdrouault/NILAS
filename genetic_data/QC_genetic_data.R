
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggpmisc)
library(agricolae)
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
phenotype_path <- paste0(root_alternative, "NILAS/main_folder/drouault/data/cleaned_data/blues/blues_day_dta_04Avr11/") ### need to change the path to access either the blues relatedd to gdd_dta_em or to gdd_dta 
fig_consistency <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/constistency_replicates/")
fig_genes_introgression <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/genotypic_data/effect_genes_introgressionNB/")
fig_phenotype <- paste0(root_alternative, "NILAS/main_folder/drouault/paper/Figures/phenotypic_data/")
file_path <- paste0(root_alternative, "/NILAS/main_folder/drouault/paper/Files/")
genetic_data_save_path <- paste0(root_alternative, "/NILAS/main_folder/drouault/data/genotypic_data/NILs/")

###### ======= Download the file 

data <-  read.delim(paste0(NILAS_genotypes_path, "nilas_compiledsummary_02Jul24.txt"), sep = "\t") %>% #nilas_compiledsummary_29Apr24.txt
  dplyr::mutate(IL_nb = str_extract(Sample, "i[0-9]{2}"), 
                qtl = str_extract(Sample, "q[0-9]{1}")) 

complementary_lines <- data %>%
  dplyr::mutate(check_C = grepl( "C", Sample) ) %>%
  dplyr::filter(check_C == T) %>%
  tidyr::unite(ID, c(Scheme, IL_nb, qtl), remove = F)


genotypic_data <- data %>%
  dplyr::mutate(check_C = grepl( "C", Sample) ) %>%
  dplyr::filter(check_C == F) %>%
  dplyr::mutate(entryid = gsub("s[0-9]p[0-9]r[0-9]", "", Sample)) %>%
  dplyr::mutate(biol_rep = str_extract(Sample, "p[0-9]"), 
                tech_rep = str_extract(Sample, "r[0-9]") )


missing_lines <- genotypic_data %>%
  dplyr::select(entryid) %>%
  dplyr::distinct() %>%
  dplyr::mutate(dp = sub(".*g(\\d+).*", "\\1", entryid), # Extract first number after 'g'
                rp = sub(".*g\\d+(\\d+).*", "\\1", entryid), # Extract second number after 'g'
                ID = sub(".*i(\\d+).*", "\\1", entryid), 
                ID = as.numeric(ID), 
                diff = ID - lag(ID)) #%>%# Extract two numbers after 'i') %>%
  dplyr::filter(diff == 2 | diff == -10)
  
  
###### ======= Checks overlapping PS and FL genes 

genes_files <- read.delim(paste0(haplot_imput_data, "ZmFTloci_agpv4_tempmod2.txt"), sep = "\t")

new_category_genes <- genes_files %>% #### remove duplicated genes between FL and PR 
  dplyr::mutate(dummi = ifelse(locus_type == "ft_gwa" | locus_type == "pr_gwa", "1", "2" )) %>%
  dplyr::mutate(locus_type = ifelse(dummi == 1 & agpv4_start == lead(agpv4_start), "ft_pr_gwa", locus_type)) %>%
  dplyr::mutate(across(locus_type:dummi, ~ifelse(agpv4_start == lag(agpv4_start), NA , .))) %>%
  dplyr::filter(!is.na(dummi)) %>%
  dplyr::select(-dummi)
 
# summary_shared_position <- genes_files %>%
#   dplyr::group_by(agpv4_start, agpv4_stop ) %>%
#   dplyr::summarise(count = n()) %>%
#   dplyr::filter(count >1) ### same start and stop used more than one time ,
# 
# shared_position <- genes_files %>%
#   dplyr::filter(agpv4_start %in% summary_shared_position$agpv4_start)

# write.csv(new_category_genes, paste0(haplot_imput_data, date_today, "ZmFTloci_agpv4_tempmod2.txt"), row.names = F)


###### ======= Checks before QC

#### Check that i have at least one biological rep per ILs
toto <- data.frame(unique(genotypic_data [,c("Sample", "entryid", "IL_nb")]))
obj <- with(toto, table(entryid, IL_nb))

check_quality1 <- genotypic_data %>%
  dplyr::filter(Scheme == "CML277 x 2369") %>%
  dplyr::select(Sample, Scheme, DP_BG:RP_FG, IL_nb) %>%
  tidyr::pivot_longer(cols = DP_BG:RP_FG,
                      names_to = c("categories"),
                      values_to = "count") %>%
ggplot() +
  geom_boxplot(aes(x= count , y = categories)) +
  geom_jitter(aes(x= count , y = categories,  color = Sample)) +
  facet_wrap(~IL_nb) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1), 
        legend.position = "none"
)

check_quality1

###### ======= Check consistancy of the data using technical replicates 

####### Create a matrix for pairwise comparison of the tech rep 

list_entryid <- unique(genotypic_data$entryid)
matrix_tech_rep <- NULL

for(i in 1:length(list_entryid)){
  rep <- NULL
  no_tech_rep <- NULL
   # i <- 14
 x <- subset(genotypic_data, entryid == list_entryid[[i]])
  
 summary_x <- x %>%
    dplyr::group_by(biol_rep) %>%
    dplyr::summarise(distinct_tech_rep = n_distinct(tech_rep)) %>%
    dplyr::filter(distinct_tech_rep >1) %>%
   dplyr::ungroup()
 
no_tech_rep <- x %>%
    dplyr::filter(!(biol_rep %in% unique(summary_x$biol_rep))) %>%
      dplyr::select(entryid, DP_BG:RP_FG, biol_rep ) %>%
      dplyr::mutate(across(c(DP_BG:RP_FG), ~  NA, .names = "{col}_diff"), 
                    test2 = "no_rep", 
                    nb_tech_rep = NA) %>%
      dplyr::select(entryid, DP_BG_diff:RP_FG_diff, biol_rep, test2, nb_tech_rep )

if(nrow(summary_x)>0){
  
  nb <- unique(summary_x$distinct_tech_rep)
  rep <- x %>%
    dplyr::filter(biol_rep %in% unique(summary_x$biol_rep)) %>% 
    dplyr::group_by(biol_rep) %>% 
    dplyr::summarise(test = combn(tech_rep,2,paste, simplify = FALSE), 
                     across(c(DP_BG:RP_FG), ~ combn(.,2,diff), .names = "{col}_diff"))%>%
    dplyr::mutate(test2 = unlist(lapply(test, function(x)paste(x,collapse="-")))) %>%
    dplyr::select(-test) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(entryid = unique(x$entryid), 
                  nb_tech_rep = summary_x$distinct_tech_rep[which(summary_x$biol_rep == biol_rep)])%>%
    dplyr::select(entryid, DP_BG_diff:RP_FG_diff, biol_rep, test2, nb_tech_rep ) %>%
    dplyr::ungroup()
  
  no_tech_rep <- no_tech_rep %>%
    dplyr::full_join(rep)
  }
  matrix_tech_rep <- rbind.data.frame(matrix_tech_rep, no_tech_rep)
}


matrix_tech_rep <- matrix_tech_rep %>% ### transform everything in abs value 
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~(abs(.))))

count_tech_rep <- matrix_tech_rep %>% ### Count the number of tech rep 
  dplyr::filter(test2 != "no_rep") %>%
  dplyr::select(entryid, nb_tech_rep)%>%
  dplyr::distinct() %>%
dplyr::summarise(sum = sum(nb_tech_rep, na.rm =T))


########## Plots 

#### Plot overall consistency 

consistency_techrep <- matrix_tech_rep %>%
  dplyr::filter(test2 != "no_rep") %>%
  dplyr::select(entryid,test2, DP_BG_diff:RP_FG_diff)%>%
tidyr::pivot_longer(cols = DP_BG_diff:RP_FG_diff,
                    names_to = c("categories"),
                    values_to = "count of pair-wise comparison") %>%
  dplyr::mutate(categories = as.factor(categories)) %>%
  ggplot() +
  facet_wrap(~categories) +
  geom_histogram(aes(x= count)) + #stat="identity", aes(fill = categories)
  theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        legend.position = "none"
  )+ xlab('divergence measured in point of percentage')

# ggsave(paste0(fig_consistency, "consistency_techrep","_", date_today,".pdf"), consistency_techrep, width = 8, height =6)


#### Plot % of lines divergent for the selected variables

summary_tech1 <-  matrix_tech_rep %>%
  dplyr::filter(test2 != "no_rep") %>%
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~  100 - .)) %>%
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. == 100), .names = "{col}_1")) %>% # here count the number of tech that are similar # tech rep that differs sum(. != 0)
 # dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 3), .names = "{col}_supp3")) %>% # here count the number of tech rep that differs and for which the diff is greater than 3
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 95  & . < 100), .names = "{col}_0.99_0.95")) %>% # here count the number of tech rep that differs and for which the diff between value > 0 and 5%
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 90 & . < 95), .names = "{col}_0.95_0.90")) %>% # here count the number of tech rep that differs and for which the diff is greater value > 0 and 10%
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 80 & . < 90), .names = "{col}_0.90_0.80")) %>%
  dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 0 & . < 80), .names = "{col}_0.80_0")) %>%
  dplyr::select(-c(entryid, biol_rep, test2, nb_tech_rep, DP_BG_diff:RP_FG_diff)) %>%
  dplyr::distinct() %>%
  #dplyr::mutate(across(DP_BG_diff_count:`RP_FG_diff_supp3`, ~ sum((. /38)*100), .names = "{col}_percent")) %>% ### 38 is the number of comparisons is tot number of tech obtained above
  dplyr::mutate(total_comparison = 38)


### à refiare avec summary_tech1
# check_percentage_unconsistency_tech <- summary_tech %>%
#   dplyr::select(DP_BG_diff_count_percent:total_comparison) %>%
#   tidyr::pivot_longer(cols = DP_BG_diff_count_percent:total_comparison,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(across(categories, ~gsub("_count", "", .) )) %>%
#   dplyr::mutate(across(categories, ~gsub("_percent", "", .) )) %>%
#   dplyr::mutate(across(categories, ~gsub("_diff", "", .) ))
#   
# list2 <- c(check_percentage_unconsistency_tech$categories)


# 
# check_percentage_unconsistency_tech_plot <- check_percentage_unconsistency_tech %>%
#   dplyr::mutate(categories = as.factor(categories), 
#                 categories = factor(categories, levels = list2)) %>%
#   dplyr::filter(categories != "total_comparison") %>%
#   ggplot(aes(x= categories, y = count)) +
#   geom_bar(stat="identity", color = "black", fill = "lightblue", size = 1) + #aes(fill = categories)
#   theme_classic() + 
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         legend.position = "none", 
#         panel.grid.major.y = element_line()
#   ) +
#   ylab("% of technical replicates \n that are similar")
# 
# check_percentage_unconsistency_tech_plot

# ggsave(paste0(fig_consistency, "percentage_unconsistancy_techrep","_", date_today,".pdf"), check_percentage_unconsistency_tech_plot, width = 8, height = 6)

matrix_tech_rep_pivoted <- matrix_tech_rep %>%
  tidyr::pivot_longer(DP_BG_diff:RP_FG_diff, 
                      names_to = c("categories"),
                      values_to = "count")

consistency_tech_plot <- matrix_tech_rep_pivoted %>%
  ggplot(aes(x= categories, y = count)) +
  geom_boxplot(color = "black", fill = "lightblue", size = 1) + #aes(fill = categories)
  geom_point(stat="identity", color = "black", fill = "lightblue", size = 1) + #aes(fill = categories)
  theme_classic() + 
  theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.y = element_text(size=20),
        legend.position = "none", 
        panel.grid.major.y = element_line()
  ) +
  ylab("value of divergence of technical replicates")


##### Summary table for overleaf 

names(summary_tech1) <- gsub("diff_", "", names(summary_tech1))

Table_overleaf_tech <- summary_tech1 %>% #summary_tech
  dplyr::select(DP_BG_1:RP_FG_0.80_0) %>% #DP_BG_diff_count:RP_FG_diff_supp3_percent
  tidyr::pivot_longer(cols = DP_BG_1:RP_FG_0.80_0, #DP_BG_diff_count:RP_FG_diff_supp3_percent
                        names_to = c("variable", ".value"),
                        names_pattern =  "(.*G)_(.*)$") %>%
    # dplyr::rename(total_differing = diff_count, sup_3per =  diff_supp3, total_percent = diff_count_percent, sup_3per_percent = diff_supp3_percent) %>%
  # dplyr::rename(total_similar = diff_count, sup_3per =  diff_supp3, total_percent = diff_count_percent, sup_3per_percent = diff_supp3_percent) %>%
  dplyr::mutate(variable = dplyr::recode(variable, DP_FG = "Donor parent in foreground (DP_FG)", DP_BG = "Donor parent in background (DP_BG)" ,
                                  RP_FG = "Recurrent parent in foreground (RP_FG)", RP_BG = "Recurrent parent in background (RP_BG)",
                                  Het_BG = "Heterosigozity background (Het_BG)", Het_FG = "Heterosigozity foreground (Het_FG)", NA_BG = "missing values background (NA_BG)", 
                                  NA_FG = "missing values foreground (NA_FG)"))
               # across(total_percent:sup_3per_percent, ~round(., 2))) 


#write.csv(Table_overleaf_tech, paste0(save_path, date_today, "Summary_consistency_techrep.csv"), row.names = F)
  
  tech_rep_summary <-matrix_tech_rep %>%
    dplyr::filter(test2 != "no_rep") %>%
    dplyr::mutate(tot = 38) #got from count_tech_rep
  
  # tech_rep_constitency <- matrix_tech_rep %>%
  #   dplyr::filter(test2 != "no_rep") %>%
  #   dplyr::group_by(entryid, biol_rep) %>%
  #   dplyr::mutate(mean = mean(Het_BG_diff, na.rm = T)) %>%
  #   dplyr::mutate(constistency = ifelse(mean < 0.5, "consistent", "unconsistent"))  
  # 

###### ======= Check consistancy of the data using biological replicates 

  #### Chose one tech rep per biological sample 
  
  tech_rep_summary <- tech_rep_summary %>%
    tidyr::unite(ID1, c(entryid, biol_rep), sep = "", remove = F) 
  
  subset <- genotypic_data %>%
  dplyr::filter(entryid %in% unique(tech_rep_summary$entryid)) 

  select_tech_rep <- genotypic_data %>%
    tidyr::unite(ID1, c(entryid, biol_rep), sep = "", remove = F) %>%
    dplyr::filter(ID1 %in% unique(tech_rep_summary$ID1)) %>% ### take the entryid and biol rep with tech rep 
    dplyr::group_by(entryid, biol_rep) %>%
    dplyr::mutate(Mean_Het_BG = mean(Het_BG, na.rm = T)) %>% ### get the mean of het per biol_rep
    dplyr::group_by(entryid, biol_rep, tech_rep) %>%
    dplyr::mutate(tech_diff = (abs(Het_BG - Mean_Het_BG))) %>% ### the diff between the biol rep mean Het and the tech rep het 
    dplyr::ungroup()%>%
    dplyr::group_by(entryid, biol_rep) %>%
    dplyr::mutate(min = min(tech_diff))%>% ### get the min value between the tech rep and the biol rep mean Het 
    dplyr::summarise(best_tech_rep = tech_rep[which(tech_diff == min )], 
                     value_Het_BG = Het_BG[which(tech_diff == min )]) %>% #### select the tech that has this value 
    dplyr::group_by(entryid, biol_rep) %>%
    dplyr::filter(value_Het_BG == min(value_Het_BG)) %>%
    tidyr::unite(ID, c(entryid, biol_rep, best_tech_rep), sep = "",remove = F)

    ### check if within the tech rep that have the same value for Het, they have different value for others variables 
    subset <- subset %>%
      tidyr::unite(ID, c(entryid, biol_rep, tech_rep), sep = "",remove = F) %>%
    dplyr::filter(ID %in% unique(select_tech_rep$ID)) %>%
      dplyr::full_join(select_tech_rep) ### all values the same across the others variable 

    ### when there is two tech remaining per biol rep and they have the same values for all the variable select the first tech rep 
    
    selected_tech_rep <- select_tech_rep %>% 
      group_by(entryid, biol_rep) %>%
      dplyr::mutate(across(ID:value_Het_BG, ~ dplyr::first(.)
      )) %>%
      dplyr::distinct()

  ###### Create a matrix for pairise comparison of biol rep 
    
    matrix_biol_rep <- NULL
  for(i in 1:length(list_entryid)){
    print(i)
    rep <- NULL
    no_tech_rep <- NULL
    i <- 19
    x <- subset(genotypic_data, entryid == list_entryid[[i]])
    
    summary_x <- x %>%
      dplyr::group_by(biol_rep) %>%
      dplyr::summarise(distinct_tech_rep = n_distinct(tech_rep)) %>%
      dplyr::ungroup()
    
    if(nrow(summary_x)>1){
      
      data <- NULL
      
      x <- x %>%
        tidyr::unite(ID, c(entryid, biol_rep, tech_rep), sep = "",remove = F)
      
      ### quand il y a pls rep tech dans un rep biol choisir celui qui est dans selected_tech_rep
      x<- x %>% split(.$biol_rep) %>% 
        map(function(x) if (nrow(x) > 1) x %>% filter(x$ID %in% unique(selected_tech_rep$ID)) else x) %>% bind_rows()
      
      nb <- unique(summary_x$distinct_tech_rep)
      rep <- x %>%
        dplyr::summarise(test = combn(biol_rep,2,paste, simplify = FALSE), 
                         across(c(DP_BG:RP_FG), ~ combn(.,2,diff), .names = "{col}_diff"))%>%
        dplyr::mutate(test2 = unlist(lapply(test, function(x)paste(x,collapse="-")))) %>%
        dplyr::select(-test) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(entryid = unique(x$entryid), 
                      nb_biol_rep = nrow(summary_x))%>%
        dplyr::select( entryid, DP_BG_diff:RP_FG_diff, test2, nb_biol_rep ) %>%
        dplyr::ungroup()
      
      data  <- rep
    } else{
      
      data <- NULL 
      ### When there is no biol-rep 
      no_biol_rep <- x %>%
        dplyr::select(entryid, DP_BG:RP_FG, biol_rep ) %>%
        dplyr::mutate(across(c(DP_BG:RP_FG), ~  NA, .names = "{col}_diff"), 
                      test2 = "no_rep_biol", 
                      nb_biol_rep = NA) %>%
        dplyr::select( entryid, DP_BG_diff:RP_FG_diff,test2, nb_biol_rep)
      
      data <- no_biol_rep
    }
    
    matrix_biol_rep <- rbind.data.frame(matrix_biol_rep, data)
    
  }

    
    matrix_biol_rep <- matrix_biol_rep %>% ### transform everything in abs value 
      dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~(abs(.))))
    
    count_biol_rep <- matrix_biol_rep %>% ### count number of biol rep 
      dplyr::filter(test2 != "no_rep_biol") %>%
      dplyr::select(entryid, nb_biol_rep)%>%
      dplyr::distinct() %>%
      dplyr::summarise(sum = sum(nb_biol_rep, na.rm =T))
    
   biol_rep_summary <-matrix_biol_rep %>%
      dplyr::filter(test2 != "no_rep_biol") %>%
      dplyr::mutate(tot = 148) #got from count_biol_rep
   
  length(biol_rep_summary$test2)
   
######### Plots 
   
#### Plot overall consistency 
    
  consistency_biol_rep <- matrix_biol_rep %>%
      dplyr::filter(test2 != "no_rep_biol") %>%
      dplyr::select(entryid,DP_BG_diff:RP_FG_diff)%>%
      tidyr::pivot_longer(cols = DP_BG_diff:RP_FG_diff,
                          names_to = c("categories"),
                          values_to = "count") %>%
      dplyr::mutate(categories = as.factor(categories)) %>%
      ggplot() +
      facet_wrap(~categories) +
      geom_histogram(aes(x= count)) + #stat="identity", aes(fill = categories)
      theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
            axis.title.x = element_text(size=20),
            axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
            axis.title.y = element_text(size=20),
            strip.text.x = element_text(size = 20),
            strip.text.y = element_text(size = 20),
            legend.position = "none"
      )+ xlab('absolute count')
   
   # ggsave(paste0(fig_consistency, "consistency_biolrep","_", date_today,".pdf"), consistency_biol_rep, width = 8, height =6)
   

   #### Plot % of lines divergent for the selected variables
   
   summary_biol1 <-  matrix_biol_rep %>% #summary_biol
     dplyr::filter(test2 != "no_rep_biol") %>%
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~  100 - .)) %>%
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. == 100), .names = "{col}_1")) %>% # here count the number of biol that are 100% similar 
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 95  & . < 100), .names = "{col}_0.99_0.95")) %>% # here count the number of tech biol that differs and for which the diff between value > 0 and 5% <
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 90 & . < 95), .names = "{col}_0.95_0.90")) %>% # here count the number of tech rep that differs and for which the diff is greater value > 5% and 10% <
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 80 & . < 90), .names = "{col}_0.90_0.80")) %>%
     dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 0 & . < 80), .names = "{col}_0.80_0")) %>%
     dplyr::select(-c(entryid, test2, nb_biol_rep, DP_BG_diff:RP_FG_diff)) %>%
     dplyr::distinct() %>%
     # dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. == 0), .names = "{col}_count")) %>% # here count the number of biol that are similar # biol rep that differs sum(. != 0)
     # dplyr::mutate(across(DP_BG_diff:RP_FG_diff, ~ sum(. > 3), .names = "{col}_supp3")) %>% # here count the number of biol for which the difference in content is greater than 3
     # dplyr::select(-c(entryid, test2, nb_biol_rep, DP_BG_diff:RP_FG_diff)) %>%
     # dplyr::mutate(across(DP_BG_diff_count:`RP_FG_diff_supp3`, ~ sum((. /92)*100), .names = "{col}_percent")) %>% ### 38 is the number of comparisons is tot number of tech obtained above
     dplyr::mutate(total_comparison = 148)

   #### à refaire avec summary_biol1
   # check_percentage_unconsistency_biol <- summary_biol %>%
   #   dplyr::select(DP_BG_diff_count_percent:total_comparison) %>%
   #   tidyr::pivot_longer(cols = DP_BG_diff_count_percent:total_comparison,
   #                       names_to = c("categories"),
   #                       values_to = "count") %>%
   #   dplyr::distinct() %>%
   #   dplyr::mutate(across(categories, ~gsub("_count", "", .) )) %>%
   #   dplyr::mutate(across(categories, ~gsub("_percent", "", .) )) %>%
   #   dplyr::mutate(across(categories, ~gsub("_diff", "", .) ))
   
   # list2 <- c(check_percentage_unconsistency_biol$categories)
   # 
   # check_percentage_unconsistency_biol_plot <- check_percentage_unconsistency_biol %>%
   #   dplyr::mutate(categories = as.factor(categories), 
   #                 categories = factor(categories, levels = list2)) %>%
   #   dplyr::filter(categories != "total_comparison") %>%
   #   ggplot(aes(x= categories, y = count)) +
   #   geom_bar(stat="identity", color = "black", fill = "lightblue", size = 1 ) + #, aes(fill = categories)
   #   theme_classic() + 
   #   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
   #         axis.title.x = element_text(size=20),
   #         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
   #         axis.title.y = element_text(size=20),
   #         legend.position = "none", 
   #         panel.grid.major.y = element_line()
   #   ) +
   #   ylab("% of biological replicates \n that are similar")
   
   check_percentage_unconsistency_biol_plot
   
   # ggsave(paste0(fig_consistency, "percentage_unconsistancy_biolrep","_", date_today,".pdf"), check_percentage_unconsistency_biol_plot, width = 8, height =6)

##### Table for overleaf 
   
   names(summary_biol1) <- gsub("diff_", "", names(summary_biol1))
   
   
   
  Table_overleaf_biol <-  summary_biol1 %>% #summary_biol
     dplyr::select(DP_BG_1:RP_FG_0.80_0) %>% #DP_BG_diff_count:RP_FG_diff_supp3_percent
     tidyr::pivot_longer(cols = DP_BG_1:RP_FG_0.80_0,
                         names_to = c("variable", ".value"),
                         names_pattern =  "(.*G)_(.*)$") %>%
     # dplyr::rename(total_differing = diff_count, sup_3per =  diff_supp3, total_percent = diff_count_percent, sup_3per_percent = diff_supp3_percent) %>%
     # dplyr::rename(total_similar = diff_count, sup_3per =  diff_supp3, total_percent = diff_count_percent, sup_3per_percent = diff_supp3_percent) %>%
     dplyr::mutate(variable = dplyr::recode(variable, DP_FG = "Donor parent in foreground (DP_FG)", DP_BG = "Donor parent in background (DP_BG)" ,
                                            RP_FG = "Recurrent parent in foreground (RP_FG)", RP_BG = "Recurrent parent in background (RP_BG)",
                                            Het_BG = "Heterosigozity background (Het_BG)", Het_FG = "Heterosigozity foreground (Het_FG)", NA_BG = "missing values background (NA_BG)", 
                                            NA_FG = "missing values foreground (NA_FG)"))
                   # across(total_percent:sup_3per_percent, ~round(., 2))) 
   
  # write.csv(Table_overleaf_biol, paste0(save_path, date_today, "Summary_consistency_biolrep.csv"), row.names = F)
   
###### ======= Select the biol rep and biol samples that should be kept for further analysis 
   
 ##### First selection of one biol rep representing one line

 ## Combine the selected tech rep with all the others lines ID that don't have tech rep
   
   all_without_tech_rep <- matrix_tech_rep %>% ### biol rep without tech rep
     dplyr::filter(test2 == "no_rep") %>%
     tidyr::unite(ID1, c(entryid, biol_rep), sep = "", remove = F)
   
  list_ID_with_genotype_data <- genotypic_data %>% #### line without tech rep 
    tidyr::unite(ID1, c(entryid, biol_rep), sep = "", remove = F) %>%
    dplyr::filter(ID1 %in% unique(all_without_tech_rep$ID1)) %>%
    tidyr::unite(ID, c(entryid, biol_rep, tech_rep), sep = "", remove = F) %>%
    dplyr::select(ID) %>%
    dplyr::distinct()
  
  list <- selected_tech_rep%>% dplyr::ungroup() %>% dplyr::select(ID) %>% distinct()
     
  list_all_ID <-  rbind(list_ID_with_genotype_data,list)
  
   x<- genotypic_data %>% ### filter on selected tech and line without rep, remains line with several biol rep 
     tidyr::unite(ID, c(entryid:tech_rep), sep = "", remove = F) %>%
     tidyr::unite(splitter, c(entryid:biol_rep), sep = "", remove = F) %>%
     split(.$splitter) %>% 
     map(function(x) if (nrow(x) > 1) x %>% filter(x$ID %in% unique(list_all_ID$ID)) else x) %>% bind_rows()

   selected_biol_rep <- x %>% #### select by filtering by variables for which the lines diverge the most and that are the most important to conserve 
     split(.$entryid) %>% 
     map(function(x) if (nrow(x) > 1) x %>%filter((DP_FG == max(DP_FG))) %>%
     filter(Het_BG == min(Het_BG))  %>% filter((DP_BG == min(DP_BG))) %>% 
    filter(Het_FG == min(Het_FG)) else x) %>% bind_rows()
   
   selected_biol_rep <- selected_biol_rep %>% #### make sure we have removed all the biol rep and select only per line
     split(.$entryid) %>% 
     map(function(x) if (nrow(x) > 1) x %>% filter(row_number()==1) else x) %>% bind_rows() 
   
   ##### Second step filter lines with too high heterozygosity and introgression outside of the ZMPR zone 
   
   summary_too_much_intro <- selected_biol_rep %>%
     dplyr::group_by(Sample) %>%
     dplyr::summarise(across(c(Nb_ZmPR1:Nb_ZmPR4), ~sum( . != 0))) %>% 
     dplyr::mutate(total = Nb_ZmPR1+ Nb_ZmPR2 + Nb_ZmPR3 + Nb_ZmPR4) %>% ### get the total number of introgressed ZPRM region per line 
     dplyr::filter(total > 1) %>% ### when there are more than 1 non-introgressed ZMPR region per lines, it's means that at least 2 ZMPR regions are introgressed so we should discard the line, cause impossible to determinate which ZMPR introgressed region causes the effect
     dplyr::select(Sample)
   
   ### Filter for heterozygosity 
   filtered_selected_biol_rep <-  selected_biol_rep %>%
     # dplyr::filter(!(Sample %in% summary_too_much_intro$Sample)) %>% #####!!!!!!!! put back this command to filter the lnes with several introgression 
     dplyr::filter(!grepl("q2g21i12s1p3r1",Sample) & !grepl("q3g31i01s2p1r1",Sample)
                   & !grepl("q2g12i04s1p3r1",Sample) & 
                     !grepl("q2g52i02s1p1r1",Sample)) ### plots from MacH 
   
   selected_biol_rep %>%
     dplyr::filter(grepl("q2g21i12s1p3r1",Sample) | grepl("q3g31i01s2p1r1",Sample)
                   | grepl("q2g12i04s1p3r1",Sample) |
                      grepl("q2g52i02s1p1r1",Sample)) ### plots from MacH 
   
   # write.csv(filtered_selected_biol_rep, paste0(genetic_data_save_path, date_today, "Selected_NILs_after_filtration_for_heterozigosity.csv"), row.names = F)
   
  
   ### Filter for several introgression 
   filtered_selected_biol_rep_all <-  filtered_selected_biol_rep %>%
     dplyr::filter(!(Sample %in% summary_too_much_intro$Sample))
   
   summary_qtl_insertion <- filtered_selected_biol_rep_all %>% # genotypic_data %>%
     dplyr::rowwise() %>%
     dplyr::mutate(sum_gene_pr  = sum(Nb_pr_qtl, Nb_pr_gwa, na.rm = T),
                   sum_gene_fl  = sum(Nb_ft_gwa, Nb_ft_cg, na.rm = T), 
                   sum_introgression = sum(Nb_ZmPR1:Nb_ZmPR4, na.rm = T),
                   several_intro = list(ifelse(across(Nb_ZmPR1:Nb_ZmPR4, ~. > 1), "T", "F"))) %>%
     dplyr::mutate(qtl = gsub("q", "", qtl)) %>%
     as.data.frame() %>%
     dplyr::select(-several_intro)
   
   # write.csv(summary_qtl_insertion, paste0(genetic_data_save_path, date_today, "Selected_NILs_after_filtration.csv"), row.names = F)

   
   # summary_too_much_intro <- filtered_selected_biol_rep_all %>%
   #   dplyr::group_by(Sample) %>%
   #   dplyr::summarise(across(c(Nb_ZmPR1:Nb_ZmPR4), ~sum( . != 0))) %>% 
   #   dplyr::mutate(total = Nb_ZmPR1+ Nb_ZmPR2 + Nb_ZmPR3 + Nb_ZmPR4) %>% ### get the total number of introgressed ZPRM region per line 
   #   dplyr::filter(total > 1) #%>% ### when there are more than 1 non-introgressed ZMPR region per lines, it's means that at least 2 ZMPR regions are introgressed so we should discard the line, cause impossible to determinate which ZMPR introgressed region causes the effect
   # # dplyr::select(Sample)
   
   # write.csv(summary_too_much_intro, paste0(genetic_data_save_path, date_today, "lines_double_triple_quadruple_intro.csv"), row.names = F)
   


###### ======= Assemble phenotypic and genotypic data 
   
### Here attention to which dataset we take either filtered_selected_biol_rep only filtered for heterozygosity or summary_qtl_insertion which is filtered for high het and several intro in targeted region  


### Get the haplotype blues values

blues_dta_em <- read.csv(paste0(phenotype_path, "haplotype_blues_20240723.csv"), header = T) %>%
  dplyr::mutate(qtl = as.character(qtl))

### Get the haplotype PRN-FTP values

# phenotypic_data <- read.csv(paste0(phenotype_path, "lines_slope_dta_em_LD_SD_lmregression_parents_20240723.csv"), header = T, sep = ",") %>% # haplotypes_effect_gdd_dta_em_env_20240508.csv
  # dplyr::mutate(qtl = as.character(qtl)) #%>%
  #dplyr::rename(indi.cof = indi.cof_corrected, intercept = intercept_corrected) 

mixed_model_PRNTFTP <- paste0(root_alternative, "NILAS/main_folder/drouault/data/cleaned_data/from_mixed_model/")


NILs_slope <- read.csv(paste0(mixed_model_PRNTFTP, "NILs_slope_20240928.csv"), header = T) %>%
  dplyr::rename(indi.cof = Estimate, entryid2 = NILs, Se_slope = Std..Error, tvalue_slope = t.value, pval_slope = Pr...t.., BH_slope = BH)

NILs_int <- read.csv(paste0(mixed_model_PRNTFTP, "NILs_int_20240928.csv"), header = T) %>%
  dplyr::rename(intercept = Estimate, entryid2 = NILs, Se_int = Std..Error, tvalue_int = t.value, pval_int = Pr...t.., BH_int = BH)

phenotypic_data <- NILs_slope %>%
  dplyr::full_join(NILs_int) %>%
  dplyr::mutate(rp = as.numeric(gsub(".*g[0-9]|i.*", "", entryid2)),
                dp = as.numeric(gsub(".*g|[0-9]i.*", "", entryid2)), 
                qtl = as.numeric(gsub(".*q|g.*", "", entryid2))) 


phenotype_genotype <- filtered_selected_biol_rep %>% #summary_qtl_insertion (file with also the dual, triple, quadruple intro removed)
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(phenotypic_data, by = c("entryid2")) %>%
  dplyr::filter(entryid2 %in% unique(filtered_selected_biol_rep$entryid)) %>% #summary_qtl_insertion
  dplyr::select(-qtl.x ) %>%
  dplyr::rename(qtl = qtl.y)

# write.csv(phenotype_genotype, paste0(phenotype_path, date_today, "PRNFT_genotype_filtered_NILs.csv"), row.names = F)

# sensitivity_genotype <- summary_qtl_insertion %>%
#   dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
#   dplyr::full_join(sensitivity_phenotype_data, by = c("entryid2") ) %>%
#   dplyr::filter(entryid2 %in% unique(summary_qtl_insertion$entryid))

# sensitivity_phenotype_data <- read.csv(paste0(phenotype_path, "haplotype_slope_BVP_dta_em_20240508.csv"), header = T, sep = ",") %>%
#   dplyr::mutate(qtl = as.character(qtl))

### make supplmentary file 

file_supp <- phenotype_genotype %>%
  dplyr::select(Sample, Scheme, ZmPR, qtl, entryid,  biol_rep, tech_rep,  indi.cof, intercept) %>% #Mean_SD, Mean_LD, SD_daylength, LD_daylength, env_SD, env_LD
  dplyr::rename(#mean_gdd_dta_em_SD = Mean_SD, mean_gdd_dta_em_LD = Mean_LD,  
                #mean_daylength_SD = SD_daylength, mean_daylength_LD = LD_daylength, 
                slope = indi.cof) #nb_env_SD = env_SD, nb_env_LD = env_LD)

# write.csv(file_supp, paste0(file_path, date_today, "Summary_LD_SD_ILs.csv"), row.names = F)

blues_dta_em_filtered <- filtered_selected_biol_rep  %>% #summary_qtl_insertion
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(blues_dta_em) %>%
  dplyr::filter(entryid2 %in% unique(filtered_selected_biol_rep$entryid)) %>% #summary_qtl_insertion
  dplyr::select(loc, year, env, LAT, LON, Sample, Scheme, entryid, ZmPR, dp, rp, dp_name,  rp_name, 
                biol_rep, tech_rep, date_planted, date_emergence, GDD_em, bvp_ti, 
                day_dta_em, Stnd_Error, blue_date_dta,  GDD_dta, GDD_dta_em, date_ti_em, doy_ti, daylength_TI) %>%
  dplyr::rename(lat = LAT, lon = LON, 
                blue_gdd_emergence = GDD_em,
                blue_day_dta_em = day_dta_em, se_day_dta_em = Stnd_Error, blue_gdd_dta = GDD_dta, blue_gdd_dta_em = GDD_dta_em, 
                date_bvpti = date_ti_em, doy_bvpti=  doy_ti, daylength_bvpti = daylength_TI)

# write.csv(blues_dta_em_filtered, paste0(file_path, date_today, "Summary_BLUEs_ILs.csv"), row.names = F)

#### ========= Plot association number of detected genes for PS and FL BVP and slope

# cols <- phenotype_genotype  %>%  #sensitivity_genotype
#   dplyr::select(DP_BG:RP_FG, Nb_ZmPR) 
# 
# list_cols <- dput(colnames(cols))
# 
# pval_slope <- list()
# pval_bvp <- list()
# data1 <- NULL
# 
# for(i in 1:length(list_cols)){
#   print(i)
#     i <- 1
#    a <- list_cols[[i]]
#   
#   anova <- phenotype_genotype  %>%  #sensitivity_genotype
#     dplyr::select(a, indi.cof, intercept)
#   
#   model_slope <- lm(indi.cof ~ anova[,a], data = anova )
#   pval_slope[[i]] <- round((summary(model_slope))$coefficients[2,4], 3)
#   model_BVP <- lm(intercept ~ anova[,a], data = anova )
#   pval_bvp[[i]] <- round((summary(model_BVP))$coefficients[2,4], 3)  
#   data <- cbind.data.frame(a, pval_bvp[[i]], pval_slope[[i]])
#   data1 <- rbind.data.frame(data, data1)
#   
# }

anova_slope <- phenotype_genotype  %>%  #sensitivity_genotype
  dplyr::rowwise() %>%
  dplyr::mutate(sum_gene_pr  = sum(Nb_pr_qtl, Nb_pr_gwa, na.rm = T),
                sum_gene_fl  = sum(Nb_ft_gwa, Nb_ft_cg, na.rm = T)) %>%
  dplyr::select(sum_gene_pr, sum_gene_fl, Nb_ft_pr_gwa, indi.cof, intercept, qtl, rp) %>%
  dplyr::mutate(sum_genes_indi.cof = sum_gene_pr + Nb_ft_pr_gwa, 
                sum_genes_BVP = sum_gene_fl + Nb_ft_pr_gwa) 
  
  


### Test significante association with the slope 
model1 <- lm(indi.cof ~ sum_genes_indi.cof, data = anova_slope)
summary(model1)
model2 <- lm(indi.cof ~ sum_gene_fl, data = anova_slope)
summary(model2)

### Test significante association with the BVP
model3 <- lm(intercept ~ sum_genes_BVP, data = anova_slope)
summary(model3)
model4 <- lm(intercept ~ sum_gene_pr, data = anova_slope)
summary(model4)


corr_pr_gene_PS <- anova_slope  %>% #sensitivity_genotype
  dplyr::mutate(qtl = as.factor(qtl),
                rp = as.factor(rp)) %>%
  dplyr::filter(!is.na(rp)) %>%
  dplyr::mutate(sum_pr_prft_genes =  Nb_ft_pr_gwa + sum_gene_pr ) %>%
  ggplot(aes(x= sum_pr_prft_genes, y = indi.cof)) + #slope
  geom_point() + #aes(color = rp) + 
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
   stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2")), size = 5) +
  ylab("Photoperiod sensitivity") +
  xlab("Number of detected alleles associated with \n photoperiod sensitivity")
  # geom_text(aes(x= 200, y = 40, label =round(pearson_coeff_slope_pr$estimate, 4)), size = 6, check_overlap = TRUE)
corr_pr_gene_PS

### FL genes and slope 


corr_fl_gene_PS <- anova_slope  %>% #sensitivity_genotype
  dplyr::mutate(qtl = as.factor(qtl),
                rp = as.factor(rp)) %>%
  dplyr::filter(!is.na(rp)) %>%
  ggplot(aes(x= sum_gene_fl, y = indi.cof)) +
  geom_point() + #aes(color = rp)
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
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2")), size = 5) +
  ylab("Photoperiod sensitivity") +
  xlab(quote(bold(atop("Number of detected alleles associated with", bold("flowering time" ~ bolditalic("per se"))))))
# geom_text(aes(x= 200, y = 40, label = round(pearson_coeff_slope_fl$estimate, 4)), size = 6, check_overlap = TRUE)

ggarrange(corr_pr_gene_PS, corr_fl_gene_PS , 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)  
  

ggsave(paste0(fig_genes_introgression, "effect_genes_introgressed_on_slope","_", date_today,".pdf"), ggarrange(corr_pr_gene_PS, corr_fl_gene_PS ,
                                                                                                        # labels = c("A", "B"),
                                                                                                        ncol = 2, nrow = 1), width = 14, height = 10)
### PS genes and BVP

anova_intercept <- phenotype_genotype %>%
  dplyr::select(sum_gene_pr, sum_gene_fl, intercept) 
model2 <- lm(intercept ~ sum_gene_fl + sum_gene_pr, data = anova_intercept)
summary(model2)

# pearson_coeff_BVP_pr <-cor.test(pearson_test$sum_gene_pr, pearson_test$intercept, method=c("pearson"))

corr_pr_gene_BVP <- anova_slope  %>% #sensitivity_genotype
  dplyr::mutate(qtl = as.factor(qtl),
                rp = as.factor(rp)) %>%
  dplyr::filter(!is.na(rp)) %>%
  ggplot(aes(x= sum_gene_pr, y = intercept)) +
  geom_point() + #aes(color = rp)
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
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2")), size = 5) +
  ylab(quote(bold("flowering time" ~ bolditalic("per se")))) +
  xlab("Number of detected alleles associated with \n photoperiod sensitivity")
  # geom_text(aes(x= 200, y = 40, label = round(pearson_coeff_BVP_pr$estimate, 4)), size = 6, check_overlap = TRUE)

corr_pr_gene_BVP

### FL genes and BVP  

# pearson_coeff_BVP_fl <-cor.test(pearson_test$sum_gene_fl, pearson_test$intercept, method=c("pearson"))

corr_fl_gene_BVP <- anova_slope %>% #sensitivity_genotype
  dplyr::mutate(qtl = as.factor(qtl),
                rp = as.factor(rp)) %>%
  dplyr::filter(!is.na(rp)) %>%
  dplyr::mutate(sum_pr_prft_genes =  Nb_ft_pr_gwa + sum_gene_fl ) %>%
  ggplot(aes(x= sum_pr_prft_genes, y = intercept)) +
  geom_point() + #aes(color = rp)
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
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2")), size = 5) +
  ylab(quote(bold("flowering time" ~ bolditalic("per se")))) +
  xlab(quote(bold(atop("Number of detected alleles associated with", bold("flowering time" ~ bolditalic("per se"))))))
  # geom_text(aes(x= 200, y = 40, label = round(pearson_coeff_BVP_fl$estimate, 4)), size = 6, check_overlap = TRUE)

corr_fl_gene_BVP

ggarrange(corr_pr_gene_BVP, corr_fl_gene_BVP, 
          labels = c("C", "D"),
          ncol = 2, nrow = 1)

ggsave(paste0(fig_genes_introgression, "effect_genes_introgressed_on_BVP","_", date_today,".pdf"), ggarrange(corr_pr_gene_BVP, corr_fl_gene_BVP,
ncol = 2, nrow = 1), width = 14, height = 10)

### all the plots 

ggarrange(corr_pr_gene_PS, corr_fl_gene_PS , 
          corr_pr_gene_BVP, corr_fl_gene_BVP, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

ggsave(paste0(fig_genes_introgression, "effect_genes_introgressed_BVP_SLOPE","_", date_today,".pdf"), ggarrange(corr_pr_gene_PS, corr_fl_gene_PS , 
                                                                                                                corr_pr_gene_BVP, corr_fl_gene_BVP, labels = c("A", "B", "C", "D"),
                                                                                                                ncol = 2, nrow = 2)
       , width = 22, height = 18)
                                                                                                             

######### ========== Nb of lines with several introgression in ZmPR zone 

 selected_biol_rep <- selected_biol_rep %>% ### work with the non filtered dataset 
  dplyr::rowwise() %>%
  dplyr::mutate(sum_gene_pr  = sum(Nb_pr_qtl, Nb_pr_gwa, na.rm = T),
                sum_gene_fl  = sum(Nb_ft_gwa, Nb_ft_cg, na.rm = T), 
                sum_introgression = sum(across(Nb_ZmPR1:Nb_ZmPR4), na.rm = T))
 
 anova_sum_introgression <- phenotype_genotype %>%
   # dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
   # dplyr::full_join(blues_dta_em)  %>%
   dplyr::mutate(across(Nb_ZmPR1:Nb_ZmPR4, ~ ifelse(. != 0, 1, 0))) %>%
 dplyr::rowwise() %>%
   dplyr::mutate(sum_introgression = sum(across(Nb_ZmPR1:Nb_ZmPR4), na.rm = T)) %>%
   dplyr::select(sum_introgression, intercept, indi.cof) 
 model_anova_intercept_sum_introgression <- lm(intercept ~ sum_introgression, data = anova_sum_introgression)
 summary(model_anova_intercept_sum_introgression)
 
 model_anova_slope_sum_introgression <- lm(indi.cof ~ sum_introgression, data = anova_sum_introgression)
 summary(model_anova_slope_sum_introgression)
 
ILs_several_intro_ZmPR <- selected_biol_rep %>%
   dplyr::filter(sum_introgression >1)  

nrow(ILs_several_intro_ZmPR)

t_test_results <- selected_biol_rep %>% 
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample),
                qtl = gsub("q", "", qtl)) %>%
  dplyr::full_join(phenotype_genotype) %>% #sensitivity_phenotype_data
  dplyr::filter(!is.na(qtl)) %>%
  dplyr::mutate(sum_introgression = as.factor(sum_introgression), 
                sum_introgression = factor(sum_introgression, levels =  c("0", "1", "2", "3", "4" ,"5", "6", "7", "8", "9", "10", "11")))%>%
  dplyr::filter(sum_introgression != "6" & sum_introgression != "7"&  sum_introgression != "8" &  sum_introgression != "11" ) %>%
  dplyr::group_by(sum_introgression) %>% 
  do(                               
    t.test(.$indi.cof, mu = 0) %>%  # compare mean to of rp 1  
      broom::tidy()                        # tidy into tibble
  )

t_test_results$Spearman_p_sign <- "ns"
t_test_results$Spearman_p_sign[t_test_results$p.value < 0.05] <- "*"
t_test_results$Spearman_p_sign[t_test_results$p.value < 0.01] <- "**"
t_test_results$Spearman_p_sign[t_test_results$p.value < 0.001] <- "***"
t_test_results$Spearman_p_sign[t_test_results$p.value < 0.0001] <- "****"


introgression_nb_slope <- selected_biol_rep %>%
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(sensitivity_phenotype_data, by = c("entryid2") ) %>%
  dplyr::filter(!is.na(qtl.y)) %>%
  dplyr::mutate(sum_introgression = as.factor(sum_introgression), 
                sum_introgression = factor(sum_introgression, levels = c("0", "1", "2", "3", "4" ,"5", "6", "7", "8", "9"))) %>%
  dplyr::full_join(t_test_results) %>%
  ggplot(aes(x= sum_introgression, y = slope)) +
  geom_boxplot() +
  geom_text(data= t_test_results,  aes(y = 40, label = Spearman_p_sign), color = "red", position = position_dodge(width = .8), size = 6) + 
  geom_jitter(aes(color = qtl.y)) + 
  theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=14),
        legend.title  = element_text(size=14),
        legend.position = "right"
  ) + ylab("slope")
introgression_nb_slope


# ggsave(paste0(fig_genes_introgression, "effect_nbZmPR_intro_on_slope","_", date_today,".pdf"), introgression_nb_slope, width = 14, height = 10)



#create a data set to put verticale lines in fig 
rp <- c("1", "1", "1","1", "2", "2", "2", "2") #"2369", "2369", "2369", "2369","LH123hT", "LH123hT", "LH123hT", "LH123hT"
qtl <- c("1", "2", "3", "4", "1", "2", "3", "4")
Int <- c(0, 0, 0, 0, 0, 0, 0, 0)
category <-c("black", "black", "black", "black", "black", "black", "black", "black")
data <- cbind.data.frame(rp, qtl, Int, category)

dataInt <- phenotype_genotype %>%
  dplyr::mutate(rp = as.character(rp), qtl = as.character(qtl)) %>% #ifelse(rp == 1, "2369", "LH123hT")
  dplyr::group_by(rp, qtl) %>%
  dplyr::summarize(Int = mean(indi.cof, na.rm = T)) %>%
  dplyr::mutate(category = "#F17CB0") %>%
  dplyr::full_join(data) %>%
  dplyr::rename(gb = rp) %>%
  as.data.frame()

sensitivity_genotype <- phenotype_genotype%>%
  dplyr::mutate(p_val_slope_adjusted = stats::p.adjust(p_val_slope, method = "BH")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(binary = ifelse(p_val_slope_adjusted < 0.2, "S", "N"), 
                binary = as.factor(binary)) 

### make the plot of the slope per rp/qtl
distribution_plot <-  sensitivity_genotype%>%
  dplyr::mutate(qtl = as.character(qtl),
                rp = as.character(rp), 
                binary  = as.character(binary))%>%
  ggplot(aes(x = indi.cof)) + #fill = binary
  geom_histogram(aes(x = indi.cof,  y = ..density.., fill = binary), # color = ifelse(binary == "s")
                 colour = 1) +
  scale_fill_manual(values = c( "#3399FF", "grey"))+
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25) +
  facet_grid(qtl~rp, scales= "fixed") +
  xlab("slope of sensitivity") +
  # geom_vline(data = dataInt,  mapping = aes(xintercept= Int), linetype="dashed", colour = dataInt$category, linewidth=2) +
  # geom_text(aes(x = 25, y = 0.08, label=paste("coef var", coef_var, sep =" ")), color = "black", size = 8, check_overlap = T) +
  # geom_text(aes(x = 25, y = 0.06, label=paste("IQR", IQR, sep =" ")), color = "black", size = 8, check_overlap = T) +
  
  # scale_x_continuous(breaks=c(-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5,0,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5), 
  # expand = c(0.05, 0)) + 
  labs(title = '') 
distribution_plot <- distribution_plot+ theme_classic(base_size = 25) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) 

ggsave(paste0(plot_path, "distribution_slope_PS_dta_em","_", date_today,".pdf"), distribution_plot, width = 18, height = 15)




##### end 

###### ======= Check consistancy of the data using technical replicates 

# 
# fin <- NULL
# data_consistancy <- do.call(rbind, lapply(split(complementary_lines,as.factor(complementary_lines$ID)), function(x) {
#   # x <- split(complementary_lines,as.factor(complementary_lines$ID))[[2]]
#   # x <- complementary_lines %>% filter (ID == "CML10 x LH123Ht_i02_q1")
#   print(unique(x$ID))
#   print(nrow(x))
#   abs_DP_FG <- abs(max(x$DP_FG) - min(x$DP_FG))
#   abs_DP_BG <- abs(max(x$DP_BG) - min(x$DP_BG))
#   abs_RP_FG <- abs(max(x$RP_FG) - min(x$RP_FG))
#   abs_RP_BG <- abs(max(x$RP_BG) - min(x$RP_BG))
#   abs_Het_BG <- abs(max(x$Het_BG) - min(x$Het_BG))
#   abs_Het_FG  <- abs(max(x$Het_FG ) - min(x$Het_FG ))
#   data <- cbind.data.frame(abs_DP_FG, abs_DP_BG, abs_RP_FG, abs_RP_BG, abs_Het_BG, abs_Het_FG)
#   data <- cbind(unique(x$ID),data)
#   colnames(data)[1] <- "ID" 
#   fin <- rbind.data.frame(fin,  data)
#   return(fin)
# }))
# 
# check_consistency_tech <- data_consistancy %>%
#   tidyr::pivot_longer(cols = abs_DP_BG:abs_Het_FG,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   tidyr::separate(ID, c("Scheme", "IL_nb", "qtl"), sep = "_", remove = F) %>%
#   ggplot() +
#   geom_histogram(aes(x= count)) +
#   facet_grid(IL_nb~categories) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20), 
#         strip.text.y = element_text(size = 20), 
#         legend.position = "none"
#   ) 
# 
# ggsave(paste0(fig_consistency, "consistancy_techrep","_", date_today,".pdf"), check_consistency_tech, width = 18, height = 14)
# 
# 
# ##### Summary content tech rep 
# 
# summary_tech <-  data_consistancy %>%
#   dplyr::mutate(count = n_distinct(ID)) %>% ### get total number of line 
#   dplyr::mutate(across(abs_DP_FG:abs_Het_FG, ~ sum(. != 0), .names = "{col}_count")) %>%
#   dplyr::mutate(across(abs_DP_FG:abs_Het_FG, ~ sum(. > 3), .names = "{col}_supp3")) %>%
#   dplyr::select(-c(ID:abs_Het_FG )) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(across(abs_DP_FG_count:`abs_Het_FG_supp3`, ~ sum((. /count)*100), .names = "{col}_percent")) 
# 
# list <- c("count", "abs_DP_FG_count", "abs_DP_BG_count", "abs_RP_FG_count", 
#           "abs_RP_BG_count", "abs_Het_BG_count", "abs_Het_FG_count", "abs_DP_FG_supp3", 
#           "abs_DP_BG_supp3", "abs_RP_FG_supp3", "abs_RP_BG_supp3", "abs_Het_BG_supp3", 
#           "abs_Het_FG_supp3") 
# 
# list02 <- c("abs_DP_FG_count_percent", "abs_DP_BG_count_percent", 
#             "abs_RP_FG_count_percent", "abs_RP_BG_count_percent", "abs_Het_BG_count_percent", 
#             "abs_Het_FG_count_percent", "abs_DP_FG_supp3_percent", "abs_DP_BG_supp3_percent", 
#             "abs_RP_FG_supp3_percent", "abs_RP_BG_supp3_percent", "abs_Het_BG_supp3_percent", 
#             "abs_Het_FG_supp3_percent")
# 
# #### Plot % of lines divergent for the selected variables 
# 
# check_percentage_unconsistency_tech <- summary_tech %>%
#   dplyr::select(abs_DP_FG_count_percent:abs_Het_FG_supp3_percent) %>%
#   tidyr::pivot_longer(cols = abs_DP_FG_count_percent:abs_Het_FG_supp3_percent,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(categories = as.factor(categories), 
#                 categories = factor(categories, levels = list02)) %>%
#   ggplot(aes(x= categories, y = count)) +
#   geom_bar(stat="identity", aes(fill = categories)) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         legend.position = "none"
#   ) +
#   ylab("percentageof technical replicates that differs")
# check_percentage_unconsistency_tech
# 
# ggsave(paste0(fig_consistency, "percentage_unconsistancy_techrep","_", date_today,".pdf"), check_percentage_unconsistency_tech, width = 10, height = 12)
# 
# 
# #### Plot number of lines divergent for the selected variables 
# 
# check_unconsistancy_tech <- summary_tech %>%
#   dplyr::select(count:abs_Het_FG_supp3) %>%
#   tidyr::pivot_longer(cols = count:abs_Het_FG_supp3,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(categories = as.factor(categories), 
#                 categories = factor(categories, levels = list)) %>%
#   ggplot(aes(x= categories, y = count)) +
#   geom_bar(stat="identity", aes(fill = categories)) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20), 
#         strip.text.y = element_text(size = 20), 
#         legend.position = "none"
#   ) 
# check_unconsistancy_tech
# 
# ggsave(paste0(fig_consistency, "values_unconsistancy_techrep","_", date_today,".pdf"), check_unconsistancy_tech, width = 10, height = 12)
# 
# 

# max_table <- NULL
# for(i in 1:length(matrix_tech_rep)){
#   print(i)
#   
#   consistency_table <- NULL
#    i <- 1
#   x <- subset(matrix_tech_rep, entryid == "NILASq1g12i02" ) # list_entryid[[i]])
#   
#   summary_x <- x %>%
#     dplyr::filter(test2 != "no_rep") %>%
#     dplyr::summarise(distinct_biol_rep = n_distinct(biol_rep)) %>%
#     dplyr::ungroup()
#   
#   summary_tot_biorep <- x %>%
#     dplyr::group_by(biol_rep) %>%
#     dplyr::summarise(distinct_biol_rep = n_distinct(biol_rep))
#   
#   #### several biological samples with tech rep but several biological rep in total 
#   
#   
#   if(unique(summary_x$distinct_biol_rep)> 1 & nrow(summary_tot_biorep) >1){
#     bio_tech_rep_name <- NULL
#     Y <- NULL
#     
#     bio_tech_rep_name <- x %>%
#       dplyr::filter(test2 != "no_rep") 
#     
#     unique(bio_tech_rep_name$biol_rep)
#     
#     y <- subset(genotypic_data, entryid == "NILASq1g12i02" 
#                 # list_entryid[[i]] 
#                 & biol_rep %in% unique(bio_tech_rep_name$biol_rep))
#     
#     best_sample_tech_rep <-y %>%
#       dplyr::group_by(biol_rep) %>%
#       dplyr::mutate(Mean_Het_BG = mean(Het_BG, na.rm = T))%>%
#       dplyr::ungroup() %>%
#       dplyr::rowwise() %>%
#       dplyr::mutate(tech_sample = (abs(Het_BG - Mean_Het_BG))) %>%
#       dplyr::mutate(min = min(abs(Het_BG - Mean_Het_BG))) %>%
#       dplyr::summarise(best_sample = Sample[which(tech_sample == min )], 
#                        value_Het_BG = Het_BG[which(tech_sample == min )])
#     
#     value_second_biorep <- subset(genotypic_data, entryid == list_entryid[[i]] & biol_rep != unique(bio_tech_rep_name$biol_rep))%>%
#       dplyr::select(Sample, Het_BG) %>%
#       dplyr::mutate(Source = "biol_rep") 
#     
# 
#     consistency_table <- best_sample_tech_rep%>%
#       dplyr::rename(Sample = best_sample, Het_BG = value_Het_BG) %>%
#       dplyr::mutate(Source = "tech_rep") %>%
#       dplyr::full_join(value_second_biorep) %>%
#       dplyr::mutate(consitency= abs(Het_BG - Het_BG[Source == 'biol_rep']), 
#     consitency = mean(consitency, na.rm =T))
#     
#     }
#   
#   #### 1 single biological sample with tech rep but several biological rep in total 
#   
#   if(summary_x$distinct_biol_rep == 1 & nrow(summary_tot_biorep) >1){
#     bio_tech_rep_name <- NULL
#     Y <- NULL
# 
#     bio_tech_rep_name <- x %>%
#       dplyr::filter(test2 != "no_rep") 
#     
#     unique(bio_tech_rep_name$biol_rep)
#     
#     y <- subset(genotypic_data, entryid == list_entryid[[i]] & biol_rep == unique(bio_tech_rep_name$biol_rep))
#     
#      best_sample_tech_rep <-y %>%
#        dplyr::mutate(Mean_Het_BG = mean(Het_BG, na.rm = T))%>%
#        dplyr::rowwise() %>%
#        dplyr::mutate(tech_sample = (abs(Het_BG - Mean_Het_BG))) %>%
#        dplyr::ungroup() %>%
#        dplyr::mutate(min = min(abs(Het_BG - Mean_Het_BG))) %>%
#        dplyr::summarise(best_sample = Sample[which(tech_sample == min )], 
#                         value_Het_BG = Het_BG[which(tech_sample == min )])
#      
#      value_second_biorep <- subset(genotypic_data, entryid == list_entryid[[i]] & biol_rep != unique(bio_tech_rep_name$biol_rep))%>%
#        dplyr::select(Sample, Het_BG) %>%
#        dplyr::mutate(Source = "biol_rep") 
#      
#      consitency <- abs(best_sample_tech_rep$value_Het_BG - value_second_biorep$Het_BG)
#      
# 
#      consistency_table <- best_sample_tech_rep%>%
#        dplyr::rename(Sample = best_sample, Het_BG = value_Het_BG) %>%
#        dplyr::mutate(Source = "tech_rep") %>%
#        dplyr::full_join(value_second_biorep) %>%
#        dplyr::mutate(consitency= consitency, 
#                      consitency = mean(consitency, na.rm =T))
#      
#   }
#  
#   #### no tech rep but several biological rep in total 
#   
#   if(nrow(summary_x) == 0 & nrow(summary_tot_biorep) > 1){
# 
#     consistency_table <- subset(genotypic_data, entryid == list_entryid[[i]])%>%
#       dplyr::select(Sample, Het_BG) %>%
#       dplyr::mutate(Source = "biol_rep") %>%
#       dplyr::mutate(consitency= abs(Het_BG - lag(Het_BG)), 
#                     consitency = mean(consitency, na.rm = T))
#   }
#   if(nrow(summary_x) == 0 & nrow(summary_tot_biorep) == 1){
#     consistency_table <- subset(genotypic_data, entryid == list_entryid[[i]])%>%
#       dplyr::select(Sample, Het_BG) %>%
#       dplyr::mutate(Source = "no_rep") 
#   }
#   
#   #### no tech rep and no biological rep in total 
#   
#   
#   max_table <- rbind.data.frame(consistency_table, max_table)
# 
# }
########

# fin <- NULL
# data_consistancy_biologicalrep <- do.call(rbind, lapply(split(genotypic_data,as.factor(genotypic_data$ID)), function(x) {
#   # x <- split(genotypic_data,as.factor(genotypic_data$ID))[[2]]
#   # x <- genotypic_data %>% filter (ID == "CML10 x LH123Ht_i02_q1")
#   print(unique(x$ID))
#   print(nrow(x))
#   abs_DP_FG <- abs(max(x$DP_FG) - min(x$DP_FG))
#   abs_DP_BG <- abs(max(x$DP_BG) - min(x$DP_BG))
#   abs_RP_FG <- abs(max(x$RP_FG) - min(x$RP_FG))
#   abs_RP_BG <- abs(max(x$RP_BG) - min(x$RP_BG))
#   abs_Het_BG <- abs(max(x$Het_BG) - min(x$Het_BG))
#   abs_Het_FG  <- abs(max(x$Het_FG ) - min(x$Het_FG ))
#   data <- cbind.data.frame(abs_DP_FG, abs_DP_BG, abs_RP_FG, abs_RP_BG, abs_Het_BG, abs_Het_FG)
#   data <- cbind(unique(x$ID),data)
#   colnames(data)[1] <- "ID" 
#   fin <- rbind.data.frame(fin,  data)
#   return(fin)
# }))
# 
# check_consistancy_biol <- data_consistancy_biologicalrep %>%
#   tidyr::pivot_longer(cols = abs_DP_BG:abs_Het_FG,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   tidyr::separate(ID, c("Scheme", "IL_nb", "qtl"), sep = "_", remove = F) %>%
#   ggplot() +
#   geom_histogram(aes(x= count)) +
#   facet_grid(IL_nb~categories) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         strip.text.x = element_text(size = 20), 
#         strip.text.y = element_text(size = 20), 
#         legend.position = "none"
#   ) 
# 
# check_consistancy_biol
# 
# ggsave(paste0(fig_consistency, "consistancy_biorep","_", date_today,".pdf"), check_consistancy_biol, width = 18, height = 20)
# 
# 
# 
# #### Plot number of lines divergent for the selected variables 
# 
# summary <-  data_consistancy_biologicalrep %>%
#   dplyr::mutate(count = n_distinct(ID)) %>% ### get total number of line 
#   dplyr::mutate(across(abs_DP_FG:abs_Het_FG, ~ sum(. != 0), .names = "{col}_count")) %>%
#   dplyr::mutate(across(abs_DP_FG:abs_Het_FG, ~ sum(. > 3), .names = "{col}_supp3")) %>%
#   dplyr::select(-c(ID:abs_Het_FG )) %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(across(abs_DP_FG_count:`abs_Het_FG_supp3`, ~ sum((. /count)*100), .names = "{col}_percent")) 
# 
# #### Plot % of lines divergent for the selected variables 
# 
# check_percentage_unconsistancy_biol <- summary %>%
#   dplyr::select(abs_DP_FG_count_percent:abs_Het_FG_supp3_percent) %>%
#   tidyr::pivot_longer(cols = abs_DP_FG_count_percent:abs_Het_FG_supp3_percent,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(categories = as.factor(categories), 
#                 categories = factor(categories, levels = list02)) %>%
#   ggplot(aes(x= categories, y = count)) +
#   geom_bar(stat="identity", aes(fill = categories)) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         legend.position = "none"
#   ) +
#   ylab("percentageof biological replicates that differs")
# check_percentage_unconsistancy_biol
# 
# ggsave(paste0(fig_consistency, "percentage_unconsistancy_biorep","_", date_today,".pdf"), check_percentage_unconsistancy_biol, width = 10, height = 12)
# 
# 
# #### Plot number of lines divergent for the selected variables 
# 
# check_unconsistancy_biol <- summary %>%
#   dplyr::select(count:abs_Het_FG_supp3) %>%
#   tidyr::pivot_longer(cols = count:abs_Het_FG_supp3,
#                       names_to = c("categories"),
#                       values_to = "count") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(categories = as.factor(categories), 
#                 categories = factor(categories, levels = list)) %>%
#   ggplot(aes(x= categories, y = count)) +
#   geom_bar(stat="identity", aes(fill = categories)) +
#   theme(axis.text.x = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.x = element_text(size=20),
#         axis.text.y = element_text(size=20, angle=90, hjust=1, vjust=1),
#         axis.title.y = element_text(size=20),
#         legend.position = "none"
#   ) 
# check_unconsistancy_biol
# 
# ggsave(paste0(fig_consistency, "values_unconsistancy_biorep","_", date_today,".pdf"), check_unconsistancy_biol, width = 10, height = 12)
# 
# 
# ##### Summary table overleaf
# 
# Table_overleaf <- summary %>%
#   dplyr::select(abs_DP_FG_count:abs_Het_FG_supp3_percent) %>%
#   tidyr::pivot_longer(cols = abs_DP_FG_count:abs_Het_FG_supp3_percent,
#                       names_to = c("variable", ".value"),
#                       names_pattern =  "(.*G)_(.*)$") %>%
#   dplyr::rename(total = count, sup_3per =  supp3, total_percent = count_percent, sup_3per_percent = supp3_percent) %>%
#   dplyr::mutate(variable = dplyr::recode(variable, abs_DP_FG = "Donor parent in foreground (DP_FG)", abs_DP_BG = "Donor parent in background (DP_BG)" ,
#                                          abs_RP_FG = "Recurrent parent in foreground (RP_FG)", abs_RP_BG = "Recurrent parent in background (RP_BG)", 
#                                          abs_Het_BG = "Heterosigozity background (Het_BG)", abs_Het_FG = "Heterosigozity foreground (Het_FG)"))
# 
# 
# 
# write.csv(Table_overleaf, paste0(save_path, date_today, "Summary_consistency_biorep.csv"), row.names = F)
# 
# #### List inconsistent biological rep 
# 
# inconsistent_list <- data_consistancy_biologicalrep %>%
#   dplyr::mutate(across(abs_DP_FG:abs_Het_FG, ~ ifelse(. != 0, "inconsistent", "consistent"), .names = "{col}_check")) %>%
#   dplyr::filter(abs_DP_FG_check== "inconsistent" | abs_DP_BG_check == "inconsistent"| abs_RP_FG_check == "inconsistent"| 
#                   abs_RP_BG_check == "inconsistent" | abs_Het_BG_check== "inconsistent"|  abs_Het_FG_check == "inconsistent")
# 
# length(inconsistent_list)
# 
# ###### ======= Select the best lines among biological replicates 
# 
# ## Criteria 
# # max_RP_BG
# # max_DP_FG
# # min_DP_BG
# # min_RP_FG
# # min_Het_BG
# # min_Het_FG
# 
# fin <- NULL
# best_lines_QC <- do.call(rbind, lapply(split(genotypic_data,as.factor(genotypic_data$ID)), function(x) {
#   x <- split(genotypic_data,as.factor(genotypic_data$ID))[[2]]
#   # x <- genotypic_data %>% filter (ID == "CML10 x LH123Ht_i02_q1")
#   print(unique(x$ID))
#   max_DP_FG <- x$Sample[which( x$DP_FG == max(x$DP_FG))]
#   min_DP_BG <- x$Sample[which( x$DP_BG == min(x$DP_BG))]
#   min_RP_FG <- x$Sample[which( x$RP_FG == min(x$RP_FG))]
#   max_RP_BG <- x$Sample[which( x$RP_BG == max(x$RP_BG))]
#   min_Het_BG <-  x$Sample[which( x$Het_BG == min(x$Het_BG))]
#   min_Het_FG <- x$Sample[which( x$Het_FG == min(x$Het_FG))]
#   n <- max(length(max_DP_FG), length(min_DP_BG), length(min_RP_FG), length(max_RP_BG), length(min_Het_BG), length(min_Het_FG))
#   data <- cbind.data.frame(max_DP_FG[1:n], min_DP_BG[1:n], min_RP_FG[1:n], max_RP_BG[1:n], min_Het_BG[1:n], min_Het_FG[1:n])
#   result <- sapply(unique(x$Sample), function(value) {
#     sum(data == value, na.rm =T)}) ### get the line that is repeated the most so fullfill the most criteria 
#   result <- as.data.frame(result)
#   result <- tibble::rownames_to_column(result, "Sample")
#   x <- x %>%
#     dplyr::filter(Sample %in% result$Sample[which(result$result == max(result$result))])
#   fin <- rbind.data.frame(fin, x)
#   return(fin)
# }))
# 
# 
# 
# 

### get the position of the genes in the chr 

directory <- paste("~/Dropbox-SE\ Dropbox/allshare/stage_eloise/NILAS_haploimputer/Outputs/", sep="")

ZmPR_list <- 1:4
rp_list <- c("2369", "LH123Ht")
dp_list <- c("CML10", "CML258", "CML277", "CML373", "CML341", "Tzi8", "Tzi9")

list.data <- NULL
list_name <- NULL
for (l in 1:length(ZmPR_list)) {
  for (r in 1:length(rp_list)) {
    for (d in 1:length(dp_list)) {
      
      # l <- 1
      # r <- 1
      # d <- 1
      
      ZmPR <- ZmPR_list[l] 
      RP <- rp_list[r]
      DP <- dp_list[d]
      i <- paste0(l,r,d)
      
      group <- paste(RP,"Cx",DP,"C",sep="")
      file1 <- paste0(RP,"Cx",DP,"CZmPR",ZmPR, "Summarydet2.xlsx")
      
      list.data[[i]] <- read_xlsx(paste0(directory,group, "/" , file1))
      
      name <- paste(RP,DP,"ZmPR",ZmPR, sep = "_")
      list_name <- cbind(list_name, name)
    }}} 


for (a in 1:length(list_name)) {
  assign(list_name[[a]], list.data[[1]], .GlobalEnv)
}

dataset <- NULL
for (i in 1:length(list.data)) {
  data <- as.data.frame(list.data[[i]])
  to_add <- rep(1, length.out = nrow(data))
  data$count <- to_add
  summed_data <- aggregate(count ~ Sample + CHROM, data = data,  # label
                           sum )
  unique_rows <- row.names(unique(summed_data[c("Sample","CHROM")]))
  unique.df <- summed_data[unique_rows,]  
  dataset <- rbind.data.frame(dataset, unique.df)
  return(dataset)
}

nb_gene_chr_pheno <- dataset %>%
  dplyr::mutate(entryid2 = gsub("s[0-9]{1}p[0-9]{1}r[0-9]{1}", "", Sample)) %>%
  dplyr::full_join(sensitivity_phenotype_data, by = c("entryid2") ) 


