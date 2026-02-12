##### Epistasis versus additivitu scenario 

U and B two loci, mm maize alleles, tt teosinte alleles

Umm_Bmm <- 18 a= 14.9 k = 1.46
Umm_Bmt <- 54.6
Umm_Btt <- 47.8
Umt_Bmm <- 40.9
Umt_Bmt <- 47.6
Umt_Btt <- 83.6
Utt_Bmm <- 61.1
Utt_Bmt <- 66.5
Utt_Btt <- 101.7 

list_alleles_QTL1 <- c("Umm", "Umm", "Umm", "Umt", "Umt", "Umt", "Utt", "Utt", "Utt")
list_alleles_QTL2 <- c("Bmm", "Bmt", "Btt", "Bmm", "Bmt", "Btt", "Bmm", "Bmt", "Btt")
list_alleles_combination <- c("Umm_Bmm", "Umm_Bmt", "Umm_Btt", "Umt_Bmm", "Umt_Bmt", "Umt_Btt", "Utt_Bmm", "Utt_Bmt", "Utt_Btt")
list_phenotypic_values <- c(18, 54.6, 47.8, 40.9, 47.9, 83.6, 61.1, 66.5, 101.7)
alleles_combination_ID <- rep(1:3, length.out = length(list_phenotypic_values))


data <- cbind.data.frame(list_alleles_QTL1, list_alleles_QTL2, list_alleles_combination, alleles_combination_ID, list_phenotypic_values)

data <- data %>%
  dplyr::group_by(list_alleles_QTL2) %>%
  dplyr::mutate(a = ((data$list_phenotypic_values[data$alleles_combination_ID == 1]- data$list_phenotypic_values[data$alleles_combination_ID == 3])/2))%>%
  dplyr::mutate(k = ((data$list_phenotypic_values[data$alleles_combination_ID == 2]- data$list_phenotypic_values[data$alleles_combination_ID == 3])/a)-1)%>%
  dplyr::ungroup()
