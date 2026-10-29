## load packages
library(patchwork)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(data.table)
library(stringr)   

## Date of the day
date_today <- Sys.Date()
date_today <- gsub("-", "", date_today) 

## empty list for storing plots
Q <- list()

##Input
args<-commandArgs(TRUE)

path=args[1]

dataset_dir = "dp6q30dr0.40maf0.01"
# map_file = "nilas_f4_dp6q30dr0.40maf0.01.alleleInfo"
# gtype_file = "nilas_f4_dp6q30dr0.40maf0.01.GT.FORMAT.final"
##set working directory
directory <- paste0("~/Documents/Drouault_etal_paper2/01_paperdev/ATLAS/paper/Stage_Eloise/Data/", dataset_dir, "/outputs/")
# List all directories in the current working directory




###### ======= Q[[4]], ZmPR4, dp 3 & 4 


ZmPR_list <- 1:4
rp_list <- c("2369", "LH123Ht")
dp_list <-  c("CML10", "CML258", "CML277", "CML373", "CML341", "Tzi8", "Tzi9")
L <- list()
data <- NULL

for (l in 1:length(ZmPR_list)) {
  for (r in 1:length(rp_list)) {
    for (d in 1:length(dp_list)) {
      # # 
            # l <- 4
      #       r <- 1
      #       d <- 1
      
      ZmPR <- ZmPR_list[l]
      RP <- rp_list[r]
      DP <- dp_list[d]
      
      #if ((!(ZmPR==2 & RP=="2369" & DP=="CML258")) & (!(ZmPR==3 & RP=="LH123Ht" & DP=="Tzi8"))){  
      
      #-------------------  Read/prep data
      
      # 'A':'DP'
      # 'A2':'DP2'
      # 'G':'RP'
      # '1':'RP2'
      # '0':'NA'
      # 'C':'Het'
      
      # RP2 is called during encoding and DP2 is called at imputation.
      # In general, I thought there was more DP2 than RP2 but going through the CML10xLH123Ht, I do see a good deal of RP2 alleles. 
      # You can see DP2 in Tzi8xLH123Ht- ZmPR3 
      
      #####
      group <- paste(RP,"Cx",DP,"C",sep="")
      file1 <- paste(RP,"Cx",DP,"C_ZmPR",ZmPR,"_Geno_Sum.txt",sep="")
      file2 <- paste(RP,"Cx",DP,"C_Coordinates.txt",sep="")
      file3 <- paste(RP,"Cx",DP,"C_ZmPR",ZmPR,"_Fig.txt",sep="")
      file4 <- paste(RP,"Cx",DP,"C_ZmPR",ZmPR,"_SUM.txt",sep="")
      #####
      
      #### Create names for expected set of lines
      if (RP=="2369") {
        RP_no <- 1
      } else if (RP=="LH123Ht") {
        RP_no <- 2
      } else 
        print("Oh no, RP not found!!!")
      
      if (DP=="CML10") {
        DP_no <- 1
      } else if (DP=="CML258") {
        DP_no <- 2
      } else if (DP=="CML277") {
        DP_no <- 3
      } else if (DP=="CML341") {
        DP_no <- 4
      } else if (DP=="CML373") {
        DP_no <- 5
      } else if (DP=="Tzi8") {
        DP_no <- 6
      } else if (DP=="Tzi9") {
        DP_no <- 7
      } else
        print("Oh no, DP not found!!!")
      
      line_no <- seq(1,12,1)
      Sample <- paste("NILASq",ZmPR,"g",DP_no,RP_no,"i", formatC(line_no, width=2, flag="0"), sep="")
      
      ####
      
      #### read and compile data from haploimputer
      
      summ <- read.table(
        paste(directory, group, "/", file1, sep=""), 
        sep='\t', header=T)
      summ$Complement <- grepl("C", summ$Sample, fixed=FALSE)
      
      summ <- summ %>%
        dplyr::rowwise() %>%
        dplyr::mutate("Recurrent_Parent"=(RP_BG +RP_FG), 
                      "Heterozygous"= (Het_BG +Het_FG),
                      "Missing"= (NA_BG + NA_FG)) %>%
        dplyr::select(-c(RP_BG, RP_FG, Het_BG, Het_FG, NA_BG, NA_FG))
      
      # summ <- cbind.data.frame(summ[c("Sample","Scheme","ZmPR")],
      #                          "Recurrent_Parent"=rowSums(summ[,c("RP_BG","RP_FG")]),
      #                          summ[,c("DP_BG","DP_FG")],
      #                          "Heterozygous"=rowSums(summ[,c("Het_BG","Het_FG")]),
      #                          "Missing"=rowSums(summ[,c("NA_BG","NA_FG")]),
      #                          "Complement"=summ$Complement)
      #filter complements
      summ <- summ[summ$Complement != "TRUE", ]
      summ <- summ[,-which(names(summ)=="Complement")]
      #filter extra reps
      #summ$plant <- as.numeric(sub(".*p(.+)r.", "\\1", summ$Sample))
      #summ$rep <- as.numeric(sub(".*r", "", summ$Sample))
      summ$name <- sub("s.*", "", summ$Sample)
      
      summ %>% group_by(name) %>% filter(DP_FG == max(DP_FG)) %>% filter(Missing == min(Missing)) %>% filter(Heterozygous == min(Heterozygous)) -> summ
      sample_subset <- summ$Sample
      
      summ$identical <- (c(FALSE, summ$name) == c(summ$name, NA))[1:length(summ$name)]
      summ <- summ[summ$identical != "TRUE", ]
      summ <- summ[,-which(names(summ)=="name")]
      summ <- summ[,-which(names(summ)=="identical")]
      
      summ$Sample <- sub("s.*", "", summ$Sample)
      ## add missing samples
      summ <- merge(as.data.frame(Sample), summ, by="Sample", all.x=T)
      summ$Sample <- as.factor(summ$Sample)
      summ.melt <- summ %>%
        dplyr::select(-X) %>%
        tidyr::pivot_longer(c(-"Sample", -"Scheme", -"ZmPR"), 
                            names_to = "variable", 
                            values_to = "values") %>%
        dplyr::mutate(variable = as.factor(variable))
      #melt(summ, id.vars=c("Sample", "Scheme", "ZmPR"))
      levels(summ.melt$variable)[levels(summ.melt$variable)=="Recurrent_Parent"] <- "RP"
      levels(summ.melt$variable)[levels(summ.melt$variable)=="DP_FG"] <- "DP-T"
      levels(summ.melt$variable)[levels(summ.melt$variable)=="DP_BG"] <- "DP-B"
      levels(summ.melt$variable)[levels(summ.melt$variable)=="Heterozygous"] <- "Het"
      levels(summ.melt$variable)[levels(summ.melt$variable)=="Missing"] <- "NA"
      #re-order levels of factor variable (for plotting)
      summ.melt$variable <- factor(summ.melt$variable, levels=c("RP", "DP-B", "DP-T", "Het", "NA"))
      #re-order levels of factor sample (for plotting)
      summ.melt$Sample <- factor(summ.melt$Sample, levels=rev(levels(summ.melt$Sample)))
      
      
      #read map
      map <- read.table(
        paste(directory, group, "/", file2, sep=""), 
        sep='\t', header=F)
      
      #read genotypes
      ti <- read.table(
        paste(directory, group, "/", file3, sep=""), 
        sep='\t', header=T)
      #rename column
      colnames(ti)[1] <- "Sample"
      #filter extra reps
      # ti$plant <- as.numeric(sub(".*p(.+)r.", "\\1", ti$Sample))
      # ti$rep <- as.numeric(sub(".*r", "", ti$Sample))
      ti <- ti[which(ti$Sample %in% sample_subset),]
      #ti$name <- sub("s.*", "", ti$Sample)
      #ti$identical <- (c(FALSE, ti$name) == c(ti$name, NA))[1:length(ti$name)]
      #ti <- ti[ti$identical != "TRUE", ]
      #ti <- ti[,-which(names(ti)=="name")]
      #ti <- ti[,-which(names(ti)=="identical")]
      
      ti$Sample <- sub("s.*", "", ti$Sample)
      ## add missing samples
      ti <- merge(as.data.frame(Sample), ti, by="Sample", all.x=T)
      ti.IDs <- ti$Sample
      ti.IDs <- sub(".*i", "", ti.IDs)
      ti.IDs <- as.numeric(ti.IDs)
      ti$Sample <- as.factor(ti$Sample)
      #add numeric id index per line
      ti <- cbind.data.frame(id=ti.IDs, ti)
      
      #get markers from map
      ti.markers <- colnames(ti)[3:ncol(ti)] 
      ti.markers <- sub("X", "", ti.markers)
      map <- map[which(map$V1 %in% ti.markers),]
      #compute end
      map <- cbind.data.frame(map, i.end=c(map$V3[-1],map$V3[nrow(map)]))
      colnames(map) <- c("marker", "chr", "i.start", "i.stop")
      #compute pos, starting from 0
      # map <- cbind.data.frame(map, diffpos=c(0,diff(map$V3)))
      # map <- cbind.data.frame(map, relpos=cumsum(map$diffpos))
      
      #rename columns
      colnames(ti) <- c("id", "Sample", map$marker)
      #gather data into table format
      ti.gather <- tidyr::gather(ti, "marker", "geno", 3:ncol(ti))
      #change object type 
      ti.gather$marker <- as.numeric(ti.gather$marker)
      #sort 
      ti.gather <- dplyr::arrange(ti.gather, id, marker)
      #add start and stop positions from map
      ti.gather <- cbind.data.frame(ti.gather, 
                                    i.start=rep(map$i.start, nrow(map)*length(ti.IDs)),
                                    i.end=rep(map$i.stop, nrow(map)*length(ti.IDs)))
      #re-order levels of factor sample (for plotting)
      ti.gather$Sample <- factor(ti.gather$Sample, levels=rev(levels(ti.gather$Sample)))
      #set geno as factor
      ti.gather$geno <- as.factor(ti.gather$geno)
      #manually define levels for geno
      #level will show in plot even if it is not in data
      ## 'A':'DP'
      ## 'A2':'DP2'
      ## 'G':'RP'
      ## '1':'RP2'
      ## '0':'NA'
      ## 'C':'Het'
      levels(ti.gather$geno)[levels(ti.gather$geno)=="G"] <- "RP"
      levels(ti.gather$geno)[levels(ti.gather$geno)=="1"] <- "RP"
      levels(ti.gather$geno)[levels(ti.gather$geno)=="A"] <- "DP"
      levels(ti.gather$geno)[levels(ti.gather$geno)=="A2"] <- "DP"
      levels(ti.gather$geno)[levels(ti.gather$geno)=="C"] <- "Het"
      levels(ti.gather$geno)[levels(ti.gather$geno)=="0"] <- "NA"
      
      #read transitions table
      breakpoints <- read.table(
        paste(directory, group, "/", file4, sep=""), 
        sep='\t', header=T)
      
      #filter extra reps
      breakpoints <- breakpoints[which(breakpoints$Sample %in% sample_subset),]
      
      # breakpoints$plant <- as.numeric(sub(".*p(.+)r.", "\\1", breakpoints$Sample))
      # breakpoints$rep <- as.numeric(sub(".*r", "", breakpoints$Sample))
      #sameSamples <- cbind.data.frame(Sample=breakpoints$Sample, 
      #                                sameSamples=(c(FALSE, breakpoints$Sample) == c(breakpoints$Sample, NA))[1:length(breakpoints$Sample)])
      #sameSamples <- sameSamples[sameSamples$sameSamples != "TRUE", ]
      #sameSamples$name <- sub("s.*", "", sameSamples$Sample)
      #sameSamples$identical <- (c(FALSE, sameSamples$name) == c(sameSamples$name, NA))[1:length(sameSamples$name)]
      #sameSamples <- sameSamples[sameSamples$identical != "TRUE", ]
      #breakpoints <- breakpoints[which(breakpoints$Sample %in% sameSamples$Sample),]
      
      #edit sample names
      breakpoints$Sample <- sub("s.*", "", breakpoints$Sample)
      ## add missing samples
      breakpoints <- merge(as.data.frame(Sample), breakpoints, by="Sample", all.x=T)
      
      #get donor introgressions only (including heterozugous segments)
      breakpoints.introg <- subset(breakpoints, Hap_Geno=="DP" | Hap_Geno=="Het")
      breakpoints.introg <- subset(breakpoints, Hap_Geno=="DP")
      
      #get breakpoint resolution
      breakpoints.introg.res <- rbind.data.frame(cbind.data.frame("Sample"=breakpoints.introg$Sample,
                                                                  "Hap_Group"=breakpoints.introg$Hap_Geno,
                                                                  "Length"=breakpoints.introg$X5_BP_Density),
                                                 cbind.data.frame("Sample"=breakpoints.introg$Sample,
                                                                  "Hap_Group"=breakpoints.introg$Hap_Geno,
                                                                  "Length"=breakpoints.introg$X3_BP_Density))
      ## add missing samples
      breakpoints.introg.res <- merge(as.data.frame(Sample), breakpoints.introg.res, by="Sample", all.x=T)
      
      #re-order levels of factor sample (for plotting)
      breakpoints.introg.res$Sample <- as.factor(breakpoints.introg.res$Sample)
      breakpoints.introg.res$Sample <- factor(breakpoints.introg.res$Sample, levels=rev(levels(breakpoints.introg.res$Sample)))
      breakpoints.introg.counts <- as.data.frame(table(breakpoints.introg$Sample)) 
      #breakpoints.introg.counts$Freq[which(is.na(summ$Scheme))] <- NA
      for (s in which(is.na(summ$Scheme))){
        breakpoints.introg.counts <- add_row(breakpoints.introg.counts, Var1 = paste("NILASq",ZmPR,"g",DP_no,RP_no,"i", formatC(s, width=2, flag="0") , sep=""), Freq = NA,.before = s)
      }
      colnames(breakpoints.introg.counts) <- c("Sample", "No_Introgressions")
      breakpoints.introg.counts$Sample <- as.factor(breakpoints.introg.counts$Sample)
      breakpoints.introg.counts$Sample <- factor(breakpoints.introg.counts$Sample, levels=rev(levels(breakpoints.introg.counts$Sample)))
      
      #-------------------  Plot introgression
      
      #set colors for geno levels
      cols.1 <- c("NA"="#F0E442", "DP-T"="#0072B2", "DP-B"="#56B4E9", "Het"="#CC79A7", "RP"="#999999")
      
      summ.melt <- cbind.data.frame(ID=as.numeric(sub(".*i", "", summ.melt$Sample)), summ.melt)
      summ.melt$Scheme <- as.factor(summ.melt$Scheme)
      
      n <- paste(ZmPR, RP, DP, sep = "_")
      
      data <- rbind.data.frame(ti.gather, data)
      
      ### PLot of the genoic coordinates 
      L[[n]] <- ggplot() +
        geom_rect(data=ti.gather, 
                  aes(xmin=i.start, xmax=i.end, ymin=rev(id)-1,
                      ymax=rev(id), fill=geno)) +
        geom_segment(aes(x = rep(min(map$i.start),13), y = c(0:12), xend = rep(max(map$i.stop),13), yend = c(0:12)),
                     colour="black") +
        geom_vline(xintercept = 94430850, color = "red") + 
        geom_vline(xintercept = 94433495, color = "red") + 
        ## ZmCCT10 AGPV4 
        # geom_vline(xintercept=94435768) + 
        ggtitle(n) +
        xlab("Target Interval\n(quantile positions, bp)") +
        ylab("") + 
        scale_fill_manual(name="",
                          values = cols.1,
                          breaks = c("NA", "Het", "DP", "RP"),
                          limits = c("NA", "Het", "DP", "RP"),
                          labels = c("NA", "Het", "DP", "RP")) +
        scale_y_continuous(breaks = c(0.5:11.5), labels=rev(ti$sample), expand = c(0.1, -1.125, 0.1, -1.15)) +
        scale_x_continuous(limits = c(min(map$i.start), max(map$i.stop)), 
                           breaks = as.vector(quantile(c(min(map$i.start)):c(max(map$i.stop))))) +
        
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
              ,panel.grid.major.y = element_blank()
              ,panel.border = element_rect(colour="grey", linewidth = 0.5, fill=NA)
              ,panel.background = element_blank()
              ,plot.title = element_text(face="bold", vjust=0.375)
              ,axis.title.x = element_text(size=14, vjust=1.5)
              ,axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)
              ,axis.title.y = element_blank()
              ,axis.text.y = element_blank()
              #top, right, bottom, left
              ,plot.margin=unit(c(1.33,0.5,0.06,0.5), "cm")
        )
      
    }
  }
  # return(data)
}


write.csv(data,  paste0(file.path(directory),"/","_", date_today,  "_","mapping_data_NILAS",".csv"), row.names = T)
