## load packages
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
#library(bbc_style)
library(gridExtra)


dataset_dir = "dp6q30dr0.40maf0.01"
# map_file = "nilas_f4_dp6q30dr0.40maf0.01.alleleInfo"
# gtype_file = "nilas_f4_dp6q30dr0.40maf0.01.GT.FORMAT.final"
##set working directory
directory <- paste0("~/Documents/Drouault_etal_paper2/01_paperdev/ATLAS/paper/Stage_Eloise/Data/", dataset_dir, "/outputs/")
# List all directories in the current working directory




## loop across groups
#### !!!! Isssue -- geno sum has no NA
## "2369 x CML373 ZmPR3" 
## "LH123Ht x CML373 ZmPR2" 
## "LH123Ht x CML341 ZmPR3"  

ZmPR_list <- 1:4
rp_list <- c("2369", "LH123Ht")
dp_list <- c("CML10", "CML258", "CML277", "CML373", "CML341", "Tzi8", "Tzi9")

# ZmPR_list <- 1:4
# rp_list <- c("2369")
# dp_list <- c("CML10")

# ZmPR_list <- 2
# rp_list <- c("2369")
# dp_list <- c("CML258")


for (l in 1:length(ZmPR_list)) {
  for (r in 1:length(rp_list)) {
    for (d in 1:length(dp_list)) {
      
      l <- 4
      r <- 1
      d <- 2
      
      ZmPR <- ZmPR_list[l]
      RP <- rp_list[r]
      DP <- dp_list[d]


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



#### read and compile data from haploimputer

summ <- read.table(
  paste(directory, group, "/", file1, sep=""), 
  sep='\t', header=T)
summ$Complement <- grepl("C", summ$Sample, fixed=FALSE)

summ <- cbind.data.frame(summ[c("Sample","Scheme","ZmPR")],
                         "Recurrent_Parent"=rowSums(summ[,c("RP_BG","RP_FG")]),
                         summ[,c("DP_BG","DP_FG")],
                         "Heterozygous"=rowSums(summ[,c("Het_BG","Het_FG")]),
                         "Missing"=rowSums(summ[,c("NA_BG","NA_FG")]),
                         "Complement"=summ$Complement)
#filter complements
summ <- summ[summ$Complement != "TRUE", ]
summ <- summ[,-which(names(summ)=="Complement")]
#filter extra reps
#summ$plant <- as.numeric(sub(".*p(.+)r.", "\\1", summ$Sample))
#summ$rep <- as.numeric(sub(".*r", "", summ$Sample))
summ$name <- sub("s.*", "", summ$Sample)

summ %>% group_by(name) %>% filter(DP_FG == max(DP_FG)) %>% filter(Missing == min(Missing)) %>% filter(Heterozygous == min(Heterozygous)) -> summ
sample_subset <- summ$Sample


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

    }}}

Q <- ggplot() +
  geom_rect(data=ti.gather, 
            aes(xmin=i.start, xmax=i.end, ymin=rev(id)-1,
                ymax=rev(id), fill=geno)) +
  geom_segment(aes(x = rep(min(map$i.start),13), y = c(0:12), xend = rep(max(map$i.stop),13), yend = c(0:12)),
               colour="black") + 
  ggtitle("d") +
  xlab("Target Interval\n(quantile positions, bp)") +
  ylab("") + 
  # scale_fill_manual(name="",
  #                   values = cols.4,
  #                   breaks = c("NA", "Het", "DP", "RP"),
  #                   limits = c("NA", "Het", "DP", "RP"),
  #                   labels = c("NA", "Het", "DP", "RP")) +
  scale_y_continuous(breaks = c(0.5:11.5), labels=rev(unique(ti.gather$id)), expand = c(0.1, -1.125, 0.1, -1.15)) +
  scale_x_continuous(limits = c(min(map$i.start), max(map$i.stop)), 
                     breaks = as.vector(quantile(c(min(map$i.start)):c(max(map$i.stop))))) +
  
  theme(plot.background = element_blank()
        # ,legend.position = "none"
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


### First defines how much the introgression have DP in common 

id.start.n <- 114449759
id.end.n <- 130678308
id.start.compare <- 126709380
id.end.compare <- 130678308


calculation_percentage_in_common <- function(id.start.n , id.end.n, id.start.compare, id.end.compare){
  if(id.start.compare < id.start.n & id.start.n < id.end.compare & id.end.n > id.end.compare){
    shared = (id.end.compare - id.start.n) / (id.end.n - id.start.n)
  }
  
  if(id.start.n < id.start.compare  &  id.end.n > id.end.compare){
    shared = (id.end.compare - id.start.compare) / (id.end.n - id.start.n)
  }
  if(id.start.n < id.start.compare  &  id.end.n == id.end.compare){
    shared = (id.end.compare - id.start.compare) / (id.end.n - id.start.n)
  }
  if ( id.start.n < id.start.compare &  id.end.n < id.end.compare) {
    shared   = ( id.end.n - id.start.compare ) / (id.end.n - id.start.n)}

  if ( id.start.n < id.start.compare &  id.end.n < id.end.compare) {
    shared   = ( id.end.n - id.start.compare ) / (id.end.n - id.start.n)}
  
  if ( id.start.n == id.start.compare &  id.end.n > id.end.compare) {
    shared   = (id.end.compare- id.start.compare)/(id.end.n - id.start.compare)}
  if ( id.start.n == id.start.compare &  id.end.n < id.end.compare) {
    shared   = 100}
  
  
  
  if ( id.start.n > id.start.compare & id.end.n < id.end.compare) {
   shared = 100}
    if (id.start.n > id.end.compare){
      shared = 0
    }
  if (id.end.n < id.start.compare){
    shared = 0
  }
  return(shared)
}

ti.gather_cropped <- ti.gather %>%
  dplyr::filter(geno == "DP") %>%
  dplyr::group_by(Sample, id) %>%
  dplyr::mutate(min.i.start = min(i.start , na.rm = T),
                max.i.end = max(i.end , na.rm = T)) %>%
  dplyr::select(id, Sample, geno, min.i.start, max.i.end) %>%
  dplyr::distinct()


 splitted_common_result <- lapply(1:12, function(x) { 
    # x <- 1
    group_lines  <-  ti.gather_cropped %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Sample, id) %>%
      dplyr::mutate(!!paste0("percentage_dp_common_with_line_", x) := ifelse(
        id == x,
        100,
        calculation_percentage_in_common(
          min.i.start,
          max.i.end,
              id.start.compare = .$min.i.start[which(.$id == x)],
              id.end.compare = .$max.i.end[which(.$id == x)])))
  })
 
 library(purrr)
 
 result <- reduce(
   splitted_common_result,
   ~ left_join(.x, .y, by = c("id", "Sample", "geno", "min.i.start", "max.i.end")),
   .init = ti.gather_cropped
 )
 
 head(result)
### Then associates the phenotypes to the dp 

func = function(x) {
  return(x)
}
num = func(8)
num
