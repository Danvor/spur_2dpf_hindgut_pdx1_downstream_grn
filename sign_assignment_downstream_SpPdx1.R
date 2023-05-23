Sys.setenv(LANG = "en")
set.seed(255)
libraries <- c("plyr", "dplyr", "tibble")
lapply(libraries, library, character.only = TRUE)
setwd("~/Work/Collaborations/downstream_pdx1")
#clear variables
rm(list=ls())
#load the tf/target and the gene1 differential expression targets within blastopore as dataframes
#make sure the tables imported have appropriate column_names
targets <- read.delim("./sorted_SpPdx1_affected_WHLs_foldchange_Sp48hpf.txt", header=TRUE) #filtered DE analysis results file, containing Sp-Pdx1 MO differentially expressed hindgut gene IDs and their fold change.
tf_target <- read.delim("./sorted_gut_peak_hindgut_Sp_48hpf_tfWHL_gWHL_TFs_SpPdx_positive.txt", header=TRUE) #interaction file from TF footprinting and filtering to contain only IDs from cells of interest (Sp-Pdx1 co-expressed in Hindgut (1)).
#get only those genes that are affected by gene1 directly
tti1 <- subset(tf_target, effector_gene == "WHL22.169409")
#Get all genes that are affected by direct gene1 targets, and then find genes directly affected by those
x <- 2
repeat {
  assign(paste0('tti', x), subset(tf_target, subset = effector_gene %in% get(paste0('tti', x-1))[['target_gene']]))
  if (isTRUE(all.equal(get(paste0('tti',x)), get(paste0('tti',x-1))))) {break}
  x = x + 1
}

#Get the differential expression information for every affected gene, by each effector gene
da1 <- left_join(tti1, targets, by = "target_gene")
y <- 2


#also make copies of resulting differential expression info tables in order to trace back effectors and targets between different dataframes
repeat {
  assign(paste0('ind', y-1), anti_join(get(paste0('tti',y)), get(paste0('tti', y-1))))
  assign(paste0('da', y), left_join(get(paste0('ind', y-1)), targets, by = "target_gene"))
  assign(paste0('da', y, 'h'), get(paste0('da', y))%>% 
           rename(
             output_gene = target_gene,
             dynamic = fold_change
           ) %>% 
           mutate(
             target_gene = effector_gene
           )%>% 
           rename(
             input_gene = effector_gene
           ))
  y = y + 1
}

hdf1 <- da1 %>% rename(output_gene = target_gene,
                       dynamic = fold_change,
                       input_gene = effector_gene) %>% add_column(fold_change = -1, .before = 1) #any negative foldchange to account for MO injection
z <- 2
repeat {
  assign(paste0('hdf', z), join(get(paste0('da', z - 1)), 
                                get(paste0('da', z, 'h')), by = "target_gene", type = "left") 
         %>% mutate(dynamic=replace(dynamic, is.na(fold_change), NA)))
  z = z + 1
}

#To add dynamic correctly, the iteration number added as a separate column
lst1 <- mget(ls(pattern = '^hdf\\d+$'))
lst1 <- head(lst1, -1)#this just removes the last table from iterations which is empty.
v1 <- as.numeric(sub("hdf", "", names(lst1))) -1
lst1 <- Map(cbind, lst1, iteration = v1)
all_int <- lapply(lst1, subset, select = c(fold_change, input_gene, output_gene, dynamic, iteration)) %>% 
  bind_rows() %>% 
  subset(subset = !is.na(input_gene) & !is.na(output_gene))


ai1 <- all_int[,c(1:2,5)]
ai2 <- all_int[,3:4]
ai2 %>% mutate(dynamic=ifelse(is.na(dynamic)==T,0,dynamic)) %>%
  group_by(output_gene) %>% summarise(dynamic = sum(as.numeric(dynamic))) -> ai2_new

join <- left_join(ai1,ai2_new,by=c("input_gene"="output_gene"))
join %>% mutate(fold_change=ifelse(dynamic==0,NA,fold_change)) %>%
  mutate(output_gene=all_int$output_gene) %>%
  relocate(fold_change,input_gene,output_gene,dynamic,iteration) %>%
  mutate(dynamic=all_int$dynamic) -> result

all_int <- result

all_int$fold_change[all_int$input_gene == "WHL22.169409"] <- "-1"
all_int$dynamic[is.na(all_int$fold_change)] <-0
all_int$sign <- "neutral"
all_int$sign[all_int$fold_change<0 & all_int$iteration %% 2 == 0 & all_int$dynamic>0] <- "negative"
all_int$sign[all_int$fold_change<0 & all_int$iteration %% 2 != 0 & all_int$dynamic>0] <- "positive"#odd iteration the sign stays
all_int$sign[all_int$fold_change<0 & all_int$iteration %% 2 == 0 & all_int$dynamic<0] <- "positive"
all_int$sign[all_int$fold_change<0 & all_int$iteration %% 2 != 0 & all_int$dynamic<0] <- "negative"#odd iteration the sign stays

all_int$sign[all_int$fold_change>0 & all_int$iteration %% 2 == 0 & all_int$dynamic>0] <- "positive"
all_int$sign[all_int$fold_change>0 & all_int$iteration %% 2 != 0 & all_int$dynamic>0] <- "positive" #iteration does not matter here but stays for code consistency
all_int$sign[all_int$fold_change>0 & all_int$iteration %% 2 == 0 & all_int$dynamic<0] <- "negative"
all_int$sign[all_int$fold_change>0 & all_int$iteration %% 2 != 0 & all_int$dynamic<0] <- "negative" #iteration does not matter here but stays for code consistency

grn_draw <- all_int[, c("input_gene", "output_gene", "sign")]
write.table(all_int, file = "Sp48hpf_gut_peaks_SpPdx_expressing_interaction_table_output_SpPdx_MO.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
