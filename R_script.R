#Loading required libraries for the analysis.
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gghighlight)
library(ggrepel)
library(Biostrings)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)

#Creating a processing workflow...
process_sample <- function(sample_name) {
  file_path <- paste0("path/to/the/directory", sample_name, "_rep1_RBD_counts.tsv")
  data <- read.delim(file_path)

  data_sel <- data %>% select(2, 3, 14, 16, 18, 20)
  data_sel$total_hits <- rowSums(data_sel[,3:6])
  data_sel[, c("A.perc", "T.perc", "G.perc", "C.perc")] <- data_sel[, c("A", "T", "G", "C")] / data_sel$total_hits * 100
  data_sel <- data_sel %>% select(1, 2, 8, 9, 10, 11)
  colnames(data_sel) <- c("position", "fa_ref", "A", "T", "G", "C")

  CS <- data_sel[, 3:6] %>% rowwise() %>% mutate(row_max = names(.)[which.max(c_across(everything()))])
  data_sel$cons_seq <- CS$row_max
  data_sel$nuc_mutation <- data_sel$fa_ref != data_sel$cons_seq
  data_sel$sig_det <- data_sel$A > 10 & data_sel$T > 10 | data_sel$A > 10 & data_sel$G > 10 | 
                      data_sel$A > 10 & data_sel$C > 10 | data_sel$T > 10 & data_sel$G > 10 | 
                      data_sel$T > 10 & data_sel$C > 10 | data_sel$G > 10 & data_sel$C > 10
  data_sel$nuc_mutation <- data_sel$nuc_mutation == TRUE | data_sel$sig_det == TRUE

  columns_to_check <- c("A", "T", "G", "C")
  df <- data_sel[,3:6] %>% rowwise() %>%
    mutate(second_max_index = ifelse(sum(c_across(all_of(columns_to_check)) > 10) >= 1, order(-c_across(all_of(columns_to_check)))[2], 1)) %>%
    mutate(second_max_column = ifelse(second_max_index > 1, names(.)[second_max_index], names(.)[1]))
  data_sel$sus_seq <- df$second_max_column

  data_sel$alt_seq <- ifelse(data_sel$fa_ref != data_sel$cons_seq, data_sel$cons_seq, data_sel$sus_seq)
  data_sel$alt_seq <- ifelse(data_sel$sig_det == TRUE, data_sel$alt_seq, data_sel$cons_seq)

  # Generate DNA string
  reference <- paste(data_sel$fa_ref, collapse = "")
  reference <- DNAString(reference)
  reference_aa <- Biostrings::translate(reference) |> as.character() |> s2c()
  data_sel$aa_ref <- rep(reference_aa, each=3)

  variable <- paste(data_sel$alt_seq, collapse = "")
  variable <- DNAString(variable)
  variable_aa <- Biostrings::translate(variable) |> as.character() |> s2c()
  data_sel$aa_seq <- rep(variable_aa, each=3)

  data_sel$position <- as.factor(data_sel$position)
  data_sel_wide <- melt(data_sel, id.vars = c("position", "fa_ref", "alt_seq", "sig_det", "nuc_mutation", "cons_seq", "sus_seq", "aa_ref", "aa_seq"))
  data_sel_wide$nuc_mutation <- data_sel_wide$fa_ref != data_sel_wide$variable & data_sel_wide$value > 10
  data_sel_wide$non_syn_mutation <- data_sel_wide$nuc_mutation & (data_sel_wide$aa_ref != data_sel_wide$aa_seq)
  data_sel_wide$annotate <- paste(data_sel_wide$position, ":", data_sel_wide$fa_ref, ">", data_sel_wide$variable)

  # Create and return final object
  return(list(
    name = sample_name,
    data_sel_wide = data_sel_wide,
    fasta = t(data_sel %>% select(alt_seq)) |> `rownames<-`(sample_name)
  ))
}

#Applying the workflow on multiple samples.
sample_list <- c("CL00", "CL01", "CL02", "CL03", "CL04", "CL05", "CL06", "CL07", "CL08", "CL09", "CL10", "CL11", "CL12", "CL13", "CL14", "CL15", "CL16", "CL17", "CL18", "CL19", "CL20", "CL21", "CL22", "CL23", "CL24", "CL25", "CL26", "CL27", "CL28", "CL29", "CL30", "CL31", "CL32", "CL33", "CL34", "CL35", "CL36", "CL37", "CL38")
results_list <- lapply(sample_list, process_sample)






























