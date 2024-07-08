### Code of Codon preference of Y by PUS7
library(Rsamtools)
library(rtracklayer)
library(BSgenome)
library(dplyr)
library(ggplot2)

# load your dataset: with Annotation2, position of Y, strand, kmer
filterD1 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D1a.csv", header = TRUE)
# filterD2 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D2a.csv", header = TRUE)
# filterD3 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D3a.csv", header = TRUE)
# filterD4 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D4a.csv", header = TRUE)
# filterD5 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D5a.csv", header = TRUE)
# filterD6 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D6a.csv", header = TRUE)
# filterD7 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D7a.csv", header = TRUE)
# filterD8 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D8a.csv", header = TRUE)
# filterD9 <- read.csv(file = "/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S5_TPM_Y_plot_5mer/CSV_file_new/filtered001_mm40_D9a.csv", header = TRUE)
head(filterD1)

#load the gff files with: Annotation2, start.point, end.point, 
genes.df2<-read.csv("/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S3_p_value/Step1_PileUp/Reference_files/UTRs/3_UTR/UTR3_final_new_output_file")
head(genes.df2)

#run the fllowing function, condon_preference, to get the codon preference and the position of Y in the codon
#sample ="filterD1"
condon_preference<- function(filterD1,sample){
  head(genes.df2$Annotation)
  head(filterD1$Annotation2)
  combine<-left_join(filterD1, genes.df2, by=c("Annotation2"="Annotation"))
  combine_select <- select(combine,Annotation2, position, start.point, end.point, strand.x, kmer)
  
  #only take these rows in the coding region
  combine_select  <- combine_select[!grepl("three_prime", combine_select$Annotation2), ]
  combine_select  <- combine_select[!grepl("five_prime", combine_select$Annotation2), ]  
  combine_select  <- combine_select[!grepl("intron", combine_select$Annotation2), ]  
  head(combine_select$Annotation2) 
  
  #generate the position of Y in codon first:
  combine_select$position_Y <- ifelse(combine_select$strand.x == "+",
                                      (combine_select$position - combine_select$start.point) %% 3 + 1,
                                      (combine_select$end.point - combine_select$position) %% 3 + 1)
  
  #generate the codon based in the position_Y and kmer
  combine_select$codon<- NA
  combine_select$codon[combine_select$position_Y == 1] <- substr(combine_select$kmer[combine_select$position_Y == 1], 3, 5)
  combine_select$codon[combine_select$position_Y == 2] <- substr(combine_select$kmer[combine_select$position_Y == 2], 2, 4)
  combine_select$codon[combine_select$position_Y == 3] <- substr(combine_select$kmer[combine_select$position_Y == 3], 1, 3)
  
  #change all the T to U
  combine_select$codon <- gsub("T", "U", combine_select$codon)
  head(combine_select$codon )
  head(combine_select$position_Y)
  
  #make the table showing the codon and Y position
  count_table <- table(combine_select$codon)
  print(count_table)
  count_table <- table(combine_select$codon, combine_select$position_Y)
  count_table
  count_df <- as.data.frame.table(count_table)
  colnames(count_df) <- c("codon", "position_Y", "count")
  print(count_df)
  count_df$position_Y <- factor(count_df$position_Y, levels = c("1", "2", "3"))
  

  # plot the result
  head(count_df)
  ggplot(count_df, aes(x = codon, y = count, fill = position_Y)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Codon", y = "Count", fill = "Position_Y") +
    ggtitle("Counts of Codons by Position_Y") +
    theme_minimal() +
    theme_classic(base_size = 14) +    scale_fill_brewer(palette="Blues") +# Use a base font size appropriate for publication
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin=margin(t=10)),
      axis.title.y = element_text(margin=margin(r=10)),
      legend.title = element_text(),
      legend.text = element_text(),
      legend.key = element_blank(),  # Removes the background for legends keys
      legend.box.background = element_rect(color="black", fill=NA, linewidth=1),  # Adds a black border around the legend
      legend.position = "right",  # Adjust legend position as needed
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill="white", color=NA),
      panel.border = element_rect(color="black", fill=NA, linewidth=1),
      plot.background = element_rect(fill = "white", color = NA),  # Ensure plot background is white
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),  # Adjust margins around the plot
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16))  
  #save the csv and plot into your path
  write.csv(count_df, paste0("/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S25_codon_preference/CSVs/condon_preference_CDS_",sample,".csv"))
  ggsave(paste0("/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/Nanopore_file/Third_try_240430/S25_codon_preference/plots/condon_preference_CDS_",sample,".png"), width =6, height = 4, dpi = 300)
}

# Iterate over group names
names_to_iterate <- seq(1,9)
for (name in names_to_iterate) {
  # Generate label dynamically
  label <- paste0("filterD", name)
  condon_preference(get(label), label)
}
