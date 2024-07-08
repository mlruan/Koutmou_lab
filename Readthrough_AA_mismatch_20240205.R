setwd("/nfs/turbo/lsa-kkoutmou/migratedData/Private_Folders/mlruan/")
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("stringdist")
library(stringdist)



#Group1: AA mismatch
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#Step1: load the file with Y info, including: gene, Start, End, position, sequence_lengths... (only pick these sites with over 40% U-to-C mismatch level)
psi_sites1 <- read.csv("Protemics/result/step3f/list_ psi_sites .csv",header = TRUE)
psi_sites2 <- read.csv("Protemics/result/step3f/list_ psi_sites_old .csv",header = TRUE)
psi_sites1$Protein_ID
distinct_merged_data<-F1_merge_and_process_data(psi_sites1,psi_sites2) 
distinct_merged_data
colnames(distinct_merged_data)
distinct_merged_data$Protein_ID
#save the protein_ID sp|P31374,sp|P27825.... for step2 to get the protien sequence
writeLines(paste(distinct_merged_data$Protein_ID, collapse = ","), paste("Protemics/result/step3f/combinded_new_old_stopcpdon_all.txt"))

#########################################################################################################################
#########################################################################################################################
#Step2: get the protein sequences: downloaded from uniprot by providing the protein_ID sp|P31374,sp|P27825.... and load here 
protein_mismatch_sequences <- readLines("Protemics/result_list_mismatch/idmapping_2023_12_22.fasta")
head(protein_mismatch_sequences)
# ">sp|P00812|ARGI_YEAST Arginase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=CAR1 PE=1 SV=1"
# [2] "METGPHYNYYKNRELSIVLAPFSGGQGKLGVEKGPKYMLKHGLQTSIEDLGWSTELEPSM"                                                       
# [3] "DEAQFVGKLKMEKDSTTGGSSVMIDGVKAKRADLVGEATKLVYNSVSKVVQANRFPLTLG"                                                       
# [4] "GDHSIAIGTVSAVLDKYPDAGLLWIDAHADINTIESTPSGNLHGCPVSFLMGLNKDVPHC"                                                       
# [5] "PESLKWVPGNLSPKKIAYIGLRDVDAGEKKILKDLGIAAFSMYHVDKYGINAVIEMAMKA"                                                       
# [6] "VHPETNGEGPIMCSYDVDGVDPLYIPATGTPVRGGLTLREGLFLVERLAESGNLIALDVV"   
#change the fasta to a dataframe with "sequence" and "gene_name"
fasta_df<- F2_parse_protein_sequences(protein_mismatch_sequences)
head(fasta_df)
colnames(fasta_df)
# [1] "name"      "sequence"  "gene_name"

#########################################################################################################################
#########################################################################################################################
#Step3: add the fasta_d$gene_name to distinct_merged_data$gene_name, and calculate potential AA mismatch position: AA_sequence_loc
#AA_sequence_loc <- ceiling((filtered_psi_sites_coding_mm40_merge_sequence$position - filtered_psi_sites_coding_mm40_merge_sequence$Start + 1) / 3)
Merge_gene_name <- F3_merge_and_filter_data(distinct_merged_data, fasta_df)
head(Merge_gene_name )
colnames(Merge_gene_name)

#########################################################################################################################
#########################################################################################################################
#Step4:generate the modified sequence according to the position and start
unique_modified_df <- F4_generate_modified_sequences(Merge_gene_name)
head(unique_modified_df)
nrow(unique_modified_df)##6367 list

#########################################################################################################################
#########################################################################################################################
#Step5: check whether any 2 psi sites are close (within 25 AA), if so, rerun the generate_function to get more mismatched sequences 20*20 posibilities for mismatched sequences with 2 close Y sites
close_Y_sites<-F5_process_unique_modified_df(unique_modified_df)
unique_modified_df_close_Y_sites<- F4_generate_modified_sequences(close_Y_sites)
#Protein.Name labedled with the second mismatched AA position

#########################################################################################################################
#########################################################################################################################
#Step6: rbind the 2 datasets: 1 AA mismatch and 2 AA mismatch
unique_modified_df_rebind<-rbind(unique_modified_df, unique_modified_df_close_Y_sites)
unique_modified_df_rebind <- unique_modified_df_rebind[!duplicated(unique_modified_df_rebind$sequence), ]
nrow(unique_modified_df_rebind) #7108 potential sequences with mismatched amino acids

#########################################################################################################################
#########################################################################################################################
#Step7: make sure the AA mismatch is within the protein CDS
unique_modified_df_rebind <- subset(unique_modified_df_rebind, AA_sequence_loc>0 & AA_sequence_loc<= sequence_lengths+2)
write.csv(unique_modified_df_rebind,  "Protemics/result_list_mismatch/protein_prediction_mismatch_mm40_old_and_new_6908.csv", row.names = FALSE)





#Group2: Stop codon readthrough
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
##############################
#Step1: load the Y sites within the stop codon
Stop_codon_annotation1<- read.csv("Protemics/result/step3e/stop_codon_modfied_LTV1_chrXI_178170_ psi_sites .csv", header = TRUE)
Stop_codon_annotation2<- read.csv("Protemics/result/step3e/stop_codon_modfied_ psi_sites_old .csv", header = TRUE)
list <- c(Stop_codon_annotation2$Annotation, "LTV1") #list <- c(Stop_codon_annotation2$Annotation, Stop_codon_annotation1$Annotation) 
#pick these rows in AA mismatch dataframe: unique_modified_df_rebind
filtered_rows <- unique_modified_df_rebind[unique_modified_df_rebind$Annotation.x %in% list, ]
distinct_rows <- filtered_rows %>% distinct(gene,position, .keep_all = TRUE)

#########################################################################################################################
#########################################################################################################################
#Step2: load the reference genome
S288C_sequences <- readDNAStringSet("Nanopore_file/Yeast_NES_wt_hs/second_alignment_toS288C/S288C_reference_genome/S288C_reference_sequence_R64_4_1_20230823_chrmt.fa")
Stanford_sequences2<- readDNAStringSet("Nanopore_file/Yeast_NES_wt_hs/20230505-nanopore/BY4741_Stanford_2014_JRIS00000000.fsa")
combined_sequences <- c(S288C_sequences, Stanford_sequences2)
current_names <- names(combined_sequences)
new_names <- sub(" .*", "", current_names)
names(combined_sequences) <- new_names

#########################################################################################################################
#########################################################################################################################
#Step3: translate these sequence to amino acids for these between start and stop codon+300, stop codon will be label as "*"
result_df <- data.frame()
for (i in 1:nrow(distinct_rows)) {
  chromosome <- distinct_rows$chr[i]
  start_position <- distinct_rows$Start[i]
  #include 300 bp sequence after the stop codon
  end_position <- distinct_rows$position[i]+300
  annotation_value <- distinct_rows$Annotation.x[i]
  
  sub_sequence <- subseq(combined_sequences[[chromosome]], start_position, end_position)
  #translate the DNA sequence to the amino acid
  protein_sequence<-translate(sub_sequence)
  protein_sequence_AA <- as.character(protein_sequence)
  # Create a new row with relevant information
  result_row <- data.frame(
    Annotation = annotation_value,
    Chromosome = chromosome,
    Start = start_position,
    End = end_position,
    protein_sequence_AA = protein_sequence_AA
  )
  # Append the row to the result dataframe
  result_df <- rbind(result_df, result_row)
}
result_df$length<-result_df$End-result_df$Start
result_df <- result_df[order(-result_df$length), ]
result_df <- result_df[!duplicated(result_df$Annotation), ]
head(result_df)

#########################################################################################################################
#########################################################################################################################
#step4: seperate the protein _sequence_AA into several readthrough_fragment depends on the "*", and name them as "protein_sequence_AA1", "protein_sequence_AA2"...
result_df <- separate(result_df, protein_sequence_AA, into = c("protein_sequence_AA1", "protein_sequence_AA2", "protein_sequence_AA3", "protein_sequence_AA4"), sep = "\\*", fill = "warn")
result_df <- merge(result_df, fasta_df[, c("gene_name", "sequence")], by.x = "Annotation", by.y = "gene_name", all.x = TRUE)
#test whether the protein sequence is right by using is_sequence
result_df <- mutate(result_df, is_sequence = protein_sequence_AA1 == sequence)
head(result_df)
#Step5: generate a readthough column us AA2 or AA3
result_df$readthrough <- ifelse((result_df$protein_sequence_AA2!=  "" ), result_df$protein_sequence_AA2, result_df$protein_sequence_AA3)

#########################################################################################################################
#########################################################################################################################
#Step6: generate the datasets with Y in stopcodon from the AA mismatch dataframe (we generated earlier in Group1: with 20 different potential AA in the Y incoportated codon)
Stopcodon_df <- unique_modified_df_rebind[
  unique_modified_df_rebind$Annotation.x %in% result_df$Annotation &
    unique_modified_df_rebind$position >= unique_modified_df_rebind$End - 2 &
    unique_modified_df_rebind$position <= unique_modified_df_rebind$End,
]
head(Stopcodon_df )

#########################################################################################################################
#########################################################################################################################
#Step7: add the readthruogh sequence from result_df to the sequence in Stopcodon_df
for (i in seq_len(nrow(result_df))) {
  Stopcodon_df <- Stopcodon_df %>%
    mutate(
      sequence = ifelse(
        Annotation.x == result_df$Annotation[i] ,
        paste0(sequence, result_df$readthrough[i]),
        sequence
      )
    )
}
head(Stopcodon_df )

#########################################################################################################################
#########################################################################################################################
#Step8: paste "stopcodon" to the Protein.Name to easy find these catogoris in the future
Stopcodon_df$Protein.Name <- ifelse(
  Stopcodon_df$Annotation.x %in% result_df$Annotation,
  paste(Stopcodon_df$Protein.Name, "stopcodon", sep = "_"),
  Stopcodon_df$Protein.Name
)
write.csv(Stopcodon_df,  "result_list_mismatch/protein_prediction_mismatch_mm40_old_new_7919_stop_codon.csv", row.names = FALSE)







#Group3: combine 2 groups togeher and export fasta files for DIA analysis
########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#rbind the AA mismatch and stopcodon reahthroug together
unique_modified_df_rebind2 <- rbind(unique_modified_df_rebind, Stopcodon_df)
write.csv(unique_modified_df_rebind2,  "result_list_mismatch/protein_prediction_mismatch_mm40_old_new_7919_stop_codon.csv", row.names = FALSE)

###save the dataset into a fasta: >Protein.Name, sequence (this file can be used for the DIA analysis to see whether any Y caused mismatch and stopcodon readthrough are detected)
fasta_file_path <- "result_list_mismatch/Minli20231223_protein_mismatch_mm40_all_stopcodon_7919.fasta"
# Open the file for writing
file_conn <- file(fasta_file_path, "w")
for (i in 1:nrow(unique_modified_df_rebind2)) {
  gene <- unique_modified_df_rebind2$gene[i]
  Annotation.x <- unique_modified_df_rebind2$Annotation.x[i]
  protein_name <- unique_modified_df_rebind2$Protein.Name[i]
  sequence <- unique_modified_df_rebind2$sequence[i]
  writeLines(paste0(">", protein_name,'_',gene, '_', Annotation.x), file_conn)
  writeLines(sequence, file_conn)
}
close(file_conn)



