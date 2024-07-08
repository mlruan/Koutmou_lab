
F1_merge_and_process_data <- function(data1, data2) {
  # Remove 'Annotation' column from data1
  data1 <- data1[, !names(data1) %in% "Annotation"]
  
  # Rename 'Start' and 'End' columns in data2
  data2$Start <- data2$X.2
  data2$End <- data2$X.3
  
  # Copy 'Annotation.x' column from data1 to 'Annotation' column
  data1$Annotation <- data1$Annotation.x
  
  # Copy 'Systematic.Name' column from data2 to 'gene' column
  data2$gene <- data2$Systematic.Name
  
  # Find common columns between data1 and data2
  common_columns <- intersect(names(data1), names(data2))
  
  # Merge data1 and data2 based on common columns
  merged_data <- rbind(data1[, common_columns], data2[, common_columns])
  
  # Print the 'gene' column of merged data
  print(merged_data$gene)
  
  # Return the merged data
  return(merged_data)
}



F2_parse_protein_sequences <- function(protein_mismatch_sequences) {
  # Count the number of lines starting with ">"
  num_sequences <- sum(grepl("^>", protein_mismatch_sequences))
  
  # Initialize empty vectors to store names and sequences
  names <- character()
  sequences <- character()
  
  # Initialize variables to store current sequence and name
  current_sequence <- character()
  current_name <- ""
  
  # Iterate through each line of protein_mismatch_sequences
  for (line in protein_mismatch_sequences) {
    if (substr(line, 1, 1) == ">") {
      # Save the previous sequence and start a new one
      if (length(current_sequence) > 0) {
        names <- c(names, current_name)
        sequences <- c(sequences, paste(current_sequence, collapse = ""))
      }
      
      # Extract the name of the sequence
      current_name <- substr(line, 2, nchar(line))
      
      # Reset the current sequence
      current_sequence <- character()
    } else {
      # Add the line to the current sequence
      current_sequence <- c(current_sequence, line)
    }
  }
  
  # Add the last sequence
  names <- c(names, current_name)
  sequences <- c(sequences, paste(current_sequence, collapse = ""))
  
  # Create a data frame with names and sequences
  fasta_df <- data.frame(name = names, sequence = sequences)
  
  # Remove duplicated entries
  fasta_df <- fasta_df[!duplicated(fasta_df), ]
  
  # Extract gene names from sequence names
  fasta_df$gene_name <- str_extract(fasta_df$name, "(?<=GN=)[^\\s]+")
  
  # Return the resulting data frame
  return(fasta_df)
}


F3_merge_and_filter_data <- function(filtered_psi_sites_coding_mm40, fasta_df) {
  filtered_psi_sites_coding_mm40_merge_sequence <- left_join(filtered_psi_sites_coding_mm40, fasta_df, by = c("Annotation" = "gene_name")) 
  
  filtered_psi_sites_coding_mm40_merge_sequence$Start <- as.numeric(filtered_psi_sites_coding_mm40_merge_sequence$Start)
  filtered_psi_sites_coding_mm40_merge_sequence$AA_sequence_loc <- ceiling((filtered_psi_sites_coding_mm40_merge_sequence$position - filtered_psi_sites_coding_mm40_merge_sequence$Start + 1) / 3)
  
  filtered_psi_sites_coding_mm40_merge_sequence$Protein.Name <- paste(filtered_psi_sites_coding_mm40_merge_sequence$Protein_ID, filtered_psi_sites_coding_mm40_merge_sequence$Annotation, filtered_psi_sites_coding_mm40_merge_sequence$AA_sequence_loc, sep = "|")
  filtered_psi_sites_coding_mm40_merge_sequence$Annotation.x<- filtered_psi_sites_coding_mm40_merge_sequence$Annotation
  
  new_dataframe <- filtered_psi_sites_coding_mm40_merge_sequence %>%
    select(gene, Annotation.x, chr, Start, End, position, mm_IVT, kmer, sequence, AA_sequence_loc, Protein.Name)
  
  new_dataframe$sequence_lengths <- nchar(new_dataframe$sequence)
  
  df <- distinct(new_dataframe)
  
  return(df)
}



F4_generate_modified_sequences <- function(df) {
  generate_random_amino_acid <- function() {
    amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    return(sample(amino_acids, 1))
  }
  
  modify_sequence <- function(sequence, position, new_amino_acid) {
    modified_sequence <- strsplit(sequence, '')[[1]]
    modified_sequence[position] <- new_amino_acid
    
    return(paste(modified_sequence, collapse = ''))
  }
  
  modified_df <- data.frame(stringsAsFactors = FALSE)
  
  for (i in 1:nrow(df)) {
    original_sequence <- df$sequence[i]
    original_position <- df$AA_sequence_loc[i]
    
    original_row <- data.frame(
      gene = df$gene[i],
      Annotation.x = df$Annotation.x[i],
      chr = df$chr[i],
      Start = df$Start[i],
      End = df$End[i],
      position = df$position[i],
      mm_IVT = df$mm_IVT[i],
      kmer = df$kmer[i],
      sequence = original_sequence,
      AA_sequence_loc = original_position,
      Protein.Name = df$Protein.Name[i],
      sequence_lengths = df$sequence_lengths[i]
    )
    
    modified_df <- rbind(modified_df, original_row)
    
    for (new_amino_acid in c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")) {
      modified_sequence <- modify_sequence(original_sequence, original_position, new_amino_acid)
      
      new_row <- data.frame(
        Protein.Name = paste0(df$Protein.Name[i], "_modified_", new_amino_acid),
        Annotation.x = df$Annotation.x[i],
        chr = df$chr[i],
        Start = df$Start[i],
        End = df$End[i],
        position = df$position[i],
        mm_IVT = df$mm_IVT[i],
        kmer = df$kmer[i],
        sequence = modified_sequence,
        AA_sequence_loc = original_position,
        gene = df$gene[i],
        sequence_lengths = df$sequence_lengths[i]
      )
      
      modified_df <- rbind(modified_df, new_row)
    }
  }
  
  unique_modified_df <- modified_df[!duplicated(modified_df$sequence), ]
  return(unique_modified_df)
}

F5_process_unique_modified_df <- function(unique_modified_df) {
  result <- unique_modified_df %>%
    group_by(gene) %>%
    mutate(has_difference_lt_25 = any(abs(diff(AA_sequence_loc)) > 0 & abs(diff(AA_sequence_loc)) < 25))
  
  true <- result %>%
    filter(has_difference_lt_25 == TRUE) %>%
    group_by(gene) %>%
    mutate(
      orig_AA_sequence_loc = AA_sequence_loc,
      AA_sequence_loc = ifelse(AA_sequence_loc == max(AA_sequence_loc), min(AA_sequence_loc), max(AA_sequence_loc))
    )
  
  return(true)
}
