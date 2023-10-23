# Functions for converting file from different soueces
# 

#' Convert proteome discoverer ouput to R
#'
#' @param pd_psm_file A tab-delimited txt file from protoeme discoverer's 
#' PSMs output. It should contains Abundance columns as well as "Average Reporter S/N".
#' A example file pd_example_export_PSMs.txt is available in data folder. 
#' @param pd_protein_file A tab-delimited txt file from protoeme discoverer's 
#' Proteins output. The file is optional, only used for getting gene name for each protein. 
#' It must contains "Gene.Symbol" and "Accession" columns
#' @param delimiter The delimiter for the files. For tsv, use "\t"; for csv, use ",".
#' @return A r dataframe of the input file, or a string start with "Error:" containing 
#' the error message. The csv file correspondings to the df will also be save in the folder. 
#' @examples
#' df <- convertPdFile("data/pd_example_export_PSMs.txt", "data/pd_example_export_Proteins.txt")
#' @export
convertPdFile <- function(pd_psm_file, pd_protein_file = "", delimiter = "\t") {
  if (!file.exists(pd_psm_file)) {
    return("Error: file not exist.")
  }
  outputFile = paste(pd_psm_file, "_converted.csv", sep="")
  
  
  # read protein file to get gene name if needed. 
  #read csv not allow special characters in name
  protein_gene_map <- list()
  if (file.exists(pd_protein_file)) {
    protein_df <- read.csv(pd_protein_file,sep=delimiter, stringsAsFactors=FALSE)
    protein_gene_map <- setNames(protein_df[,"Gene.Symbol"],protein_df[,"Accession"])
    rm(protein_df)
  }
  
  # Export to CSV without row names
  #write.csv(df,file=outputFile, row.names=FALSE)
  
  # Read file line by line to handle large files
  input = file(pd_psm_file, "r")
  output = file(outputFile,"w")
  headerLine = readLines(input, n = 1)
  headerLine <- gsub("\"","",headerLine) # remove \"
  if ( length(headerLine) == 0 ) {
    return("Error: file is empty")
  }
  headers <- strsplit(headerLine, delimiter)[[1]]
  
  # get index 
  # some version of pd uses "Average Reporter S/N" while some uses "Average Reporter SN"
  avgSN_index <- grep("Average Reporter S", headers, fixed=TRUE)
  
  # two possible matches: Master Protein Accessions and	Protein Accessions
  protein_index <- grep("Protein Accessions", headers, fixed=TRUE)
  seq_index <- which(headers == "Annotated Sequence")
  ptm_index <- which(headers == "Modifications")
  plex_index <- which(headers == "File ID")
  # some version of pd uses "Abundance: 126" while some uses "Abundance 126"
  abd_index <- which(grepl("Abundance", headers, fixed=TRUE))
  
  if (length(avgSN_index) == 0 ) { 
    return("Error: missing Average Reporter S/N column")
  }
  if (length(protein_index) == 0 ) { 
    return("Error: missing Master Protein Accessions or	Protein Accessions column")
  }
  if (length(seq_index) == 0 ) { 
    return("Error: missing File ID column")
  }
  if (length(plex_index) == 0 ) { 
    return("Error: missing Average Reporter S/N column")
  }
  if (length(abd_index) == 0 ) { 
    return("Error: missing Abundance columns")
  }
  
  fixedColumns = c(
    "PA.Gene.Symbol",
    "Protein.ID",
    "Peptide",
    "Plex"
  )
  
  # channels looks like "126"  "127n" "127c" etc
  channels <- tolower(gsub("Abundance:? ", "", headers[abd_index]))
  print(channels)
  intensityChannels <- paste("X", channels, ".Adjusted.Intensity", sep="")
  snChannels <- paste("X", channels, ".Sn" , sep="")
    
  newHeader <- c(fixedColumns, intensityChannels, snChannels)

  writeLines(paste(newHeader, collapse=","), output)
  
  while ( TRUE ) {
    row = readLines(input, n = 1)
    if ( length(row) == 0 ) {
      break
    }
    
    rowData <- strsplit(row, delimiter)[[1]]
    rowData <- gsub("\"","",rowData) # remove \"
    
    # skip if average SN is NA or 0
    if (is.na(rowData[avgSN_index]) | rowData[avgSN_index] == "" | as.numeric(rowData[avgSN_index]) == 0) {
      next
    }
    
    
    protein <- rowData[protein_index[1]]
    gene <- ""

    # MSTR-135: fix parsing bug
    if (protein %in% names(protein_gene_map)) {
      gene <- protein_gene_map[[protein]]
    }
    
    cleanedPeptide <- gsub("\\[", "", rowData[seq_index[1]])
    cleanedPeptide <- gsub("\\]", "", cleanedPeptide)
    newFixedColumns <- c(
      gene,  
      protein,
      paste(cleanedPeptide,rowData[ptm_index[1]]),
      rowData[plex_index[1]]
    )
    
    intensityChannelsData<- as.numeric(rowData[abd_index])
    
    #convert NA to 0 for noise calculation
    intensityChannelsData_noNA <- intensityChannelsData
    intensityChannelsData_noNA[is.na(intensityChannelsData_noNA)] <-0 
    noise <- mean(intensityChannelsData_noNA) / as.numeric(rowData[avgSN_index])
    
    if (noise > 0) {
      snChannelsData <- intensityChannelsData_noNA / noise
    } else {
      snChannelsData <- rep(0, length(rowData[abd_index]))
    }
    
    
    newRow <- c(newFixedColumns, intensityChannelsData_noNA, snChannelsData)
    writeLines(paste(newRow, collapse=","), output)
    #print(line)
  }
  
  close(input)
  close(output)
  
  # read as r dataframe in memory
  DF <- read.csv(file = outputFile, stringsAsFactors = FALSE)
  return(DF)
}


