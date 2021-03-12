#!/usr/bin/env Rscript

# List of required packages
packages <- c('httr', 'jsonlite', 'xml2', 'stringr', 'optparse', 'Biostrings')

# Check if the required packages are installed; install if not
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")

library(httr)
library(jsonlite)
library(xml2)
library(stringr)
library(optparse)
suppressMessages(library(Biostrings))

# Set up optional argument flags
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="input file (csv format)", metavar="filename"),
  make_option(c("-d", "--distance"), type="numeric", default=150,
              help="flanking sequence length", metavar="number")
);

# Parse arguments into opt object
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in the input file
SNP_list <- read.csv(paste(opt$infile), header = T, sep = ",")

# Test if the input file is missing any data
input_test <- which(is.na(SNP_list), arr.ind=TRUE)

if (length(input_test) != 0){
  message("Input file is missing data! Located at index:")
  print(input_test)
  stop()
}

# Set up the main extract_flank function
extract_flank <- function(chr, SNP_pos, distance, allele){
  
  # Specify ensembl server
  server <- "http://rest.ensembl.org"
  

  ## First server query; get the flanking sequences from Ensembl
  
  # Calculate start and end positions
  start <- SNP_pos - distance
  end <- SNP_pos + distance
  
  # Set up sequence server request params
  seq_query <- paste0("/sequence/region/horse/", chr, ":", start, "..", end, ":1?")
  
  # Execute sequence server request
  seq_r <- GET(paste(server, seq_query, sep = ""), content_type("text/plain"))
  
  # Extract sequence from result object
  stop_for_status(seq_r)
  seq <- content(seq_r)

  
  ## Second server query; get info on all SNPs within the query region
  
  # Set up variant server request params
  variant_query <- paste0("/overlap/region/horse/", chr, ":", start, "-", end, "?feature=variation")
  
  # Execute variant server request
  variant_r <- GET(paste(server, variant_query, sep = ""), content_type("application/json"))
  
  # Extract variant list from result object
  stop_for_status(variant_r)
  variant <- content(variant_r)
  
  # Extract SNP alleles within the selected region and mark with their corresponding IUPAC code
  if(length(variant) != 0){
    for (i in 1:length(variant)){
      
      # Extract alleles and merge into a single string, then convert to IUPAC code
      al <- str_c(variant[[i]]$alleles[[1]], variant[[i]]$alleles[[2]])
      IUPAC_al <- mergeIUPACLetters(al)
      
      # Substitute allele into SNP position, in IUPAC code form
      str_sub(seq, (variant[[i]][["start"]] - start + 1), (variant[[i]][["start"]] - start + 1)) <- IUPAC_al
      
    }
  }

  
  ## Third server request to get the query SNP RS ID
  
  # Set up server request params
  RS_query <- paste0("/overlap/region/horse/", chr, ":", SNP_pos, "-", SNP_pos, "?feature=variation")
  
  # Execute server request
  RS_r <- GET(paste(server, RS_query, sep = ""), content_type("application/json"))
  
  # Extract from result object
  stop_for_status(RS_r)
  RS <- content(RS_r)
  
  # If there is a known SNP at the query position, then note its RS ID
  if(length(RS) != 0){
    SNP_ID <- RS[[1]]$id
  } else {
    SNP_ID <- 'N/A'
  }
  

  ## Final section; format the output sequence with IUPAC codes marking SNPs
  
  # Define str_comp function, which checks if two strings contain the same letters disregarding character order
  # Disappointing that stringr doesn't contain this function
  str_comp <- function(a, b){
    a <- str_sort(strsplit(a, "")[[1]])
    b <- str_sort(strsplit(b, "")[[1]])
    identical(a, b)
  }
  
  # Mark query position with an S
  str_sub(seq, distance + 1, distance + 1) <- "S"
  
  # Get supplied allele in desired format
  if (grepl("/", allele) == TRUE){
    allele <- str_remove(allele, "/")
  }
  
  # If there is a known SNP at the query position, then insert it into the sequence in IUPAC code format
  if(length(RS) != 0){
    
    # Get the alleles at the query position according to Ensembl
    RS_allele <- str_c(RS[[1]]$alleles[[1]], RS[[1]]$alleles[[2]])
    
    # Check if supplied alleles and Ensembl alleles at query position match
    if (str_comp(allele, RS_allele) == FALSE){
      message(paste0("Supplied allele does not match Ensembl allele! Using Ensembl allele ", RS[[1]]$alleles[[1]], "/", RS[[1]]$alleles[[2]], " instead."))
    }
    
    # Sub in the alleles to the query SNP position, in IUPAC form
    seq <- gsub("S", mergeIUPACLetters(RS_allele), seq)
    
  } else {
    message("Important: there is no SNP documented on the Ensembl database at the specified position.")
    seq <- gsub("S", mergeIUPACLetters(allele), seq)
  }
  
  # Set up final results object
  results <- cbind(seq, SNP_ID)

}


## This section runs the extract_flank function on supplied SNPs in the input file

# Open counter
counter = 0

# Loop over input file, execute the function for each row
for (i in 1:nrow(SNP_list)){
  counter = counter + 1
  message(paste0('Getting flanking sequences for SNP No. ', counter, '...'))
  SNP_list[i,4:5] <- extract_flank(chr = SNP_list[i,1], SNP_pos = SNP_list[i,2], distance = 150, allele = SNP_list[i,3])
  message('Done.')
}

# Write output file
write.csv(SNP_list, paste0("output_", opt$infile ))

message("All SNP flanking sequences extracted!")
message(paste0("Output file is: output_", opt$infile))
