#!/usr/bin/env Rscript
if (!require("optparse")) install.packages("optparse", quiet = TRUE)
suppressWarnings(library(optparse, warn.conflicts = FALSE, quietly = TRUE,))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input csv filename", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output filename (FASTA and CSV)", metavar="character"),
  make_option(c("-k", "--kmer"), type="integer", default=9L, 
              help="kmer size", metavar="character"),
  make_option(c("-c", "--conservation"), type="character", default="CCS",  
              help="conservation level: completely conserved(CCS); both completely conserved and highly conserved(HCS) [Default: CCS]", 
              metavar="character"),
  make_option(c("-p", "--pct"), type="integer", default=NULL, 
              help="threshold for conseration level to use instead of standard HCS (90%) and CCS (100%) thresholds", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 3){
  print_help(opt_parser)
  stop("Input and output files are required. Conservation level set to CCS by default.", call.=FALSE)
}


# ---- Main script ----

# load libs
libs <- c("dplyr", "tidyr", "stringr")
for (lib in libs) {
  suppressWarnings(library(lib, character.only = TRUE,
                           warn.conflicts = FALSE, quietly = TRUE))
}


# read csv file
proteins <- read.csv(opt$input)

# HCS or CCS conservation level
conservation <- opt$conservation
# threshold HCS / CCS
if (is.null(opt$pct)) {
  threshold <- ifelse(conservation == "CCS", 100, 90)
} else {
  threshold <- opt$pct
}


# filter whole dataset by index.incidence (HCS/CCS)
df <- proteins %>% 
  filter(index.incidence >= threshold)

# stop if no peptides were found
if (nrow(df) == 0) {
  stop("No sequences with given conservative level were found", call.=FALSE)
} 

# kmer size
kmer <- opt$kmer


# ---- 1. Protein Sequences ----

# remember whole sequence of the protein
proteins_seq <- proteins %>% 
  select(proteinName, indexSequence) %>% 
  group_by(proteinName) %>% 
  summarise(seq = paste0(str_sub(indexSequence, 1, 1), collapse = "")) %>% 
  spread(key = proteinName, value = seq)

# add missing last amino acids
for (protein in unique(proteins$proteinName)) {
  proteins_seq[[protein]] <- 
    paste0(proteins_seq[[protein]], 
           proteins %>% 
           filter(proteinName == protein) %>% 
             select(indexSequence) %>% 
             dplyr::slice(n()) %>% 
             as.character() %>% str_sub(2))
}




# ---- 2. Dataset manipulations ----


# first split by protein name
# then separately proceed with each table (each protein independently)
# then rbind

csv_df <- bind_rows(
  
  lapply(split(df, df$proteinName), function(df_x) {
    
    # remember protein name
    prot_name <- df_x[1, "proteinName", drop = TRUE]
    
    # create table of all indexes falling into peptides (from filtered df)
    # cols: indexes of all amino acids
    # rows: peptides
    # value: TRUE for all amino acids from peptide on their indexes
    # example: row VKRP (from 22 to 25) - cols 22-25 will have TRUE
    index <- 
      bind_rows(
        apply(df_x, 1, function(row) {
          start_pos <- as.numeric(row["position"])
          data.frame(matrix(data = TRUE,
                            nrow = 1, ncol = kmer,
                            dimnames = list(row["indexSequence"],
                                            seq(start_pos, start_pos + kmer - 1))))
        })
        
      ) %>% 
      # if for any index there is amino acid from conservative peptide - set TRUE
      apply(2, any) %>% names() %>% as.data.frame() %>% 
      # get rid of "X" in front of indexes
      separate(col = 1, into = c("x", "index"), sep = 1) %>% 
      mutate(index = as.integer(index)) %>% 
      select(index)
    
    # set start and end indexes for each consecutive index series
    # calculated as difference between previous and next elements
    # diff = 1 means that sequence still continuing
    # diff > 1 means that there was a gap -> next sequence started
    # example: 4-5-6-10-11-12 // diff = 1-1-4-1-1
    # sequences: 4-6 (ind.1-3), 10-12 (ind.4-6)
    start_ind <- index$index[c(1, which(diff(index$index) > 1) + 1)]
    end_ind <- index$index[c(which(diff(index$index) > 1), length(index$index))]
    
    # format columns for output
    index_df <- data.frame(
      n = seq(length(start_ind)),
      start = start_ind,
      end = end_ind
    ) %>% 
      mutate(
        !!conservation := sprintf("%s_%s_%i", conservation, prot_name, n),
        Position = sprintf("%i-%i", start, end),
        Sequence = str_sub(proteins_seq[[prot_name]], start, end)
      ) %>% 
      select(-c(n, start, end))
    
  })
  
)


# create df to store info for fasta file
fasta_df <- do.call(rbind, lapply(seq(nrow(csv_df)), function(i) {
  csv_df[i, ] %>% 
    select(-Position) %>% 
    mutate(!!conservation := paste0(">", get(conservation))) %>% 
    t()
  }))



# ---- 3. Write tables ----

write.table(csv_df, file = sprintf("%s.csv", opt$output),
            sep = ",", row.names = FALSE, quote = FALSE)

write.table(fasta_df, file = sprintf("%s.fasta", opt$output),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

                    
                    
                    
                    
                    
                    
