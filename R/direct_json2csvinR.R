###################################
# Date: 29/05/2022                
# Convert DiMA JSON to CSV file   
# for input option 2 where user directly provide a JSON file      
# a separate json2csv() function is needed as the original json2csv() function directly rename the output csv file based on the json filepath       
# and this is not applicable in this case as the filepath provided via server will be in hex code (? not sure what is this called but not human readable)
# so we include another argument here which is the file name
###################################

#to extract:
#results - low_support, entropy, distinct_variants_incidence
#variants - position, count

# require data transformation
# indexPeptide - extract peptide with Index label
# multiIndex -  TRUE / FALSE
# host - user defined

#install library
#install.packages("RJSONIO")
library(RJSONIO)
library(dplyr)
library(tidyr)
library(plyr)
library(mlr3misc)
library(stringr)
#user-defined value: host
#host<-"human"

#infile<-"C:\\Users\\pendy\\Downloads\\DiMA_output_2022-05-20 (1)\\ns1_n_1.json"
#infile<-"C:\\desktop\\Research_BVU\\protocol\\Rshiny\\DiveR_myVersion_workingOnDIMA\\www\\NS3_9mer.json"
#host<-'human'
#json2csvinR(infile, "human")
#write("\r\n", file = infile, append = TRUE, sep = "\n")
direct_json2csvinR <-function(infile, host, proteinName, outfilename){
  #read JSON file
  con <- file(infile)
  jsonfile <- readLines(con)
  close(con)
  json_data <- jsonlite::fromJSON(paste(jsonfile, collapse = "\n"))
  #print(json_data)
  data_flatten <- as.data.frame(json_data) %>% 
    tidyr::unnest()
  # data transformation  
  motifs_incidence <- aggregate(data_flatten$incidence, list(data_flatten$results.position,data_flatten$motif_short), FUN=sum) %>%
    spread(Group.2,x) %>%                           # transpose the rows (motif-long & incidence) to columns (index, major, minor, unique)
    mutate_if(is.numeric, ~(replace_na(., 0)))      # replace NAN with 0
  #rename column 'Group.1' to 'results.position'
  colnames(motifs_incidence)[colnames(motifs_incidence) == 'Group.1'] <- "results.position" 
  
  #sum the number of Index motif found in each position (if > 1 => multiIndex == TRUE)
  multiIndex<-data_flatten%>%
    group_by(results.position) %>%
    dplyr::summarize(multiIndex = sum(motif_short=="I")) %>%
    as.data.frame()
  
  #merge multiIndex to motifs_incidence df
  motifs_incidence <-right_join(motifs_incidence,multiIndex, by='results.position')%>%
    distinct()
  
  #replace multiIndex with boolean: x > 1 = TRUE and vice versa
  motifs_incidence$multiIndex <- ifelse(motifs_incidence$multiIndex>1, TRUE,FALSE)
  
  #extract the first encountered index kmer info if > 1 index is encountered for each position (rarely happens)
  index_data<-subset(data_flatten, motif_short == "I") %>% #extract rows that are of index
    group_by(results.position) %>% #group them based on position
    slice(1) %>% #take the first index encountered per position
    as.data.frame() #return the data in dataframe
  
  #replace low_support of NAN to FALSE
  index_data$results.low_support[is.na(index_data$results.low_support)] <- FALSE 
  
  #combine both the index and variant motif information to motifs
  motifs<-right_join(motifs_incidence,index_data[c("sample_name","results.support","results.low_support","results.entropy","results.distinct_variants_incidence","results.position","results.total_variance","sequence")],by='results.position')%>%
    distinct()
  
  #assign host
  motifs['host'] <- host
  
  #rename columns
  colnames(motifs)<-c('position','index.incidence','major.incidence','minor.incidence','unique.incidence','multiIndex','proteinName','count','lowSupport','entropy','distinctVariant.incidence','totalVariants.incidence','indexSequence','host')
  #reorder the columns
  motifs<-motifs[,c(7,1,8,9,10,13,2,3,4,5,12,11,6,14)]
  #assign protein name
  motifs['proteinName']<-proteinName
  
  #write to csv file
  write.table(motifs, sep=",", row.names = FALSE , file = outfilename)

}

