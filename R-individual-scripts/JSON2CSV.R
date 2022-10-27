###################################
# Date: 23/07/2022                #
# Convert DiMA JSON to CSV file   #
# DiMA version: 4.1.1             #
###################################

library(RJSONIO)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

#Users are required to modify the input file and host name
infile<-"json file path"
host<-'human'

write("\r\n", file = infile, append = TRUE, sep = "\n")
#read JSON file
print("new json2csv")
con <- file(infile)
jsonfile <- readLines(con)
close(con)
json_data <- jsonlite::fromJSON(paste(jsonfile, collapse = "\n"))
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
motifs<-right_join(motifs_incidence,index_data[c('query_name',"results.support","results.low_support","results.entropy","results.distinct_variants_incidence","results.position","results.total_variants_incidence","sequence",'highest_entropy.position','highest_entropy.entropy','average_entropy')],by='results.position')%>%
      distinct()
  
#assign host
motifs['host'] <- host
  
#rename columns
colnames(motifs)<-c('position','index.incidence','major.incidence','minor.incidence','unique.incidence','multiIndex','proteinName','count','lowSupport','entropy','distinctVariant.incidence','totalVariants.incidence','indexSequence','highestEntropy.position','highestEntropy','averageEntropy','host')
#reorder the columns
motifs<-motifs[,c(7,1,8,9,10,13,2,3,4,5,12,11,6,17,14,15,16)]
  
write.table(motifs, sep=",", row.names = FALSE , file = paste0(sub(".json","",infile),".csv"))
