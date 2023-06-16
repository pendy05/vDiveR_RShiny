plot_wp<-function(data,wordsize){
  colnames(data) <- c('region','count')
  world_map <- map_data("world")
  gg <- ggplot(world_map, aes(x = long, y = lat, group = group)) + geom_polygon(fill="lightgray", colour = "#888888")
  pathogens.map <- left_join(data, world_map, by = "region")
  gg <- gg + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
             scale_fill_gradient(low = "#FFFFFF", high = "#E63F00", name = 'Number of sequences') +
             theme(plot.background = element_rect(fill = "transparent", colour = NA),
                   panel.border = element_blank(), panel.grid = element_blank(),
                   axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                   legend.position = "right") +
             theme_classic(base_size = wordsize)
  gg
}

plot_tm<-function(data, wordsize, scale){

    gg <- ggplot(data, aes(Date, Total)) + geom_col(position = "stack", color='black') + 
               ylab('Number of protein sequence records') +
               scale_x_date(labels = date_format("%b-%Y"), date_breaks = "1 months") +
               theme_classic(base_size = wordsize) + 
               theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
                     axis.title.x = element_blank())
  if(scale == 'log'){
    gg <- gg + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)))
  }
  gg
}
metadataExtraction <- function(file_path, source){
  if(source == 'NCBI'){
    meta <- extract_from_NCBI(file_path)
  }
  if(source == 'GISAID'){
    meta <- extract_from_GISAID(file_path)
  }
  return(meta)
}
extract_from_NCBI <- function(file_path){
  IDs <- c()
  lines <- readLines(file_path, warn=FALSE)
  for(line in lines){
    if(grepl('>', line)){
      ID <- strsplit(line, ' ')[[1]][1]
      ID <- substr(ID,2,nchar(ID))
      IDs <- c(IDs, ID)
    }
  }
  countrys <- c(); dates <- c(); dropsample <- c(); keepsample <- c()
  for(ID in IDs){
    if(substr(ID,1,3) == 'pdb'){
      dropsample <- c(dropsample, ID)
      next
    }
    keepsample <- c(keepsample, ID)
    search_result <- entrez_search(db = "protein", term = ID, retmax = 1)
    accession <- search_result$ids[[1]]
    info <- entrez_fetch(db = "protein", id = accession, rettype = "gb", retmode = "text")
    info <- strsplit(info, '\n')[[1]]
    idx1 <- grep('/country=',  info); idx2 <- grep('/collection_date=',  info)
    info1 <- info[idx1];               info2 <- info[idx2]
    country <- strsplit(info1, '\\"')[[1]][2]
    if(grepl(':', country)){country <- strsplit(country,':')[[1]][1]}
    date <- strsplit(info2, '\\"')[[1]][2]
    countrys <- c(countrys, country); dates <- c(dates, date)
  }
  tmp <- data.frame('ID' = keepsample, 'country' = countrys, 'date' = dates)

  for(i in 1:nrow(tmp)){
    if(!is.na(as.Date(tmp$date[i],format='%Y-%m-%d'))){
      tt <- 1
    } else if(!is.na(as.Date(tmp$date[i],format='%d-%B-%Y'))){
      tmp$date[i] <- as.character(as.Date(tmp$date[i],format="%d-%B-%Y"))
    } else {
      dropsample <- c(dropsample, tmp$ID[i])
    }
  }
  wraminfo <- paste(c(dropsample, 'did not provide a clear date, so it is excluded.'), collapse = " ")
  warning(wraminfo)
  tmp <- tmp[! tmp$ID %in% dropsample, ]
  colnames(tmp) <- c("Accession_ID","Country","Date")
  return(tmp)
}
extract_from_GISAID <- function(file_path){
  heads <- c()
  lines <- readLines(file_path, warn=FALSE)
  for(line in lines){
    if(grepl('>', line)){
      line <- substr(line, 2, nchar(line))
      heads <- c(heads, line)
    }
  }
  IDs <- c(); countrys <- c(); dates <- c()
  for(head in heads){
    ID <- str_extract(head, "EPI[^|]*"); IDs <- c(IDs, ID)
    country <- strsplit(head, '/', fixed=T)[[1]][2]; countrys <- c(countrys, country)
    date <- str_extract(head, "[0-9]{4}-[0-9]{2}-[0-9]{2}"); dates <- c(dates, date)
  }
  return(data.frame('Accession_ID' = IDs, 'Country' = countrys, 'Date' = dates))
}

refineCounty <- function(metatable){
  metatable$Country[metatable$Country == "DRC"] = "Democratic Republic of the Congo"
  metatable$Country[metatable$Country == "NewCaledonia"] = "New Caledonia"
  metatable$Country[metatable$Country == "Northern Ireland"] = "New Caledonia"
  metatable$Country[metatable$Country %in% c("England","Scotland","Wales")] = "UK"
  metatable$Country[metatable$Country %in% c("Shangahi", "Xinjiang","Sichuan", "Guangdong","Shannxi", "Chongqing", "Inner_Mongolia","Shenzhen",
                                             "Fujian", "Inner Mongolia", "Tianjing", "Hebei","Jiangsu", "Shandong", "Zhejiang",
                                             "Liaoning","Shanxi", "Henan", "Chongqin", "Yunnan", "Beijing","Heilongjiang",
                                             "Hunan", "Guangxi","Ningxia","Jilin","Tibet","Hainan", "Macao", 
                                             "Jiangxi","Qinghai", "Hubei", "Gansu", "Anhui" ,"Guizhou" )] = "China"
  return(metatable)
}
