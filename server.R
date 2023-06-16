#set maximum upload file size
options(shiny.maxRequestSize=3000*1024^2)
#https://www.shinyapps.io/admin/#/application/5536925/logs
library(ggplot2)
library(gridExtra) #tutorial: https://ggplot2.tidyverse.org/reference/facet_grid.html 
library(plyr)
library(dplyr)
library(tidyr)
library(ggtext)
library(ggpubr)
library(grid)
library(cowplot)
library(RJSONIO)
library(stringr)
library(seqinr)
library(purrr)
library("tools")  
library("gghalves")
library(rjson)
library(readr)
library(zip)
library(shinyjs)
library(jsonlite)
library(glue)
library(shinyThings) # devtools::install_github("gadenbuie/shinyThings")
library(maptools)
library(lubridate)
library(scales)
library(rentrez)
#library(NGLVieweR)
#library(bio3d)
#library(Biostrings)
#library(NGLVieweR)
#library(ggmsa) #devtools::install_github("YuLab-SMU/ggmsa")
               #devtools::install_github("hrbrmstr/ggalt", ref = "noproj")  
#reticulate::use_virtualenv("python_env", required = TRUE)

#individual functions
#reset input and output settings
resetInput_to_initialState <-function(output){
  shinyjs::reset("proteinOrder")
  shinyjs::reset("line_dot_size")
  shinyjs::reset("wordsize")
  shinyjs::reset("supportLimit")
  shinyjs::reset("kmerlength")
  shinyjs::reset("filetype")
  shinyjs::reset("host")
  shinyjs::reset("proteinNames")
  shinyjs::reset("hostname")
  shinyjs::reset("proteinNames_secondHost")
  shinyjs::reset("hostname_secondHost")
  shinyjs::reset("MSAfile")
  shinyjs::reset("MSAfile_secondHost")
  shinyjs::reset("Metafile")
  shinyjs::reset("Metafasta")
  shinyjs::reset("Meta")
  shinyjs::reset("inmetafilename")
  shinyjs::reset("inmetafasta")
  shinyjs::disable(id="downloadDiMA")

  #clear output
  output$alert <- renderUI({})
  output$alertSample <- renderUI({})
  #output$protein_selection <- renderUI({})
  output$plot_worldmap<- renderPlot({})
  output$countrytable<- renderDataTable({})
  output$plot_time<- renderPlot({})
  output$timetable<- renderDataTable({})
  output$plot1<- renderPlot({})
  output$plotEntropy<- renderPlot({})
  output$plot2<- renderPlot({})
  output$plot3<- renderPlot({})
  output$plot4<- renderPlot({})
  output$plot7<- renderPlot({})
  output$entropyTable<-renderDataTable({})
  output$plot7_seqs<-renderDataTable({})
}

#download sample data
downloadSampleData<-function(output){
  output$downloadSampleData <- downloadHandler(
    filename <- function() {
        paste("vDiveR_sample_input_", Sys.Date(), ".zip", sep = "")
    },
    
    content <- function(file) {
      file.copy("www/vDiveR_sample_input.zip", file)
    },
    contentType = "application/zip"
  # filename = function() {
  #   paste("vDiveR_sample_input_", Sys.Date(), ".zip", sep = "")
  # },
  # content = function(file) { 
  #   print('in temp dir')
  #   #To create temporary directory
  #   temp_directory_sample <- file.path(tempdir(), as.integer(Sys.time()))
  #   dir.create(temp_directory_sample)
  #   
  #   csv_dataset <- read.csv("www/DiMA_HCV.csv")
  #   MSA_dataset_Core <- seqinr::read.fasta("www/Core_mafft.fasta")
  #   csv_dataset_Core <- read.csv("www/core_9mer.csv")
  #   JSON_dataset_Core<- rjson::fromJSON(file = "www/core_9mer.json")
  #   MSA_dataset_NS3 <- seqinr::read.fasta("www/NS3_mafft.fasta")
  #   csv_dataset_NS3 <- read.csv("www/NS3_9mer.csv")
  #   JSON_dataset_NS3<- rjson::fromJSON(file = "www/NS3_9mer.json")
  #   csv_dataset_metadata <- read.csv("www/oneID_GISAID.csv")
  # 
  #   #write sample dataset in CSV, FA & JSON formats into temp directory
  #   write.csv(csv_dataset,paste0(temp_directory_sample,"/HCV_DiMA.csv"))
  #   seqinr::write.fasta(MSA_dataset_Core,names=names(MSA_dataset_Core),file.out = paste0(temp_directory_sample,"/HCV_aligned_Core.fasta"))
  #   write.csv(csv_dataset_Core,paste0(temp_directory_sample,"/HCV_Core.csv"))
  #   write_json(JSON_dataset_Core,paste0(temp_directory_sample,"/HCV_Core.json"))
  #   seqinr::write.fasta(MSA_dataset_NS3,names=names(MSA_dataset_NS3),file.out = paste0(temp_directory_sample,"/HCV_aligned_NS3.fasta"))
  #   write.csv(csv_dataset_NS3,paste0(temp_directory_sample,"/HCV_NS3.csv"))
  #   write_json(JSON_dataset_NS3,paste0(temp_directory_sample,"/HCV_NS3.json"))
  #   write.csv(csv_dataset_metadata,paste0(temp_directory_sample,"/metadata.csv"))
  #   zip::zip(zipfile = file,files = dir(temp_directory_sample), root = temp_directory_sample)
  # },
  # contentType = "application/zip"
)
}

#plot generators
generate_worldmap <-function(input, output, plot_worldmap){
  output$plot_worldmap <- renderPlot({
    plot_worldmap()  
  })
  
  output$plot_worldmap_download <- downloadHandler(
    filename = function() { paste("plot_world_map", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot_worldmap(), width=input$width_wm, height=input$height_wm,unit="in", device = "jpg", dpi=input$dpi_wm)
    }
  )
  
}

generate_timeplot <-function(input, output, plot_time){
  output$plot_time <- renderPlot({
    plot_time()  
  })
  
  output$plot_time_download <- downloadHandler(
    filename = function() { paste("plot_time", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot_time(), width=input$width_tm, height=input$height_tm,unit="in", device = "jpg", dpi=input$dpi_tm)
    }
  )
  
}

generate_plot1<-function(input, output, plot1){
  output$plot1 <- renderPlot({
    plot1()  
  })
  
  output$plot1_download <- downloadHandler(
    filename = function() { paste("plot_entropy_incidence", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot1(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
    }
  )
  
}

generate_plotEntropy<-function(input, output, plotEntropy){
  output$plotEntropy <- renderPlot({
    plotEntropy()  
  })
  
  output$plotEntropy_download <- downloadHandler(
    filename = function() { paste("plot_entropy", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plotEntropy(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
    }
  )
  
}

generate_plot2<-function(input, output, plot2){
  output$plot2<- renderPlot({
    plot2()
  })
  
  output$plot2_download <- downloadHandler(
    filename = function() { paste("plot_relationship_entropy_total_variants", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi2)
    }
  )
  
  output$info_plot2 <- renderText({
    HTML(paste0("Total variants (%) = ", input$plot2_click$x, "<br>",em("k"),"-mer entropy (bits) = ",input$plot2_click$y))
  })
}

generate_plot3<-function(input, output, plot3){
  output$plot3<- renderPlot({
    plot3()
  })
  
  output$plot3_download <- downloadHandler(
    filename = function() { paste("plot_dynamics_diversity_motifs_proteome", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot3(), width=input$width3, height=input$height3,unit="in", device = "jpg", dpi=input$dpi3,bg='white')
    }
  )
}

generate_plot4<-function(input, output, plot4){
  output$plot4<- renderPlot({
    plot4()
  })
  
  output$plot4_download <- downloadHandler(
    filename = function() { paste("plot_dynamics_diversity_motifs_proteins", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot4(), width=input$width4, height=input$height4,unit="in", device = "jpg", dpi=input$dpi4, bg='white')
    }
  )
}

generate_plot7<-function(input, output, plot7){
  output$plot7<- renderPlot({
    plot7()
  })
  
  output$plot7_download <- downloadHandler(
    filename = function() { paste("plot_conservationLevels_protein", '.jpg', sep='') },
    content = function(file) {
      ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
  })
}


generate_entropyTable<-function(data, output, proteinName){
  #get position of min entropy, min, max of entropy and total variant
  entropyTable <- data %>% 
    dplyr::group_by(proteinName) %>%
    dplyr::summarise(
      Position = gsub("((?:\\d+,){2}\\d+),", "\\1,\n", paste0(position[which(entropy == min(entropy))], collapse = ",")),
      minEntropy = format(round(min(entropy),digits=2),nsmall=2),
      maxEntropy = format(round(max(entropy), digits = 2),nsmall=2),
      minTotalVariants = format(round(min(totalVariants.incidence),digits=2),nsmall=2),
      maxTotalVariants = format(round(max(totalVariants.incidence), digits = 2),nsmall=2)
    )
  #rename table df
  names(entropyTable)<- c("Protein Name","Position (Minimum Entropy)","Minimum Entropy","Maximum Entropy","Minimum Total Variants (%)","Maximum Total Variants (%)")
  
  output$entropyTable <- renderDataTable(
    entropyTable,
    width = "100%",
    options = list(scrollX=TRUE, scrollCollapse=TRUE)
  )
}

generate_CCS_HCS_table<-function(input, output, data){
  # for now not splitted by hosts
  output$plot7_seqs <- renderDataTable({
     seqConcatenation(input_file=data.frame(data), kmer=input$kmerlength, 
                     threshold_pct = as.numeric(input$conserv_percent),
                     conservation=input$conserv_lvl)[[input$table_type]]
  })
  
  output$conservSeq_download <- downloadHandler(
    filename =  function() {paste0(input$conserv_lvl, ".", input$table_type)},
    content = function(fname) {
      df <- seqConcatenation(input_file=data.frame(data), kmer=input$kmerlength, 
                             threshold_pct = as.numeric(input$conserv_percent),
                             conservation=input$conserv_lvl)[[input$table_type]]
      write.table(df, file = fname, col.names = ifelse(input$table_type == "csv", TRUE, FALSE), sep = ",", row.names = FALSE, quote = FALSE)
    }
  )
}

#server main function
server <- function(input, output,session) {
  # initial state of downloadDiMA button is disabled
  shinyjs::disable(id="downloadDiMA")
  shinyjs::disable(selector = '.nav-tabs a[data-value="Second Host"')

  #=====================================================#
  #                Sample Input Format                  #
  #          display a table of sample input            #
  #=====================================================#
  output$mainDataSample <- DT::renderDT({
    mainDataSample <- read.csv("www/DiMA_HCV.csv", header = T, stringsAsFactors = F)
    mainDataSample
  })

  
  #if "two hosts" is selected, allow user to input data at the second tab
  observeEvent(input$host, {
    if(input$host == 2){
      shinyjs::enable(selector = '.nav-tabs a[data-value="Second Host"')
    }else{
      shinyjs::disable(selector = '.nav-tabs a[data-value="Second Host"')
    }
  })
  
  #reset every input and output to initial state
  observeEvent(input$reset, {
    resetInput_to_initialState(output)
  })
  
  #reset every input and output to initial state
  observeEvent(input$reset1, {
    resetInput_to_initialState(output)
  })
  
  #redirect user to input tab when they click on start button
  #https://stackoverflow.com/questions/32971921/navigate-to-particular-sidebar-menu-item-in-shinydashboard
  observeEvent(input$start, {
    newtab <- switch(input$tabs,
                "description" = "inputdata_description",
                "MetaData" = "inputdata_description",
                "plotEntropy" = "inputdata_description",
                "plot1" = "inputdata_description",
                "plot2" = "inputdata_description",
                "plot3" = "inputdata_description",
                "plot4" = "inputdata_description",
                "plot7" = "inputdata_description",
                "helppage" = "inputdata_description"
    )
    updateTabItems(session, "tabs", newtab)
  })
  
  #------------------------------------------------------------------------------#
  #                   Print Input File Names to DiveR's UI                       #
  #   -To let users input protein names followed the sequence of input files     #
  #------------------------------------------------------------------------------#
  output$infilename<-renderText({
    if (is.null(input$MSAfile)){
      return(NULL)
    }else{
      paste("File(s):",toString(input$MSAfile$name),sep=" ")
    }
  })
  
  output$infilename_secondHost<-renderText({
    if (is.null(input$MSAfile_secondHost)){
      return(NULL)
    }else{
      paste("File(s):",toString(input$MSAfile_secondHost$name),sep=" ")
    }
  })
  #=====================================================#
  #                   MetaData                          #
  #         show the World map and time plot            #
  #=====================================================#
  output$metademo <- DT::renderDT({
    demoMeta <- read.csv("www/oneID_GISAID.csv", header = T, stringsAsFactors = F)
    demoMeta
  })
  output$inmetafilename <- renderText({
    if (is.null(input$Metafile)){
      return(NULL)
    }else{
      paste("File:",toString(input$Metafile$name),sep=" ")
    }
  })
  output$inmetafasta <- renderText({
    if (is.null(input$Metafasta)){
      return(NULL)
    }else{
      paste("File:",toString(input$Metafasta$name),sep=" ")
    }
  })
  observeEvent(input$submitMeta1, {
    shinyjs::addClass(id = "submitmeta1", class = "loading dots")
    Meta <- reactive({
      req(input$Metafile)
      filepath <- input$Metafile$datapath
      Meta <- read.csv(filepath, header = T, stringsAsFactors = F)
      Meta$Country[Meta$Country == "DRC"] = "Democratic Republic of the Congo"
      Meta$Country[Meta$Country == "NewCaledonia"] = "New Caledonia"
      Meta$Country[Meta$Country == "Northern Ireland"] = "New Caledonia"
      Meta$Country[Meta$Country %in% c("England","Scotland","Wales")] = "UK"
      Meta
    })
    WorldmapInFo <- reactive({
      req(Meta())
      meta <- Meta()
      countrylist <- meta$Country
      countrylist <- data.frame(table(countrylist))
      colnames(countrylist) <- c('Country','Number of sequences')
      countrylist
    })
    TimeInFo <- reactive({
      req(Meta())
      temporal <- Meta()
      temporal$count <- 1
      temporal$Date <- as.Date(temporal$Date, format = "%d/%m/%Y")
      temporal <- aggregate(temporal$count, by=list(temporal$Date), sum)
      colnames(temporal) <- c('Date', 'Total')
      temporal
    })
    plot_worldmap <- reactive({
      req(WorldmapInFo())
      plot_wp(WorldmapInFo(), input$wordsize)
    })
    generate_worldmap(input,output,plot_worldmap)
    output$countrytable = DT::renderDataTable({
      req(WorldmapInFo())
      WorldmapInFo()
    })
    output$table_worldmap_download <- downloadHandler(
      filename = function() {paste("CountryInfo", '.csv', sep='')},
      content = function(file) {write.csv(WorldmapInFo(), file, quote = F)}
    )
    plot_time <- reactive({
      req(TimeInFo())
      plot_tm(TimeInFo(), input$wordsize, input$time_scale)
    })
    generate_timeplot(input,output,plot_time)
    output$timetable = DT::renderDataTable({
      req(TimeInFo())
      TimeInFo()
    })
    output$table_time_download <- downloadHandler(
      filename = function() {paste("TimeInfo", '.csv', sep='')},
      content = function(file) {write.csv(TimeInFo(), file, quote = F)}
    )
    shinyjs::removeClass(id = "submitmeta1", class = "loading dots")
  })
  
  observeEvent(input$submitMeta2, {
    shinyjs::addClass(id = "submitmeta2", class = "loading dots")
    Meta <- reactive({
      req(input$Metafasta)
      filepath <- input$Metafasta$datapath
      Meta <- metadataExtraction(filepath, input$MetafastaSource)
      Meta <- refineCounty(Meta)
      Meta
    })
    output$metademoSee <- DT::renderDT({
      req(input$Metafasta)
      Meta()
    })
    WorldmapInFo <- reactive({
      req(Meta())
      meta <- Meta()
      countrylist <- meta$Country
      countrylist <- data.frame(table(countrylist))
      colnames(countrylist) <- c('Country','Number of sequences')
      countrylist
    })
    TimeInFo <- reactive({
      req(Meta())
      temporal <- Meta()
      temporal$count <- 1
      temporal$Date <- as.Date(temporal$Date)
      temporal <- aggregate(temporal$count, by=list(temporal$Date), sum)
      colnames(temporal) <- c('time', 'sum_count')
      temporal
    })
    plot_worldmap <- reactive({
      req(WorldmapInFo())
      plot_wp(WorldmapInFo(), input$wordsize)
    })
    generate_worldmap(input,output,plot_worldmap)
    output$countrytable = DT::renderDataTable({
      req(WorldmapInFo())
      WorldmapInFo()
    })
    output$table_worldmap_download <- downloadHandler(
      filename = function() {paste("CountryInfo", '.csv', sep='')},
      content = function(file) {write.csv(WorldmapInFo(), file, quote = F)}
    )
    plot_time <- reactive({
      req(TimeInFo())
      plot_tm(TimeInFo(), input$wordsize, input$time_scale)
    })
    generate_timeplot(input,output,plot_time)
    output$timetable = DT::renderDataTable({
      req(TimeInFo())
      TimeInFo()
    })
    output$table_time_download <- downloadHandler(
      filename = function() {paste("TimeInfo", '.csv', sep='')},
      content = function(file) {write.csv(TimeInFo(), file, quote = F)}
    )
    shinyjs::removeClass(id = "submitmeta2", class = "loading dots")
  })
  #To create directory
  temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
  dir.create(temp_directory)
  #Reactive value to store the path of all the files in the directory
  mylist <- reactiveValues()
  mylist$files <- NULL
  
  #----------------------------------------------------#
  #     Trigger the Diver Run (DiMA & Plotting)        #
  #         -By Clicking on Submit Button              #
  #----------------------------------------------------#
  observeEvent(input$submitDiMA, ignoreInit=TRUE,{
    req(input$MSAfile)
    shinyjs::addClass(id = "submitAnimate", class = "loading dots")
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize)
    )
    #proceed if the input csv file is provided
    if (is.null(input$MSAfile)){
      return(NULL)
    }
    #default host name: Unknown Host 1
    if (input$hostname == ""){
      hostname <-"Unknown Host 1"
    }else{
      hostname<-input$hostname
    }
    
    #default host name: Unknown Host 2
    if (input$hostname_secondHost == ""){
      hostname_secondHost <-"Unknown Host 2"
    }else{
      hostname_secondHost<-input$hostname_secondHost
    }
    
    #-------------------------------------------------------------------#
    #          Assign protein name to each input file (using lapply)    #
    #   -Depends on the order of protein names provided by users        #
    #-------------------------------------------------------------------#
    
    #get filepath(s)
    filepath <- input$MSAfile$datapath
    #get protein name(s)
    if (input$proteinNames == ""){
      proteinName <- list(rep("Unknown",each = length(filepath)))
    }else{
      proteinName <- unlist(lapply(strsplit(input$proteinNames,','),trimws))
         
      if (length(proteinName) < length(filepath)){ #if number of protein names != number of files, then assign NA to that protein
        proteinName <- append(proteinName,rep("Unknown",each = length(filepath)-length(proteinName)))
      }
      # }else if (length(proteinName) > length(filepath)){
      #   #do something
      # }
    }
    
    #---------------------#
    #     Second Host     #
    #---------------------#
    
    filepath_secondHost <- input$MSAfile_secondHost$datapath

    if (input$proteinNames_secondHost != ""){
      proteinName_secondHost <- unlist(lapply(strsplit(input$proteinNames_secondHost,','),trimws))
      
      if (length(proteinName_secondHost) < length(filepath_secondHost)){ #if number of protein names != number of files, then assign NA to that protein
        proteinName_secondHost <- append(proteinName_secondHost,rep("Unknown",each = length(filepath_secondHost)-length(proteinName_secondHost)))
      }# }else if (length(proteinName_secondHost) > length(filepath_secondHost)){
      #   #do something
      # }
    }
    

    #--------------------------------------------------------------#
    #                Check the Input Type of Files                 #
    #    1. MSA Files - Run DiMA                                   #
    #    2. DiMA JSON Output Files - convert to CSV files          #
    #    3. DiMA CSV-converted Output Files                        #
    #  All 1,2,3 data in CSV files will be stored in 1 final CSV   #
    #  file that act as the source of subsequent plotting          #
    #--------------------------------------------------------------#
    #detect the input data file type
    if(input$filetype == 1){ #if input is MSA
      outfile<-""
      csvfilelist<-c()
      if (input$inputtype == 1){
          inputtype <- "protein"
      }else{
          inputtype <- "nucleotide"
      }
      #run DiMA
      for (i in 1:length(filepath)){        
        outfile<- paste0(proteinName[[i]][1],"_",hostname,"_",i,".json",sep="")

        #print('before dima-cli', paste0(temp_directory, "/",outfile))
        system(paste("python_env/Scripts/dima-cli.exe -i", filepath[i], "-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-q",proteinName[[i]][1], "-l",input$kmerlength, "-a",inputtype))

        print('after dima-cli')
        #https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
        #write("\r\n", file = paste0("./",outfile), append = TRUE, sep = "\n")
        print('before unnest')
        json2csvinR_unnest(paste0(temp_directory, "/",outfile),hostname, proteinName[[i]][1])
        print('after unnest')

        #store the DiMA csv output names into a list (for further concatenation into one file)
        csvfile<-paste0(temp_directory, "/",proteinName[[i]][1],"_",hostname,"_",i,".csv",sep="")
        csvfilelist <- append(csvfilelist, csvfile)
      }
      
      #------------------------------------------------------------------------------------------------------#
      #  Placing DiMA JSON and CSV output files into a zipped folder for download purpose                    #
      #  1. Store all CSV output in a final csv file named "DiMA.csv" for further plotting                   #
      #  2. Place the DiMA JSON and CSV output as well as "DiMA.csv" into "DiMAoutput/" folder and zip it    #
      #------------------------------------------------------------------------------------------------------#
      
      data <-  read.csv(csvfilelist[1])
      #Saving the first DiMA csv output file in the filepath array in the temp zipped folder
      
      #reading each file within the range and append them to create one file
      for (f in csvfilelist[-1]){
        df <- read.csv(f)      # read the file
        #Saving the following DiMA csv output file in the filepath array in the temp zipped folder
        data <- rbind(data, df)    # append the current file
      }
      
      
      #---------------------#
      #     Second Host     #
      #---------------------#
      if(input$host == 2){
        outfile_secondHost<-""
        csvfilelist_secondHost<-c()
        #run DiMA
        for (i in 1:length(filepath_secondHost)){
          outfile<- paste0(proteinName_secondHost[[i]][1],"_",hostname_secondHost,"_",i,".json",sep="")
          system(paste("python_env/Scripts/dima-cli.exe -i", filepath_secondHost[i], "-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-q",proteinName_secondHost[[i]][1], "-l",input$kmerlength, "-a",inputtype))
          #py_run_string("from dima import Dima")
          #dima_input<- paste0("sequences=r'",filepath_secondHost[i],"',kmer_length=",input$kmerlength,",support_threshold=",input$supportLimit,",query_name=r'",proteinName_secondHost[[i]][1],"',alphabet=r'",inputtype,"'")
          #dima_input<- gsub("/", "\\", dima_input)
          #py_run_string(glue("results = Dima({dima_input}).run()"))
          #py_run_string(glue("results = Dima(sequences=\"{filepath[i]}\", kmer_length={input$kmerlength},support_threshold={input$supportLimit}, query_name=\"{proteinName[[i]][1]}\").run()"))
          #py_run_string(glue("jsonFile = open(r'{temp_directory}/{outfile}', 'w')"))
          #py_run_string(glue("jsonFile.write(str(results))"))
          #py_run_string("jsonFile.close()")
          #system(paste("dima-cli.exe -i ",filepath_secondHost[i]," -o ",paste0(temp_directory, "/",outfile)," -s ",input$supportLimit, " -p ",proteinName_secondHost[[i]][1], " -l ",input$kmerlength))
          
          #system(paste("./python_env/Scripts/dima-cli.exe -i",filepath_secondHost[i],"-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-p",proteinName_secondHost[[i]][1], "-l",input$kmerlength))
          #system2(command="./python_env/Scripts/dima-cli.exe",args = c("-i",filepath_secondHost[i],"-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-p",proteinName_secondHost[[i]][1], "-l",input$kmerlength))
          #system2(command="./python_env/Scripts/dima-cli.exe",args = c("-i",filepath[i],"-o",outfile,"-s",input$supportLimit, "-p",proteinName[[i]][1], "-l",input$kmerlength))
          #print(input$supportLimit, proteinName[[i]][1], input$kmerlength)
          #append "\n" to the end of file to solve the issue of 'incomplete final line'
          #https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
          write("\r\n", file = outfile, append = TRUE, sep = "\n")
          json2csvinR_unnest(paste0(temp_directory, "/",outfile),hostname_secondHost,proteinName_secondHost[[i]][1])
          #json2csvinR(paste0(temp_directory, "/",outfile),hostname)
          #json2csvinR(outfile,hostname)
          #store the DiMA csv output names into a list (for further concatenation into one file)
          csvfile<-paste0(temp_directory, "/",proteinName_secondHost[[i]][1],"_",hostname_secondHost,"_",i,".csv",sep="")
          #csvfile<-paste0(proteinName[[i]][1],"_",i,".csv",sep="")
          
          csvfilelist_secondHost <- append(csvfilelist_secondHost, csvfile)
          #outfile<- paste0(proteinName[[i]][1],"_",i,".json",sep="")
        }
        
        data_secondHost <-  read.csv(csvfilelist_secondHost[1])
        #reading each file within the range and append them to create one file
        for (f in csvfilelist_secondHost[-1]){
          df <- read.csv(f)      # read the file
          #Saving the following DiMA csv output file in the filepath array in the temp zipped folder
          data_secondHost <- rbind(data_secondHost, df)    # append the current file
        }
        
        data<- rbind(data,data_secondHost)
      }
      
      #write to a final csv file "DiMAoutput.csv", consists of all the submitted proteins in the temp zipped folder
      write.table(data, sep=",", row.names = FALSE , file = paste0(temp_directory,"/DiMA.csv"))
      #Store all the path of the files in the directory in the reactive value
      mylist$files <- list.files(temp_directory,"*.*")
      
    }else if (input$filetype == 2){ #if the data is DiMA json output
      csvfilelist<-c()
      #convert DiMA output from JSON to CSV
      for (i in 1:length(filepath)){
        #----------------------NOTE (2/5/2022)-------------------#
        #"OUTFILE" SHOULD directly be the name of the input files user provided
        #"proteinName" are needed?
        csvfile<-paste0(strsplit(input$MSAfile$name[i], ".json")[[1]][1],".csv",sep="")
        direct_json2csvinR(filepath[i],hostname,proteinName[[i]][1], paste0(temp_directory,"/",csvfile))

        #store the DiMA csv output names into a list (for further concatenation into one file)
        csvfilelist <- append(csvfilelist, paste0(temp_directory,"/",csvfile))
      }
      data <-  read.csv(csvfilelist[1])
      
      #------------------------------------------------------------------------------------------------------#
      #  Placing DiMA JSON and CSV output files into a zipped folder for download purpose                    #
      #  1. Store all CSV output in a final csv file named "DiMA.csv" for further plotting                   #
      #  2. Place the DiMA JSON and CSV output as well as "DiMA.csv" into "DiMAoutput/" folder and zip it    #
      #------------------------------------------------------------------------------------------------------#
      
      #reading each file within the range and append them to create one file
      for (f in csvfilelist[-1]){
        df <- read.csv(f)      # read the file
        data <- rbind(data, df)    # append the current file
      }
      
      #---------------------#
      #     Second Host     #
      #---------------------#
      if(input$host == 2){
        csvfilelist_secondHost<-c()
        #run DiMA
        for (i in 1:length(filepath_secondHost)){
          csvfile<-paste0(strsplit(input$MSAfile_secondHost$name[i], ".json")[[1]][1],".csv",sep="")
          direct_json2csvinR(filepath[i],hostname_secondHost,proteinName_secondHost[[i]][1], paste0(temp_directory,"/",csvfile))
          #store the DiMA csv output names into a list (for further concatenation into one file)
          csvfilelist_secondHost <- append(csvfilelist_secondHost, paste0(temp_directory,"/",csvfile))
        }
        
        data_secondHost <-  read.csv(csvfilelist_secondHost[1])
        
        #reading each file within the range and append them to create one file
        for (f in csvfilelist_secondHost[-1]){
          df <- read.csv(f)      # read the file
          data_secondHost <- rbind(data_secondHost, df)    # append the current file
        }
        
        data<- rbind(data,data_secondHost)
      }
      
      #write to a final csv file "DiMAoutput.csv", consists of all the submitted proteins in the temp zipped folder
      write.table(data, sep=",", row.names = FALSE , file = paste0(temp_directory,"/DiMA.csv"))
      
    }else if (input$filetype == 3){ #if the data is DiMA csv output
      #read the input file
      if(input$host == 1){ #1 host
        data <- filepath %>%
          lapply(read_csv) %>%
          bind_rows
      }else{ #2 hosts
        filepath<-append(filepath,filepath_secondHost)
        print('filepath')
        print(filepath)
        data <- filepath %>%
          lapply(read_csv) %>%
          bind_rows
      }
      print('csv data')
      print(data)
    }
    
    #Alert results are ready
    output$alert <- renderUI({
      h5("DiMA Output is ready! Click on other tabs to visualize the diversity dynamics of viral sequences!", style = "color:green") 
    })
    shinyjs::enable(id="downloadDiMA")
    
    df<-data
    
    #determine number of host
    if (input$host == 1){
      #single host
      if (!"host" %in% colnames(df)){
        print("Host column is not detected!")
      }
      #absent of user input for protein order; follow the protein order provided in csv file
      if (input$proteinOrder ==""){
        #get protein name and number of positions
        a<-table(df$proteinName)
        proteinName<-as.vector(names(a))
        position<-as.vector(a)
        #store protein as factor
        df$size_f = factor(df$proteinName,levels = proteinName) #ALERT! df$level -> df$size_f
        #scale the nonamer positions
        scales_x<-mapply(function(x,y){
          x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
        }, proteinName,position)
        
      }else{
        if (!"host" %in% colnames(df)){
          print("Host column is not detected!")
        }
        #order the proteins based on user input
        #level<-strsplit(input$proteinOrder, ',')[[1]]
        level<-unlist(lapply(strsplit(input$proteinOrder,','),trimws))
        position<-c()
        for (i in level){
          position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
        }
        #store protein as factor
        df$size_f = factor(df$proteinName,levels = level) #ALERT! df$level -> df$size_f
        #scale the nonamer positions
        scales_x<-mapply(function(x,y){
          x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
        }, level,position)
      }
      
    }else{ # multihost
      if (!"host" %in% colnames(df)){
        print("Host column is not detected!")
      }
      #categorise data based on host
      df$host<- factor(df$host)
      print('host')
      print(unique(df$host))
      #count the aa length for each proteins (each host is expected to have same number of proteins with same length)
      df_sub<-df[df$host==unique(df$host[1]),]
      
      #absent of user input for protein order; follow the protein order provided in csv file
      if (input$proteinOrder ==""){
        #get protein name and number of positions
        a<-table(df_sub$proteinName)
        proteinName<-as.vector(names(a))
        position<-as.vector(a)
        #store protein as factor
        df$size_f = factor(df$proteinName,levels = proteinName)
        #scale the nonamer positions
        scales_x<-mapply(function(x,y){
          x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
        }, proteinName,position)
        
      }else{
        
        #order the proteins based on user input
        level<-unlist(lapply(strsplit(input$proteinOrder,','),trimws))
        position<-c()
        for (i in level){
          position<-append(position,table(df_sub$proteinName)[names(table(df_sub$proteinName)) == i])
        }
        #store protein as factor
        df$size_f = factor(df$proteinName,levels = level)
        #scale the nonamer positions
        scales_x<-mapply(function(x,y){
          x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
        }, level,position)
        
      }
    }

    group_names<-c("Index","Major","Minor","Unique","Total variants","Nonatypes")    
    
    #--------------------Table Output----------------------#
    generate_entropyTable(data, output, proteinName)
    
    #----------------------Plotting-----------------------#
    
    #Tab 1: entropy incidence plot
    plot1<-reactive({
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder, input$kmerlength)
    })

    generate_plot1(input,output,plot1)


    #Additional plot: Entropy plot
    plotEntropy<-reactive({
      plot_entropy(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder, input$kmerlength)
    })

    generate_plotEntropy(input, output, plotEntropy)

    #Tab 2: correlation plot
    plot2<-reactive({
      plot_correlation(df,input$line_dot_size,input$wordsize,input$host)
    })

    generate_plot2(input, output, plot2)
    
    plot3<-reactive({
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot3<- renderPlot({
      plot3()
    })
    
    output$plot3_download <- downloadHandler(
      filename = function() { paste("plot_dynamics_proteome", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot3(), width=input$width3, height=input$height3,unit="in", device = "jpg", dpi=input$dpi3)
      }
    )
    
    
    plot4<-reactive({
      plot_dynamics_protein(data,input$line_dot_size,input$wordsize,input$host,input$proteinOrder)
    })
    
    output$plot4<- renderPlot({
      plot4()
    })
    
    output$plot4_download <- downloadHandler(
      filename = function() { paste("plot_dynamics_proteins", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot4(), width=input$width4, height=input$height4,unit="in", device = "jpg", dpi=input$dpi4)
      }
    )
    
    plot7<-reactive({
      data<- data.frame(data)
      data<-data%>%mutate(ConservationLevel = case_when(
        data$index.incidence == 100 ~ "Completely conserved (CC)",
        data$index.incidence >= 90 ~ "Highly conserved (HC)",
        data$index.incidence >= 20 ~ "Mixed variable (MV)",
        data$index.incidence >= 10  ~ "Highly diverse (HD)",
        data$index.incidence < 10 ~ "Extremely diverse (ED)"
        
      ))
      if (input$host == 1){
        #single host
        plot_conservationLevel(data,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)
      }else{#multihost
        data$host = factor(data$host)
        #split the data into multiple subsets (if multiple hosts detected)
        plot7_list<-split(data,data$host)
        plot7_multihost<-lapply(plot7_list,plot_conservationLevel,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)
        
        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$host))))
        
      }
      
    })

    generate_plot7(input, output, plot7)
    output$conserv_threshold_box <- renderUI({
      textInput(inputId="conserv_percent", label=NULL, 
                value = ifelse(input$conserv_lvl == "CCS", 100, 90))
    })
    generate_CCS_HCS_table(input, output, data)
    
    output$downloadDiMA <- downloadHandler(
      filename = function(){
        paste("DiMA_output_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        zip::zip(zipfile = file, 
                 files = dir(temp_directory),
                 root = temp_directory)
      },
      
      contentType = "application/zip"
    )
    shinyjs::removeClass(id = "submitAnimate", class = "loading dots")
    
  })
  

  output$plot3_hosts<-renderUI({
    if (input$host==1){
      fluidRow(
        box(width=6,
            plotOutput("plot3", height = 650)),
        box(id="container",style = "background-color:#ecf0f5;",width=4,solidHeader = TRUE,style.display = "none"),
        box(width=2,
            title="Download Option", status = "primary", solidHeader = TRUE,
            numericInput(inputId="height3", label="Height (inch):", value = 8.0),
            numericInput(inputId="width3", label="Width (inch):", value = 8.0),
            numericInput(inputId="dpi3", label="DPI:", value = 500),
            downloadButton('plot3_download')
        )
      ) 
    }else if (input$host ==2){
      fluidRow(
        box(width=10,
            plotOutput("plot3", height = 650)),
        box(width=2,
            title="Download Option", status = "primary", solidHeader = TRUE,
            numericInput(inputId="height3", label="Height (inch):", value = 8.0),
            numericInput(inputId="width3", label="Width (inch):", value = 8.0),
            numericInput(inputId="dpi3", label="DPI:", value = 500),
            downloadButton('plot3_download')
        )
      )
    }
  })
  
  output$plot3_description<-renderUI({
    if (input$host == 1){
      box(width=6,title="Description", status = "primary", solidHeader = TRUE,
          collapsible = TRUE,
          HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>).  \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and distinct variant motifs."))
    }else if (input$host==2){
      box(width=10,title="Description", status = "primary", solidHeader = TRUE,
          collapsible = TRUE,
          HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>). \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and distinct variant motifs."))
    }
    
  })

  output$plot7_hosts<-renderUI({
    if(input$host==1){
      plotOutput(outputId = "plot7",height = 650)
    }else {
      plotOutput(outputId = "plot7",height = 1000)
    }    
  })
  
  downloadSampleData(output)
  
  observeEvent(input$samplesubmit,ignoreInit=TRUE,{
    
    shinyjs::addClass(id = "UpdateAnimate", class = "loading dots")
    data<-read.csv("www/DiMA_HCV.csv")
    df <- data.frame(data)
 
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize)
    )
    #absent of user input for protein order; follow the protein order provided in csv file
    if (input$proteinOrder ==""){
      #get protein name and number of positions
      a<-table(df$proteinName)
      proteinName<-as.vector(names(a))
      position<-as.vector(a)
      #store protein as factor
      df$size_f = factor(df$proteinName,levels = proteinName) #ALERT! df$level -> df$size_f
      #scale the nonamer positions
      scales_x<-mapply(function(x,y){
        x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
      }, proteinName,position)
      
    }else{
      #order the proteins based on user input
      level<-unlist(lapply(strsplit(s,','),trimws))
      position<-c()
      for (i in level){
        position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
      }
      #store protein as factor
      df$size_f = factor(df$proteinName,levels = level) #ALERT! df$level -> df$size_f
      #scale the nonamer positions
      scales_x<-mapply(function(x,y){
        x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
      }, level,position)
    }

    #--------------------------------------#
    #               Output                 #
    #--------------------------------------#

    #Tab 1: Entropy and incidence of total variants for each aligned <i>k</i>-mer positions of a viral protein(s)
    plot1<-reactive({
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder,input$kmerlength)
    })

    generate_plot1(input,output,plot1)
    generate_entropyTable(data, output, proteinName)

    #Additional plot: Entropy plot
    plotEntropy<-reactive({
      plot_entropy(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder, input$kmerlength)
    })

    generate_plotEntropy(input, output, plotEntropy)
    

    #Tab 2: Relationship between entropy and total variants for <i>k</i>-mer positions of the viral protein(s)
    plot2<-reactive({
      plot_correlation(df,input$line_dot_size,input$wordsize,input$host)
    })

    generate_plot2(input, output, plot2)
    
    #Tab 3: Dynamics of diversity motifs of viral proteome
    plot3<-reactive({
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,input$host)
    })
    
    generate_plot3(input, output, plot3)
    
    #Tab 4: Dynamics of diversity motifs (Protein)
    plot4<-reactive({
      plot_dynamics_protein(data,input$line_dot_size,input$wordsize,input$host,input$proteinOrder)
    })
    
    generate_plot4(input, output, plot4)
    
    #Tab 7: Conservation levels of viral <i>k</i>-mer positions for each individual protein
    plot7<-reactive({
      data<- data.frame(data)
      data<-data%>%mutate(ConservationLevel = case_when(
        data$index.incidence == 100 ~ "Completely conserved (CC)",
        data$index.incidence >= 90 ~ "Highly conserved (HC)",
        data$index.incidence >= 20 ~ "Mixed variable (MV)",
        data$index.incidence >= 10  ~ "Highly diverse (HD)",
        data$index.incidence < 10 ~ "Extremely diverse (ED)"
        
      ))

       plot_conservationLevel(data,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)  
    })
    
    generate_plot7(input, output, plot7)
    output$conserv_threshold_box <- renderUI({
      textInput(inputId="conserv_percent", label=NULL, 
                value = ifelse(input$conserv_lvl == "CCS", 100, 90))
    })
    generate_CCS_HCS_table(input, output, data)
    
    shinyjs::removeClass(id = "UpdateAnimate", class = "loading dots")
    #Alert results are ready
    output$alertSample <- renderUI({
      h5("Click on other tabs for sample run visualization!", style = "color:white")
    })
  })
  
}
