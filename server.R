#set maximum upload file size
options(shiny.maxRequestSize=3000*1024^2)
#https://www.shinyapps.io/admin/#/application/5536925/logs
library(ggplot2)
library(gridExtra) #tutorial: https://ggplot2.tidyverse.org/reference/facet_grid.html 
library(facetscales) #https://stackoverflow.com/a/54074323/13970350
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
library(reticulate)
library(glue)
library(shinyThings) # devtools::install_github("gadenbuie/shinyThings")
#py_run_string("from dima import Dima")
#Sys.setenv(RETICULATE_PYTHON = "python_env/Scripts/python.exe")
#reticulate::use_virtualenv("./python_env", required = TRUE)

# #uncomment codes from line 26 to 28 if you would like to run DiveR locally (prerequiste: a Python virtual environment is needed; refer README for more instructions)
# virtualenv_create(envname = "python_env", python= "python3")
# virtualenv_install("python_env", packages = c('pandas','numpy','dima-cli==4.1.1'))
#reticulate::use_virtualenv("python_env", required = TRUE)

#server side
server <- function(input, output,session) {
  # initial state of downloadDiMA button is disabled
  shinyjs::disable(id="downloadDiMA")
  shinyjs::disable(selector = '.nav-tabs a[data-value="Second Host"')
  
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
    #reset input
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
    shinyjs::disable(id="downloadDiMA")
    
    #clear output
    output$alert <- renderUI({})
    output$alertSample <- renderUI({})
    output$plot1<- renderPlot({})
    output$plot2<- renderPlot({})
    output$plot3<- renderPlot({})
    #output$plot3_hosts<-renderUI({})
    output$plot4<- renderPlot({})
    #output$plot7_hosts<-renderUI({})
    output$plot7<- renderPlot({})
    #output$infilename<-renderText({NULL })
    #output$infilename_secondHost<-renderText({ NULL})
  })
  
  #reset every input and output to initial state
  observeEvent(input$reset1, {
    #reset input
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
    shinyjs::disable(id="downloadDiMA")
    
    #clear output
    output$alert <- renderUI({})
    output$alertSample <- renderUI({})
    output$plot1<- renderPlot({})
    output$plot2<- renderPlot({})
    output$plot3<- renderPlot({})
    #output$plot3_hosts<-renderUI({})
    output$plot4<- renderPlot({})
    #output$plot7_hosts<-renderUI({})
    output$plot7<- renderPlot({})
    #output$infilename<-renderText({NULL })
    #output$infilename_secondHost<-renderText({ NULL})
  })
  
  #change tab to input tab when user clicks on start button
  #https://stackoverflow.com/questions/32971921/navigate-to-particular-sidebar-menu-item-in-shinydashboard
  observeEvent(input$start, {
    newtab <- switch(input$tabs,
                     "description" = "inputdata_description",
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
  observeEvent(input$submitDiMA,ignoreInit=TRUE,{
    req(input$MSAfile)
    shinyjs::addClass(id = "submitAnimate", class = "loading dots")
    
    #proceed if the input csv file is provided
    if (is.null(input$MSAfile)){
      return(NULL)
    }
    #default host name: unknown
    if (input$hostname == ""){
      hostname <-"Unknown"
    }else{
      hostname<-input$hostname
    }
    
    #default host name: unknown
    if (input$hostname_secondHost == ""){
      hostname_secondHost <-"Unknown"
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
        print("different len")
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
        print("different len")
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
          print("input$inputtype == 1")
          inputtype <- "protein"
      }else{
          inputtype <- "nucleotide"
      }
      #run DiMA
      for (i in 1:length(filepath)){
        print(proteinName[[i]][1])
        
        outfile<- paste0(proteinName[[i]][1],"_",hostname,"_",i,".json",sep="")

        system(paste("python_env/Scripts/dima-cli.exe -i", filepath[i], "-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-q",proteinName[[i]][1], "-l",input$kmerlength, "-a",inputtype))

        #https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
        write("\r\n", file = paste0("./",outfile), append = TRUE, sep = "\n")
        json2csvinR_unnest(paste0(temp_directory, "/",outfile),hostname, proteinName[[i]][1])

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
        print("secondHost")
        outfile_secondHost<-""
        csvfilelist_secondHost<-c()
        #run DiMA
        for (i in 1:length(filepath_secondHost)){
          outfile<- paste0(proteinName_secondHost[[i]][1],"_",hostname_secondHost,"_",i,".json",sep="")
          print(outfile)
          py_run_string("from dima import Dima")
          dima_input<- paste0("sequences=r'",filepath_secondHost[i],"',kmer_length=",input$kmerlength,",support_threshold=",input$supportLimit)
          #dima_input<- gsub("/", "\\", dima_input)
          print(dima_input)
          print("results")
          py_run_string(glue("results = Dima({dima_input}).run()"))
          #py_run_string(glue("results = Dima(sequences=\"{filepath[i]}\", kmer_length={input$kmerlength},support_threshold={input$supportLimit}, query_name=\"{proteinName[[i]][1]}\").run()"))
          print("jsonfile")
          py_run_string(glue("jsonFile = open(r'{temp_directory}/{outfile}', 'w')"))
          print("jsonfile writing...")
          py_run_string(glue("jsonFile.write(str(results))"))
          print("jsonfile closing...")
          py_run_string("jsonFile.close()")
          #system("chmod +x python_env/Scripts/dima-cli.exe dima-cli.exe")
          #system("ls -lh python_env/Scripts/dima-cli.exe")
          #system("ls -lh")
          #system("file ./vim")
          #system(paste("dima-cli.exe -i ",filepath_secondHost[i]," -o ",paste0(temp_directory, "/",outfile)," -s ",input$supportLimit, " -p ",proteinName_secondHost[[i]][1], " -l ",input$kmerlength))
          
          #system(paste("./python_env/Scripts/dima-cli.exe -i",filepath_secondHost[i],"-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-p",proteinName_secondHost[[i]][1], "-l",input$kmerlength))
          #system2(command="./python_env/Scripts/dima-cli.exe",args = c("-i",filepath_secondHost[i],"-o",paste0(temp_directory, "/",outfile),"-s",input$supportLimit, "-p",proteinName_secondHost[[i]][1], "-l",input$kmerlength))
          #system2(command="./python_env/Scripts/dima-cli.exe",args = c("-i",filepath[i],"-o",outfile,"-s",input$supportLimit, "-p",proteinName[[i]][1], "-l",input$kmerlength))
          #print(input$supportLimit, proteinName[[i]][1], input$kmerlength)
          #append "\n" to the end of file to solve the issue of 'incomplete final line'
          #https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
          write("\r\n", file = outfile, append = TRUE, sep = "\n")
          print("json2csv...")
          json2csvinR_unnest(paste0(temp_directory, "/",outfile),hostname_secondHost,proteinName_secondHost[[i]][1])
          #json2csvinR(paste0(temp_directory, "/",outfile),hostname)
          #json2csvinR(outfile,hostname)
          print("json2csv ends...")
          #store the DiMA csv output names into a list (for further concatenation into one file)
          csvfile<-paste0(temp_directory, "/",proteinName_secondHost[[i]][1],"_",hostname_secondHost,"_",i,".csv",sep="")
          #csvfile<-paste0(proteinName[[i]][1],"_",i,".csv",sep="")
          
          csvfilelist_secondHost <- append(csvfilelist_secondHost, csvfile)
          #outfile<- paste0(proteinName[[i]][1],"_",i,".json",sep="")
          #system2(command="cat", args=c(paste0(proteinName[[i]][1],"_",i,".csv",sep=""),">>", "dima.csv"))
        }
        
        data_secondHost <-  read.csv(csvfilelist_secondHost[1])
        #Saving the first DiMA csv output file in the filepath array in the temp zipped folder
        #write.csv(data, file = paste0(session$token, "/", proteinName[[1]][1],"_",i,".csv",sep=""))
        
        #reading each file within the range and append them to create one file
        for (f in csvfilelist_secondHost[-1]){
          #print("f: ",f)
          df <- read.csv(f)      # read the file
          #Saving the following DiMA csv output file in the filepath array in the temp zipped folder
          #write.csv(data, file = paste0(session$token, "/", f, sep=""))
          data_secondHost <- rbind(data_secondHost, df)    # append the current file
        }
        
        data<- rbind(data,data_secondHost)
      }
      
      #write to a final csv file "DiMAoutput.csv", consists of all the submitted proteins in the temp zipped folder
      #write.table(data, sep=",", row.names = FALSE , file = paste0(session$token,"/DiMA.csv"))
      write.table(data, sep=",", row.names = FALSE , file = paste0(temp_directory,"/DiMA.csv"))
      
      #Store all the path of the files in the directory in the reactive value
      mylist$files <- list.files(temp_directory,"*.*")
      #print("mylist$files: ",mylist$files)
      
      
      
    }else if (input$filetype == 2){ #if the data is DiMA json output
      csvfilelist<-c()
      #convert DiMA output from JSON to CSV
      for (i in 1:length(filepath)){
        #print(i)
        print("input file type 2")
        print(input$MSAfile$name[i])
        
        #----------------------NOTE (2/5/2022)-------------------#
        #"OUTFILE" SHOULD directly be the name of the input files user provided
        #"proteinName" are needed?
        
        #outfile<- paste0(proteinName[[i]][1],"_",i,".json",sep="") 
        csvfile<-paste0(strsplit(input$MSAfile$name[i], ".json")[[1]][1],".csv",sep="")
        print("csv file")
        print(csvfile)
        print("json2csv...input type 2")
        direct_json2csvinR(filepath[i],hostname,proteinName[[i]][1], paste0(temp_directory,"/",csvfile))
        print("json2csv ends...")
        #store the DiMA csv output names into a list (for further concatenation into one file)
        
        csvfilelist <- append(csvfilelist, paste0(temp_directory,"/",csvfile))
        #outfile<- paste0(proteinName[[i]][1],"_",i,".json",sep="")
        #system2(command="cat", args=c(paste0(proteinName[[i]][1],"_",i,".csv",sep=""),">>", "dima.csv"))
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
        print("option 2 secondHost")
        csvfilelist_secondHost<-c()
        #run DiMA
        #showModal(modalDialog("Running DiMA...", footer=NULL))
        for (i in 1:length(filepath_secondHost)){
          print(proteinName_secondHost[[i]][1])
          csvfile<-paste0(strsplit(input$MSAfile_secondHost$name[i], ".json")[[1]][1],".csv",sep="")
          print("csv file")
          print(csvfile)
          print("json2csv...input type 2")
          direct_json2csvinR(filepath[i],hostname_secondHost,proteinName_secondHost[[i]][1], paste0(temp_directory,"/",csvfile))
          print("json2csv ends...")
          #store the DiMA csv output names into a list (for further concatenation into one file)
          csvfilelist_secondHost <- append(csvfilelist_secondHost, paste0(temp_directory,"/",csvfile))
        }
        
        data_secondHost <-  read.csv(csvfilelist_secondHost[1])
        #Saving the first DiMA csv output file in the filepath array in the temp zipped folder
        #write.csv(data, file = paste0(session$token, "/", proteinName[[1]][1],"_",i,".csv",sep=""))
        
        #reading each file within the range and append them to create one file
        for (f in csvfilelist_secondHost[-1]){
          #print("f: ",f)
          df <- read.csv(f)      # read the file
          #Saving the following DiMA csv output file in the filepath array in the temp zipped folder
          #write.csv(data, file = paste0(session$token, "/", f, sep=""))
          data_secondHost <- rbind(data_secondHost, df)    # append the current file
        }
        
        data<- rbind(data,data_secondHost)
      }
      
      #write to a final csv file "DiMAoutput.csv", consists of all the submitted proteins in the temp zipped folder
      write.table(data, sep=",", row.names = FALSE , file = paste0(temp_directory,"/DiMA.csv"))
      
      #Store all the path of the files in the directory in the reactive value
      #mylist$files <- list.files(session$token,"*.*")
      #print(mylist$files)
      
    }else if (input$filetype == 3){ #if the data is DiMA csv output
      #read the input file
      if(input$host == 1){ #1 host
        data <- filepath %>%
          lapply(read_csv) %>%
          bind_rows
      }else{ #2 hosts
        filepath<-append(filepath,filepath_secondHost)
        data <- filepath %>%
          lapply(read_csv) %>%
          bind_rows
      }
      
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
    
    
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize))
    
    
    #--------------------Table Output----------------------#
    #get position of min entropy, min, max of entropy and total variant
    outputTable <- data %>% 
      dplyr::group_by(proteinName) %>%
      dplyr::summarise(
        Position = gsub("((?:\\d+,){2}\\d+),", "\\1,\n", paste0(position[which(entropy == min(entropy))], collapse = ",")),
        minEntropy = format(round(min(entropy),digits=2),nsmall=2),
        maxEntropy = format(round(max(entropy), digits = 2),nsmall=2),
        minTotalVariants = format(round(min(totalVariants.incidence),digits=2),nsmall=2),
        maxTotalVariants = format(round(max(totalVariants.incidence), digits = 2),nsmall=2)
      )
    
    print(outputTable)
    #rename table df
    names(outputTable)<- c("Protein Name","Position (Minimum Entropy)","Minimum Entropy","Maximum Entropy","Minimum Total Variants (%)","Maximum Total Variants (%)")
    
    
    output$table <- renderDataTable(
      outputTable#,
      #width = "100%"
    )
    
    #----------------------Plotting-----------------------#
    
    output$text <- renderText({
      print(input$proteinOrder)
    })

    plot1<-reactive({
      req(df)
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder, input$kmerlength)
    })
    
    output$plot1 <- renderPlot({
      plot1()
      
    })
    
    output$plot1_download <- downloadHandler(
      filename = function() { paste("plot1", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot1(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
      }
    )
    
    output$info_plot1 <- renderText({
      paste0("Position = ", input$plot1_click$x,"\nEntropy (bits) = ",input$plot1_click$y)
    })
    
    plot2<-reactive({
      print("plot 2 reactive func")
      print(unique(df$proteinName))
      plot_correlation(df,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot2<- renderPlot({
      plot2()
    })
    
    output$plot2_download <- downloadHandler(
      filename = function() { paste("plot2", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi2)
      }
    )
    
    output$info_plot2 <- renderText({
      paste0("Total variants (%) = ", input$plot2_click$x,"\nNonamer entropy (bits) = ",input$plot2_click$y)
    })
    
    output$plot3_description<-renderUI({
      if (input$host ==1){
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
                    for detailed definition</i>).  \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and distinct variant motifs."))
      }
      
    })
    
    plot3<-reactive({
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot3<- renderPlot({
      plot3()
    })
    
    output$plot3_hosts<-renderUI({
      if (input$host==1){
        fluidRow(
          box(width=6,
              plotOutput("plot3", height = 650))
          ,
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
              plotOutput("plot3", height = 650))
          ,
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
    
    
    output$plot3_download <- downloadHandler(
      filename = function() { paste("plot3", '.jpg', sep='') },
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
      filename = function() { paste("plot4", '.jpg', sep='') },
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
        print("plot 7 - 2 host")
        plot7_multihost<-lapply(plot7_list,plot_conservationLevel,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)
        
        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$host))))
        
      }
      
    })
    
    output$plot7_hosts<-renderUI({
      if(input$host==1){
        plotOutput(outputId = "plot7",height = 650)
      }else{
        plotOutput(outputId = "plot7",height = 1000)
      }
      
    })
    
    output$plot7<- renderPlot({
      plot7()
    })
    
    output$plot7_download <- downloadHandler(
      filename = function() { paste("plot7", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
      })
    
    
    
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
            plotOutput("plot3", height = 650))
        ,
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
            plotOutput("plot3", height = 650))
        ,
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
    if (input$host ==1){
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
    plotOutput(outputId = "plot7",height = 650)
    
  })
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("DiMA_sample_output_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) { 
      #   write.csv(CSV_dataset, file)
      #   write.fasta(sequences = MSA_dataset, names = names(MSA_dataset))
      #   write(JSON_dataset)
      # }
      print('donwloading')
      #To create temporary directory
      temp_directory_sample <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory_sample)
      
      csv_dataset <- read.csv("www/DiMA_HCV.csv")
      MSA_dataset_Core <- read.fasta("www/Core_mafft.fasta")
      csv_dataset_Core <- read.csv("www/core_9mer.csv")
      JSON_dataset_Core<- rjson::fromJSON(file = "www/core_9mer.json")
      MSA_dataset_NS3 <- read.fasta("www/NS3_mafft.fasta")
      csv_dataset_NS3 <- read.csv("www/NS3_9mer.csv")
      JSON_dataset_NS3<- rjson::fromJSON(file = "www/NS3_9mer.json")
      print('donwloading.....')
      #write sample dataset in CSV, FA & JSON formats into temp directory
      write.csv(csv_dataset,paste0(temp_directory_sample,"/HCV_DiMA.csv"))
      write.fasta(MSA_dataset_Core,names=names(MSA_dataset_Core),file.out = paste0(temp_directory_sample,"/HCV_aligned_Core.fasta"))
      write.csv(csv_dataset_Core,paste0(temp_directory_sample,"/HCV_Core.csv"))
      write_json(JSON_dataset_Core,paste0(temp_directory_sample,"/HCV_Core.json"))
      write.fasta(MSA_dataset_NS3,names=names(MSA_dataset_NS3),file.out = paste0(temp_directory_sample,"/HCV_aligned_NS3.fasta"))
      write.csv(csv_dataset_NS3,paste0(temp_directory_sample,"/HCV_NS3.csv"))
      write_json(JSON_dataset_NS3,paste0(temp_directory_sample,"/HCV_NS3.json"))
      print('donwloading..............')
      zip::zip(zipfile = file, 
               files = dir(temp_directory_sample),
               root = temp_directory_sample)
    },
    
    contentType = "application/zip"
  )
  
  #wait for user's input
  observeEvent(input$submit,ignoreInit=TRUE,{
    req(input$csvfile)
    #proceed if the input csv file is provided
    if (is.null(input$csvfile)){
      return(NULL)
    }
    #read the csv file and store into dataframe
    data <- read.csv(input$csvfile$datapath)
    df <- data.frame(data)
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
        df$size_f = factor(df$proteinName,levels = proteinName) 
        #scale the nonamer positions
        scales_x<-mapply(function(x,y){
          x = scale_x_continuous(limits = c(0,y),breaks = seq(0,y,50))
        }, proteinName,position)
        
      }else{
        if (!"host" %in% colnames(df)){
          print("Host column is not detected!")
        }
        #order the proteins based on user input
        level<-unlist(lapply(strsplit(input$proteinOrder,','),trimws))
        position<-c()
        for (i in level){
          position<-append(position,table(df$proteinName)[names(table(df$proteinName)) == i])
        }
        #store protein as factor
        df$size_f = factor(df$proteinName,levels = level)
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
    
    
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize))
    
    
    #--------------------Table Output----------------------#
    #get position of min entropy, min, max of entropy and total variant
    outputTable <- data %>% 
      dplyr::group_by(proteinName) %>%
      dplyr::summarise(
        Position = paste0(position[which(entropy == min(entropy))], collapse = ","),
        minEntropy = format(round(min(entropy),digits=2),nsmall=2),
        maxEntropy = format(round(max(entropy), digits = 2),nsmall=2),
        minTotalVariants = format(round(min(totalVariants.incidence),digits=2),nsmall=2),
        maxTotalVariants = format(round(max(totalVariants.incidence), digits = 2),nsmall=2)
      )
    
    print(outputTable)
    #rename table df
    names(outputTable)<- c("Protein Name","Position (Minimum Entropy)","Minimum Entropy (%)","Maximum Entropy (%)","Minimum Total Variants (%)","Maximum Total Variants (%)")
    
    
    output$table <- renderTable(
      outputTable,
      width = "100%"
    )
    
    #----------------------Plotting-----------------------#
    
    output$text <- renderText({
      print(input$proteinOrder)
    })
    
    plot1<-reactive({
      req(df)
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder, input$kmerlength)
    })
    
    output$plot1 <- renderPlot({
      plot1()
      
    })
    
    output$plot1_download <- downloadHandler(
      filename = function() { paste("plot1", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot1(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
      }
    )
    
    plot2<-reactive({
      plot_correlation(df,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot2<- renderPlot({
      plot2()
    })
    
    output$plot2_download <- downloadHandler(
      filename = function() { paste("plot2", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi2)
      }
    )
    
    output$info_plot2 <- renderUI({
      HTML(paste0("Total variants (%) = ", input$plot2_click$x,"<br>",em("k"),"-mer entropy (bits) = ",input$plot2_click$y))
    })
    
    output$plot3_hosts<-renderUI({
      if (input$host==1){
        fluidRow(
          box(width=6,
              plotOutput("plot3", height = 650))
          ,
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
              plotOutput("plot3", height = 650))
          ,
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
      if (input$host ==1){
        box(width=6,title="Description", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>). \
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
    
    plot3<-reactive({
      print("plot3")
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot3<- renderPlot({
      plot3()
    })
    
    output$plot3_download <- downloadHandler(
      filename = function() { paste("plot3", '.jpg', sep='') },
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
      filename = function() { paste("plot4", '.jpg', sep='') },
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
        plot_plot7conservation(data,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)
      }else{#multihost
        data$host = factor(data$host)
        #split the data into multiple subsets (if multiple hosts detected)
        plot7_list<-split(data,data$host)
        plot7_multihost<-lapply(plot7_list,plot_plot7conservation,input$line_dot_size, input$wordsize, input$host, input$proteinOrder,input$conservationLabel)
        
        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$host))))
        
      }
      
    })
    
    output$plot7_hosts<-renderUI({
      if(input$host==1){
        plotOutput(outputId = "plot7",height = 650)
      }else{
        plotOutput(outputId = "plot7",height = 1000)
      }
      
    })
    
    output$plot7<- renderPlot({
      plot7()
    })
    
    output$plot7_download <- downloadHandler(
      filename = function() { paste("plot7", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
      })
    
    
  })
  
  
  observeEvent(input$samplesubmit,ignoreInit=TRUE,{
    shinyjs::addClass(id = "UpdateAnimate", class = "loading dots")
    data<-read.csv("www/DiMA_HCV.csv")
    df <- data.frame(data)
    
    output$plot3_hosts<-renderUI({
      fluidRow(
        box(width=6,
            plotOutput("plot3", height = 650))
        ,
        box(id="container",style = "background-color:#ecf0f5;",width=4,solidHeader = TRUE,style.display = "none"),
        box(width=2,
            title="Download Option", status = "primary", solidHeader = TRUE,
            numericInput(inputId="height3", label="Height (inch):", value = 8.0),
            numericInput(inputId="width3", label="Width (inch):", value = 8.0),
            numericInput(inputId="dpi3", label="DPI:", value = 500),
            downloadButton('plot3_download')
        )
      )
    })
    
    output$plot3_description<-renderUI({
      
      box(width=6,title="Description", status = "primary", solidHeader = TRUE,
          collapsible = TRUE,
          HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>).  \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and distinct variant motifs."))
      
    })
    #determine number of host
    host = 1
    if (host == 1){
      #single host
      
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
      
    }
    
    
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize))
    
    #--------------------Table Output----------------------#
    #get position of min entropy, min, max of entropy and total variant
    outputTable <- data %>% 
      dplyr::group_by(proteinName) %>%
      dplyr::summarise(
        Position = paste0(position[which(entropy == min(entropy))], collapse = ","),
        minEntropy = format(round(min(entropy),digits=2),nsmall=2),
        maxEntropy = format(round(max(entropy), digits = 2),nsmall=2),
        minTotalVariants = format(round(min(totalVariants.incidence),digits=2),nsmall=2),
        maxTotalVariants = format(round(max(totalVariants.incidence), digits = 2),nsmall=2)
      )
    
    print(outputTable)
    #rename table df
    names(outputTable)<- c("Protein Name","Position (Minimum Entropy)","Minimum Entropy (%)","Maximum Entropy (%)","Minimum Total Variants (%)","Maximum Total Variants (%)")
    
    
    output$table <- renderTable(
      outputTable,
      width = "100%"
    )
    
    #----------------------Plotting-----------------------#
    
    output$text <- renderText({
      print(input$proteinOrder)
    })
    
    
    output$text <- renderText({
      print(input$proteinOrder)
    })
    
    plot1<-reactive({
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder,input$kmerlength)
    })
    
    
    output$plot1 <- renderPlot({
      plot1()
      
    })
    
    output$plot1_download <- downloadHandler(
      filename = function() { paste("plot1", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot1(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
      }
    )
    
    output$info_plot1 <- renderText({
      paste0("Position = ", input$plot1_click$x,"\nEntropy (bits) = ",input$plot1_click$y)
    })
    
    plot2<-reactive({
      plot_correlation(df,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot2<- renderPlot({
      plot2()
    })
    
    output$plot2_download <- downloadHandler(
      filename = function() { paste("plot2", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi)
      }
    )
    
    output$info_plot2 <- renderUI({
      HTML(paste0("Total variants (%) = ", input$plot2_click$x, "<br>",em("k"),"-mer entropy (bits) = ",input$plot2_click$y))
    })
    
    plot3<-reactive({
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,host)
    })
    
    
    output$plot3<- renderPlot({
      plot3()
    })
    
    
    
    output$plot3_download <- downloadHandler(
      filename = function() { paste("plot3", '.jpg', sep='') },
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
      filename = function() { paste("plot4", '.jpg', sep='') },
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
      
      plot_plot7<- function(data){
        #add word 'protein' in front of each protein name
        data$proteinName<-paste("Protein",data$proteinName)
        #create data for proteome bar "All" from existing data
        data1<-data
        data1$proteinName <- "All"
        data1$level <- "All"
        #set up the order of proteins in plot from left to right
        if (input$proteinOrder ==""){ #follow the default order in csv file
          level<-c("All",unique(data$proteinName))
        }else{ #order the proteins based on user input
          proteins<-paste("Protein",unlist(lapply(strsplit(input$proteinOrder,','),trimws)))
          level<-c("All",proteins)
        }
        #determine the protein order
        data$level = factor(data$proteinName, levels=level)
        #combine proteome bar with protein bars
        data<-rbind(data1,data)
        #determine the protein order
        data$level = factor(data$proteinName, levels=level)
        
        #calculation for total and percentage of conservation levels for each protein
        #sum up the total positions for each conservation level of proteins
        plot7_data<-ddply(data,.(proteinName,ConservationLevel),nrow)
        names(plot7_data)[3]<-"Total"
        
        C_level<- c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)")
        #check the presence of conservation level: insert value 0 if it is absent
        if (input$conservationLabel == 1){ #full label
          #check the presence of conservation level: insert value 0 if it is absent
          for ( conservation in C_level){ #conservation level
            for (name in level){ #proteinName
              if (!(conservation %in% plot7_data[plot7_data$proteinName==name,]$ConservationLevel)){
                plot7_data<-rbind(plot7_data,c(name,conservation,0))
              }}}
        }
        #sort the dataframe
        plot7_data[order(plot7_data$proteinName),]
        plot7_data$Total<- as.integer(plot7_data$Total)
        #get the percentage of each conservation level for each protein
        plot7_data<-ddply(plot7_data,.(proteinName),transform, percent=Total/sum(Total)*100)
        
        #gather the protein label in multicolor
        plot7_data<-plot7_data%>%mutate(Label = case_when(
          plot7_data$ConservationLevel == "Completely conserved (CC)" ~ paste0(sprintf("<span style =
    'color:#000000;'>CC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f (%.1f %%) </span>",plot7_data$Total, round(plot7_data$percent,1))),
          
        ))
        
        #set conservation level in specific order (CC,HC,MV,HD,ED)
        plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
        
        #combine all conservation level labels into one for each protein
        Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
        #get number of protein for labelling
        nProtein<-nrow(Proteinlabel)
        #plotting
        plot7<- ggplot(data, aes(x=level,y=index.incidence))+
          # gghalves
          geom_half_boxplot(outlier.shape = NA) +
          geom_half_point(aes(col = ConservationLevel), side = "r", 
                          position = position_jitter(width = 0, height=-0.7), size=(input$line_dot_size),alpha=0.7) +
          ylim(0,105) +
          labs(x=NULL, y="Index incidence (%)\n", fill="Conservation level")+
          theme_classic(base_size = input$wordsize)+
          theme(
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.position = 'bottom',
            plot.margin = unit(c(5, 1, 1, 1), "lines"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)
          )+
          scale_colour_manual('Conservation Level',
                              breaks = c("Completely conserved (CC)",
                                         "Highly conserved (HC)",
                                         "Mixed variable (MV)",
                                         "Highly diverse (HD)",
                                         "Extremely diverse (ED)"),
                              values = c("Completely conserved (CC)"="black",
                                         "Highly conserved (HC)"="#0057d1",
                                         "Mixed variable (MV)"="#02d57f",
                                         "Highly diverse (HD)"="#8722ff", 
                                         "Extremely diverse (ED)"="#ff617d")) +  
          geom_richtext(data = Proteinlabel, 
                        aes(x=proteinName,label = Label, y=c(rep(105,nProtein)),
                            label.size=0, label.color="transparent"),
                        position = position_dodge(width=0.1), 
                        size=((input$wordsize/2)-2), color="black", hjust=0, angle=90) + 
          guides(color = guide_legend(override.aes = list(size = 2), nrow=2))+
          coord_cartesian(clip = "off")+
          ggtitle(unique(data$host))
        
        # plot7<-ggplot(data) +
        #   geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
        #   geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
        #               size=input$line_dot_size)+
        #   labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
        #   scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
        #   theme_classic(base_size = input$wordsize)+
        #   theme(
        #     legend.key = element_rect(fill = "transparent", colour = "transparent"),
        #     legend.position="bottom"
        #   )+
        #   scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
        #                       values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
        #   geom_richtext(data = Proteinlabel, aes(x=proteinName,label = Label, y=c(rep(105,nProtein)), label.size=0, label.color="transparent"),
        #                 position = position_dodge(width=0.1),size=((input$wordsize/2)-2),color="black", fill="white",hjust=0,angle=90) + 
        #   guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
        #   theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
        #   coord_cartesian(clip = "off") #allow ggtext outside of the plot
        plot7
      }
      
      #single host
      if (input$host == 1){
        plot_plot7(data)
      }else{#multihost
        data$host = factor(data$host)
        #split the data into multiple subsets (if multiple hosts detected)
        plot7_list<-split(data,data$host)
        plot7_multihost<-lapply(plot7_list,plot_plot7)
        
        #create spacing between multihost plots
        theme = theme(plot.margin = unit(c(2.5,1.0,0.1,0.5), "cm"))
        do.call("grid.arrange", c(grobs=lapply(plot7_multihost,"+",theme), nrow = length(unique(data$host))))
        
      }
    })
    
    output$plot7<- renderPlot({
      plot7()
    })
    
    # for now not splitted by hosts
    output$plot7_seqs <- renderDataTable({
      seqConcatenation(input_file=data.frame(data), kmer=input$kmerlength, conservation=input$conserv_lvl)[[input$table_type]]
    })
    
    output$plot7_download <- downloadHandler(
      filename = function() { paste("plot7", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
      })
    
    output$conservSeq_download <- downloadHandler(
      filename =  function() {paste0(input$conserv_lvl, ".", input$table_type)},
      content = function(fname) {
        df <- seqConcatenation(input_file=data.frame(data), kmer=input$kmerlength, conservation=input$conserv_lvl)[[input$table_type]]
        write.table(df, file = fname, col.names = ifelse(input$table_type == "csv", TRUE, FALSE),
                    sep = ",", row.names = FALSE, quote = FALSE)
      }
    )
    
    shinyjs::removeClass(id = "UpdateAnimate", class = "loading dots")
    #Alert results are ready
    output$alertSample <- renderUI({
      h5("Click on other tabs for sample run visualization!", style = "color:white")
    })
  })
  
}
