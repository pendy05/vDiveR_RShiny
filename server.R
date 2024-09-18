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
library(vDiveR)

source("functions/helpers.R")


#server main function
server <- function(input, output,session) {
    # initial state of downloadDiMA button is disabled
    shinyjs::disable(id="downloadDiMA")
    shinyjs::disable(selector = '.nav-tabs a[data-value="Second Host"')
    
    output$footer_wording <- renderText({
        currentYear <- as.numeric(format(Sys.Date(), "%Y"))
        paste("&copy;2021-", currentYear, ", Tok et al. All Rights Reserved.<br>")
  })

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

    observeEvent(input$resetMeta1,{
        resetMetaDataInput(output)
    })

    observeEvent(input$resetMeta2,{
        resetMetaDataInput(output)
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
        shinyjs::addClass(id = "submitMeta1", class = "loading dots")
        Meta <- reactive({
            req(input$Metafile)
            filepath <- input$Metafile$datapath
            Meta <- read.csv(filepath, header = T, stringsAsFactors = F)
            Meta
        })
        WorldmapInfoDf <- reactive({
            req(Meta())
            vDiveR::plot_worldmap(Meta(), input$wordsize)$df
        })
        WorldmapInfoPlot <- reactive({
            req(Meta())
            vDiveR::plot_worldmap(Meta(), input$wordsize)$plot
        })
        TimeInfoDf <- reactive({
            req(Meta())
            date_format <- detect_date_format(Meta()$date)
            vDiveR::plot_time(metadata = Meta(), 
                              base_size = input$wordsize,
                              date_format = date_format,
                              scale = input$time_scale)$df
                              
        })
        TimeInfoPlot <- reactive({
            req(Meta())
            date_format <- detect_date_format(Meta()$date)
            tick_interval <- get_plot_time_tick_interval(Meta()$date)
            vDiveR::plot_time(metadata = Meta(), 
                              base_size = input$wordsize,
                              date_format = date_format,
                              scale = input$time_scale, 
                              date_break = paste(tick_interval, "month"))$plot
                                
        })
        plot_worldmap <- reactive({
            req(WorldmapInfoPlot())
            WorldmapInfoPlot()
            
            
        })
        generate_worldmap(input,output,plot_worldmap)
        output$countrytable = DT::renderDataTable({
            req(WorldmapInfoDf())
            
            WorldmapInfoDf()
        })
        output$table_worldmap_download <- downloadHandler(
            filename = function() {paste("CountryInfo", '.csv', sep='')},
            content = function(file) {
                data <- WorldmapInfoDf()
                if (!is.null(data)) {
                    write.csv(data, file, quote = F)
                } else {
                    stop("No data available for download")
                }}
                )
        plot_time <- reactive({
            req(TimeInfoPlot())
            TimeInfoPlot()
                    
        })
        generate_timeplot(input,output,plot_time)
        output$timetable = DT::renderDataTable({
            req(TimeInfoDf())
            tryCatch({
                TimeInfoDf()
            }, error = function(e) {
                # Log the error in the console
                print(paste("Error in plot_time$df:", e$message))
                # Return NULL or a user-friendly message
                showNotification("Error generating plot time data.", type = "error")
                NULL
            })
            
            
        })
        output$table_time_download <- downloadHandler(
            filename = function() {paste("TimeInfo", '.csv', sep='')},
            content = function(file) {
                data <- TimeInfoDf()
                if (!is.null(data)) {
                    write.csv(data, file, quote = F)
                } else {
                    stop("No data available for download")
                }}
        )
        shinyjs::removeClass(id = "submitMeta1", class = "loading dots")
    })

    observeEvent(input$submitMeta2, {
        shinyjs::addClass(id = "submitMeta2", class = "loading dots")
        Meta <- reactive({
            req(input$Metafasta)
            filepath <- input$Metafasta$datapath
            Meta <- vDiveR::metadata_extraction(filepath, input$MetafastaSource)
            Meta
        })
        output$metademoSee <- DT::renderDT({
            req(input$Metafasta)
            Meta()
        })
        WorldmapInfoDf <- reactive({
            req(Meta())
            vDiveR::plot_worldmap(Meta(), input$wordsize)$df
        })
        WorldmapInfoPlot <- reactive({
            req(Meta())
            vDiveR::plot_worldmap(Meta(), input$wordsize)$plot
        })
        TimeInfoDf <- reactive({
            req(Meta())
            date_format <- detect_date_format(Meta()$date)
            vDiveR::plot_time(metadata = Meta(), 
                              base_size = input$wordsize,
                              date_format = date_format,
                              scale = input$time_scale)$df
        })
        TimeInfoPlot <- reactive({
            req(Meta())
            date_format <- detect_date_format(Meta()$date)
            tick_interval <- get_plot_time_tick_interval(Meta()$date)
            paste(tick_interval, " month")
            vDiveR::plot_time(metadata = Meta(), 
                              base_size = input$wordsize,
                              date_format = date_format,
                              scale = input$time_scale,
                              date_break = paste(tick_interval, "month"))$plot
                                
        })
        plot_worldmap <- reactive({
            req(WorldmapInfoPlot())
            WorldmapInfoPlot()
        })
        generate_worldmap(input,output,plot_worldmap)
        output$countrytable = DT::renderDataTable({
            req(WorldmapInfoDf())
            WorldmapInfoDf()
        })
        output$table_worldmap_download <- downloadHandler(
            filename = function() {paste("CountryInfo", '.csv', sep='')},
            content = function(file) {
                data <- WorldmapInfoDf()
                if (!is.null(data)) {
                    write.csv(data, file, quote = F)
                } else {
                    stop("No data available for download")
                }}
        )
        plot_time <- reactive({
            req(TimeInfoPlot())
            TimeInfoPlot()
        })
        generate_timeplot(input,output,plot_time)
        output$timetable = DT::renderDataTable({
            req(TimeInfoDf())
            TimeInfoDf()
        })
        output$table_time_download <- downloadHandler(
            filename = function() {paste("TimeInfo", '.csv', sep='')},
            content = function(file) {
                data <- TimeInfoDf()
                if (!is.null(data)) {
                    write.csv(data, file, quote = F)
                } else {
                    stop("No data available for download")
                }}
        )
        shinyjs::removeClass(id = "submitMeta2", class = "loading dots")
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

        shinyjs::addClass(id = "submitDiMA", class = "loading dots")
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

            if (!file.exists('venv/Scripts/dima-cli.exe')){
                stop("dima-cli.exe is not found at the specified path 'venv/Scripts/'. Please check the path and try again.")
            }

            #run DiMA
            for (i in 1:length(filepath)){
                outfile<- paste0(proteinName[[i]][1],"_",hostname,"_",i,".json",sep="")

                system(paste("venv/Scripts/dima-cli.exe -i", filepath[i], "-o", file.path(temp_directory,outfile),"-s",input$supportLimit, "-q",proteinName[[i]][1], "-l",input$kmerlength, "-a",inputtype))
 
                json_file <- file.path(temp_directory, outfile)
                csv_file <- json_file %>% str_replace(".json",".csv")
                
                json_data <- fromJSON(json_file)
                dima_df <- vDiveR::json2csv(json_data, hostname_secondHost, proteinName_secondHost[[i]][1])
                write.table(dima_df, sep=",", row.names = FALSE , file = csv_file)

                #store the DiMA csv output names into a list (for further concatenation into one file)
                csvfilelist <- append(csvfilelist, csv_file)
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
                    system(paste("venv/Scripts/dima-cli.exe -i", filepath_secondHost[i], "-o", file.path(temp_directory,outfile),"-s",input$supportLimit, "-q",proteinName_secondHost[[i]][1], "-l",input$kmerlength, "-a",inputtype))
                    
                    #https://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
                    write("\r\n", file = outfile, append = TRUE, sep = "\n")

                    json_file <- file.path(temp_directory, outfile)
                    csv_file <- json_file %>% str_replace(".json",".csv")
                    
                    json_data <- fromJSON(json_file)
                    dima_df <- vDiveR::json2csv(json_data, hostname_secondHost, proteinName_secondHost[[i]][1])
                    write.table(dima_df, sep=",", row.names = FALSE , file = csv_file)
                    
                    #store the DiMA csv output names into a list (for further concatenation into one file)
                    csvfile<-file.path(temp_directory, paste0(proteinName_secondHost[[i]][1],"_",hostname_secondHost,"_",i,".csv",sep=""))

                    csvfilelist_secondHost <- append(csvfilelist_secondHost, csv_file)
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
            write.table(data, sep=",", row.names = FALSE , file = file.path(temp_directory,"DiMA.csv"))
            #Store all the path of the files in the directory in the reactive value
            mylist$files <- list.files(temp_directory,"*.*")

        }else if (input$filetype == 2){ #if the data is DiMA json output
            csvfilelist<-c()
            #convert DiMA output from JSON to CSV
            for (i in 1:length(filepath)){
                #----------------------NOTE (2/5/2022)-------------------#
                #"OUTFILE" SHOULD directly be the name of the input files user provided
                #"proteinName" are needed?
                csvfile<-file.path(temp_directory, paste0(strsplit(input$MSAfile$name[i], ".json")[[1]][1],".csv",sep=""))
                json_data<-fromJSON(filepath[i])
                dima_df<-vDiveR::json2csv(json_data, hostname, proteinName[[i]][1])
                write.table(dima_df, sep=",", row.names = FALSE , file = csvfile)

                #store the DiMA csv output names into a list (for further concatenation into one file)
                csvfilelist <- append(csvfilelist, csvfile)
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
                    csvfile<-file.path(temp_directory, paste0(strsplit(input$MSAfile_secondHost$name[i], ".json")[[1]][1],".csv",sep=""))
                    
                    json_data<-fromJSON(filepath_secondHost[i])
                    dima_df<-vDiveR::json2csv(json_data, hostname_secondHost, proteinName_secondHost[[i]][1])
                    write.table(dima_df, sep=",", row.names = FALSE , file = csvfile)
                    
                    #store the DiMA csv output names into a list (for further concatenation into one file)
                    csvfilelist_secondHost <- append(csvfilelist_secondHost, csvfile)
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
            write.table(data, sep=",", row.names = FALSE , file = file.path(temp_directory,"DiMA.csv"))

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
        plot1 <- reactive({
            vDiveR::plot_entropy(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                # scales_x,
                protein_order = input$proteinOrder,
                kmer_size = input$kmerlength,
                all = T
            )
        })

        generate_plot1(input,output,plot1)


        #Additional plot: Entropy plot
        plotEntropy<-reactive({
            vDiveR::plot_entropy(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                # scales_x,
                protein_order = input$proteinOrder,
                kmer_size = input$kmerlength,
                all = F
            )
        })

        generate_plotEntropy(input, output, plotEntropy)

        #Tab 2: correlation plot
        plot2<-reactive({
            vDiveR::plot_correlation(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host
            )
        })

        generate_plot2(input, output, plot2)

        plot3<-reactive({
            vDiveR::plot_dynamics_proteome(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host
            )
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


        plot4 <- reactive({
            vDiveR::plot_dynamics_protein(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                protein_order = input$proteinOrder
            )
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
            vDiveR::plot_conservationLevel(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                protein_order = input$proteinOrder,
                label_size = (input$wordsize / 2) - 2,
                conservation_label = input$conservationLabel
            )


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
        shinyjs::removeClass(id = "submitDiMA", class = "loading dots")

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

        shinyjs::addClass(id = "samplesubmit", class = "loading dots")
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
            vDiveR::plot_entropy(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                # scales_x,
                protein_order = input$proteinOrder,
                kmer_size = input$kmerlength,
                all = T
            )
        })

        generate_plot1(input,output,plot1)
        generate_entropyTable(data, output, proteinName)

        #Additional plot: Entropy plot
        plotEntropy<-reactive({
            vDiveR::plot_entropy(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                # scales_x,
                protein_order = input$proteinOrder,
                kmer_size = input$kmerlength,
                all = F
            )
        })

        generate_plotEntropy(input, output, plotEntropy)


        #Tab 2: Relationship between entropy and total variants for <i>k</i>-mer positions of the viral protein(s)
        plot2<-reactive({
            vDiveR::plot_correlation(
                df,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host
            )
        })

        generate_plot2(input, output, plot2)

        #Tab 3: Dynamics of diversity motifs of viral proteome
        plot3<-reactive({
            vDiveR::plot_dynamics_proteome(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host
            )
        })

        generate_plot3(input, output, plot3)

        #Tab 4: Dynamics of diversity motifs (Protein)
        plot4<-reactive({
            vDiveR::plot_dynamics_protein(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                protein_order = input$proteinOrder
            )
        })

        generate_plot4(input, output, plot4)

        #Tab 7: Conservation levels of viral <i>k</i>-mer positions for each individual protein
        plot7<-reactive({
            vDiveR::plot_conservationLevel(
                data,
                line_dot_size = input$line_dot_size,
                base_size = input$wordsize,
                host = input$host,
                protein_order = input$proteinOrder,
                label_size = (input$wordsize / 2) - 2,
                conservation_label = input$conservationLabel
            )
        })

        generate_plot7(input, output, plot7)
        output$conserv_threshold_box <- renderUI({
            textInput(inputId="conserv_percent", label=NULL,
                      value = ifelse(input$conserv_lvl == "CCS", 100, 90))
        })
        generate_CCS_HCS_table(input, output, data)

        shinyjs::removeClass(id = "samplesubmit", class = "loading dots")
        #Alert results are ready
        output$alertSample <- renderUI({
            h5("Click on other tabs for sample run visualization!", style = "color:white")
        })
    })

}
