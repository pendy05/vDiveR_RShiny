library(ggplot2)
library(gridExtra) #tutorial: https://ggplot2.tidyverse.org/reference/facet_grid.html 
library(facetscales) #https://stackoverflow.com/a/54074323/13970350
library(plyr)
library(dplyr)
library(ggtext)
library(ggpubr)
library(grid)



#server side
server <- function(input, output,session) {
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
                    for detailed definition</i>). Nonatypes are defined as distinct nonamers for a given <i>k</i>-mer position. \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and <i>k</i>-merTypes motifs."))
    }else if (input$host==2){
      box(width=10,title="Description", status = "primary", solidHeader = TRUE,
          collapsible = TRUE,
          HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>). Nonatypes are defined as distinct nonamers for a given <i>k</i>-mer position. \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and <i>k</i>-merTypes motifs."))
    }
    
  })
  
  output$plot7_hosts<-renderUI({
    plotOutput(outputId = "plot7",height = 650)
    
  })
  dataset <- read.csv("www/sample.csv")

  output$downloadData <- downloadHandler(
    filename = function() {
      myfile <- "sample.csv"
      myfile
    },
    content = function(file) { 
      write.csv(dataset, file)
    })
  
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
    
    #start from here blm check
    
    group_names<-c("Index","Major","Minor","Unique","Total variants","Nonatypes")
    
    
    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize))

    output$text <- renderText({
      print(input$proteinOrder)
    })
    
    plot1<-reactive({
      req(df)
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder)
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
    
    output$info_plot2 <- renderText({
      paste0("Total variants (%) = ", input$plot2_click$x,"\nNonamer entropy (bits) = ",input$plot2_click$y)
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
                    for detailed definition</i>). Nonatypes are defined as distinct nonamers for a given <i>k</i>-mer position. \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and <i>k</i>-merTypes motifs."))
      }else if (input$host==2){
        box(width=10,title="Description", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            HTML("<i>k</i>-mer (peptide sequences of <i>k</i> bases) are classified into four different motifs, \
                    namely index, major, minor and unique, based on their incidences (<i>please refer <b>Project Description</b> section\
                    for detailed definition</i>). Nonatypes are defined as distinct nonamers for a given <i>k</i>-mer position. \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and <i>k</i>-merTypes motifs."))
      }
      
    })
    
    plot3<-reactive({
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
    data<-read.csv("www/sample.csv")
    df <- data.frame(data)
    #determine number of host
    if (input$host == 1){
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
      
    }else{ # multihost
      
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
        level<-unlist(lapply(strsplit(s,','),trimws))
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

    MY_THEME<-theme(
      axis.title.x = element_text(size = input$wordsize),
      axis.text.x = element_text(size = input$wordsize),
      axis.title.y = element_text(size = input$wordsize))

    output$text <- renderText({
      print(input$proteinOrder)
    })
    
    plot1<-reactive({
      plot_entropy_incidence(df,input$line_dot_size,input$wordsize,input$host,scales_x,input$proteinOrder)
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
        ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi)
      }
    )
    
    output$info_plot2 <- renderText({
      paste0("Total variants (%) = ", input$plot2_click$x,"\nNonamer entropy (bits) = ",input$plot2_click$y)
    })
    
    plot3<-reactive({
      plot_dynamics_proteome(data,input$line_dot_size,input$wordsize,input$host)
    })
    
    output$plot3<- renderPlot({
      plot3()
    })
    
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
                    for detailed definition</i>). Nonatypes are defined as distinct nonamers for a given <i>k</i>-mer position. \
                    The diversity of the position was depicted by the decline of the index incidences (black),\
                     the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and <i>k</i>-merTypes motifs."))
      
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
    'color:#000000;'>CC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Highly conserved (HC)" ~ paste0(sprintf("<span style =
    'color:#0057d1;'>HC: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Mixed variable (MV)" ~ paste0(sprintf("<span style =
    'color:#02d57f;'>MV: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Highly diverse (HD)" ~ paste0(sprintf("<span style =
    'color:#A022FF;'>HD: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
          plot7_data$ConservationLevel == "Extremely diverse (ED)" ~ paste0(sprintf("<span style =
    'color:#ff617d;'>ED: %.0f; %.1f %% </span>",plot7_data$Total, round(plot7_data$percent,1))),
          
        ))
        
        #set conservation level in specific order (CC,HC,MV,HD,ED)
        plot7_data<-plot7_data[order(factor(plot7_data$ConservationLevel, levels=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"))),]
        
        #combine all conservation level labels into one for each protein
        Proteinlabel<- aggregate(Label~proteinName, plot7_data, paste, collapse="<br>")
        #get number of protein for labelling
        nProtein<-nrow(Proteinlabel)
        #plotting
        plot7<-ggplot(data) +
          geom_boxplot(aes(x=level,y=index.incidence),outlier.shape=NA,width=0.5)+ 
          geom_jitter(aes(x=level,y=index.incidence,col=ConservationLevel),position = position_jitter(width = .15, height=-0.7),
                      size=input$line_dot_size)+
          labs(x=NULL,y="Index Incidence (%)\n",fill="Conservation level")+ #, title = unique(data$host)
          scale_y_continuous(breaks = c(0,25,50,75,100),labels=c("0","25","50","75","100"),limits = c(0,110))+
          theme_classic(base_size = input$wordsize)+
          theme(
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.position="bottom"
          )+
          scale_colour_manual('Conservation Level',breaks=c("Completely conserved (CC)","Highly conserved (HC)","Mixed variable (MV)","Highly diverse (HD)","Extremely diverse (ED)"),
                              values = c("Completely conserved (CC)"="black","Highly conserved (HC)"="#0057d1","Mixed variable (MV)"="#02d57f","Highly diverse (HD)"="#8722ff", "Extremely diverse (ED)"="#ff617d")) +  
          geom_richtext(data = Proteinlabel, aes(x=proteinName,label = Label, y=c(rep(105,nProtein)), label.size=0, label.color="transparent"),
                        position = position_dodge(width=0.1),size=((input$wordsize/2)-2),color="black", fill="white",hjust=0,angle=90) + 
          guides(color = guide_legend(override.aes = list(size = 2),nrow=2))+
          theme(plot.margin = unit(c(5, 1, 1, 1), "lines"),axis.text.x = element_text(angle = 55, vjust = 0.5, hjust=0.5)) +
          coord_cartesian(clip = "off") #allow ggtext outside of the plot
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
    
    output$plot7_download <- downloadHandler(
      filename = function() { paste("plot7", '.jpg', sep='') },
      content = function(file) {
        ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
      })
    
  })
  
}
