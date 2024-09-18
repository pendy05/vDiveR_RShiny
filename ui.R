library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT) #https://stackoverflow.com/questions/32149487/controlling-table-width-in-shiny-datatableoutput
library(shinyjs)
library(shinycssloaders)

shinyThings::radioSwitchButtons_default_style(selected_background = "#265071")  #00589a

#https://stackoverflow.com/questions/31703241/activate-tabpanel-from-another-tabpanel

sideBar<-dashboardSidebar(
  width = 320,
  sidebarMenu(
    id="tabs",
    menuItem("Project Description", tabName = "description", icon = icon("r-project")),
    menuItem("Input Data Description", tabName = "inputdata_description", icon = icon("info-circle")),
    menuItem("Metadata", tabName = "MetaData", icon = icon("th")),
    menuItem("Entropy", tabName = "plotEntropy", icon = icon("area-chart")),
    menuItem("Entropy And Incidence Of Total Variants", tabName = "plot1", icon = icon("area-chart")),
    menuItem("Correlation Of Entropy", tabName = "plot2", icon = icon("chart-line")),
    menuItem("Dynamics Of Diversity Motifs (Proteome)", tabName = "plot3", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Dynamics Of Diversity Motifs (Protein)", tabName = "plot4", icon = icon("area-chart")),
    menuItem("Distribution Of Conservation Levels", tabName = "plot7", icon = icon("bar-chart")),
    menuItem("Help Page", tabName = "helppage", icon = icon("question")),
    br(),
    
    # Horizontal line ----
    tags$hr(),
    
    #Protein Names in Order
    textInput(inputId = "proteinOrder", 
              label="Protein Names in Order (Plotting)", 
              placeholder="Core,NS3"),
    
    sliderInput(inputId = "line_dot_size",
                label = "Line and Dot Size:",
                min = 2,
                max = 15,
                value = 2,
                step=0.5),
    
    #word size
    sliderInput(inputId = "wordsize",
                label = "Font Size:",
                min = 8,
                max = 20,
                value = 11,
                step=0.5),
    
    # Horizontal line ----
    tags$hr(),
    div(style="display:inline-block;",actionButton("start","Start")),
    div(style="display:inline-block;",
      actionButton("samplesubmit","Load Sample Dataset",icon("paper-plane", id="samplesubmit", class=""), onclick="www/HCV_protein.csv")),
    tags$br(),
    div(style="display:inline-block; margin: 6px 5px 6px 15px;color: #000000;",downloadButton("downloadSampleData", "Download Sample",style="color: #000000")),
    div(style="display:inline-block",actionButton("reset1","Clear")),
    #Alert results are ready
    div(style="margin: 6px 5px 6px 15px;",uiOutput("alertSample"))
  )
)

body<-## Body content
  dashboardBody(
    #initiate the usage of shinyjs
    shinyjs::useShinyjs(),
    tags$script(src = "www/functions.js"),
    includeCSS("www/styles.css"),
    
    tags$title("vDiveR"), 
    tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass"> vDiveR: Viral Protein Diversity Dynamics Visualization in R </span>\');
      })
     ')),
    tags$head(
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "vDiveR_logo.png")),
    tabItems( 
      tabItem(
        tabName = "description",
        h2("Dissecting the dynamics of viral protein sequence diversity"),
        fluidRow(
          box(
            title = "Project Description", 
            width = 12, 
            status = "primary", 
            solidHeader = TRUE,
            HTML('
              Sequence diversity, driven by various evolutionary forces, challenges the design of interventions against viruses.
              The publicly available tool, Diversity Motif Analyser (DiMA; <a href="https://github.com/PU-SDS/DiMA">DiMA</a>),
              was developed to explore viral sequence diversity dynamics. vDiveR, a DiMA wrapper, is a web application
              that eases the visualization of diversity motifs (index, major, minor, unique variants) to elucidate underlying sequence dynamics.
            '),
            div(
              img(src = 'glossary_bold.jpg', height = 'auto', width = '60%', align = 'center'),
              style = "text-align: center;"
            ),
            tags$ol(
              tags$li(HTML("Diversity motif: a term referring to the index or any variants (major, minor, unique), which elucidate sequence diversity patterns.")),
              tags$ol(
                HTML("
                  (i) Index: sequence with highest incidence at a given k-mer position<br>
                  (ii) Major: predominant sequence(s) among the variants<br>
                  (iii) Minor: distinct sequences that occur more than once but less than major variants<br>
                  (iv) Unique: sequences that occur only once<br>
                ")
              ),
              tags$li("Total variants: sequences that vary from the index, including major, minor, and unique variants."),
              tags$li(HTML("Distinct variant incidence: frequency of the distinct k-mer variant"))
            )
          )
        )
      ),
      tabItem(
        tabName = "inputdata_description",
        fluidRow(
          box(
            title = "Input Data", width = 12, status = "primary", solidHeader = TRUE,
            sidebarLayout(
              position = 'right',
              sidebarPanel(
                width = 3,
                numericInput("supportLimit", "Minimum Support Threshold", value = 30, min = 0, step = 1),
                numericInput("kmerlength", HTML("<i>k</i>-mer length"), value = 9, min = 0, step = 1),
                splitLayout(
                  radioButtons(
                    "filetype", 
                    HTML("Aligned Sequence / DiMA Output File Format  <span style='color:red'>*</span>"),
                    choices = list(
                      "FASTA (.fasta/.fas/.fa/.faa/.fnn/.fna)" = 1, 
                      "DiMA (.json)" = 2, 
                      "DiMA (.csv)" = 3
                    ),
                    selected = 1
                  )
                )
              ),
              mainPanel(
                width = 9,
                tabsetPanel(
                  id = "hostSelection_input",
                  tabPanel(
                    "Description", tags$br(),
                    HTML('
                      vDiveR requires aligned sequence file(s) or DiMA output file(s) as inputs, 
                      which vDiveR converts and concatenates into a single CSV file for subsequent data visualization. 
                      Each file represents one viral protein. Supported formats: FASTA, JSON, CSV. Users define parameters such as 
                      <i>k</i>-mer size, support threshold, and host/protein names. 
                      <b>Please assign files of the same host under one tab.</b> Additional plotting parameters can be adjusted, such as 
                      font size, line size, and order of protein names.
                    '),
                    tags$br(), tags$br(),
                    div(
                      style = "display: flex;",
                      div(
                        style = "flex: 1;",
                        radioButtons(
                          "host", 
                          HTML("Number of hosts"), 
                          choices = list("One Host" = 1, "Two Hosts" = 2),
                          selected = 1
                        )
                      ),
                      div(
                        style = "flex: 1;",
                        radioButtons(
                          "inputtype", 
                          HTML("Input File Type"), 
                          choices = list("Protein" = 1, "Nucleotide" = 2),
                          selected = 1
                        )
                      )
                    )
                  ),
                  tabPanel(
                    "First Host", tags$br(),
                    fileInput(
                      "MSAfile", 
                      HTML("Aligned Sequences / DiMA Output File(s)"), 
                      accept = c(".fa", ".faa", ".fasta", ".fas", ".json", ".csv"), 
                      placeholder = "alignedNS1.fa, alignedCore.fa", 
                      multiple = TRUE
                    ),
                    uiOutput("infilename"), tags$br(),
                    splitLayout(
                      textInput(
                        "proteinNames", 
                        HTML("Protein Name(s) in Ascending Order <span style='color:red'>*</span>"), 
                        placeholder = "Core, NS3"
                      ),
                      textInput(
                        "hostname", 
                        HTML("Host Name <span style='color:red'>*</span>"), 
                        placeholder = "Human"
                      )
                    )
                  ),
                  tabPanel(
                    "Second Host", tags$br(),
                    fileInput(
                      "MSAfile_secondHost", 
                      HTML("Aligned Sequences / DiMA Output File(s)"), 
                      accept = c(".fa", ".faa", ".fasta", ".fas", ".json", ".csv"), 
                      placeholder = "alignedNS1.fa, alignedCore.fa", 
                      multiple = TRUE
                    ),
                    uiOutput("infilename_secondHost"), tags$br(),
                    splitLayout(
                      textInput(
                        "proteinNames_secondHost", 
                        HTML("Protein Name(s) in Ascending Order <span style='color:red'>*</span>"), 
                        placeholder = "Core, NS3"
                      ),
                      textInput(
                        "hostname_secondHost", 
                        HTML("Host Name <span style='color:red'>*</span>"), 
                        placeholder = "Bat"
                      )
                    )
                  )
                )
              )
            ),
            tags$br(),
            tags$hr(style = "border-width: 1px; border-color: #265071; margin: 0.1em;"),
            uiOutput("alert"),
            div(style = "display:inline-block;", 
              actionButton("submitDiMA", "Submit", icon("paper-plane", id="submitDiMA_icon"))),
            div(style = "display:inline-block;", downloadButton("downloadDiMA", "Download")),
            div(style = "display:inline-block;", actionButton("reset", "Clear")),
            tags$br()
          )
        ),
        fluidRow(
          box(
            title = "DiMA JSON-Converted CSV Output Format", width = 12, status = "primary", solidHeader = TRUE,
            div(DT::dataTableOutput('mainDataSample'), style = "overflow-x: scroll; display: block;"),
            HTML("
              <br><br>
              <ol>
                <li>proteinName: name of the protein</li>
                <li>position: starting position of the aligned, overlapping <i>k</i>-mer window</li>
                <li>count: number of <i>k</i>-mer sequences at the given position</li>
                <li>lowSupport: positions with sequence counts below the support threshold (TRUE)</li>
                <li>entropy: variability at the <i>k</i>-mer position (0 is conserved)</li>
                <li>indexSequence: the predominant sequence (index motif)</li>
                <li>index.incidence: percentage of index sequences</li>
                <li>major.incidence: percentage of major sequences</li>
                <li>minor.incidence: percentage of minor sequences</li>
                <li>unique.incidence: percentage of unique sequences (singletons)</li>
                <li>totalVariants.incidence: percentage of all variant sequences</li>
                <li>distinctVariant.incidence: incidence of distinct <i>k</i>-mer peptides</li>
                <li>multiIndex: more than one index sequence present</li>
                <li>host: species name of the organism</li>
                <li>highestEntropy.position: position with highest entropy</li>
                <li>highestEntropy: highest entropy value</li>
                <li>averageEntropy: average entropy across all positions</li>
              </ol>
            ")
          )
        )
      ),
      tabItem(tabName = "MetaData",
              h2("Sequence Metadata"),
              fluidRow(
                box(title = "Metadata Input Format", width =12,status = "primary", solidHeader = TRUE,
                    div(DT::dataTableOutput("metademo")),
                    HTML("<br><ol>
                          <li>Accession_ID : ID of the protein sequences</li>\
                          <li>Country : country of isolation for the respective sequence based on sequence database source (<i>e.g.</i> GISAID)</li>\
                          <li>Year : date of the sequence was discovered, format: YYYY/MM/DD </li>\
                          </ol> "))
                ),
              fluidRow(box(width=9, title = "Data input in csv format", status = "primary" ,solidHeader = T,
                           fileInput(inputId = "Metafile",label = HTML("Input your metadata with csv format :"), accept = c(".csv"), placeholder = "metadata.csv", multiple = T),
                           textOutput("inmetafilename"),
                           actionButton("submitMeta1","Get metadata",icon("paper-plane", id="submitMeta1"), style="display:inline-block"),
                           actionButton("resetMeta1","Clear", style="display:inline-block")

              )),
              fluidRow(box(width=9, title = "Data input in fasta format", status = "primary" ,solidHeader = T,
                           fileInput(inputId = "Metafasta",label = HTML("input your fasta file for metadata extracting:"), accept = c(".fa",".faa",".fasta",".fas",".json",".JSON"), placeholder = "alignedNS1.fa, alignedCore.fa", multiple = T),
                           radioButtons("MetafastaSource", "fasta from :", c("NCBI"='NCBI', "GISAID EpiCoV"="GISAID"), selected="NCBI"),
                           textOutput("inmetafasta"),
                           actionButton("submitMeta2","Get metadata",icon("paper-plane", id="submitMeta2"), style="display:inline-block"),
                            actionButton("resetMeta2","Clear", style="display:inline-block"),
                           div(DT::dataTableOutput("metademoSee"), style="margin-top:5%;")
              )),
              fluidRow(
                box(
                  width=10, title = "Country", status = "primary" ,solidHeader = T,
                  h4("Geographical Location"),
                  plotOutput("plot_worldmap", height = 500),
                  h4("Geographical Location (Table)"),
                  DT::dataTableOutput("countrytable")),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height_wm", label="Height (inch):", value = 4.0),
                  numericInput(inputId="width_wm", label="Width (inch):", value = 8.0),
                  numericInput(inputId="dpi_wm", label="DPI:", value = 500),
                  downloadButton('plot_worldmap_download', label = 'Download figure'),
                  HTML("<br><br>"),
                  downloadButton('table_worldmap_download', label = 'Download table')
                   )),
              fluidRow(
                box(
                  width=10, title = "Time", status = "primary" ,solidHeader = T,
                  shinyThings::radioSwitchButtons(inputId = 'time_scale',label = "Display y-scale", choices = c("count", "log"), selected = "count"),
                  h4("Date"),
                  plotOutput("plot_time", height = 500),
                  h4("Date (Table)"),
                  DT::dataTableOutput("timetable")),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height_tm", label="Height (inch):", value = 4.0),
                  numericInput(inputId="width_tm", label="Width (inch):", value = 8.0),
                  numericInput(inputId="dpi_tm", label="DPI:", value = 500),
                  downloadButton('plot_time_download', label = 'Download figure'),
                  HTML("<br><br>"),
                  downloadButton('table_time_download', label = 'Download table')
                ))
      ),
      tabItem(tabName = "plot1",
              HTML("<h2>Entropy and incidence of total variants for each aligned <i>k</i>-mer positions of a viral protein(s)</h2>"),
              fluidRow(
                box(
                  width=10,
                  plotOutput("plot1", height = 650)),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height", label="Height (inch):", value = 7.0),
                  numericInput(inputId="width", label="Width (inch):", value = 25.5),
                  numericInput(inputId="dpi", label="DPI:", value = 500),
                  downloadButton('plot1_download')
                )
                
              ),
              fluidRow(box(width = 10,title="Description", status = "primary", solidHeader = TRUE,
                           collapsible = TRUE,
                           HTML("Entropy (black) and incidence of total variants (pink) were measured for each aligned <i>k</i>-mer (<i>k</i> number of amino acids) \
                    position (1-<i>k</i>, 2-(<i>k</i>+1), etc.) of the proteins. The entropy values indicate the level of variability at the corresponding \
                    <i>k</i>-mer positions, with zero representing completely conserved positions (total variants incidence of 0%). Benchmark reference for \
                    values for entropy (black dotted line) and total variants (pink dotted line) are provided. For both individual protein and across proteome, \
                    the minimum entropy value is zero while the maximum entropy value at y-axis is 100. \
                    <i>k</i>-mers with zero entropy are highlighted in light yellow."))),
              fluidRow(box(width=10,title="Entropy Table", status = "primary", solidHeader = TRUE,
              p(style="text-align: right;","*values rounded to 2 decimal places"),
              column(12, align="center", dataTableOutput("entropyTable", width="100%"))
              )
              )
      ),
       tabItem(tabName = "plotEntropy",
              HTML("<h2>Entropy for each aligned <i>k</i>-mer positions of a viral protein(s)</h2>"),
              fluidRow(
                box(
                  width=10,
                  plotOutput("plotEntropy", height = 650)),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height", label="Height (inch):", value = 7.0),
                  numericInput(inputId="width", label="Width (inch):", value = 25.5),
                  numericInput(inputId="dpi", label="DPI:", value = 500),
                  downloadButton('plotEntropy_download')
                )
                
              ),
              fluidRow(box(width = 10,title="Description", status = "primary", solidHeader = TRUE,
                           collapsible = TRUE,
                           HTML("Entropy (black) was measured for each aligned <i>k</i>-mer (<i>k</i> number of amino acids) \
                    position (1-<i>k</i>, 2-(<i>k</i>+1), etc.) of the proteins. The entropy values indicate the level of variability at the corresponding \
                    <i>k</i>-mer positions, with zero representing completely conserved positions. Benchmark reference for \
                    values for entropy (black dotted line) is provided. For both individual protein and across proteome, \
                    the minimum entropy value is zero while the maximum entropy value at y-axis is 100. \
                    <i>k</i>-mers with zero entropy are highlighted in light yellow.")))
      ),
      
      # Second plot content
      tabItem(tabName = "plot2",
              HTML("<h2>Relationship between entropy and total variants for <i>k</i>-mer positions of the viral protein(s)</h2>"),
              fluidRow(
                box(
                  width=10,
                  plotOutput("plot2", height = 650, click = "plot2_click"),
                  uiOutput("info_plot2")),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height2", label="Height (inch):", value = 8.0),
                  numericInput(inputId="width2", label="Width (inch):", value = 8.0),
                  numericInput(inputId="dpi2", label="DPI:", value = 500),
                  downloadButton('plot2_download')
                )
              ),
              fluidRow(box(width = 10,title="Description", status = "primary", solidHeader = TRUE,
                           collapsible = TRUE,
                           HTML("At y-axis, the minimum entropy value is zero while the maximum entropy value is obtained by rounding the highest entropy encountered up to integer.")))
              
      ),
      tabItem(tabName = "plot3",
              h2("Dynamics of diversity motifs of viral proteome"),
              uiOutput("plot3_hosts"),
              fluidRow(uiOutput("plot3_description"))),
      
      tabItem(tabName = "plot4",
              h2("Dynamics of diversity motifs (Protein)"),
              fluidRow(
                box(width=10,
                    plotOutput("plot4", height = 650)),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height4", label="Height (inch):", value = 8.0),
                  numericInput(inputId="width4", label="Width (inch):", value = 8.0),
                  numericInput(inputId="dpi4", label="DPI:", value = 500),
                  downloadButton('plot4_download')
                )
              ),
              
              fluidRow(box(width=10,title="Description", status = "primary", solidHeader = TRUE,
                           collapsible = TRUE,
                           HTML("<i>k</i>-mers (<i>k</i> number of amino acids) are classified into four \
                different motifs, namely index, major, minor and unique, based on their incidences \
                (please refer <b>Project Description</b> section for detailed definition). The diversity of the position was depicted by the decline of \
                the index incidences (black) and the increase of total variant incidences (pink).")))
      ),
      
      tabItem(tabName = "plot7",
              h2(HTML("Conservation levels of viral <i>k</i>-mer positions for each individual protein")),
              fluidRow(
                box(
                  width=10,
                  uiOutput("plot7_hosts")),
                column(2,style="padding-right: 0px;padding-left: 0px;",
                       box(
                         width=12,title="Download Option", status = "primary", solidHeader = TRUE,
                         numericInput(inputId="height7", label="Height (inch):", value = 8.0),
                         numericInput(inputId="width7", label="Width (inch):", value = 8.0),
                         numericInput(inputId="dpi7", label="DPI:", value = 500),
                         downloadButton('plot7_download')
                       ),
                       box(width=12,title="Conservation Level Label",solidHeader = TRUE, status="primary",
                           p("You may choose to have all or only the conservation level labels with values shown in the plot:"),
                           br(),
                           radioButtons("conservationLabel",label=NULL,
                                        choices = list("Full version" = 1, "Simplified version" = 0), 
                                        selected = 1))
                )),
              fluidRow(box(width=10,title="Description", status = "primary", solidHeader = TRUE,
                           collapsible = TRUE,
                           HTML("The <i>k</i>-mer positions of the proteome and the individual proteins were classified as completely \
                    conserved (black) (index incidence = 100%), highly conserved (blue) \
                    (90% <= index incidence < 100%), mixed variable (green) (20% <= index incidence < 90%),\
                    highly diverse (purple) (10% <= index incidence < 20%) and extremely diverse (pink) \
                    (index incidence < 10%)."))),
              
              fluidRow(box(width=10, height=15, style = "overflow-y: scroll;",
                           title="HCS/CCS Sequences", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                           column(2, div(shinyThings::radioSwitchButtons("conserv_lvl", choices = c("CCS", "HCS"), selected = "HCS"),
                                         style="float:left;margin-bottom: 10px; margin-top: 7px;")),
                           column(2, div("Conservation threshold",
                                         style="margin-bottom: 10px; margin-top: 5px;"),
                                  align = "right"),
                           column(2, div(uiOutput("conserv_threshold_box"),
                                         style="margin-bottom: 10px; margin-top: 7px;"),
                                  align = "left"),
                           column(6, div(shinyThings::radioSwitchButtons("table_type", choices = c("csv", "fasta"), selected = "csv"),
                                         style="float:right;margin-bottom: 10px; margin-top: 7px;")),
                           DT::dataTableOutput("plot7_seqs"),
                           downloadButton('conservSeq_download')))
      ),
      tabItem(tabName = "helppage",
              h2("Help Page"),
              fluidRow(box(width=12,title="Contact",solidHeader = TRUE, status="primary",
                           HTML('For technical assistance or bug report, please reach us out via GitHub (<a href=\'https://github.com/pendy05/vDiveR\'>https://github.com/pendy05/vDiveR</a>). For the general correspondence, please email Dr. Asif M. Khan (<a href=\'asif.khan@udst.edu.qa\'>asif.khan@udst.edu.qa</a>).')                           
              )),
              fluidRow(box(width=12,title="Frequently Asked Questions (FAQs)",solidHeader = TRUE,status="primary",
                           HTML("1. What can I do if the elements in the plot appear to be overlapping each other due to the displayed plot size? <br>\
                       You may want to increase the height and/or width of the plot offered in the download option, based on your need and download the plot.<br><br>\
                       2. Where can I get the source code of these R plots if I would like to modify the code based on my need? <br>\
                       You may visit this GitHub repository (<a href = 'https://github.com/pendy05/vDiveR'>https://github.com/pendy05/vDiveR</a>) to get the corresponding source codes.<br><br>
                       3. What is the maximum image size (in inches) that can be downloaded?<br>Maximum 50 (H) x50 (W) inches to prevent the common error of specifying dimensions in pixels encountered in R ggsave() function. <br><br>\
                       4. I encountered 'Error in x$clone: attempt to apply non-function' in plot 'entropy and incidence of total variants' when I submit files for two hosts. Other plots work fine. Why does this happen?<br>\
                       vDiveR expects the proteins with same protein name have same length (number of positions) across both the hosts to carry out the comparison plot.")
              )
              )
              
      )
      
    )
  )

#client side
ui <- dashboardPage(title = "vDiveR",
  dashboardHeader(title = NULL, tags$li(class="dropdown",tags$a(href='https://github.com/pendy05/vDiveR',
                                                                tags$img(src='GitHub-lightLogo.png',height = "18px")
  ))),
  sideBar,
  body,
  footer = dashboardFooter(
    left = htmlOutput("footer_wording"),
    right = ""
  )
  
)
