library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT) #https://stackoverflow.com/questions/32149487/controlling-table-width-in-shiny-datatableoutput
library(shinyjs)
library(shinycssloaders)
shinyThings::radioSwitchButtons_default_style(selected_background = "#265071")  #00589a

#https://stackoverflow.com/questions/31703241/activate-tabpanel-from-another-tabpanel
css <- "
.nav li a.disabled {
background-color: #dddddd !important;
color: #474747 !important;
cursor: not-allowed !important;
border-color: #dddddd !important;
}

/*
          box status
      */
      
      .box{box-shadow: 0 0px 0px rgb(0 0 0 / 10%)}
      
      /*change primary status properties*/
      .box.box-solid.box-primary{
        border-bottom-color:#265071;
        border-left-color:#265071;
        border-right-color:#265071;
        border-top-color:#265071;
      }
      
      /*change properties of primary status header*/
      .box.box-solid.box-primary>.box-header {
      color:#fff;
       background:#265071
      }
      
      //change box header properties
      .box-header.with-border {
        border-bottom: 0px solid #f4f4f4;
      }
    
      .btn-primary{
        background:#265071
      }
    
      /*
          navigation bar modification
      */
      
      /*change logo background color*/
      .skin-blue .main-header .logo {
        background-color: #265071;
      }
      
      /*change header background color*/
      .skin-blue .main-header .navbar {
        background-color: #265071;
      }
      
      /*change header font properties*/
      .myClass { 
        font-size: 25px;
        line-height: 50px;
        text-align: center;
        font-family: Arial;
        padding: 0 95px;
        overflow: hidden;
        color: white;
      }
      
      /*
          modification on elements in all pages
      */
      
      /*anchor tag right-aligned*/
      .a {
        text-align: right;
      }
      
      .wrapper {
      height: auto  !important; 
      position:relative; 
      overflow-x:hidden; 
      overflow-y:auto}
      .box box-solid{
      
      }
      
      /*set table width to 50%*/
      .tablewidth{
        width: 50%;
      }
      
      /*set container properties*/
      #container{
        background-color:#ecf0f5;
      }
      
      /*set font family to Arial*/
      * {font-family: \"Arial\"};
"

jscode <-"
shinyjs.disableTab = function(name) {
  var tab = $('.nav-tabs li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

sideBar<-dashboardSidebar(
  width = 320,
  sidebarMenu(
    id="tabs",
    menuItem("Project Description", tabName = "description", icon = icon("r-project")),
    menuItem("Input Data Description", tabName = "inputdata_description", icon = icon("info-circle")),
    menuItem("Entropy and incidence of total variants", tabName = "plot1", icon = icon("area-chart")),
    menuItem("Correlation of entropy", tabName = "plot2", icon = icon("chart-line")),
    menuItem("Dynamics of diversity motifs (Proteome)", tabName = "plot3", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Dynamics of diversity motifs (Protein)", tabName = "plot4", icon = icon("area-chart")),
    menuItem("Distribution of conservation levels", tabName = "plot7", icon = icon("bar-chart")),
    menuItem("Help Page", tabName = "helppage", icon = icon("question")),
    br(),
    
    # Horizontal line ----
    tags$hr(),
    
    #Protein Names in Order
    textInput(inputId = "proteinOrder", 
              label="Protein Names in Order (Plotting)", 
              placeholder="Core,NS3"),
    #tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
    #line / dot size
    tags$style(type = "text/css", 
               ".irs-grid-text:nth-child(n) {color: white}",
               ".irs-grid-pol:nth-of-type(n) {background: white}"),
    
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
    div(style="display:inline-block",actionButton("samplesubmit","Load Sample Dataset",icon("samplesubmit", id="UpdateAnimate", class=""), onclick="www/HCV_protein.csv")),
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
    shinyjs::extendShinyjs(text = jscode, functions = c("disableTab","enableTab")),
    shinyjs::inlineCSS(css),
    tags$title("DiveR"),
    tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "DiveR_logo.png"),
    tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass"> DiveR: Diversity dynamics Visualization in R </span>\');
      })
     ')),
    tags$head(tags$style(type="text/css", '
            .loading {
                display: inline-block;
                overflow: hidden;
                height: 1.3em;
                margin-top: -0.3em;
                line-height: 1.5em;
                vertical-align: text-bottom;
                box-sizing: border-box;
            }
            .loading.dots::after {
                text-rendering: geometricPrecision;
                content: "⠋\\A⠙\\A⠹\\A⠸\\A⠼\\A⠴\\A⠦\\A⠧\\A⠇\\A⠏";
                animation: spin10 1s steps(10) infinite;
                animation-duration: 1s;
                animation-timing-function: steps(10);
                animation-delay: 0s;
                animation-iteration-count: infinite;
                animation-direction: normal;
                animation-fill-mode: none;
                animation-play-state: running;
                animation-name: spin10;
            }
            .loading::after {
                display: inline-table;
                white-space: pre;
                text-align: left;
            }
            @keyframes spin10 { to { transform: translateY(-15.0em); } }
            ')),
    tags$script("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';"), #hide sidebar toggle
    tabItems( 
      tabItem(tabName = "description",
              h2("Dissecting the dynamics of viral protein sequence diversity"),
              fluidRow(
                box(title="Project Description",width = 12,status = "primary", solidHeader = TRUE,
                    HTML('Sequence diversity, as a result of various evolutionary forces, \
                challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. \
                 The publicly available tool, Diversity Motif Analyser\
                (DiMA; <a href=\'https://github.com/PU-SDS/DiMA\'>https://github.com/PU-SDS/DiMA</a>) was developed to facilitate the dissection of sequence diversity dynamics for viruses. \
                Herein, we present DiveR, a DiMA wrapper implemented as a web-based application \
                                  to ease the visualization of outputs from DiMA. DiveR allows visualization of the diversity motifs\
                         (index and its variants – major, minor and unique) for elucidation of the underlying inherent dynamics.\
                         '),
                    div(img(src='glossary_bold.jpg' ,height='auto',width="60%", align = "center"), style="text-align: center;"),
                    tags$ol(
                      tags$li(HTML("Diversity motif: a signature term used to refer to the index or any of its variants (major, minor, unique), as well as distinct variant. \
                            Collectively these terms elucidate the inherent patterns of sequence diversity.")), 
                      tags$ol(
                        HTML("(i) Index: the sequence with the highest incidence at a given <i>k</i>-mer position in a protein alignment <br>\
                        (ii) Major: the predominant sequence(s) amongst the variants <br>\
                        (iii) Minor: distinct sequences with frequency lesser than the major variant, but occur more than once <br>\ 
                             (iv) Unique: distinct sequences that occur only once <br>")), 
                      tags$li("Total variants: sequences that are variant to the index; comprises the major variant, minor variants and unique variants"), 
                      tags$li(HTML("distinct variant incidence: incidence of the distinct <i>k</i>-mer variant"))
                    ))
              )
      ),
      tabItem(tabName = "inputdata_description",
              fluidRow(box(title = "Input Data", width =12,status = "primary", solidHeader = TRUE,
                           sidebarLayout(position = 'right',
                                         sidebarPanel(width=3 ,
                                                      numericInput(inputId = "supportLimit", label = "Minimum Support Threshold", value = 30, min = 0, step =1),
                                                      numericInput(inputId = "kmerlength", label = HTML("<i>k</i>-mer length"), value = 9, min = 0, step =1),
                                                      splitLayout(
                                                        radioButtons(inputId="filetype", label = HTML("Aligned Sequence / DiMA Output File Format  <span style='color:red'>*</span>"),
                                                                     choices = list("FASTA (.fasta/.fas/.fa/.faa/.fnn/.fna)" = 1, "DiMA (.json)" = 2, "DiMA (.csv)"=3), 
                                                                     selected = 1)
                                                      )),
                                         mainPanel(width=9,
                                                   tabsetPanel(id="hostSelection_input",
                                                               
                                                               tabPanel("Description",tags$br(),
                                                                        HTML('DiveR requires either aligned sequence file(s) or DiMA output file(s) as input file(s), \
                                                               where DiveR will convert and concatenate them (the inputs) into a single CSV file. This CSV file will act as the source for subsequent data visualisation. \
                                                               Each file is treated as one viral protein. Currently, DiveR accepts FASTA or JSON/CSV \
                                                               files generated using multiple sequence alignment (MSA) tools and DiMA, respectively. <br><br>\
                                                               Parameters such as host number selection (one or two hosts), <i>k</i>-mer size, support threshold, host name, and protein name are defined by the user. So, <b>please assign files of same host under one tab</b>. Users can also manipulate additional plotting parameters: \
                                                               order of protein name, font, line, and dot size.'),
                                                                        tags$br(),tags$br(),
                                                                        div(
                                                                            style = "display: flex;",
                                                                            div(
                                                                                style = "flex: 1;",
                                                                                radioButtons(inputId="host", label = HTML("Number of host"),
                                                                                     choices = list("One Host" = 1, "Two Host" = 2),
                                                                                     selected = 1)
                                                                                           ),
                                                                            div(
                                                                                style = "flex: 1;",
                                                                                radioButtons(inputId="inputtype", label = HTML("Input File Type"),
                                                                                             choices = list("Protein" = 1, "Nucleotide" = 2),
                                                                                             selected = 1)
                                                                            )
                                                                        )
                                                               ),
                                                               tabPanel("First Host",tags$br(),fileInput(inputId = "MSAfile",label = HTML("Aligned Sequences / DiMA Output File(s)"), accept = c(".fa",".faa",".fasta",".fas",".json",".JSON",".csv"), placeholder = "alignedNS1.fa,alignedCore.fa", multiple = TRUE),
                                                                        uiOutput("infilename"),tags$br(),
                                                                        splitLayout(
                                                                          #Protein Names in Order
                                                                          textInput(inputId = "proteinNames", label=HTML("Protein Name(s) in Ascending Order <span style='color:red'>*</span>"), placeholder="Core, NS3"),
                                                                          textInput(inputId = "hostname", label=HTML("Host Name <span style='color:red'>*</span>"), placeholder="Human")
                                                                        ),),
                                                               tabPanel("Second Host",tags$br(),fileInput(inputId = "MSAfile_secondHost",label = HTML("Aligned Sequences / DiMA Output File(s)"), accept = c(".fa",".faa",".fasta",".fas",".json",".csv"), placeholder = "alignedNS1.fa,alignedCore.fa", multiple = TRUE),
                                                                        uiOutput("infilename_secondHost"),tags$br(),
                                                                        splitLayout(
                                                                          #Protein Names in Order
                                                                          textInput(inputId = "proteinNames_secondHost", label=HTML("Protein Name(s) in Ascending Order <span style='color:red'>*</span>"), placeholder="Core, NS3"),
                                                                          textInput(inputId = "hostname_secondHost", label=HTML("Host Name <span style='color:red'>*</span>"), placeholder="Bat")
                                                                        ))
                                                   ))
                           ),
                           
                           # Horizontal line ----
                           tags$br(),
                           tags$hr(style="border-width: 1px;border-color:#265071;margin: 0em 0.1em 1.5em 0.1em;"),
                           #Alert results are ready
                           uiOutput("alert"),
                           div(style="display:inline-block",actionButton("submitDiMA","Submit",icon("submitDiMA", id="submitAnimate", class=""))),
                           div(style="display:inline-block",downloadButton("downloadDiMA","Download")),
                           div(style="display:inline-block",actionButton("reset","Clear")),
                           tags$br()
              )),
              fluidRow(box(title = "DiMA JSON-Converted CSV Output Format", width =12,status = "primary", solidHeader = TRUE,
                           div(img(src='inputFileformat.JPG' ,width="95%", height='auto'), style="text-align: center;"),
                           HTML("<br><br>\
                           <ol>
                           <li>proteinName: name of the protein</li>\
                             <li>position: starting position of the aligned, overlapping <i>k</i>-mer window</li>\
                             <li>count: number of <i>k</i>-mer sequences at the given position</li>\
                             <li>lowSupport: <i>k</i>-mer position with sequences lesser than the minimum support threshold (TRUE) are considered of low support, in terms of sample size </li>\
                             <li>entropy: level of variability at the <i>k</i>-mer position, with zero representing completely conserved</li>\
                             <li>indexSequence: the predominant sequence (index motif) at the given <i>k</i>-mer position</li>\
                             <li>index.incidence: the fraction (in percentage) of the index sequences at the <i>k</i>-mer position </li>\
                             <li>major.incidence: the fraction (in percentage) of the major sequence (the predominant variant to the index) at the <i>k</i>-mer position </li>\
                            <li> minor.incidence: the fraction (in percentage) of minor sequences (of frequency lesser than the major variant, but not singletons) at the <i>k</i>-mer position</li>
                             <li>unique.incidence: the fraction (in percentage) of unique sequences (singletons, observed only once) at the <i>k</i>-mer position</li>\
                             <li> totalVariants.incidence: the fraction (in percentage) of sequences at the  <i>k</i>-mer position that are variants to the index (includes: major, minor and unique variants)</li>\
                             <li> distinctVariant.incidence: incidence of the distinct <i>k</i>-mer peptides at the <i>k</i>-mer position </li>\
                             <li> multiIndex: presence of more than one index sequence of equal incidence </li>\
                             <li> host: species name of the organism host to the virus</li>\
                             <li>highestEntropy.position: <i>k</i>-mer position that has the highest entropy value</li>\
                             <li>highestEntropy: highest entropy values observed in the studied protein</li>\
                             <li>averageEntropy: average entropy values across all the <i>k</i>-mer positions</li>\
                             </ol>
                                ")))
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
                                the minimum entropy value is zero while the maximum entropy value at y-axis is 100."))),
              fluidRow(box(width=10,title="Entropy Table", status = "primary", solidHeader = TRUE,
              p(style="text-align: right;","*values rounded to 2 decimal places"),
              column(12, align="center", dataTableOutput("entropyTable", width="100%"))
              )
              )
      ),
      
      # Second tab content
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
                           column(6, div(shinyThings::radioSwitchButtons("conserv_lvl", choices = c("CCS", "HCS"), selected = "HCS"),
                                         style="float:left;margin-bottom: 10px; margin-top: 7px;")),
                           column(6, div(shinyThings::radioSwitchButtons("table_type", choices = c("csv", "fasta"), selected = "csv"),
                                         style="float:right;margin-bottom: 10px; margin-top: 7px;")),
                           DT::dataTableOutput("plot7_seqs"),
                           downloadButton('conservSeq_download')))
      ),
      tabItem(tabName = "helppage",
              h2("Help Page"),
              fluidRow(box(width=12,title="Contact",solidHeader = TRUE, status="primary",
                           HTML('For technical assistance or bug report, please reach us out via GitHub (<a href=\'https://github.com/pendy05/DiveR\'>https://github.com/pendy05/DiveR</a>). For the general correspondence, please email Dr. Asif M. Khan (<a href=\'asif@perdanauniversity.edu.my\'>asif@perdanauniversity.edu.my</a>, <a href=\'makhan@bezmialem.edu.tr\'>makhan@bezmialem.edu.tr</a>).')                           
              )),
              fluidRow(box(width=12,title="Frequently Asked Questions (FAQs)",solidHeader = TRUE,status="primary",
                           HTML("1. What can I do if the elements in the plot appear to be overlapping each other due to the displayed plot size? <br>\
                       You may want to increase the height and/or width of the plot offered in the download option, based on your need and download the plot.<br><br>\
                       2. Where can I get the source code of these R plots if I would like to modify the code based on my need? <br>\
                       You may visit this GitHub repository (<a href = 'https://github.com/pendy05/DiveR'>https://github.com/pendy05/DiveR</a>) to get the corresponding source codes.<br><br>
                       3. What is the maximum image size (in inches) that can be downloaded?<br>Maximum 50 (H) x50 (W) inches to prevent the common error of specifying dimensions in pixels encountered in R ggsave() function. <br><br>\
                       4. I encountered 'Error in x$clone: attempt to apply non-function' in plot 'entropy and incidence of total variants' when I submit files for two hosts. Other plots work fine. Why does this happen?<br>\
                       DiveR expects the proteins with same protein name have same length (number of positions) across both the hosts to carry out the comparison plot.")
              )
              )
              
      )
      
    )
  )





#client side
ui <- dashboardPage(title = "DiveR",
  dashboardHeader(title = NULL, tags$li(class="dropdown",tags$a(href='https://github.com/pendy05/DiveR',
                                                                tags$img(src='GitHub-lightLogo.png',height = "18px")
  ))),
  sideBar,
  body,
  footer = dashboardFooter(
    left = HTML("&copy;2021-2022, Tok et al. All Rights Reserved.<br>"),
    right = ""
  )
  
)