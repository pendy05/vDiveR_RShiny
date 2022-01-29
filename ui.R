library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
sideBar<-dashboardSidebar(
  width = 350,
  sidebarMenu(
    menuItem("Project Description", tabName = "description", icon = icon("area-chart")),
    menuItem("Entropy and incidence of total variants", tabName = "plot1", icon = icon("area-chart")),
    menuItem("Correlation of entropy", tabName = "plot2", icon = icon("th")),
    menuItem("Dynamics of diversity motifs (Proteome)", tabName = "plot3", icon = icon("dashboard")),
    menuItem("Dynamics of diversity motifs (Protein)", tabName = "plot4", icon = icon("th")),
    menuItem("Distribution of conservation levels", tabName = "plot7", icon = icon("th")),
    menuItem("Help Page", tabName = "helppage", icon = icon("th")),
    br(),
    fileInput(inputId = "csvfile",label = "Input CSV File", accept = ".csv", placeholder = "virus.csv"),
    radioButtons(inputId="host", label = HTML("Host Selection"),
                 choices = list("One Host" = 1, "Two Host" = 2),# "Triple host"=3), 
                 selected = 1),
    
    # Horizontal line ----
    tags$hr(),
    
    #Protein Names in Order
    textInput(inputId = "proteinOrder", label="Protein Names in Order", placeholder="Core,NS3"),
    
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
    div(style="display:inline-block",actionButton("submit","Submit")),
    div(style="display:inline-block",actionButton("samplesubmit","Load Sample Dataset", onclick="www/HCV_protein.csv")),
    div(style="display:block; margin: 6px 5px 6px 15px;color: #000000;",downloadButton("downloadData", "Download Sample",style="color: #000000"))
  )
)

body<-## Body content
  dashboardBody(
    tags$head(tags$style(HTML(
      '
      .myClass { 
        font-size: 25px;
        line-height: 50px;
        text-align: center;
        font-family: Arial;
        padding: 0 95px;
        overflow: hidden;
        color: white;
      }
      
      .skin-blue .main-header .logo {
                              background-color: #3c8dbc;
      }
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
      
      .box-header.with-border {
        border-bottom: 0px solid #f4f4f4;
      }
      
      #container{
        background-color:#ecf0f5;
      }
      
      .box{box-shadow: 0 0px 0px rgb(0 0 0 / 10%)}
      
      * {font-family: "Arial"};
      
    '))),
    tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass"> Viral Diversity Plotting </span>\');
      })
     ')),
    tags$script("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';"),
    tabItems(
      
      tabItem(tabName = "description",
              h2("Dissecting the dynamics of viral protein sequence diversity"),
              fluidRow(
                box(title="Project Description",width = 12,status = "primary", solidHeader = TRUE,
                p('Sequence diversity, as a result of various evolutionary forces, \
                challenges the design of diagnostic, prophylactic and therapeutic interventions against viruses. \
                Diversity Motif Analyser (DiMA), a publicly available tool (', a(href = 'https://github.com/PU-SDS/DiMA', 'https://github.com/PU-SDS/DiMA', .noWS = "outside"), .noWS = c("after-begin", "before-end"),') has been developed \
                                  to facilitate the dissection of protein sequence diversity dynamics for viruses, \
                                  the understanding of which is critical to develop effective intervention strategies. \
                                  Herein, we provide an extension to DiMA in the form of an R Shiny application that \
                                  allows plotting of the diversity motifs for elucidation of the underlying inherent dynamics.'),
                    div(img(src='glossary_bold.jpg' ,height='500px', align = "center"), style="text-align: center;"),
                       HTML("1. Diversity motif: a signature term used to refer to the index or any of its variants (major, minor, unique), as well as <i>k</i>-merTypes. \
                            Collectively these terms elucidate the inherent patterns of sequence diversity.<br><br>\
                            a) Index: the sequence(s) with the highest ranked incidence <br>\
                            b) Major: the sequence(s) with the second highest ranked incidence<br>\
                            c) Minor: sequences observed more than once, but of lesser frequency than the major variant<br>\
                            d) Unique: singleton sequences, each observed only once<br>"),
                       HTML("<br>2. Total variants: sequences that are variant to the index; comprises the major variant, minor variants and unique variants<br><br>"),
                       HTML("3. <i>k</i>-merTypes: the distinct <i>k</i>-mers<br>"))
                
                
              ),
              fluidRow(box(title = "Note", width =12,status = "primary", solidHeader = TRUE,
                           HTML("You may utilize the JsonConverter Python script available in our <a href = 'https://github.com/pendy05/Protocol-of-viral-diversity-dynamics'>github</a> \
                                to convert the DiMA JSON output file to this CSV input file format needed in this R shiny app.")
                           )),
              fluidRow(box(title = "Input File Format", width =12,status = "primary", solidHeader = TRUE,
                           div(img(src='inputformat_courierNew_includedHost.jpg' ,width="1300px", height='auto'), style="text-align: center;"),
                           HTML("<br><br>\
                           1. proteinName: name of the protein<br>\
                             2. position: starting position of the aligned, overlapping <i>k</i>-mer window<br>\
                             3. count: number of <i>k</i>-mer sequences at the given position<br>\
                             4. lowSupport: <i>k</i>-mer position with less than 100 sequences (TRUE) are considered of low support, in terms of sample size <br>\
                             5. indexPeptide: the amino acid sequence of the predominant peptide (index motif) at the given <i>k</i>-mer position<br>\
                             6. index.incidence: the fraction (in percentage) of the index sequences at the <i>k</i>-mer position <br>\
                             7. major.incidence: the fraction (in percentage) of the major sequence (the predominant variant to the index) at the <i>k</i>-mer position <br>\
                             8. minor.incidence: the fraction (in percentage) of minor sequences (of frequency lesser than the major variant, but not singletons) at the <i>k</i>-mer position<br>\
                             9. unique.incidence: the fraction (in percentage) of unique sequences (singletons, observed only once) at the <i>k</i>-mer position<br>\
                             10. totalVariants.incidence: the fraction (in percentage) of sequences at the  <i>k</i>-mer position that are variants to the index (includes: major, minor and unique variants)<br>\
                             11. kmerTypes.incidence: incidence of the distinct <i>k</i>-mer peptides at the position <br>\
                             12. multiIndex: presence of more than one index sequence of equal incidence <br>\
                             13. host: species name of the organism host to the virus<br><br><br>\
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
                    position (1-k, 2-(k+1), etc.) of the proteins. The entropy values indicate the level of variability at the corresponding \
                    <i>k</i>-mer positions, with zero representing completely conserved positions (total variants incidence of 0%). Benchmark reference for \
                    values for entropy (black dotted line) and total variants (pink dotted line) are provided. ")))
      ),
      
      # Second tab content
      tabItem(tabName = "plot2",
              HTML("<h2>Relationship between entropy and total variants for <i>k</i>-mer positions of the viral protein(s)</h2>"),
              fluidRow(
                box(
                  width=10,
                  plotOutput("plot2", height = 650, click = "plot2_click"),
                  verbatimTextOutput("info_plot2")),
                box(
                  width=2,title="Download Option", status = "primary", solidHeader = TRUE,
                  numericInput(inputId="height2", label="Height (inch):", value = 8.0),
                  numericInput(inputId="width2", label="Width (inch):", value = 8.0),
                  numericInput(inputId="dpi2", label="DPI:", value = 500),
                  downloadButton('plot2_download')
                )
              )
              
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
                (please refer <b>Project Description</b> section for detailed definition). Nonatypes defines as the distinct \
                <i>k</i>-mers for a given position. The diversity of the position was depicted by the decline of \
                the index incidences (black) and the increase of total variant incidences (pink).")))
      ),
      
      tabItem(tabName = "plot7",
              h2("Conservation levels of viral nonamer positions for each individual protein"),
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
                    conserved (black) (index incidence = 100% ), highly conserved (blue) \
                    (90% <= index incidence < 100%), mixed variable (green) ( 20% <= index incidence < 90%),\
                    highly diverse (purple) (10% <= index incidence < 20%) and extremely diverse (pink) \
                    ( index incidence < 10% ) ")))
      ),
      tabItem(tabName = "helppage",
              h2("Help Page"),
              fluidRow(box(width=12,title="Contact",solidHeader = TRUE, status="primary",
              p('For technical assistance or bug report, please reach us out via github (', a(href = 'https://github.com/pendy05/Protocol-of-viral-diversity-dynamics', 'https://github.com/pendy05/Protocol-of-viral-diversity-dynamics'),').\
              For the general correspondence, please email Dr. Mohammad Asif Khan (asif@perdanauniversity.edu.my, makhan@bezmialem.edu.tr).')                           
              )),
              fluidRow(box(width=12,title="Frequently Asked Questions (FAQ)",solidHeader = TRUE,status="primary",
                       HTML("1. What can I do if the elements in the plot appear to be overlapping each other due to the displayed plot size? <br>\
                       You may want to increase the height and/or width of the plot offered in the download option, based on your need and download the plot.<br><br>\
                       2. Where can I get the source code of these R plots if I would like to modify the code based on my need? <br>\
                       You may visit this github repository ( <a href = 'https://github.com/pendy05/Protocol-of-viral-diversity-dynamics'>https://github.com/pendy05/Protocol-of-viral-diversity-dynamics</a>",") to get the source codes.<br><br>
                       3. What is the maximum image size (in inches) that can be downloaded?<br>Maximum 50x50 inches as mentioned by R ggsave() function, to prevent the common error of specifying dimensions in pixels. <br><br>")
                       )
              )
              
      )
    
  )
  )





#client side
ui <- dashboardPage(
  dashboardHeader(title = NULL, tags$li(class="dropdown",tags$a(href='https://github.com/pendy05/Protocol-of-viral-diversity-dynamics',
                                         tags$img(src='githublogo.png',height = "18px")))),
  sideBar,
   body,
  footer = dashboardFooter(
    left = HTML("&copy;2021-2022, Tok et al. All Rights Reserved.<br>"),
    right = ""
  )
  
)
