# edited by Ying on 10/10/2019
library(shiny)
#require(rCharts)
library(DT)
library(plotly)
#library(shinyBS)
#library(shinyIncubator)
#library(D3TableFilter)
library(shinyMatrix)
library(rhandsontable)
# Assign shinyUI as a function (same as in shinyServer)
bincode <- function(n)do.call(expand.grid,rep(list(0:1),n))
# Creates shinyUI, navbarPage creates the overall layout
shinyUI <- navbarPage("Rust Fight",
                      # For each of the sections a tabPanel is used to divide the interface
                      tabPanel("Inputs",
                               tags$head(
                                 tags$style(HTML("
                                                 .shiny-output-error-validation {
                                                 color: green;
                                                 }
                                                 "))
                                 ),
                               
                               #mainPanel(
                               # tabsetPanel further divides the interface into smaller
                               # subsections (here in inputs)
                               tabsetPanel(
                                 # Tried to include a picture of wheat rust, did not quite work
                                 # tabPanel("title page", img(src="wr1.png", style="float:right; padding-right:25px")),
                                 # First "subsection"
                                 tabPanel("Choose n, m, t",
                                          sidebarLayout(sidebarPanel(
                                            # creates helpful information texts for the user
                                            # Have to be further developed
                                            helpText("Choose virulence loci n to be simulated below."),
                                            helpText("(A standard value of 3 loci will result in 8 virus phenotypes)"),
                                            # Slider input with name of the input, description, and further information
                                            sliderInput('loci.n','No. of virulence loci involved in the simulation (n)',min = 1, max=8, step = 1, value=3),
                                            tags$hr(),
                                            helpText("Choose host genotypes m to be simulated below."),
                                            helpText("(In this scenario the hosts are wheat cultivars.)"),
                                            sliderInput('host.m','No. of cultivars that can be deployed (m)',min = 1,max=10,step = 1,value=4),
                                            tags$hr(),
                                            helpText("Choose seasons t to be simulated below."),
                                            helpText("(One season equals to one year in terms of planting new cultivars)"),
                                            sliderInput('seasons.t','No. of seasons to be simulated (t)',min = 1, max=35, step = 1,value=10)),
                                            # This main panel shows the ixV matrix
                                            mainPanel(
                                              helpText("n loci will have 2^n phenotype.
                                                       Here you can see all combinations of virulences and avirulences"),
                                              uiOutput("ixV"), width = 15))
                                          ),
                                 
                                 tabPanel("Initial rel. path. phenotype freq.",
                                          sidebarLayout(
                                            sidebarPanel(
                                              radioButtons(
                                                "distr",
                                                label = "Choose distribution type of the initial pathogen phenotype frequencies",
                                                c("Exponential" = "exp",
                                                  "Uniform" = "unif",
                                                  "Quadratic exponential" = "qexp",
                                                  "Double exponential" = "exp2"
                                                ))
                                            ),
                                            mainPanel(
                                              # With h4 you can create a larger text
                                              # h4("Initial relative pathogen phenotype frequency"),
                                              # h4("This is the initial vector of [...]"),
                                              helpText("Initial relative pathogen phenotype frequency."),
                                              helpText("This is the vector of the initial distribution of cultivar phenotypes."),
                                              uiOutput("fi"))),
                                          fluidRow(column(width = 4)), column(width = 6, offset = 2)),
                                 tabPanel("Fitness costs",
                                          sidebarLayout(
                                            sidebarPanel(
                                              # Here you create the input
                                              uiOutput("fitnesscostInputUi"),
                                              helpText("Please choose numbers from the interval [0,1].")
                                              
                                            ),
                                            mainPanel(
                                              # This output indicates how many numbers are left to input/delete
                                              verbatimTextOutput("fitness.vector.count")
                                            )
                                          )),
                                 
                                 tabPanel("Aggressiveness",
                                          sidebarLayout(
                                            sidebarPanel(
                                              uiOutput("ai.vector.inputUi"),
                                              helpText("Please choose numbers from the interval [0,1].
                                                       A higher value indicates a higher aggressiveness for the respective virulence.")
                                              ),
                                            mainPanel(
                                              verbatimTextOutput("ai.vector.count")
                                            )
                                            )),
                                 
                                 tabPanel("Partial resistance vector",
                                          sidebarLayout(
                                            sidebarPanel(
                                              uiOutput("pj.vector.inputUi"),
                                              helpText("Please choose numbers from the interval [0,1].")
                                            ),
                                            mainPanel(
                                              verbatimTextOutput("pj.vector.count")
                                            )
                                          )),
                                 
                                 
                                 tabPanel("Resistance by cultivar",
                                          helpText("This matrix contains information about ..."),
                                          uiOutput("jxr.matrix.inputUi")),
                                 tabPanel("Cultivar area fraction by season",
                                          helpText("This matrix contains the cultivar area fraction by season.
                                                   Sum of row entries has to add up to 1."),
                                          uiOutput("sjxt.matrix.inputUi")),
                                 
                                 tabPanel("Relative fitness by pathogen phenotype and cultivar",
                                          helpText("This matrix contains information about ..."),
                                          uiOutput("ujxi"), width = 15)
                                 )
                                 ),
                      
                      
                      
                      # Next section
                      # Could be switched with the section above
                      tabPanel('File Input',
                               sidebarLayout(
                                 sidebarPanel(
                                   helpText('Input an excel file with parameters in different sheets'),
                                   radioButtons('readFileInputs',"Input Options", c('Use Shiny input','Upload from PC'), selected = 'Use Shiny input'),
                                   conditionalPanel( condition = "input.readFileInputs == 'Upload from PC'",
                                                     fileInput('file1', 'Choose file to upload',
                                                               accept = c(
                                                                 '.xlsx')))
                                 ) ,
                                 # Empty main panel, maybe show the inputs, or some calculations
                                 mainPanel(
                                   tableOutput('sample.parameter'),
                                   DT::dataTableOutput("fileParameter")
                                   #tableOutput("mlgSummaryTable"),
                                 )
                               )
                      ),
                      
                      # Last section
                      
                      tabPanel("Output",
                               sidebarPanel(
                                 # Creates the download button
                                 h6('Download a pdf report here'),
                                 downloadButton('report'), width = 2
                               ),
                               mainPanel(
                                 tabsetPanel(
                                   tabPanel("Cultivar area fraction resp. mean fitness",
                                            sidebarLayout(
                                              sidebarPanel(
                                                helpText("This plot describes the area fractions of each host genotype (cultivar) planted in the model.
                                                         The dashed black line corresponds to the mean fitness."), width = 12),
                                              mainPanel(plotOutput("plot1"), width = 12), position = "right")),
                                   tabPanel("Relative virulence freq. resp. mean fitness",
                                            sidebarLayout(sidebarPanel(
                                              helpText("This plot describes the relative frequencies of virulence properties.
                                                       The dashed black line corresponds to the mean fitness."), width = 12),
                                              mainPanel(plotOutput("plot2"), width = 12), position = "right")),
                                   tabPanel("Relative pathogen phenotype freq. resp. mean fitness ",
                                            sidebarLayout(
                                              sidebarPanel(
                                                helpText("This plot describes relative virulence frequencies resp. to the mean fitness.
                                                         The dashed black line corresponds to the mean fitness."), width = 12),
                                              mainPanel(plotOutput("plot3"), width = 12), position = "right"))
                                            )
                                   )
                               
                                   ))


