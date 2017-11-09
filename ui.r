
filenames <- list.files(path = "data",pattern="\\.txt$")
names(filenames) <- gsub(pattern = "\\.txt$", "", filenames)
shinyUI(fluidPage(theme = "bootstrap.css",
                  (navbarPage("ShinyMDE BETA Version",position = c("fixed-top"), fluid=TRUE,  ##position = c("fixed-top"), fluid=TRUE
                  navbarMenu("Help",
                             
                             tabPanel(
                               a("List of GPLs available",
                                 target="_blank",href="GPLs.pdf")),
                             
                             tabPanel(
                               a("Reference Manual",
                                 target="_blank", href="UserGuide.pdf")))
                  )),
                  
                  br(),
                  titlePanel( 
                    headerPanel( 
                      h3("Shiny Application for Identifying differentially expressed genes", align="center", style="bold"))),
                  
                  br(),
                  br(),
                  
                  sidebarLayout(
                    
                    sidebarPanel(
                      h5("Upload Data Files",style="bold"),
                      fileInput("files", "Choose CSV/txt processed files or raw files", multiple = "TRUE",
                                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv','.cel','.TXT','.txt')),
                      
                      selectInput('dataset',"Choose platform annotation file", c("Please select a file" ='',filenames)),
                      
                      fluidRow(
                        column(5,
                               radioButtons("radio", label = h5("Data uploaded"),
                                            choices = list("Affymetrix" = 1, "Illumina" = 3),selected = NULL))),
                        br(),
                      fluidRow(
                        column(12,
                               radioButtons("radio1", label = h5("Choose a meta-analysis method"),
                                            choices = list(
                                              "Fisher" = 1, "Stouffer" = 2,
                                              "AWFisher" = 3, "minP" = 4,  
                                              "maxP" = 5, "SR" = 6,
                                              "PR" = 7,"minP (One Sided Correction)" = 8,
                                              "maxP(One Sided Correction)" = 9,
                                              "Fisher (One Sided Correction)" = 10,
                                              "Stouffer (One Sided Correction)" = 11),
                                            selected = NULL))),
                        br(),
                        #column(7,
                               #checkboxInput("checkbox", label = "Differential Expression", value = FALSE))),
                        #br(),
                        
                      #fluidRow(
                        #h4("Affymetrix Differential Expression :"),
                        #column(6,
                               #numericInput("num", label=h5("Fold Change"), value=NULL, min = 0.5, max = 5, step = 0.5)),
                       # column(6,
                              # numericInput("num1", label=h5("P-Value"), value=NULL, min = 0.01, max = 0.05, step = 0.01))
                       # ),
                      
                     # fluidRow(
                        #column(6,
                               #textInput("text", label = h5("Make Contrasts"), value = "Enter conditions.."))
                        #),
                      #br(),
                      
                      #fluidRow(
                        #h4("Codlink Differential Expression :"),
                        #column(6,
                                #textInput("text1", label = h5("Make Contrasts"), value = "Enter conditions..")),
                        #column(6,
                               #selectizeInput("select", label = h5("Choose Annotation"),
                                              #choices = list("h10kcod.db"="h10kcod","h20kcod.db"="h20kcod","hwgcod.db"="hwgcod"),
                                              #options = list(placeholder = 'Choose DB',
                                                             #onInitialize = I('function() { this.setValue(""); }')))

                          #)),
                      
                      #fluidRow(
                    
                        #column(6,
                               #numericInput("num3", label=h5("Fold Change"), value=NULL, min = 0.5, max = 5, step = 0.5)),
                        #column(6,
                               #numericInput("num4", label=h5("P-Value"), value=NULL, min = 0.01, max = 0.05, step = 0.01))
                        #),
                      #br(),
                      
                      #fluidRow(
                        #h4("Illumina Differential Expression :"),
                        #column(6,
                               #textInput("text3", label = h5("Make Contrasts"), value = "Enter conditions..")),
                        #column(6,
                               #numericInput("num5", label=h5("Fold Change"), value=NULL, min = 0.5, max = 5, step = 0.5))            
                      #),
                      #fluidRow(
                        #column(6,
                               #numericInput("num6", label=h5("P-Value"), value=NULL, min = 0.01, max = 0.05, step = 0.01))
                      
                      #),
                        
                      #fluidRow(
                        #br(),
                      
                     column(5,
                            actionButton("Submit", label = "Submit")),
                      
                      
                      br(),
                      br(),
                      titlePanel(
                       
                        fluidRow(
                          column(5,
                                 h5("Downloads", style="bold"),
                                 downloadButton('downloadData1', 'Summary')),
                          br(),
                                                  
                          column(7,
                                 downloadButton("downloadData2",'Supplementary')))),
                          
                      
                      br(),
                      br(),
                       
                      h4("Development:",style = "color:#0000"),
                      p("The application is developed in R.",
                        "This version deals with",
                        a("Affymetrix", 
                          href = "http://en.wikipedia.org/wiki/Affymetrix"),
                        ",",
                        a("Codelink",
                          href="http://www.appliedmicroarrays.com/back_up/pdf/CodeLink_Human_Whole_Genome_Bioarray.pdf"),
                        "and",
                        a("Illumina",
                          href= "http://www.illumina.com/applications/transcriptome-analysis
                          /gene-expression-analysis/gene-expression-arrays.html"),
                        "gene expression data ."),
                      p("Each of the three data platform type has three functionalities, which cover meta-analysis of
                        processed data, processing raw data to perform meta-analysis and finally differential expression"),
                         
                      br(),
                      br(),
                      
                      img(src = "bigorb.png", height = 100, width = 100),
                      br(),
                      h5("Powered by:", 
                         span("HLS and team", colour="Yellow"))),
                    
                    
                    mainPanel(
                      tabsetPanel(
                        
                        tabPanel("File-list",dataTableOutput("file")),
                        tabPanel("Source-data", dataTableOutput("sourced")),
                        tabPanel("Annotation-data",dataTableOutput("annotation")),
                        tabPanel("Summary",dataTableOutput("final")),
                        tabPanel("Supplementary file", dataTableOutput("full"))
                        #tabPanel("Files", dataTableOutput("files")) 
                      )
                    )
                    )))