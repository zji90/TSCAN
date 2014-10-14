######################################################
##                      TSCAN                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

sidebarPanel3 <- function (...) 
{
      div(class = "span3", tags$form(class = "well", ...))
}

shinyUI(      
      pageWithSidebar(
            
            headerPanel('TSCAN: Tools for Single-Cell ANalysis'),
            
            sidebarPanel3(
                  
                  wellPanel(
                        radioButtons("MainMenu","Main Menu",
                                     list("Input dataset"="Input",
                                          "Preprocess"="Preprocess",
                                          "Cell ordering"="Ordering",
                                          "Gene differential analysis and visualization"="Changepoint",
                                          "Miscellaneous"="Miscellaneous",
                                          "About"="About")
                        )
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Input'",
                                   wellPanel(
                                         h5("Input Dataset"),
                                         fileInput('InputFile', 'Choose File'),
                                         p(actionButton("Inputreadin","Read in")),
                                         checkboxInput('Inputheader', 'Header', TRUE),
                                         radioButtons('Inputsep', 'Separator',
                                                      c('Tab'='\t',
                                                        'Space'=' ',
                                                        'Comma(csv)'=',',
                                                        'Semicolon'=';'
                                                      ),
                                                      '\t'),
                                         radioButtons('Inputquote', 'Quote',
                                                      c('None'='',
                                                        'Double Quote'='"',
                                                        'Single Quote'="'"),
                                                      '')
                                   )                                   
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Preprocess'",
                                   wellPanel(
                                         h5("Filter out genes not expressed"),
                                         textInput("Preprocessexpvalcutoff","Expression value cutoff",1),
                                         sliderInput("Preprocessexppercent","Percentage of expressed cells",min=0,max=1,value=0.3,step=0.01,format="##%"),
                                         h5("Filter out genes with little expression variation"),
                                         textInput("Preprocesscvcutoff","coefficient of variance cutoff",1),
                                         h5("Exclude cells"),
                                         uiOutput("Preprocessexcludeui"),
                                         checkboxInput("Preprocesslogtf","Take log of current data",value = T),
                                         conditionalPanel(condition="input.Preprocesslogtf==1",
                                                          radioButtons("Preprocesslogbase","Choose log base",choices = c("2","10","e")),
                                                          textInput("Preprocesslogpseudocount","Enter pseudo count added when taking log",value = 1)
                                         )
                                   )
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Ordering'",
                                   wellPanel(
                                         checkboxInput("Orderinguploadordering","Upload your own cell ordering",value=F),
                                         conditionalPanel(condition="input.Orderinguploadordering=='1'",
                                                          wellPanel(
                                                                h5("Input Dataset"),
                                                                fileInput('OrderingFile', 'Choose File'),
                                                                p(actionButton("Orderingreadin","Read in")),
                                                                checkboxInput('Orderingheader', 'Header', TRUE),
                                                                radioButtons('Orderingsep', 'Separator',
                                                                             c('Tab'='\t',
                                                                               'Space'=' ',
                                                                               'Comma(csv)'=',',
                                                                               'Semicolon'=';'
                                                                             ),
                                                                             '\t'),
                                                                radioButtons('Orderingquote', 'Quote',
                                                                             c('None'='',
                                                                               'Double Quote'='"',
                                                                               'Single Quote'="'"),
                                                                             '')
                                                          )               
                                         ),
                                         
                                         conditionalPanel(condition="input.Orderinguploadordering=='0'",                                                      
                                                          radioButtons("Orderingchoosestep","",list("Step 1: Reduce dimension"="reduction","Step 2: Calculate pseudotime"="ptime","Save results (optional)"="save")),
                                                          conditionalPanel(condition="input.Orderingchoosestep=='reduction'",
                                                                           wellPanel(
                                                                                 radioButtons("Orderingdimredmet","Choose dimension reduction method",c("Principal Component Analysis (PCA)"="PCA","Independent Component Analysis (ICA)"="ICA"))
                                                                           ),
                                                                           helpText("Warning: ICA could be extremely slow for large datasets, use with care!"),
                                                                           sliderInput("Orderingdimredncomp","Choose number of components",min = 2,max = 20,step = 1,value = 2),
                                                                           conditionalPanel(condition="input.Orderingdimredmet=='PCA'",
                                                                                            p("Automatically select optimal dimension for PCA"),
                                                                                            p(actionButton("Orderingdimredoptbut","Select")),
                                                                                            checkboxInput("Orderingshowvarianceplottf","Show Variance proportion plot",value=F)
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.Orderingchoosestep=='ptime'",
                                                                           wellPanel(
                                                                                 #radioButtons("Orderingptimechoosemethod","Choose calculation method",choices=list("Shortest cell distance"="TSCAN","Minimum spanning tree"="Monocle","Principal curve"="PC"))
                                                                                 radioButtons("Orderingptimechoosemethod","Choose calculation method",choices=list("TSCAN"="TSCAN","Monocle"="Monocle"))
                                                                           ),                                                                           
                                                                           uiOutput("Orderingptimeui"),checkboxInput("Orderingptimetrimtf","Trim branch/cell",value=F),
                                                                           uiOutput("Orderingptimetrimui"),
                                                                           checkboxInput("Orderingptimezoomintf","Zoom in plot",value=F),
                                                                           textInput("Orderingptimescale","Scale pseudotime",value=100),
                                                                           uiOutput("Orderingptimezoominui")
                                                          ),
                                                          conditionalPanel(condition="input.Orderingchoosestep=='start'&&input.Orderingptimechoosemethod=='TSCAN'",
                                                                           helpText("Use marker genesets to determine starting point. Average expression value is used for each geneset."),
                                                                           uiOutput("Orderingstartchoosemarkerui"),
                                                                           selectInput("Orderingstartchoosegenetrend","Choose geneset trend",choices = list("Monotone increasing"="increasing","Monotone decreasing"="decreasing","Not clear"="No")),
                                                                           p(actionButton("Orderingstartaddbutton","Add geneset")),
                                                                           uiOutput("Orderingstartincludegenesetui"),
                                                                           checkboxInput("Orderingstartscalegeneset","Scale gene expression levels",value=T)
                                                                           
                                                          ),
                                                          conditionalPanel(condition="input.Orderingchoosestep=='save'",
                                                                           p("Save pseudotime ordering list"),
                                                                           selectInput("Orderingsavepdatatype","Choose file type",choices = c("txt","csv")),
                                                                           p(downloadButton("Orderingsavepdata")),
                                                                           p("Save pseudotime ordering plot"),
                                                                           checkboxInput("Orderingsaveplotparatf","Change titles",value=F),
                                                                           sliderInput("Orderingsaveplotfontsize","Adjust font size",min = 1,max=50,step=1,value=12),
                                                                           uiOutput("Orderingsaveplotparaui"),
                                                                           selectInput("Orderingsaveplottype","Choose plot type",choices = c("pdf","ps")),
                                                                           textInput("Orderingsaveplotfilewidth","Enter plot width (inches)",12),
                                                                           textInput("Orderingsaveplotfileheight","Enter plot height (inches)",12),
                                                                           p(downloadButton("Orderingsaveplot"))
                                                          )
                                         )
                                   )
                                   
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Changepoint'",
                                   wellPanel(
                                         uiOutput("Changpointpdataselectui"),
                                         radioButtons("Changepointchoosemethod","",choices = list("Gene differential expression analysis"="Difexpr","View gene expression"="View")),
                                         conditionalPanel(condition="input.Changepointchoosemethod=='Difexpr'",
                                                          helpText("Likelihood ratio test of comparing GAM and constant fit models. P-values are adjusted for multiple testing using FDR."),
                                                          p(actionButton("Changepointfilterbutton","Calculate adjusted p-value")),
                                                          textInput("Changpointfilterfdrval","Select FDR cutoff",value=0.5),
                                                          radioButtons("Changpointfiltershowresultopt","",choices=list("Show all results"="all","Show filtered results"="filtering")),
                                                          wellPanel(
                                                                helpText("Save results"),
                                                                selectInput("Changepointsavepvaltabletype","Choose file type",choices = c("txt","csv")),
                                                                p(downloadButton("Changpointsavepvaltable"))
                                                          )
                                         ),
                                         conditionalPanel(condition="input.Changepointchoosemethod=='View'",
                                                          uiOutput("Changepointviewgeneselectui"),
                                                          uiOutput("Changepointviewmethodui"),
                                                          wellPanel(
                                                                helpText("Save plots"),
                                                                selectInput("Changepointviewplottype","Choose plot type",choices = c("pdf","ps")),
                                                                textInput("Changepointviewfilewidth","Enter plot width (inches)",12),
                                                                textInput("Changepointviewfileheight","Enter plot height (inches)",12),
                                                                p(downloadButton("Changepointviewsaveplot"))
                                                          )
                                         )
                                   )
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Miscellaneous'",
                                   h5("Tools for comparing different cell orderings"),
                                   helpText("This is an independent tool which does not directly depend on previous analysis"),
                                   wellPanel(
                                         radioButtons("Compareinputopt","",list("Step 1: Input cell sub-population information"="sub","Step 2: Input cell ordering information"="order")),
                                         conditionalPanel(condition="input.Compareinputopt=='sub'",
                                                          h5("Input Subpopulation Dataset"),
                                                          fileInput('ComparesubFile', 'Choose File'),
                                                          p(actionButton("Comparesubreadin","Read in")),
                                                          checkboxInput('Comparesubheader', 'Header', TRUE),
                                                          radioButtons('Comparesubsep', 'Separator',
                                                                       c('Tab'='\t',
                                                                         'Space'=' ',
                                                                         'Comma(csv)'=',',
                                                                         'Semicolon'=';'
                                                                       ),
                                                                       '\t'),
                                                          radioButtons('Comparesubquote', 'Quote',
                                                                       c('None'='',
                                                                         'Double Quote'='"',
                                                                         'Single Quote'="'"),
                                                                       '')            
                                         ),
                                         conditionalPanel(condition="input.Compareinputopt=='order'",
                                                          h5("Input Cell Ordering Dataset"),
                                                          uiOutput("compareordernameui"),
                                                          fileInput('CompareorderFile', 'Choose File'),
                                                          p(actionButton("Compareorderreadin","Read in"), actionButton("Compareorderadddata","Add Ordering Data")),
                                                          checkboxInput('Compareorderheader', 'Header', TRUE),
                                                          radioButtons('Compareordersep', 'Separator',
                                                                       c('Tab'='\t',
                                                                         'Space'=' ',
                                                                         'Comma(csv)'=',',
                                                                         'Semicolon'=';'
                                                                       ),
                                                                       '\t'),
                                                          radioButtons('Compareorderquote', 'Quote',
                                                                       c('None'='',
                                                                         'Double Quote'='"',
                                                                         'Single Quote'="'"),
                                                                       '')
                                         )
                                   )                                                                     
                  )
                  
            ),
            
            mainPanel(
                  
                  uiOutput("showbusybar"),
                  
                  conditionalPanel(condition="input.MainMenu=='Input'",
                                   checkboxInput("Inputshowinstructiontf","Show instructions",value=T),
                                   uiOutput("Inputshowinstructionui"),
                                   uiOutput("Inputshowsummaryui")
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Preprocess'",
                                   uiOutput("Preprocessstatusui")
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Ordering'",
                                   conditionalPanel(condition="input.Orderinguploadordering=='1'",
                                                    checkboxInput("Orderinguploadshowinstructiontf","Show instructions",value=T),
                                                    uiOutput("Orderinguploadshowinstructionui"),
                                                    uiOutput("Orderinguploadshowpdataui")
                                   ),
                                   conditionalPanel(condition="input.Orderinguploadordering=='0'",
                                                    conditionalPanel(condition="input.Orderingchoosestep=='reduction'",
                                                                     plotOutput("Orderingreductionshowplot",width = "800px",height = "800px"),
                                                                     plotOutput("Orderingreductionshowvariance",width = "800px",height = "800px")
                                                    ),
                                                    conditionalPanel(condition="input.Orderingchoosestep=='ptime'",
                                                                     tabsetPanel(
                                                                           tabPanel("Plot",plotOutput("Orderingptimeshowplot",width = "800px",height = "800px")),
                                                                           tabPanel("Pseudotime",dataTableOutput("Orderingptimeshowptime")),
                                                                           tabPanel("trim expression",
                                                                                    h5("List of criterion:"),
                                                                                    tableOutput("trimexprlistshowtable"),
                                                                                    h5("Trimmed cells:"),
                                                                                    textOutput("trimexprshowcelllist"),
                                                                                    h5("Gene expression heatmap:"),
                                                                                    plotOutput("trimexprshowheatmap")                                                                                    
                                                                           )
                                                                     )                                                    
                                                    ),
                                                    conditionalPanel(condition="input.Orderingchoosestep=='start'",                                                    
                                                                     uiOutput("Orderingstartmainui")                                         
                                                    ),
                                                    conditionalPanel(condition="input.Orderingchoosestep=='save'",
                                                                     tabsetPanel(
                                                                           tabPanel("Plot",plotOutput("Orderingsaveshowplot",width = "800px",height = "800px")),
                                                                           tabPanel("Pseudotime",dataTableOutput("Orderingsaveshowptime"))
                                                                     )                                                    
                                                    )
                                   )
                  ),
                  conditionalPanel(condition="input.MainMenu=='Changepoint'",
                                   uiOutput("Changpointmainui")
                  ),
                  conditionalPanel(condition="input.MainMenu=='Miscellaneous'",
                                   uiOutput("Miscshowresultsui")
                  ),
                  conditionalPanel(condition="input.MainMenu=='About'",
                                   p('SCAtool: Single-cell Analysis Tool'),
                                   p('Current Version: 0.99.0'),
                                   p('Release Date: 2014-8-20'),
                                   p('Author: Zhicheng Ji,Hongkai Ji'),
                                   p('Maintainer: Zhicheng Ji <zji4@jhu.edu>'),
                                   p(a("Visit my homepage",href="http://www.biostat.jhsph.edu/~zji4/",target="_blank")),
                                   p(a("Visit web page of our lab",href="http://www.biostat.jhsph.edu/~hji/",target="_blank"))                                   
                  )
                  
            )
            
      ))


