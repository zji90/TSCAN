######################################################
##                      TSCAN                       ##
##             Interactive User Interface           ##
##                     Server File                  ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################
library(shiny)
library(fastICA)
library(grid)
library(igraph)
library(ggplot2)
library(plyr)
library(combinat)
library(mgcv)
library(gplots)
library(mclust)

options(shiny.maxRequestSize=30000*1024^2)

shinyServer(function(input, output,session) {
      
      output$showbusybar <- renderUI({
            tagList(
                  tags$head(
                        tags$link(rel="stylesheet", type="text/css",href="style.css"),
                        tags$script(type="text/javascript", src = "busy.js")
                  ),
                  
                  div(class = "busy",  
                      p("Calculation in progress.."), 
                      img(src="ajaxloaderq.gif")
                  )
            )        
      })
      
      Maindata <- reactiveValues()
      
      ### Input ###
      observe({
            if (input$Inputreadin > 0)
                  isolate({
                        FileHandle <- input$InputFile
                        if (!is.null(FileHandle)) {
                              tmpdata <- read.table(FileHandle$datapath,header=input$Inputheader,sep=input$Inputsep,quote=input$Inputquote,stringsAsFactors=F,blank.lines.skip=TRUE)
                              Maindata$rawdata <- as.matrix(tmpdata[,-1])
                              row.names(Maindata$rawdata) <- make.names(tmpdata[,1])
                        }
                  })
      })
      
      output$Inputshowrawtable <- renderTable(head(Maindata$rawdata))
      
      output$Inputshowsummaryui <- renderUI({
            if (!is.null(Maindata$rawdata)) {
                  tagList(
                        h5(paste("The dataset contains",nrow(Maindata$rawdata),"numbers of genes and",ncol(Maindata$rawdata),"numbers of cells")),
                        h5("Head of the input file:"),
                        tableOutput("Inputshowrawtable")
                  )
            }
      })
      
      output$Inputshowinstructionui <- renderUI({
            if (input$Inputshowinstructiontf) {
                        wellPanel(
                        h5("Instructions:"),
                        p("Single cell data should be prepared in a matrix-like data format. Each row corresponds to a gene/feature and each column corresponds to a single cell."),
                        p("Notice that each row should have the same number of entries, especially for the header (first row)"),
                        p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to read in different file formats."),
                        h5("A typical example of tab-delimited file format:"),
                        p("Gene Cell1  Cell2  Cell3  Cell4"),
                        p("SOX2 0.455  0.543  0.000  2.188"),
                        p("PAT1 0.231  2.792  1.222  0.000")                      
                        )
            }
      })
      
      ### Preprocess ###
      
      observe({
            if (input$MainMenu == "Preprocess" && !is.null(Maindata$rawdata)) {
                  tmpdata <- Maindata$rawdata
                  if (input$Preprocesslogtf) {
                        if (input$Preprocesslogbase == "2") {
                              tmpdata <- log2(tmpdata+as.numeric(input$Preprocesslogpseudocount))
                        } else if (input$Preprocesslogbase == "10") {
                              tmpdata <- log10(tmpdata+as.numeric(input$Preprocesslogpseudocount))
                        } else if (input$Preprocesslogbase == "e") {
                              tmpdata <- log(tmpdata+as.numeric(input$Preprocesslogpseudocount))
                        }
                  }
                  tmpdata <- tmpdata[rowSums(tmpdata) > 0,]
                  Maindata$fullrawlogdata <- tmpdata
                  clures <- hclust(dist(tmpdata))
                  cluster <- cutree(clures,0.05*nrow(tmpdata))                  
                  aggdata <- aggregate(tmpdata,list(cluster),mean)   
                  aggdata <- aggdata[,-1]
                  tmpdata <- t(apply(aggdata,1,scale))
                  colnames(tmpdata) <- colnames(aggdata)
                  Maindata$fullprocdata <- as.matrix(tmpdata)                  
            }
      })
      
      output$preprocessshowdata <- renderDataTable(head(Maindata$fullprocdata))
      
#       observe({
#             if (input$MainMenu == "Preprocess" && !is.null(Maindata$rawdata)) {
#                   tmpdata <- Maindata$rawdata
#                   if (input$Preprocesslogtf) {
#                         if (input$Preprocesslogbase == "2") {
#                               tmpdata <- log2(tmpdata+as.numeric(input$Preprocesslogpseudocount))
#                         } else if (input$Preprocesslogbase == "10") {
#                               tmpdata <- log10(tmpdata+as.numeric(input$Preprocesslogpseudocount))
#                         } else if (input$Preprocesslogbase == "e") {
#                               tmpdata <- log(tmpdata+as.numeric(input$Preprocesslogpseudocount))
#                         }
#                   }
#                   Maindata$fullrawlogdata <- tmpdata
#                   tmpdata <- tmpdata[rowMeans(tmpdata > as.numeric(input$Preprocessexpvalcutoff)) > as.numeric(input$Preprocessexppercent),]
#                   tmprowcv <- apply(tmpdata,1,sd)/rowMeans(tmpdata)
#                   Maindata$fullprocdata <- tmpdata[tmprowcv > as.numeric(input$Preprocesscvcutoff),]
#             }
#       })
      
#       output$Preprocesshowproctable <- renderDataTable({
#             if (!is.null(Maindata$fullprocdata)) {
#                   tmpdata <- Maindata$fullprocdata[row.names(Maindata$fullprocdata) %in% input$Preprocesschoosegene,,drop=F]
#                   tmpdata <- cbind(row.names(tmpdata),tmpdata)
#                   colnames(tmpdata)[1] <- "GENE"
#                   tmpdata
#             }
#       })
#       
#       output$Preprocesshowgenesummary <- renderPrint(summary(t(Maindata$fullprocdata[row.names(Maindata$fullprocdata) %in% input$Preprocesschoosegene,,drop=F])))
#       
#       output$Preprocesshowcellsummary <- renderPrint(summary(Maindata$fullprocdata))
      
#       output$Preprocessstatusui <- renderUI({
#             if (!is.null(Maindata$procdata)) {
#                   tagList(
#                         p(h5("Filter summary")),
#                         p(paste(nrow(Maindata$fullprocdata),"genes out of",nrow(Maindata$rawdata),"genes are preserved, which is",round(nrow(Maindata$fullprocdata)/nrow(Maindata$rawdata)*100,3),"percent.")),
#                         hr(),
#                         p(h5("Choose genes to display:")),
#                         selectInput("Preprocesschoosegene","",choices = row.names(Maindata$fullprocdata),selected = head(row.names(Maindata$fullprocdata),n=2),multiple = T),                        
#                         p(h5("Expression of selected genes:")),
#                         dataTableOutput("Preprocesshowproctable"),
#                         hr(),
#                         p(h5("Gene level summary table:")),
#                         verbatimTextOutput("Preprocesshowgenesummary"),
#                         p(h5("Cell level summary table (for all retained genes after filtering):")),
#                         verbatimTextOutput("Preprocesshowcellsummary")                        
#                   )
#             }
#       })
      
      ### Cell Ordering ###
      
      observe({
            if (input$MainMenu != "Ordering")
                  updateRadioButtons(session,"Orderingchoosestep","",list("Step 1: Reduce dimension"="reduction","Step 2: Pseudo time reconstruction"="ptime","Save results (optional)"="save"),selected = "reduction")
      })
      
      #upload own cell ordering
      
      observe({
            if (input$Orderingreadin > 0)
                  isolate({
                        FileHandle <- input$OrderingFile
                        if (!is.null(FileHandle)) {
                              Maindata$uploadpdata <- read.table(FileHandle$datapath,header=input$Orderingheader,sep=input$Orderingsep,quote=input$Orderingquote,stringsAsFactors=F,blank.lines.skip=TRUE)
                        }
                  })
      })
      
      output$Orderinguploadshowpdataui <- renderUI({
            if (!is.null(Maindata$uploadpdata)) {
                  tagList(
                        h5("Uploaded Pseudotime and cell ordering:"),
                        dataTableOutput("Orderinguploadshowpdata")    
                  )
            }
      })
      
      output$Orderinguploadshowpdata <- renderDataTable({
            if (!is.null(Maindata$uploadpdata))
                  Maindata$uploadpdata
      })
      
      output$Orderinguploadshowinstructionui <- renderUI({
            if (input$Orderinguploadshowinstructiontf) {
                  tagList(
                        h5("Instructions:"),
                        p("Cell ordering data should contain exactly three columns. First column: cell name (character value). Second column: cell state (numeric value). Third column: pseudo time (numeric value)"),
                        p("Please ensure the cell names should match exactly to the cell names of the gene expression data, otherwise unpredictable errors may occur."),
                        h5("A typical example of tab-delimited file format is as follows:"),
                        p("Cell   State Pseudotime"),
                        p("Cell1  1     0  "),
                        p("Cell2  1     25  "),
                        p("Cell3  2     75  "),
                        p("Cell4  3     100 ")                        
                  )
            }
      })
      
      #Step 1: Dimension Reduction
      
      observe({
            if (input$MainMenu == "Ordering" && !is.null(Maindata$procdata)) {
                  tmpdata <- t(as.matrix(Maindata$fullprocdata))
                  if (input$Orderingdimredmet == "ICA") {
                        set.seed(12345)
                        icares <- fastICA(tmpdata,n.comp=as.numeric(input$Orderingdimredncomp),row.norm=T,method="C")
                        Maindata$fullreduceres <- t(icares$W) %*% t(tmpdata %*% icares$K)
                  } else if (input$Orderingdimredmet == "PCA") {
                        Maindata$fullreducepc <- tmppc <- prcomp(tmpdata,scale=T)
                        Maindata$fullreduceres <- t(tmpdata %*% tmppc$rotation[,1:as.numeric(input$Orderingdimredncomp)])
                  }
            }
      })
      
      observe({
            if (!is.null(input$Orderingdimredoptbut) && input$Orderingdimredoptbut > 0 ) {
                  isolate({
                        sdev <- Maindata$fullreducepc$sdev[1:20]
                        x <- 1:20
                        optpoint <- which.min(sapply(2:10, function(i) {
                              #sum(lm(sdev[1:i]~x[1:i])$residuals^2)+sum(lm(sdev[(i):20]~x[(i):20])$residuals^2)
                              x2 <- pmax(0,x-i)
                              sum(lm(sdev~x+x2)$residuals^2)
                        }))
                        updateSliderInput(session,"Orderingdimredncomp","Choose number of components",value = optpoint + 1) 
                        
                  })
            }
      })
      
      output$Orderingreductionshowplot <- renderPlot({
            if (!is.null(Maindata$fullreduceres))
                  plot(t(Maindata$fullreduceres[1:2,]),main="First two components")
      })
      
      output$Orderingreductionshowvariance <- renderPlot({
            if (input$Orderingdimredmet == "PCA" && !is.null(Maindata$fullreducepc) && input$Orderingshowvarianceplottf) {
                  plot(Maindata$fullreducepc$sdev[1:20],xlab="number of PC",ylab="Standard deviation proportion explained")
                  abline(v=as.numeric(input$Orderingdimredncomp),col="red")
            }
      })
      
      #Step 2: ptime
      
      #trim
      output$Orderingptimetrimui <- renderUI({
            if (input$Orderingptimetrimtf) {
                  wellPanel(
                        radioButtons("Orderingptimetrimmethod","",choices=list("branch/cluster"="branch","cells"="cell","expression values"="expression")),
                        conditionalPanel(condition="input.Orderingptimetrimmethod=='branch'",
                                         selectInput("Orderingptimetrimbranchselect","Select branch/cluster id on original plot",choices = sort(unique(Maindata$scapdata$State)),multiple = T)
                        ),
                        conditionalPanel(condition="input.Orderingptimetrimmethod=='cell'",
                                         selectInput("Orderingptimetrimcellselect","Select cell name",choices = colnames(Maindata$procdata),multiple = T)
                        ),
                        conditionalPanel(condition="input.Orderingptimetrimmethod=='expression'",
                                         helpText("Cells meeting following criterion simultaneously will be trimmed. Refer to 'Trim expression' tab on the main panel."),
                                         selectInput("Orderingexpressiontrimgene","Gene",choices=row.names(Maindata$rawlogdata)),
                                         radioButtons("Orderingexpressiontrimgtlt","",choices=list("greater than"="greater","smaller than"="smaller")),
                                         textInput("Orderingexpressiontrimvalue","Value",value=0),
                                         p(actionButton("Orderingexpressiontrimaddbutton","Add criteria"))                                         
                        ),
                        p(actionButton("Orderingptimetrimbutton","Trim"),actionButton("Orderingptimetrimreset","Reset"))
                  )
            }
      })
      
      observe({
            if (is.null(Maindata$trimlist)) {                  
                  Maindata$reduceres <- Maindata$fullreduceres
                  Maindata$procdata <- Maindata$fullprocdata
                  Maindata$rawlogdata <- Maindata$fullrawlogdata
            } else {
                  Maindata$reduceres <- Maindata$reduceres[,!colnames(Maindata$reduceres) %in% Maindata$trimlist]
                  Maindata$procdata <- Maindata$procdata[,!colnames(Maindata$procdata) %in% Maindata$trimlist] 
                  Maindata$rawlogdata <- Maindata$rawlogdata[,!colnames(Maindata$rawlogdata) %in% Maindata$trimlist] 
            }
      })
      
      observe({
            if (!is.null(input$Orderingexpressiontrimaddbutton) && input$Orderingexpressiontrimaddbutton > 0) {
                  isolate({
                        if (is.null(Maindata$trimexprlist)) {
                              Maindata$trimexprlist <- data.frame(gene=input$Orderingexpressiontrimgene,relationship=input$Orderingexpressiontrimgtlt,value=input$Orderingexpressiontrimvalue,stringsAsFactors = F)
                        } else {
                              Maindata$trimexprlist <- rbind(Maindata$trimexprlist,data.frame(gene=input$Orderingexpressiontrimgene,relationship=input$Orderingexpressiontrimgtlt,value=input$Orderingexpressiontrimvalue,stringsAsFactors = F))
                        }    
                        celllist <- 1:ncol(Maindata$rawlogdata)
                        for (i in 1:nrow(Maindata$trimexprlist)) {
                              tmpgene <- Maindata$trimexprlist[i,1]
                              tmprela <- Maindata$trimexprlist[i,2]
                              tmpvalue <- Maindata$trimexprlist[i,3]
                              if (tmprela == "greater") {
                                    celllist <- intersect(celllist,which(Maindata$rawlogdata[tmpgene,] > tmpvalue))
                              } else {
                                    celllist <- intersect(celllist,which(Maindata$rawlogdata[tmpgene,] < tmpvalue))
                              }                              
                        }
                        Maindata$trimexprcelllist <- celllist
                        
                  })                       
            }            
      })
      
      output$trimexprlistshowtable <- renderTable(Maindata$trimexprlist)
      
      output$trimexprshowcelllist <- renderText(colnames(Maindata$rawlogdata)[Maindata$trimexprcelllist])
      
      observe({
            if (!is.null(input$Orderingptimetrimbutton) && input$Orderingptimetrimbutton > 0) {
                  isolate({
                        if (input$Orderingptimetrimmethod == "branch") {
                              Maindata$trimlist <- Maindata$scapdata$sample_name[Maindata$scapdata$State %in% as.numeric(input$Orderingptimetrimbranchselect)]      
                        } else if (input$Orderingptimetrimmethod == "cell") {
                              Maindata$trimlist <- input$Orderingptimetrimcellselect
                        } else if (input$Orderingptimetrimmethod == "expression") {
                              Maindata$trimlist <- colnames(Maindata$rawlogdata)[Maindata$trimexprcelllist]
                        }
                  })
            }
      })
      
      observe({
            if (!is.null(input$Orderingptimetrimreset) && input$Orderingptimetrimreset > 0){
                  isolate({
                        Maindata$trimlist <- NULL
                        Maindata$trimexprlist <- NULL
                        Maindata$trimexprcelllist <- NULL
                        Maindata$reduceres <- Maindata$fullreduceres
                        Maindata$procdata <- Maindata$fullprocdata
                        Maindata$rawlogdata <- Maindata$fulllograwdata
                  })
            }
      })
      
      output$trimexprshowheatmap <- renderPlot({
            if (!is.null(Maindata$trimexprlist)) {
                  colcolorall <- rep("cyan",ncol(Maindata$rawlogdata))
                  colcolorall[Maindata$trimexprcelllist] <- "blue"
                  if (nrow(Maindata$trimexprlist) == 1) {
                        plot(Maindata$rawlogdata[Maindata$trimexprlist[,1],],col=colcolorall,lty=19,ylab="Expression value")     
                        legend("topleft",legend=c("Trimmed cells","Retained cells"),lty=19,col=c("blue","cyan"))                        
                  } else {
                        heatmap.2(Maindata$rawlogdata[Maindata$trimexprlist[,1],,drop=F],col=bluered,Colv=F,dendrogram="none",trace="none",Rowv=F,ColSideColors=colcolorall,useRaster=T,cexRow=0.7,srtRow=-45)
                        legend("bottomleft",legend=c("Trimmed cells","Retained cells"),lwd=1,col=c("blue","cyan"))                        
                  }                  
            }            
      })
      
      output$Orderingptimeui <- renderUI({                                                    
            if (input$Orderingptimechoosemethod=="Monocle") {
                  tagList(
                        sliderInput("OrderingMonoclepathnum","Choose number of paths",min=1,max=20,step=1,value=3),
                        checkboxInput("OrderingMonoclereversetf","Reverse the ordering",value=F),
                        checkboxInput("OrderingMonocleshow_tree","Show MST",value = T),
                        checkboxInput("OrderingMonocleshow_backbone","Show diameter path (backbone)",value = T),
                        selectInput("OrderingMonoclebackbone_color","Select backbone color",choices = c("black","red","blue","green","yellow")),
                        checkboxInput("OrderingMonocleshow_cell_names","Show cell names",value = T),
                        conditionalPanel(condition="input.OrderingMonocleshow_cell_names=='1'",textInput("OrderingMonoclecell_name_size","Choose the size of cell name labels",value = 3)),
                        checkboxInput("OrderingMonoclemarkertf","Use marker gene to define node size",value=F),
                        conditionalPanel(condition="input.OrderingMonoclemarkertf=='1'",helpText("Node size is defined by gene expression"),uiOutput("OrderingMonoclemarkerui")),
                        checkboxInput("OrderingMonoclerootcelltf","Choose root cell",value=F),
                        conditionalPanel(condition="input.OrderingMonoclerootcelltf=='1'",uiOutput("OrderingMonoclerootcellui")),
                        sliderInput("OrderingMonoclexcomp","Component displayed on x axis",1,as.numeric(input$Orderingdimredncomp),value=1,step=1),
                        sliderInput("OrderingMonocleycomp","Component displayed on y axis",1,as.numeric(input$Orderingdimredncomp),value=2,step=1))                  
            } else if (input$Orderingptimechoosemethod=="TSCAN") {
                  tagList(
                        sliderInput("OrderingTSCANclunum","Choose number of cell clusters",min=1,max=20,step=1,value=3),
                        p(actionButton("OrderingTSCANoptclunum","Use optimal cluster number")),
                        checkboxInput("OrderingTSCANreversetf","Reverse the ordering",value=F),
                        checkboxInput("OrderingTSCANtuneordertf","Manually tune ordering"),
                        conditionalPanel(condition="input.OrderingTSCANtuneordertf=='1'",
                                         selectInput("OrderingTSCANtuneorderchoose","",choices=list("List cluster order"="order","Use marker gene expression"="gene")),
                                         conditionalPanel("input.OrderingTSCANtuneorderchoose=='order'",textInput("OrderingTSCANtuneorder","Ordering should be seperated by comma (For example: 1,3,2)")),
                                         conditionalPanel("input.OrderingTSCANtuneorderchoose=='gene'",helpText("Clusters will be reordered so the averaged marker gene expression of clusters changes monotonically"),selectInput("OrderingTSCANtunegene","Select marker gene",choices = row.names(Maindata$fullrawlogdata)))
                                         ),
                        checkboxInput("OrderingTSCANshow_tree","Show MST",value = F),                        
                        checkboxInput("OrderingTSCANshow_cell_names","Show cell name",value=T),
                        conditionalPanel(condition="input.OrderingTSCANshow_cell_names=='1'",textInput("OrderingTSCANcell_name_size","Choose the size of cell name labels",value = 3)),
                        checkboxInput("OrderingTSCANmarkertf","Visualize marker gene",value=F),
                        conditionalPanel(condition="input.OrderingTSCANmarkertf=='1'",helpText("Node size is defined by gene expression"),uiOutput("OrderingTSCANmarkerui")),
                        sliderInput("OrderingTSCANxcomp","Component displayed on x axis",1,as.numeric(input$Orderingdimredncomp),value=1,step=1),
                        sliderInput("OrderingTSCANycomp","Component displayed on y axis",1,as.numeric(input$Orderingdimredncomp),value=2,step=1) 
                  )
            } else if (input$Orderingptimechoosemethod=="PC") {
                  tagList(
                        checkboxInput("OrderingPCchoosestartcelltf","Specify start cell",value=F),
                        conditionalPanel(condition="input.OrderingPCchoosestartcelltf=='1'",selectInput("OrderingPCchoosestartcell","Choose start cell",choices = colnames(Maindata$procdata))),
                        checkboxInput("OrderingPCchooseanchorcelltf","Specify anchor cell",value=F),
                        conditionalPanel(condition="input.OrderingPCchooseanchorcelltf=='1'",
                                         selectInput("OrderingPCchooseanchorcell","Choose anchor cell",choices = colnames(Maindata$procdata),multiple = T),
                                         textInput("OrderingPCchooseanchorweight","Choose anchor weight (1 for equal weight)",value=50)
                        ),
                        sliderInput("OrderingPCclunum","Choose number of cell clusters",min=1,max=20,step=1,value=3),
                        checkboxInput("OrderingPCshowcellname","Show cell name",value=T)
                  )
            }
      })      
      
      output$Orderingptimezoominui <- renderUI({
            if (input$Orderingptimezoomintf && !is.null(Maindata$reduceres)) {
                  if (input$Orderingptimechoosemethod=="Monocle") {
                        x = as.numeric(input$OrderingMonoclexcomp)
                        y = as.numeric(input$OrderingMonocleycomp)
                        tagList(
                              sliderInput("Orderingptimezoominxaxis","X-axis range",min=min(Maindata$fullreduceres[x,]),max=max(Maindata$fullreduceres[x,]),value=c(min(Maindata$fullreduceres[x,]),max(Maindata$fullreduceres[x,]))),
                              sliderInput("Orderingptimezoominyaxis","Y-axis range",min=min(Maindata$fullreduceres[y,]),max=max(Maindata$fullreduceres[y,]),value=c(min(Maindata$fullreduceres[y,]),max(Maindata$fullreduceres[y,])))
                        )
                  } else if (input$Orderingptimechoosemethod=="TSCAN") {
                        x = as.numeric(input$OrderingTSCANxcomp)
                        y = as.numeric(input$OrderingTSCANycomp)
                        tagList(
                              sliderInput("Orderingptimezoominxaxis","X-axis range",min=min(Maindata$fullreduceres[x,]),max=max(Maindata$fullreduceres[x,]),value=c(min(Maindata$fullreduceres[x,]),max(Maindata$fullreduceres[x,]))),
                              sliderInput("Orderingptimezoominyaxis","Y-axis range",min=min(Maindata$fullreduceres[y,]),max=max(Maindata$fullreduceres[y,]),value=c(min(Maindata$fullreduceres[y,]),max(Maindata$fullreduceres[y,])))
                        )
                  }                  
            }
      })
      
      observe({
            if (!is.null(Maindata$reduceres) && input$Orderingptimechoosemethod=="TSCAN") {
                  input$OrderingTSCANoptclunum
                  isolate({
                        optnum <- suppressWarnings(Mclust(t(Maindata$reduceres),modelNames="VVV"))$G
                        updateSliderInput(session,"OrderingTSCANclunum","Choose number of cell clusters",value=as.numeric(optnum))
                  })
            }
      })
      
      observe({
            if (!is.null(Maindata$reduceres)) {
                  if (input$Orderingptimechoosemethod=="Monocle" && !is.null(input$OrderingMonoclepathnum) && !is.null(Maindata$reduceres)) { 
                        if (!is.null(input$OrderingMonoclerootcell) && input$OrderingMonoclerootcelltf && nchar(input$OrderingMonoclerootcell) != 0) {
                              root_cell <- input$OrderingMonoclerootcell 
                        } else {
                              root_cell <- NULL      
                        }
                        num_paths <- as.numeric(input$OrderingMonoclepathnum)
                        reverse <- input$OrderingMonoclereversetf                            
                        dp <- as.matrix(dist(t(Maindata$reduceres)))
                        gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
                        dp_Monocle <- minimum.spanning.tree(gp)
                        Maindata$dp_Monocle <- dp_Monocle
                        next_node <<- 0
                        res <- pq_helper(dp_Monocle, use_weights = FALSE, root_node = root_cell)
                        cc_ordering <- extract_good_branched_ordering(res$subtree, res$root, dp, num_paths, reverse)
                        colnames(cc_ordering) <- c("sample_name","State","Pseudotime")
                        cc_ordering$Pseudotime <- cc_ordering$Pseudotime/max(cc_ordering$Pseudotime) * as.numeric(input$Orderingptimescale)
                        Maindata$scapdata <- cc_ordering
                  } else if (input$Orderingptimechoosemethod=="TSCAN" && !is.null(input$OrderingTSCANclunum)) {
                        set.seed(12345)
                        res <- suppressWarnings(Mclust(t(Maindata$reduceres),G=as.numeric(input$OrderingTSCANclunum),modelNames="VVV"))
                        clusterid <- apply(res$z,1,which.max)
                        clucenter <- matrix(0,ncol=ncol(t(Maindata$reduceres)),nrow=res$G)
                        for (cid in 1:res$G) {
                              clucenter[cid,] <- colMeans(t(Maindata$reduceres)[names(clusterid[clusterid==cid]),,drop=F])
                        }
                        dp <- as.matrix(dist(clucenter))                                            
                        gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
                        Maindata$dp_TSCAN <- MSTtree <- minimum.spanning.tree(gp)
                        Maindata$clucenter_TSCAN=clucenter
                        allsp <- shortest.paths(MSTtree)
                        longestsp <- which(allsp == max(allsp), arr.ind = T)
                        if (input$OrderingTSCANtuneordertf) { 
                              if (input$OrderingTSCANtuneorderchoose=='order') {                                    
                                    if (input$OrderingTSCANtuneorder=="") {
                                          Maindata$MSTorder_TSCAN <- MSTorder <- get.shortest.paths(MSTtree,longestsp[1,1],longestsp[1,2])$vpath[[1]]      
                                    } else {
                                          tmp <- as.numeric(strsplit(input$OrderingTSCANtuneorder,",")[[1]])
                                          tmp <- tmp[tmp <= nrow(clucenter)]
                                          Maindata$MSTorder_TSCAN <- MSTorder <- tmp            
                                    }                                    
                              } else if (input$OrderingTSCANtuneorderchoose=='gene') {
                                    if (res$G == 1) {
                                          Maindata$MSTorder_TSCAN <- MSTorder <- 1
                                    } else {
                                          meanval <- sapply(1:res$G, function(i) {
                                                mean(Maindata$rawlogdata[input$OrderingTSCANtunegene,clusterid[clusterid==i]])
                                          })
                                          Maindata$MSTorder_TSCAN <- MSTorder <- order(meanval)
                                    }                                    
                              }
                        } else {
                              Maindata$MSTorder_TSCAN <- MSTorder <- get.shortest.paths(MSTtree,longestsp[1,1],longestsp[1,2])$vpath[[1]]      
                        }   
                        
                        trimclu <- setdiff(as.vector(V(MSTtree)),MSTorder)
                        row.names(clucenter) <- paste0("clu",1:nrow(clucenter)) 
                        trimcell <- names(clusterid[clusterid %in% trimclu])
                        reduceres <- t(Maindata$reduceres)[setdiff(row.names(t(Maindata$reduceres)),trimcell),]
                        if (res$G == 1 || length(MSTorder)==1) {
                              TSCANorder <- names(sort(Maindata$reduceres[1,]))
                        } else {
                              TSCANorder <- NULL                              
                              for (i in 1:length(MSTorder)) {
                                    if (i == 1) {
                                          currentcluid <- MSTorder[i]
                                          nextcluid <- MSTorder[i + 1]
                                          currentclucenter <- clucenter[currentcluid,]
                                          nextclucenter <- clucenter[nextcluid,]
                                          difvec <- nextclucenter - currentclucenter
                                          tmppos <- reduceres[names(clusterid[clusterid==currentcluid]),] %*% difvec
                                          pos <- as.vector(tmppos)
                                          names(pos) <- row.names(tmppos)
                                          TSCANorder <- c(TSCANorder,names(sort(pos)))                  
                                    } else if (i == length(MSTorder)) {
                                          currentcluid <- MSTorder[i]
                                          lastcluid <- MSTorder[i - 1]
                                          currentclucenter <- clucenter[currentcluid,]
                                          lastclucenter <- clucenter[lastcluid,]
                                          difvec <- currentclucenter - lastclucenter
                                          tmppos <- reduceres[names(clusterid[clusterid==currentcluid]),] %*% difvec
                                          pos <- as.vector(tmppos)
                                          names(pos) <- row.names(tmppos)
                                          TSCANorder <- c(TSCANorder,names(sort(pos)))   
                                    } else {
                                          currentcluid <- MSTorder[i]
                                          nextcluid <- MSTorder[i + 1]
                                          lastcluid <- MSTorder[i - 1]
                                          currentclucenter <- clucenter[currentcluid,]
                                          nextclucenter <- clucenter[nextcluid,]
                                          lastclucenter <- clucenter[lastcluid,]
                                          clupoints <- names(clusterid[clusterid==currentcluid])
                                          distlast <- rowSums((reduceres[clupoints,]-lastclucenter)^2)
                                          distnext <- rowSums((reduceres[clupoints,]-nextclucenter)^2)
                                          lastpoints <- names(which(distlast < distnext))
                                          nextpoints <- names(which(distlast >= distnext))                                          
                                          difvec <- currentclucenter - lastclucenter
                                          tmppos <- reduceres[lastpoints,] %*% difvec
                                          pos <- as.vector(tmppos)
                                          names(pos) <- row.names(tmppos)
                                          TSCANorder <- c(TSCANorder,names(sort(pos)))                                            
                                          difvec <- nextclucenter - currentclucenter
                                          tmppos <- reduceres[nextpoints,] %*% difvec
                                          pos <- as.vector(tmppos)
                                          names(pos) <- row.names(tmppos)
                                          TSCANorder <- c(TSCANorder,names(sort(pos)))  
                                    }
                              }            
                        }
                        if (input$OrderingTSCANreversetf)
                              TSCANorder <- rev(TSCANorder)
                        datadist <- dist(t(Maindata$reduceres))
                        distmat <- as.matrix(datadist)
                        alldist <- sapply(1:(length(TSCANorder)-1), function(x) {
                              distmat[TSCANorder[x],TSCANorder[x+1]]
                        })
                        ptime <- c(0,cumsum(alldist))
                        ptime <- ptime/max(ptime) * as.numeric(input$Orderingptimescale)
                        if (res$G == 1 || length(MSTorder)==1) {
                              Maindata$scapdata <- data.frame(sample_name=TSCANorder,State=1,Pseudotime=ptime,stringsAsFactors = F)      
                        } else {
                              Maindata$scapdata <- data.frame(sample_name=TSCANorder,State=clusterid[TSCANorder],Pseudotime=ptime,stringsAsFactors = F)      
                        }                        
                  } else if (input$Orderingptimechoosemethod=="PC" && !is.null(input$OrderingPCclunum)) {
                        coord <- t(Maindata$reduceres)
                        w <- rep(1,nrow(coord))
                        if (!is.null(input$OrderingPCchooseanchorcelltf) && input$OrderingPCchooseanchorcelltf && !is.null(input$OrderingPCchooseanchorcell)) {
                              w[row.names(coord) %in% input$OrderingPCchooseanchorcell] <- as.numeric(input$OrderingPCchooseanchorweight)
                        }
                        if (!is.null(input$OrderingPCchoosestartcelltf) && input$OrderingPCchoosestartcelltf && !is.null(input$OrderingPCchoosestartcell)) {
                              lpcobj <- lpc(coord,x0=coord[row.names(coord) %in% input$OrderingPCchoosestartcell,],weights=w)      
                        } else {
                              lpcobj <- lpc(coord,weights=w)      
                        }
                        Maindata$lpcobj <- lpcobj
                        ptime <- lpc.project(lpcobj,coord)$closest.or.pi
                        ptime <- ptime/max(ptime) * as.numeric(input$Orderingptimescale)
                        state <- kmeans(t(Maindata$reduceres),centers=as.numeric(input$OrderingPCclunum))$cluster
                        Maindata$scapdata <- data.frame(sample_name=row.names(coord),State=state,Pseudotime=ptime,stringsAsFactors = F)
                  }
            }
      })
      
      output$OrderingMonoclemarkerui <- renderUI({
            selectInput("OrderingMonoclemarker","select marker gene",choices = row.names(Maindata$fullrawlogdata))
      })
      
      output$OrderingTSCANmarkerui <- renderUI({
            selectInput("OrderingTSCANmarker","select marker gene",choices = row.names(Maindata$fullrawlogdata))
      })
      
      output$OrderingMonoclerootcellui <- renderUI({
            selectInput("OrderingMonoclerootcell","select root cell",choices = colnames(Maindata$procdata))
      })
      
      Monocledrawplot <- function(xlabtext="X",ylabtext="Y",titletext="") {
            x = as.numeric(input$OrderingMonoclexcomp)
            y = as.numeric(input$OrderingMonocleycomp)
            color_by = "State"
            show_tree = input$OrderingMonocleshow_tree
            show_backbone = input$OrderingMonocleshow_backbone
            backbone_color = input$OrderingMonoclebackbone_color
            if (!is.null(input$OrderingMonoclemarkertf) && input$OrderingMonoclemarkertf) {
                  markers <- input$OrderingMonoclemarker
            } else {
                  markers <- NULL
            }
            show_cell_names = input$OrderingMonocleshow_cell_names
            cell_name_size = as.numeric(input$OrderingMonoclecell_name_size)
            
            lib_info_with_pseudo <- Maindata$scapdata
            lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
            S_matrix <- Maindata$reduceres
            ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
            colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
            ica_space_df$sample_name <- row.names(ica_space_df)
            ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
            dp_Monocle <- Maindata$dp_Monocle
            edge_list <- as.data.frame(get.edgelist(dp_Monocle))
            colnames(edge_list) <- c("source", "target")
            edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", ICA_dim_2 = "source_ICA_dim_2"))
            edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", ICA_dim_2 = "target_ICA_dim_2"))
            diam <- as.data.frame(as.vector(V(dp_Monocle)[get.diameter(dp_Monocle, weights = NA)]$name))
            colnames(diam) <- c("sample_name")
            diam <- arrange(merge(ica_space_with_state_df, diam, by.x = "sample_name", by.y = "sample_name"), Pseudotime)
            if (!is.null(markers)) {
                  markers_exprs <- data.frame(value=Maindata$rawlogdata[markers, ])
                  edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", by.y = "row.names")
                  g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, y = source_ICA_dim_2, size = log10(value + 0.1)))
            } else {
                  g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, y = source_ICA_dim_2))
            }
            if (show_tree) {
                  g <- g + geom_segment(aes_string(xend = "target_ICA_dim_1", yend = "target_ICA_dim_2", color = color_by), size = 0.3, linetype = "solid", na.rm = TRUE)
            }
            g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
            if (show_backbone) {
                  g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2), 
                                     color = I(backbone_color), size = 0.75, data = diam, 
                                     na.rm = TRUE) + geom_point(aes_string(x = "ICA_dim_1", 
                                                                           y = "ICA_dim_2", color = color_by), size = I(1.5), 
                                                                data = diam, na.rm = TRUE)
            }
            if (show_cell_names) {
                  g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
            }
            g <- g + guides(colour = guide_legend(override.aes = list(size=5))) + 
                  theme_minimal(base_size = as.numeric(input$Orderingsaveplotfontsize)) + theme(panel.border = element_blank(), axis.line = element_line()) + 
                  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
                  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
                  labs(title=titletext) + ylab(ylabtext) + xlab(xlabtext) + theme(legend.position = "top", 
                                                                                  legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
                  theme(panel.background = element_rect(fill = "white")) +
                  theme(axis.text.x = element_text(size=17,color="darkred"),
                        axis.text.y = element_text(size=17,color='black'),
                        axis.title.x = element_text(size=20,vjust=-1),
                        axis.title.y = element_text(size=20,vjust=1),
                        plot.margin=unit(c(1,1,1,1),"cm"))
            
            if (input$Orderingptimezoomintf) {
                  g <- g + coord_cartesian(xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
            } else {
                  sc1 <- 0.05 * (max(Maindata$fullreduceres[x,])-min(Maindata$fullreduceres[x,]))
                  sc2 <- 0.05 * (max(Maindata$fullreduceres[y,])-min(Maindata$fullreduceres[y,]))
                  g <- g + coord_cartesian(xlim = c(min(Maindata$fullreduceres[x,]) - sc1, max(Maindata$fullreduceres[x,])+sc1),ylim=c(min(Maindata$fullreduceres[y,])-sc2, max(Maindata$fullreduceres[y,])+sc2))                  
            }                  
            g 
      }
      
      TSCANdrawplot <- function(xlabtext="X",ylabtext="Y",titletext="") {
            x = as.numeric(input$OrderingTSCANxcomp)
            y = as.numeric(input$OrderingTSCANycomp)
            color_by = "State"
            
            show_tree = input$OrderingTSCANshow_tree
            show_cell_names = input$OrderingTSCANshow_cell_names
            cell_name_size = as.numeric(input$OrderingTSCANcell_name_size)
            if (!is.null(input$OrderingTSCANmarkertf) && input$OrderingTSCANmarkertf) {
                  markers <- input$OrderingTSCANmarker
            } else {
                  markers <- NULL
            }
            
            lib_info_with_pseudo <- Maindata$scapdata
            lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
            S_matrix <- Maindata$reduceres
            ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
            colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
            ica_space_df$sample_name <- row.names(ica_space_df)
            edge_df <- merge(ica_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
            
            if (!is.null(markers)) {
                  markers_exprs <- data.frame(markerexpr=Maindata$rawlogdata[markers, ])
                  edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", by.y = "row.names")
                  g <- ggplot(data = edge_df, aes(x = ICA_dim_1, y = ICA_dim_2,size=markerexpr))
                  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
            } else {
                  g <- ggplot(data = edge_df, aes(x = ICA_dim_1, y = ICA_dim_2))
                  g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,size=3)
            }
            
            if (show_cell_names) {
                  g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
            }            
            if (show_tree && max(Maindata$scapdata$State) > 1) {
                  clucenter <- Maindata$clucenter_TSCAN[,c(x,y)]
                  clulines <- NULL
                  for (i in 1:(length(Maindata$MSTorder_TSCAN)-1)) {
                        clulines <- rbind(clulines, c(clucenter[Maindata$MSTorder_TSCAN[i],],clucenter[Maindata$MSTorder_TSCAN[i+1],]))
                  }
                  clulines <- data.frame(x=clulines[,1],xend=clulines[,3],y=clulines[,2],yend=clulines[,4])
                  g <- g + geom_segment(aes_string(x="x",xend="xend",y="y",yend="yend",size=NULL),data=clulines,size=1)
                  clucenter <- data.frame(x=clucenter[,1],y=clucenter[,2],id=1:nrow(clucenter))
                  g <- g + geom_text(aes_string(label="id",x="x",y="y",size=NULL),data=clucenter,size=10)                  
            }            
            g <- g + guides(colour = guide_legend(override.aes = list(size=5))) + 
                  theme_minimal(base_size = as.numeric(input$Orderingsaveplotfontsize))+theme(panel.border = element_blank(), axis.line = element_line()) + 
                  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
                  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            labs(title=titletext) + ylab(ylabtext) + xlab(xlabtext) + theme(legend.position = "top", 
                                                                            legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
                  theme(panel.background = element_rect(fill = "white")) +
                  theme(axis.text.x = element_text(size=17,color="darkred"),
                        axis.text.y = element_text(size=17,color='black'),
                        axis.title.x = element_text(size=20,vjust=-1),
                        axis.title.y = element_text(size=20,vjust=1),
                        plot.margin=unit(c(1,1,1,1),"cm"))
            
            if (input$Orderingptimezoomintf) {
                  g <- g + coord_cartesian(xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
            } else {
                  sc1 <- 0.05 * (max(Maindata$fullreduceres[x,])-min(Maindata$fullreduceres[x,]))
                  sc2 <- 0.05 * (max(Maindata$fullreduceres[y,])-min(Maindata$fullreduceres[y,]))
                  g <- g + coord_cartesian(xlim = c(min(Maindata$fullreduceres[x,]) - sc1, max(Maindata$fullreduceres[x,])+sc1),ylim=c(min(Maindata$fullreduceres[y,])-sc2, max(Maindata$fullreduceres[y,])+sc2))
            }  
            g       
      }
      
      output$Orderingptimeshowplot <- renderPlot({
            if (!is.null(Maindata$reduceres)) {
                  if (input$Orderingdimredmet == "ICA") {                        
                        tmpxlabtext <- paste0("ICA_dimension_",as.numeric(input$OrderingTSCANxcomp))
                        tmpylabtext <- paste0("ICA_dimension_",as.numeric(input$OrderingTSCANycomp))
                  } else {
                        tmpxlabtext <- paste0("PCA_dimension_",as.numeric(input$OrderingTSCANxcomp))
                        tmpylabtext <- paste0("PCA_dimension_",as.numeric(input$OrderingTSCANycomp))
                  }                  
                  if (input$Orderingptimechoosemethod=="Monocle" && !is.null(Maindata$dp_Monocle)) {
                        Monocledrawplot(xlabtext=tmpxlabtext,ylabtext=tmpylabtext)
                  } else if (input$Orderingptimechoosemethod=="TSCAN" && !is.null(input$OrderingTSCANshow_cell_names)) {
                        TSCANdrawplot(xlabtext=tmpxlabtext,ylabtext=tmpylabtext)
                  } else if (input$Orderingptimechoosemethod=="PC") {
                        coord <- t(Maindata$reduceres)
                        if (input$Orderingptimezoomintf) {
                              plot(Maindata$lpcobj,datcol=Maindata$scapdata$State,datpch=19,xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
                        } else {
                              plot(Maindata$lpcobj,datcol=Maindata$scapdata$State,datpch=19)
                        }
                        if (input$OrderingPCshowcellname)
                              text(coord,row.names(coord))     
                  }
            }
      })
      
      output$Orderingptimeclustershowplot <- renderPlot({
            if (!is.null(input$OrderingTSCANmarkertf) && input$OrderingTSCANmarkertf) {
            x = as.numeric(input$OrderingTSCANxcomp)
            y = as.numeric(input$OrderingTSCANycomp)
            color_by = "State"
            
            lib_info_with_pseudo <- Maindata$scapdata
            lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
            S_matrix <- Maindata$reduceres
            pca_space_df <- data.frame(t(S_matrix[c(x, y), ]))
            colnames(pca_space_df) <- c("pca_dim_1","pca_dim_2")            
            pca_space_df$sample_name <- row.names(pca_space_df)
            edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")    
            markers_exprs <- Maindata$rawlogdata[input$OrderingTSCANmarker, ]
            edge_df$markerexpr <- markers_exprs[edge_df$sample_name]
            data <- data.frame(pca_dim_1=tapply(edge_df$pca_dim_1,edge_df$State,mean),pca_dim_2=tapply(edge_df$pca_dim_2,edge_df$State,mean),expr=scale(tapply(edge_df$markerexpr,edge_df$State,mean))[,1],state=1:length(unique(edge_df$State)))
            data$state <- as.factor(data$state)
            g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
            g <- g + geom_point(color = "white", na.rm = TRUE)
            
            clucenter <- Maindata$clucenter_TSCAN[,c(x,y)]
            clulines <- NULL
            MSTorder <- Maindata$MSTorder_TSCAN
            for (i in 1:(length(MSTorder)-1)) {
                  clulines <- rbind(clulines, c(clucenter[MSTorder[i],],clucenter[MSTorder[i+1],]))
            }
            clulines <- data.frame(x=clulines[,1],xend=clulines[,3],y=clulines[,2],yend=clulines[,4])
            g <- g + geom_segment(aes_string(x="x",xend="xend",y="y",yend="yend",size=NULL),data=clulines,size=1)
            
            g <- g + geom_point(aes(x=pca_dim_1,y=pca_dim_2,color=expr),data=data,size=10)
            clucenter <- data.frame(x=clucenter[,1],y=clucenter[,2],id=1:nrow(clucenter))
            g <- g + geom_text(aes_string(label="id",x="x",y="y",size=NULL),data=clucenter,size=10)
            g <- g + scale_color_continuous(low="white",high="red")            
            g <- g + guides(colour = guide_legend(override.aes = list(size=5))) + 
                  xlab(paste0("PCA_dimension_",x)) + ylab(paste0("PCA_dimension_",y)) +
                  theme(panel.border = element_blank(), axis.line = element_line()) + 
                  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
                  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
                  theme(legend.position = "top", legend.key.size = unit(0.3, "in"),legend.text = element_text(size = 20),legend.title=element_text(size = 20)) + theme(legend.key = element_blank()) + 
                  theme(panel.background = element_rect(fill = "white")) +
                  theme(axis.text.x = element_text(size=17,color="darkred"),
                        axis.text.y = element_text(size=17,color='black'),
                        axis.title.x = element_text(size=20,vjust=-1),
                        axis.title.y = element_text(size=20,vjust=1),
                        plot.margin=unit(c(1,1,1,1),"cm"))
            if (input$Orderingptimezoomintf) {
                  g <- g + coord_cartesian(xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
            } else {
                  sc1 <- 0.05 * (max(Maindata$fullreduceres[x,])-min(Maindata$fullreduceres[x,]))
                  sc2 <- 0.05 * (max(Maindata$fullreduceres[y,])-min(Maindata$fullreduceres[y,]))
                  g <- g + coord_cartesian(xlim = c(min(Maindata$fullreduceres[x,]) - sc1, max(Maindata$fullreduceres[x,])+sc1),ylim=c(min(Maindata$fullreduceres[y,])-sc2, max(Maindata$fullreduceres[y,])+sc2))
            }  
            g       
            }
      })

      output$Orderingptimeshowptime <- renderDataTable(Maindata$scapdata)
      
      #       #Step 3: starting point
      #       
      #       output$Orderingstartmainui <- renderUI({
      #             if (input$Orderingptimechoosemethod == "TSCAN") {
      #                   tabsetPanel(
      #                         tabPanel("Ordering", uiOutput("Orderingstartshowheatui")
      #                         ),
      #                         tabPanel("Genesets",dataTableOutput("Orderingstartsuggestshowgeneset")),
      #                         tabPanel("Pseudotime",dataTableOutput("Orderingstartsuggestshowpseudotime"))
      #                   )
      #             } else {
      #                   h5("This function is not available for minimum spanning tree approach.")
      #             }
      #       })
      #       
      #       output$Orderingstartchoosemarkerui <- renderUI({
      #             if (!is.null(Maindata$procdata))
      #                   selectInput("Orderingstartchoosemarker","select genes",choices = row.names(Maindata$rawdata),multiple = T)
      #       })
      #       
      #       observe({
      #             if (!is.null(input$Orderingstartaddbutton) && input$Orderingstartaddbutton > 0) {
      #                   isolate({
      #                         if (is.null(Maindata$fullstartgeneset)) {
      #                               Maindata$fullstartgeneset <- data.frame(gene=input$Orderingstartchoosemarker,trend=input$Orderingstartchoosegenetrend,genesetname="Geneset 1",stringsAsFactors = F)
      #                         } else {
      #                               genesetid <- max(as.numeric(sub("Geneset ","",unique(Maindata$fullstartgeneset$genesetname))))+1
      #                               Maindata$fullstartgeneset <- rbind(Maindata$fullstartgeneset,data.frame(gene=input$Orderingstartchoosemarker,trend=input$Orderingstartchoosegenetrend,genesetname=paste("Geneset",genesetid),stringsAsFactors = F))
      #                         }
      #                   })                        
      #             }          
      #       })
      #       
      #       output$Orderingstartincludegenesetui <- renderUI({
      #             selectInput("Orderingstartincludegeneset","Select geneset to include",choices = unique(Maindata$fullstartgeneset$genesetname),multiple = T,selected = unique(Maindata$fullstartgeneset$genesetname))
      #       })
      #       
      #       observe({
      #             if (!is.null(Maindata$fullstartgeneset) && !is.null(input$Orderingstartincludegeneset)) {
      #                   Maindata$startgeneset <- Maindata$fullstartgeneset[Maindata$fullstartgeneset$genesetname %in% input$Orderingstartincludegeneset,]      
      #                   increaseexpr <- NULL
      #                   increasename <- NULL
      #                   decreaseexpr <- NULL
      #                   decreasename <- NULL
      #                   for (i in unique(Maindata$startgeneset$genesetname)) {
      #                         tmp <- Maindata$startgeneset[Maindata$startgeneset$genesetname == i,]             
      #                         if (tmp[1,2] != "No") {
      #                               tmpexpr <- Maindata$rawlogdata[tmp[,1],,drop=F]
      #                               if (input$Orderingstartscalegeneset)
      #                                     tmpexpr <- t(scale(t(tmpexpr)))
      #                               tmpexpr <- colMeans(tmpexpr)
      #                               
      #                               if (tmp[1,2] == "increasing") {
      #                                     increaseexpr <- rbind(increaseexpr,tmpexpr)
      #                                     increasename <- c(increasename, i)
      #                               } else if (tmp[1,2] == "decreasing") {
      #                                     decreaseexpr <- rbind(decreaseexpr,tmpexpr)
      #                                     decreasename <- c(decreasename, i)
      #                               }
      #                         }                        
      #                   }
      #                   row.names(increaseexpr) <- increasename
      #                   row.names(decreaseexpr) <- decreasename
      #                   Maindata$increaseexpr <- increaseexpr
      #                   Maindata$decreaseexpr <- decreaseexpr                                          
      #             }
      #       })
      #       
      #       observe({
      #             if (input$Orderingchoosestep=='start' && !is.null(input$Orderingstartslider) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
      #                   pdata <- Maindata$pdata
      #                   pdata <- pdata[order(pdata$Pseudotime),]
      #                   start <- as.numeric(input$Orderingstartslider)
      #                   if (start == 1) {
      #                         order <- 1:nrow(pdata)
      #                   } else {
      #                         order <- c(start:(nrow(pdata)),1:(start-1))
      #                   }
      #                   if (input$Orderingstartflip)
      #                         order <- rev(order)
      #                   pdata <- pdata[order,]
      #                   datadist <- dist(t(Maindata$reduceres))
      #                   distmat <- as.matrix(datadist)
      #                   alldist <- sapply(1:(nrow(pdata)-1), function(x) {
      #                         distmat[pdata[x,1],pdata[x+1,1]]
      #                   })
      #                   ptime <- c(0,cumsum(alldist))
      #                   ptime <- ptime/max(ptime) * as.numeric(input$Orderingptimescale)
      #                   pdata[,3] <- ptime
      #                   Maindata$scapdata <- pdata
      #             }
      #       })
      #       
      #       observe({
      #             if (input$Orderingchoosestep=='start' && !is.null(Maindata$pdata) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
      #                   pdata <- Maindata$pdata[,-2]
      #                   pdata <- pdata[order(pdata$Pseudotime),]
      #                   if (is.null(Maindata$decreaseexpr)) {
      #                         allexpr <- Maindata$increaseexpr
      #                   } else {
      #                         allexpr <- rbind(Maindata$increaseexpr,-Maindata$decreaseexpr)     
      #                   }
      #                   allres <- NULL
      #                   for (flip in c(TRUE,FALSE)) {
      #                         for (start in 1:nrow(pdata)) {
      #                               if (start == 1) {
      #                                     order <- 1:nrow(pdata)
      #                               } else {
      #                                     order <- c(start:(nrow(pdata)),1:(start-1))
      #                               }
      #                               if (flip)
      #                                     order <- rev(order)            
      #                               cellorder <- pdata[order,1]                              
      #                               pcos <- sum(apply(allexpr[,cellorder,drop=F],1,function(expr) {
      #                                     sum(sapply(1:(length(expr)-1),function(x) {
      #                                           sum(expr[(x+1):length(expr)] - expr[x])
      #                                     }))                                    
      #                               }))
      #                               allres <- rbind(allres,c(start,flip,pcos))
      #                         }
      #                   }
      #                   Maindata$Orderingsuggestres <- allres
      #             }
      #       })
      #       
      #       output$Orderingstartshowheatmap <- renderPlot({
      #             if (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr)) {
      #                   allexpr <- rbind(Maindata$increaseexpr,Maindata$decreaseexpr)
      #                   par(mar=c(0,0,0,0))
      #                   if (nrow(allexpr) == 1) {
      #                         plot(allexpr[,Maindata$scapdata[,1]])
      #                   } else {
      #                         image(t(allexpr[,Maindata$scapdata[,1]]),axes=F,col=bluered(100))      
      #                         names <- rev(row.names(allexpr))
      #                         dim <- length(names)
      #                         pos <- seq(0,1,length.out=dim)
      #                         for (i in 1:dim) {
      #                               text(0.8,pos[i],names[i],cex=1.3)
      #                         }
      #                   }
      #             }            
      #       })
      #       
      #       output$Orderingstartshowheatui <- renderUI({
      #             if (!is.null(Maindata$procdata) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
      #                   tagList(                  
      #                         p(actionButton("Orderingstartsetoptimalbutton","Set optimal value")),
      #                         checkboxInput("Orderingstartflip","Flip the ordering"),
      #                         sliderInput("Orderingstartslider","Slide to select the starting point",min=1,max=ncol(Maindata$procdata),step=1,value=1,width='800px'),
      #                         plotOutput("Orderingstartshowheatmap",width='800px'),
      #                         p(textOutput("Orderingstartsuggestshowres"))
      #                   )      
      #             }
      #       })
      #       
      #       observe({
      #             if (!is.null(input$Orderingstartsetoptimalbutton) && input$Orderingstartsetoptimalbutton) {
      #                   isolate({
      #                         maxpos <- Maindata$Orderingsuggestres[which.max(Maindata$Orderingsuggestres[,3]),]
      #                         updateCheckboxInput(session = session,"Orderingstartflip",value=as.logical(maxpos[2]))
      #                         updateSliderInput(session = session, "Orderingstartslider",value=as.integer(maxpos[1]))
      #                   })
      #             }
      #       })
      #       
      #       output$Orderingstartsuggestshowres <- renderText({
      #             if (!is.null(input$Orderingstartsetoptimalbutton)) {
      #                   pcos <- Maindata$Orderingsuggestres[Maindata$Orderingsuggestres[,1]== as.numeric(input$Orderingstartslider) & Maindata$Orderingsuggestres[,2] == input$Orderingstartflip,3]
      #                   paste("Pseudotemporal cell ordering score:",pcos)      
      #             }
      #       })
      #       
      #       output$Orderingstartsuggestshowgeneset <- renderDataTable(Maindata$fullstartgeneset)
      #       
      #       output$Orderingstartsuggestshowpseudotime <- renderDataTable(Maindata$scapdata)
      
      #save results
      
      output$Orderingsaveplotparaui <- renderUI({
            if (input$Orderingsaveplotparatf) {
                  tagList(
                        textInput("Orderingsaveplotchangetitle","Change title",value="Pseudotime ordering plot"),
                        textInput("Orderingsaveplotchangexlab","Change x axis title",value="Component 2"),
                        textInput("Orderingsaveplotchangeylab","Change y axis title",value="Component 1")
                  )                  
            }
      })
      
      output$Orderingsavepdata <- downloadHandler(
            filename = function() { paste0("Pseudotime_time.",ifelse(input$Orderingsavepdatatype == 'txt','txt','csv')) },
            content = function(file) {
                  if (input$Orderingsavepdatatype == 'txt') {
                        write.table(Maindata$scapdata,file,row.names=F,quote=F)     
                  } else {
                        write.csv(Maindata$scapdata,file,row.names=F,quote=F)     
                  }
                  
            }
      )
      
      output$Orderingsaveshowptime <- renderDataTable(Maindata$scapdata)
      
      output$Orderingsaveshowplot <- renderPlot({
            if (input$Orderingptimechoosemethod=="Monocle" && !is.null(Maindata$dp_Monocle)) {
                  Monocledrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle)
            } else if (input$Orderingptimechoosemethod=="TSCAN" && !is.null(input$OrderingTSCANshow_cell_names)) {
                  TSCANdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle)
            }
      })
      
      output$Orderingsaveplot <- downloadHandler(
            filename = function() { paste0('Cell_ordering.',input$Orderingsaveplottype) },
            content = function(file) {
                  if (input$Orderingsaveplottype == 'pdf') {
                        pdf(file,width=as.numeric(input$Orderingsaveplotfilewidth),height=as.numeric(input$Orderingsaveplotfileheight))
                  } else if (input$Orderingsaveplottype == 'ps') {
                        postscript(file,width=as.numeric(input$Orderingsaveplotfilewidth),height=as.numeric(input$Orderingsaveplotfileheight),paper="special")
                  }
                  if (input$Orderingptimechoosemethod=="Monocle") {
                        print(Monocledrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle))
                  } else if (input$Orderingptimechoosemethod=="TSCAN") {
                        print(TSCANdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle))
                  }
                  dev.off()
            }
      )
      
      ### Difftest ###
      
      #choose uploaded pdata or sca pdata
      
      output$Difftestpdataselectui <- renderUI({
            if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata)) {
                  radioButtons("Difftestpdataselect","Select the cell ordering you want to use",choices = list("Uploaded pseudotime"="Uploaded","SCA generated pseudotime"="SCA"))     
            }            
      })
      
      observe({
            if (input$MainMenu=="Difftest") {
                  if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata) && !is.null(input$Difftestpdataselect)) {
                        if (input$Difftestpdataselect == 'Uploaded') {
                              Maindata$finalpdata <- Maindata$uploadpdata
                        } else {
                              Maindata$finalpdata <- Maindata$scapdata
                        }
                  } else if (!is.null(Maindata$uploadpdata) && is.null(Maindata$scapdata)) {
                        Maindata$finalpdata <- Maindata$uploadpdata
                  } else if (is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata)) {
                        Maindata$finalpdata <- Maindata$scapdata
                  }            
            }
      })
      
      observe({
            if (!is.null(input$Difftestbutton) && input$Difftestbutton > 0) {
                  isolate({
                        progress <- shiny::Progress$new(session, min=1, max=1)
                        on.exit(progress$close())
                        
                        progress$set(message = 'Calculation in progress')
                        
                        for (i in 1:nrow(Maindata$rawlogdata)) {
                              progress$set(value = i)
                              Sys.sleep(0.2)
                        }
                        tmpdata <- Maindata$rawlogdata[,Maindata$finalpdata[,1]]
                        ptime <- 1:ncol(tmpdata)
                        
                        Difftestpval <- apply(tmpdata,1,function(x) {                                                                                          
                              if (sum(x) == 0) {
                                    1
                              } else {
                                    model <- mgcv::gam(x~s(ptime,k=3))
                                    pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F)                  
                              }   
                        })
                        Maindata$Difftestadjpval <- p.adjust(Difftestpval,method="fdr")   
                        Maindata$Difftestcalculatestatus <- "End"
                  })
            }
      })
      
      observe({
            if (!is.null(Maindata$Difftestadjpval)) {
                  Maindata$Difftestdata <- Maindata$rawlogdata[Maindata$Difftestadjpval < as.numeric(input$Difftestfdrval),Maindata$finalpdata[,1]]
                  Maindata$Difftestresalltable <- data.frame(GENE=row.names(Maindata$rawlogdata),adjusted.p.value=Maindata$Difftestadjpval)      
                  Maindata$Difftestrestable <- data.frame(GENE=row.names(Maindata$Difftestdata),adjusted.p.value=Maindata$Difftestadjpval[Maindata$Difftestadjpval < as.numeric(input$Difftestfdrval)])      
            }
      })
      
      output$Difftestsummaryui <- renderUI({
            if (!is.null(Maindata$Difftestdata) && nrow(Maindata$Difftestdata)!=0)
                  p(paste(nrow(Maindata$Difftestdata),"genes out of",nrow(Maindata$rawlogdata),"genes are differentially expressed, which is",round(nrow(Maindata$Difftestdata)/nrow(Maindata$procdata)*100,3),"percent."))
      })      
      
      output$Difftestcalculatecomplete <- renderText({
            if (!is.null(Maindata$Difftestcalculatestatus) && Maindata$Difftestcalculatestatus == "End") {
                  "Calculation completed!"
            }
      })
      
      output$Difftestshowresult <- renderDataTable({
            if (!is.null(Maindata$Difftestdata) && nrow(Maindata$Difftestdata)!=0)
                  if (input$Difftestshowresultopt == "all") {
                        Maindata$Difftestresalltable   
                  } else {
                        Maindata$Difftestrestable
                  }                  
      })
      
      output$Difftestsavepvaltable <- downloadHandler(
            filename = function() { paste0("Test_for_differentially_expressed.",ifelse(input$Difftestsavepvaltabletype == 'txt','txt','csv')) },
            content = function(file) {
                  if (input$Difftestsavepvaltabletype == 'txt') {
                        if (input$Difftestshowresultopt == "all") {
                              write.table(Maindata$Difftestresalltable,file,row.names=F,quote=F)                              
                        } else {
                              write.table(Maindata$Difftestrestable,file,row.names=F,quote=F)     
                        }                          
                  } else {
                        if (input$Difftestshowresultopt == "all") {
                              write.csv(Maindata$Difftestresalltable,file,row.names=F,quote=F)                              
                        } else {
                              write.csv(Maindata$Difftestrestable,file,row.names=F,quote=F)     
                        }
                  }                  
            }
      )
      
      output$Difftestheatmapui <- renderUI({
            tagList({
                  helpText("All differentially expressed genes:")
                  plotOutput("Difftestheatmap")                  
            })
            
      })
      
      output$Difftestheatmap <- renderPlot({
            heatmap.2(Maindata$Difftestdata,trace="none")
      })
      
      output$Difftestsaveplot <- downloadHandler(
            filename = function() { paste0('Gene_expression.',input$Orderingsaveplottype) },
            content = function(file) {
                  if (input$Difftestplottype == 'pdf') {
                        pdf(file,width=as.numeric(input$Difftestfilewidth),height=as.numeric(input$Difftestfileheight))
                  } else if (input$Difftestplottype == 'ps') {
                        postscript(file,width=as.numeric(input$Difftestfilewidth),height=as.numeric(input$Difftestfileheight),paper="special")
                  }                  
                  heatmap.2(Maindata$Difftestdata,trace="none")
                  dev.off()
            }
      )
      
      ### Visualization ###
      
      output$Visualizationpdataselectui <- renderUI({
            if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata)) {
                  radioButtons("Visualizationpdataselect","Select the cell ordering you want to use",choices = list("Uploaded pseudotime"="Uploaded","SCA generated pseudotime"="SCA"))     
            }            
      })
      
      observe({
            if (input$MainMenu=="Visualization") {
                  if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata) && !is.null(input$Visualizationpdataselect)) {
                        if (input$Visualizationpdataselect == 'Uploaded') {
                              Maindata$finalpdata <- Maindata$uploadpdata
                        } else {
                              Maindata$finalpdata <- Maindata$scapdata
                        }
                  } else if (!is.null(Maindata$uploadpdata) && is.null(Maindata$scapdata)) {
                        Maindata$finalpdata <- Maindata$uploadpdata
                  } else if (is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata)) {
                        Maindata$finalpdata <- Maindata$scapdata
                  }            
            }
      })
      
      output$Visualizationmainui <- renderUI({
            plotOutput("Visualizationshowplot",width = "100%",height=ifelse(length(input$Visualizationgeneselect)==1,800,ifelse(!is.null(input$Visualizationmethod) && input$Visualizationmethod, 400,800*length(input$Visualizationgeneselect))))            
      })
      
      
      output$Visualizationgeneselectui <- renderUI({            
            selectInput("Visualizationgeneselect","Select gene",choices = row.names(Maindata$rawlogdata),multiple = T)
      })
      
      output$Visualizationmethodui <- renderUI({
            if (!is.null(input$Visualizationgeneselect) && length(input$Visualizationgeneselect) > 1) {
                  checkboxInput("Visualizationmethod","Display heatmap",value = F)
            }
      })
      
      Visualizationplotfunc <- function() {
            tmpdata <- Maindata$rawlogdata[input$Visualizationgeneselect,Maindata$finalpdata[,1]]
            if (is.vector(tmpdata)) {
                  ptime <- 1:length(tmpdata)
                  plot(ptime,tmpdata,xlab="Cell ordering index",ylab="Expression value",pch=19,col=Maindata$finalpdata[,2],main=input$Visualizationgeneselect)      
                  model <- mgcv::gam(tmpdata~s(ptime,k=3))
                  lines(ptime,fitted(model))                  
            } else {
                  if (!is.null(input$Visualizationmethod) && input$Visualizationmethod) {
                        heatmap.2(tmpdata,Rowv = F,Colv = F,col = bluered,symbreaks=F,dendrogram="none",trace="none",cexRow=1,srtRow=-45,lwid=c(0.2,1))          
                  } else {
                        par(mfrow=c(nrow(tmpdata),1))
                        for (i in 1:nrow(tmpdata)) {                                    
                              x <- tmpdata[i,]
                              ptime <- 1:length(x)
                              plot(1:length(x),x,xlab="Cell ordering index",ylab="Expression value",pch=19,col=Maindata$finalpdata[,2],main=row.names(tmpdata)[i])      
                              model <- mgcv::gam(x~s(ptime,k=3))
                              lines(ptime,fitted(model))                  
                        }
                  }                              
            }
      }
      
      output$Visualizationshowplot <- renderPlot({
            if (!is.null(input$Visualizationgeneselect) && input$Visualizationgeneselect[1]!="") {
                  Visualizationplotfunc()                        
            }            
      })
      
      output$Visualizationfileheightui <- renderUI({
            textInput("Visualizationfileheight","Enter plot height (inches)",value = ifelse(length(input$Visualizationgeneselect)==1,12,ifelse(!is.null(input$Visualizationmethod) && input$Visualizationmethod, 6,12*length(input$Visualizationgeneselect))))
      })
      
      output$Visualizationsaveplot <- downloadHandler(
            filename = function() { paste0('Gene_expression.',input$Visualizationplottype) },
            content = function(file) {
                  if (input$Visualizationplottype == 'pdf') {
                        pdf(file,width=as.numeric(input$Visualizationfilewidth),height=as.numeric(input$Visualizationfileheight))
                  } else if (input$Visualizationplottype == 'ps') {
                        postscript(file,width=as.numeric(input$Visualizationfilewidth),height=as.numeric(input$Visualizationfileheight),paper="special")
                  }                  
                  Visualizationplotfunc() 
                  dev.off()
            }
      )
      
      
      ###  Miscellaneous ###
      
      observe({
            if (input$MainMenu != "Ordering")
                  updateRadioButtons(session,"Orderingchoosestep","",list("Step 1: Reduce dimension"="reduction","Step 2: Pseudo time reconstruction"="ptime","Save results (optional)"="save"),selected = "reduction")
      })
      
      Miscdata <- reactiveValues()
      
      output$Miscshowresultsui <- renderUI({
            if (input$Compareinputopt=='sub') {
                  tagList(
                        checkboxInput("Miscsubshowinstructiontf","Show instructions",value=T),
                        uiOutput("Miscsubshowinstruction"),
                        dataTableOutput("Miscsubshowtable")
                  )
            } else {
                  tabsetPanel(
                        tabPanel("Current order",
                                 checkboxInput("Miscordershowinstructiontf","Show instructions",value=T),
                                 uiOutput("Miscordershowinstruction"),
                                 tableOutput("Miscordershowtable")
                        ),
                        tabPanel("Ordering scores",
                                 tableOutput("Miscordershowscores")
                        ),
                        tabPanel("All order",
                                 tableOutput("Miscordershowalltable")                                 
                        )                        
                  )
            }            
      })
      
      output$Miscsubshowtable <- renderDataTable({
            Miscdata$sub
      })
      
      output$Miscsubshowinstruction <- renderUI({
            if (input$Miscsubshowinstructiontf) {
                  tagList(
                        h5("Instructions:"),
                        p("Sub-population information should be prepared in a data.frame. The first column is the cell name and the second column is the subpopulation id."),
                        p("Sub-population id should be integers starting from 0. Larger id stands for subpopulation collected at latter time point"),
                        p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to make changes."),
                        h5("A typical example of tab-delimited file format is as follows:"),
                        p("Cell  ID"),
                        p("Cell1  0"),
                        p("Cell2  1"),
                        p("Cell3  1")                      
                  )
            }
      })
      
      observe({
            if (input$Comparesubreadin > 0)
                  isolate({
                        FileHandle <- input$ComparesubFile
                        if (!is.null(FileHandle)) {
                              Miscdata$sub <- read.table(FileHandle$datapath,header=input$Comparesubheader,sep=input$Comparesubsep,quote=input$Comparesubquote,stringsAsFactors=F,blank.lines.skip=TRUE)
                        }
                  })
      })
      
      output$Miscordershowtable <- renderTable({
            Miscdata$order
      })
      
      output$Miscordershowalltable <- renderTable({
            tmp <- sapply(Miscdata$allorder,cbind)
            colnames(tmp) <- Miscdata$allordername
            tmp
      })
      
      output$Miscordershowscores <- renderTable({
            subinfo <- Miscdata$sub[,2]
            names(subinfo) <- Miscdata$sub[,1]
            scorefunc <- function(order) {
                  scoreorder <- subinfo[order]
                  optscoreorder <- sort(scoreorder)
                  optscore <- sum(sapply(1:(length(optscoreorder)-1),function(x) {
                        sum(optscoreorder[(x+1):length(optscoreorder)] - optscoreorder[x])
                  })) 
                  sum(sapply(1:(length(scoreorder)-1),function(x) {
                        sum(scoreorder[(x+1):length(scoreorder)] - scoreorder[x])
                  })) / optscore
            }      
            data.frame(order = Miscdata$allordername, score = sapply(Miscdata$allorder,scorefunc))
            
      })
      
      output$Miscordershowinstruction <- renderUI({
            if (input$Miscordershowinstructiontf) {
                  tagList(
                        h5("Instructions:"),
                        p("order information should be prepared in a data.frame with one column. The first column is the ordered cell name. If the data has multiple columns, other columns will be omitted."),
                        p("Please make sure that the cell names agree exactly with the cell names given in the subpopulation information, otherwise unpredictable errors may occur."),
                        p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to make changes."),
                        h5("A typical example of tab-delimited file format is as follows:"),
                        p("Cell"),
                        p("Cell3"),
                        p("Cell1"),
                        p("Cell2")                      
                  )
            }
      })
      
      observe({
            if (input$Compareorderreadin > 0)
                  isolate({
                        FileHandle <- input$CompareorderFile
                        if (!is.null(FileHandle)) {
                              Miscdata$order <- read.table(FileHandle$datapath,header=input$Compareorderheader,sep=input$Compareordersep,quote=input$Compareorderquote,stringsAsFactors=F,blank.lines.skip=TRUE)
                        }
                  })
      })
      
      observe({
            if (input$Compareorderadddata && input$Compareorderadddata > 0) {
                  isolate({
                        tmporder <- as.vector(Miscdata$order[,1])
                        if (is.null(Miscdata$allorder)) {
                              Miscdata$allorder <- list(tmporder)
                              Miscdata$allordername <- input$Compareordername
                        } else {
                              Miscdata$allorder <- c(Miscdata$allorder, list(tmporder))
                              Miscdata$allordername <- c(Miscdata$allordername,input$Compareordername)
                        }                        
                  })                  
            }                        
      })
      
      output$compareordernameui <- renderUI({
            if (is.null(Miscdata$allorder)) {
                  textInput("Compareordername","Input name for order","Order_1")            
            } else {
                  textInput("Compareordername","Input name for order",paste0("Order_",length(Miscdata$allorder)+1))
            }
      })
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      ###Monocle internal functions
{
            get_next_node_id <- function () {
                  next_node <<- next_node + 1
                  return(next_node)
            }
            
            pq_helper <- function (mst, use_weights = TRUE, root_node = NULL) {
                  new_subtree <- graph.empty()
                  root_node_id <- paste("Q_", get_next_node_id(), sep = "")
                  new_subtree <- new_subtree + vertex(root_node_id, type = "Q", 
                                                      color = "black")
                  if (is.null(root_node) == FALSE) {
                        sp <- get.all.shortest.paths(mst, from = V(mst)[root_node])
                        sp_lengths <- sapply(sp$res, length)
                        target_node_idx <- which(sp_lengths == max(sp_lengths))[1]
                        diam <- V(mst)[unlist(sp$res[target_node_idx])]
                  }
                  else {
                        if (use_weights) {
                              diam <- V(mst)[get.diameter(mst)]
                        }
                        else {
                              diam <- V(mst)[get.diameter(mst, weights = NA)]
                        }
                  }
                  V(new_subtree)[root_node_id]$diam_path_len = length(diam)
                  diam_decisiveness <- igraph::degree(mst, v = diam) > 2
                  ind_nodes <- diam_decisiveness[diam_decisiveness == TRUE]
                  first_diam_path_node_idx <- head(as.vector(diam), n = 1)
                  last_diam_path_node_idx <- tail(as.vector(diam), n = 1)
                  if (sum(ind_nodes) == 0 || (igraph::degree(mst, first_diam_path_node_idx) == 
                                                    1 && igraph::degree(mst, last_diam_path_node_idx) == 
                                                    1)) {
                        ind_backbone <- diam
                  }
                  else {
                        last_bb_point <- names(tail(ind_nodes, n = 1))[[1]]
                        first_bb_point <- names(head(ind_nodes, n = 1))[[1]]
                        diam_path_names <- V(mst)[as.vector(diam)]$name
                        last_bb_point_idx <- which(diam_path_names == last_bb_point)[1]
                        first_bb_point_idx <- which(diam_path_names == first_bb_point)[1]
                        ind_backbone_idxs <- as.vector(diam)[first_bb_point_idx:last_bb_point_idx]
                        ind_backbone <- V(mst)[ind_backbone_idxs]
                  }
                  mst_no_backbone <- mst - ind_backbone
                  for (backbone_n in ind_backbone) {
                        if (igraph::degree(mst, v = backbone_n) > 2) {
                              new_p_id <- paste("P_", get_next_node_id(), sep = "")
                              new_subtree <- new_subtree + vertex(new_p_id, type = "P", 
                                                                  color = "grey")
                              new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, 
                                                                  type = "leaf", color = "white")
                              new_subtree <- new_subtree + edge(new_p_id, V(mst)[backbone_n]$name)
                              new_subtree <- new_subtree + edge(root_node_id, new_p_id)
                              nb <- graph.neighborhood(mst, 1, nodes = backbone_n)[[1]]
                              for (n_i in V(nb)) {
                                    n <- V(nb)[n_i]$name
                                    if (n %in% V(mst_no_backbone)$name) {
                                          sc <- subcomponent(mst_no_backbone, n)
                                          sg <- induced.subgraph(mst_no_backbone, sc, 
                                                                 impl = "copy_and_delete")
                                          if (ecount(sg) > 0) {
                                                sub_pq <- pq_helper(sg, use_weights)
                                                for (v in V(sub_pq$subtree)) {
                                                      new_subtree <- new_subtree + vertex(V(sub_pq$subtree)[v]$name, 
                                                                                          type = V(sub_pq$subtree)[v]$type, color = V(sub_pq$subtree)[v]$color, 
                                                                                          diam_path_len = V(sub_pq$subtree)[v]$diam_path_len)
                                                }
                                                edge_list <- get.edgelist(sub_pq$subtree)
                                                for (i in 1:nrow(edge_list)) {
                                                      new_subtree <- new_subtree + edge(V(sub_pq$subtree)[edge_list[i, 
                                                                                                                    1]]$name, V(sub_pq$subtree)[edge_list[i, 
                                                                                                                                                          2]]$name)
                                                }
                                                new_subtree <- new_subtree + edge(new_p_id, 
                                                                                  V(sub_pq$subtree)[sub_pq$root]$name)
                                          }
                                          else {
                                                new_subtree <- new_subtree + vertex(n, type = "leaf", 
                                                                                    color = "white")
                                                new_subtree <- new_subtree + edge(new_p_id, 
                                                                                  n)
                                          }
                                    }
                              }
                        }
                        else {
                              new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, 
                                                                  type = "leaf", color = "white")
                              new_subtree <- new_subtree + edge(root_node_id, V(mst)[backbone_n]$name)
                        }
                  }
                  return(list(root = root_node_id, subtree = new_subtree))
            }
            
            order_p_node <- function (q_level_list, dist_matrix) {
                  q_order_res <- permn(q_level_list, fun = order_q_node, dist_matrix)
                  all_perms <- lapply(q_order_res, function(x) {
                        x$ql
                  })
                  all_perms_weights <- unlist(lapply(q_order_res, function(x) {
                        x$wt
                  }))
                  opt_perm_idx <- head((which(all_perms_weights == min(all_perms_weights))), 
                                       1)
                  opt_perm <- all_perms[[opt_perm_idx]]
                  stopifnot(length(opt_perm) == length(q_level_list))
                  return(opt_perm)
            }
            
            order_q_node <- function (q_level_list, dist_matrix) {
                  new_subtree <- graph.empty()
                  if (length(q_level_list) == 1) {
                        return(list(ql = q_level_list, wt = 0))
                  }
                  for (i in 1:length(q_level_list)) {
                        new_subtree <- new_subtree + vertex(paste(i, "F"), type = "forward")
                        new_subtree <- new_subtree + vertex(paste(i, "R"), type = "reverse")
                  }
                  for (i in (1:(length(q_level_list) - 1))) {
                        cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], 
                                            q_level_list[[i + 1]][1]]
                        new_subtree <- new_subtree + edge(paste(i, "F"), paste(i + 
                                                                                     1, "F"), weight = cost)
                        cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], 
                                            q_level_list[[i + 1]][length(q_level_list[[i + 1]])]]
                        new_subtree <- new_subtree + edge(paste(i, "F"), paste(i + 
                                                                                     1, "R"), weight = cost)
                        cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i + 
                                                                                      1]][1]]
                        new_subtree <- new_subtree + edge(paste(i, "R"), paste(i + 
                                                                                     1, "F"), weight = cost)
                        cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i + 
                                                                                      1]][length(q_level_list[[i + 1]])]]
                        new_subtree <- new_subtree + edge(paste(i, "R"), paste(i + 
                                                                                     1, "R"), weight = cost)
                  }
                  first_fwd = V(new_subtree)[paste(1, "F")]
                  first_rev = V(new_subtree)[paste(1, "R")]
                  last_fwd = V(new_subtree)[paste(length(q_level_list), "F")]
                  last_rev = V(new_subtree)[paste(length(q_level_list), "R")]
                  FF_path <- unlist(get.shortest.paths(new_subtree, from = as.vector(first_fwd), 
                                                       to = as.vector(last_fwd), mode = "out", output = "vpath")$vpath)
                  FR_path <- unlist(get.shortest.paths(new_subtree, from = as.vector(first_fwd), 
                                                       to = as.vector(last_rev), mode = "out", output = "vpath")$vpath)
                  RF_path <- unlist(get.shortest.paths(new_subtree, from = as.vector(first_rev), 
                                                       to = as.vector(last_fwd), mode = "out", output = "vpath")$vpath)
                  RR_path <- unlist(get.shortest.paths(new_subtree, from = as.vector(first_rev), 
                                                       to = as.vector(last_rev), mode = "out", output = "vpath")$vpath)
                  FF_weight <- sum(E(new_subtree, path = FF_path)$weight)
                  FR_weight <- sum(E(new_subtree, path = FR_path)$weight)
                  RF_weight <- sum(E(new_subtree, path = RF_path)$weight)
                  RR_weight <- sum(E(new_subtree, path = RR_path)$weight)
                  paths <- list(FF_path, FR_path, RF_path, RR_path)
                  path_weights <- c(FF_weight, FR_weight, RF_weight, RR_weight)
                  opt_path_idx <- head((which(path_weights == min(path_weights))), 
                                       1)
                  opt_path <- paths[[opt_path_idx]]
                  stopifnot(length(opt_path) == length(q_level_list))
                  directions <- V(new_subtree)[opt_path]$type
                  q_levels <- list()
                  for (i in 1:length(directions)) {
                        if (directions[[i]] == "forward") {
                              q_levels[[length(q_levels) + 1]] <- q_level_list[[i]]
                        }
                        else {
                              q_levels[[length(q_levels) + 1]] <- rev(q_level_list[[i]])
                        }
                  }
                  return(list(ql = q_levels, wt = min(path_weights)))
            }
            
            assign_cell_lineage <- function (pq_tree, curr_node, assigned_state, node_states) {
                  if (V(pq_tree)[curr_node]$type == "leaf") {
                        node_states[V(pq_tree)[curr_node]$name] = assigned_state
                        return(node_states)
                  }
                  else {
                        for (child in V(pq_tree)[nei(curr_node, mode = "out")]) {
                              node_states <- assign_cell_lineage(pq_tree, child, 
                                                                 assigned_state, node_states)
                        }
                        return(node_states)
                  }
            }
            
            extract_good_ordering <- function (pq_tree, curr_node, dist_matrix) {
                  if (V(pq_tree)[curr_node]$type == "leaf") {
                        return(V(pq_tree)[curr_node]$name)
                  }
                  else if (V(pq_tree)[curr_node]$type == "P") {
                        p_level <- list()
                        for (child in V(pq_tree)[nei(curr_node, mode = "out")]) {
                              p_level[[length(p_level) + 1]] <- extract_good_ordering(pq_tree, 
                                                                                      child, dist_matrix)
                        }
                        p_level <- order_p_node(p_level, dist_matrix)
                        p_level <- unlist(p_level)
                        return(p_level)
                  }
                  else if (V(pq_tree)[curr_node]$type == "Q") {
                        q_level <- list()
                        for (child in V(pq_tree)[nei(curr_node, mode = "out")]) {
                              q_level[[length(q_level) + 1]] <- extract_good_ordering(pq_tree, 
                                                                                      child, dist_matrix)
                        }
                        q_level <- order_q_node(q_level, dist_matrix)
                        q_level <- q_level$ql
                        q_level <- unlist(q_level)
                        return(q_level)
                  }
            }
            
            weight_of_ordering <- function (ordering, dist_matrix) {
                  time_delta <- c(0)
                  curr_weight <- 0
                  ep <- 0.01
                  for (i in 2:length(ordering)) {
                        d <- dist_matrix[ordering[[i]], ordering[[i - 1]]]
                        curr_weight <- curr_weight + d + ep
                        time_delta <- c(time_delta, curr_weight)
                  }
                  return(time_delta)
            }
            
            extract_good_branched_ordering <- function (orig_pq_tree, curr_node, dist_matrix, num_branches, 
                                                        reverse_main_path = FALSE) {
                  pq_tree <- orig_pq_tree
                  branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
                  names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
                  branch_node_counts <- sort(branch_node_counts, decreasing = TRUE)
                  cell_states <- rep(NA, length(as.vector(V(pq_tree)[type == 
                                                                           "leaf"])))
                  names(cell_states) <- V(pq_tree)[type == "leaf"]$name
                  cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, 
                                                     cell_states)
                  branch_point_roots <- list()
                  branch_tree <- graph.empty()
                  for (i in 1:num_branches) {
                        branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
                        branch_id <- names(branch_node_counts)[i]
                        branch_tree <- branch_tree + vertex(branch_id)
                        parents <- V(pq_tree)[nei(names(branch_node_counts)[i], 
                                                  mode = "in")]
                        if (length(parents) > 0 && parents$type == "P") {
                              p_node_parent <- V(pq_tree)[nei(names(branch_node_counts)[i], 
                                                              mode = "in")]
                              parent_branch_id <- V(pq_tree)[nei(p_node_parent, 
                                                                 mode = "in")]$name
                              branch_tree <- branch_tree + edge(parent_branch_id, 
                                                                branch_id)
                        }
                        pq_tree[V(pq_tree)[nei(names(branch_node_counts)[i], 
                                               mode = "in")], names(branch_node_counts)[i]] <- FALSE
                  }
                  branch_pseudotimes <- list()
                  for (i in 1:length(branch_point_roots)) {
                        branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], 
                                                                 dist_matrix)
                        branch_ordering_time <- weight_of_ordering(branch_ordering, 
                                                                   dist_matrix)
                        names(branch_ordering_time) <- branch_ordering
                        branch_pseudotimes[[length(branch_pseudotimes) + 1]] = branch_ordering_time
                        names(branch_pseudotimes)[length(branch_pseudotimes)] = branch_point_roots[[i]]
                  }
                  cell_ordering_tree <- graph.empty()
                  curr_branch <- "Q_1"
                  extract_branched_ordering_helper <- function(branch_tree, 
                                                               curr_branch, cell_ordering_tree, branch_pseudotimes, 
                                                               dist_matrix, reverse_ordering = FALSE) {
                        curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
                        curr_branch_root_cell <- NA
                        for (i in 1:length(curr_branch_pseudotimes)) {
                              cell_ordering_tree <- cell_ordering_tree + vertex(names(curr_branch_pseudotimes)[i])
                              if (i > 1) {
                                    if (reverse_ordering == FALSE) {
                                          cell_ordering_tree <- cell_ordering_tree + 
                                                edge(names(curr_branch_pseudotimes)[i - 1], 
                                                     names(curr_branch_pseudotimes)[i])
                                    }
                                    else {
                                          cell_ordering_tree <- cell_ordering_tree + 
                                                edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i - 
                                                                                                                             1])
                                    }
                              }
                        }
                        if (reverse_ordering == FALSE) {
                              curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
                        }
                        else {
                              curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
                        }
                        for (child in V(branch_tree)[nei(curr_branch, mode = "out")]) {
                              child_cell_ordering_subtree <- graph.empty()
                              child_head <- names(branch_pseudotimes[[child]])[1]
                              child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]
                              curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
                              head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
                              closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == 
                                                                                      min(head_dist_to_curr))]
                              head_dist_to_anchored_branch = NA
                              branch_index_for_head <- NA
                              head_dist_to_anchored_branch <- dist_matrix[closest_to_head, 
                                                                          child_head]
                              tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
                              closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == 
                                                                                      min(tail_dist_to_curr))]
                              tail_dist_to_anchored_branch = NA
                              branch_index_for_tail <- NA
                              tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, 
                                                                          child_tail]
                              if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
                                    reverse_child <- TRUE
                              }
                              else {
                                    reverse_child <- FALSE
                              }
                              res <- extract_branched_ordering_helper(branch_tree, 
                                                                      child, child_cell_ordering_subtree, branch_pseudotimes, 
                                                                      dist_matrix, reverse_child)
                              child_cell_ordering_subtree <- res$subtree
                              child_subtree_root <- res$root
                              for (v in V(child_cell_ordering_subtree)) {
                                    cell_ordering_tree <- cell_ordering_tree + vertex(V(child_cell_ordering_subtree)[v]$name)
                              }
                              edge_list <- get.edgelist(child_cell_ordering_subtree)
                              for (i in 1:nrow(edge_list)) {
                                    cell_ordering_tree <- cell_ordering_tree + edge(V(cell_ordering_tree)[edge_list[i, 
                                                                                                                    1]]$name, V(cell_ordering_tree)[edge_list[i, 
                                                                                                                                                              2]]$name)
                              }
                              if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
                                    cell_ordering_tree <- cell_ordering_tree + edge(closest_to_tail, 
                                                                                    child_subtree_root)
                              }
                              else {
                                    cell_ordering_tree <- cell_ordering_tree + edge(closest_to_head, 
                                                                                    child_subtree_root)
                              }
                        }
                        return(list(subtree = cell_ordering_tree, root = curr_branch_root_cell, 
                                    last_cell_state = 1, last_cell_pseudotime = 0))
                  }
                  res <- extract_branched_ordering_helper(branch_tree, curr_branch, 
                                                          cell_ordering_tree, branch_pseudotimes, dist_matrix, 
                                                          reverse_main_path)
                  cell_ordering_tree <- res$subtree
                  curr_state <- 1
                  assign_cell_state_helper <- function(ordering_tree_res, curr_cell) {
                        cell_tree <- ordering_tree_res$subtree
                        V(cell_tree)[curr_cell]$cell_state = curr_state
                        children <- V(cell_tree)[nei(curr_cell, mode = "out")]
                        ordering_tree_res$subtree <- cell_tree
                        if (length(children) == 1) {
                              ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, 
                                                                            V(cell_tree)[children]$name)
                        }
                        else {
                              for (child in children) {
                                    curr_state <<- curr_state + 1
                                    ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, 
                                                                                  V(cell_tree)[child]$name)
                              }
                        }
                        return(ordering_tree_res)
                  }
                  res <- assign_cell_state_helper(res, res$root)
                  assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, 
                                                       last_pseudotime, curr_cell) {
                        cell_tree <- ordering_tree_res$subtree
                        curr_cell_pseudotime <- last_pseudotime
                        V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
                        ordering_tree_res$subtree <- cell_tree
                        children <- V(cell_tree)[nei(curr_cell, mode = "out")]
                        for (child in children) {
                              next_node <- V(cell_tree)[child]$name
                              delta_pseudotime <- dist_matrix[curr_cell, next_node]
                              ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, 
                                                                            dist_matrix, last_pseudotime + delta_pseudotime, 
                                                                            next_node)
                        }
                        return(ordering_tree_res)
                  }
                  res <- assign_pseudotime_helper(res, dist_matrix, 0, res$root)
                  cell_names <- V(res$subtree)$name
                  cell_states <- V(res$subtree)$cell_state
                  cell_pseudotime <- V(res$subtree)$pseudotime
                  ordering_df <- data.frame(sample_name = cell_names, cell_state = cell_states, pseudo_time = cell_pseudotime,stringsAsFactors = F)
                  ordering_df <- arrange(ordering_df, pseudo_time)
                  return(ordering_df)
            }
            
      }    















})