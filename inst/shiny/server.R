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
library(TSP)
library(combinat)
library(mgcv)
library(gplots)

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
                  tagList(
                        h5("Instructions:"),
                        p("Single cell data should be prepared in a matrix-like data format. Each row corresponds to a gene/feature and each column corresponds to a single cell."),
                        p("Notice that each row should have the same number of entries, especially for the header (first row)"),
                        p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to make changes."),
                        h5("A typical example of tab-delimited file format is as follows:"),
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
                  Maindata$fullrawlogdata <- tmpdata
                  tmpdata <- tmpdata[rowMeans(tmpdata > as.numeric(input$Preprocessexpvalcutoff)) > as.numeric(input$Preprocessexppercent),!colnames(tmpdata) %in% input$Preprocessexclude]
                  tmprowcv <- apply(tmpdata,1,sd)/rowMeans(tmpdata)
                  Maindata$fullprocdata <- tmpdata[tmprowcv > as.numeric(input$Preprocesscvcutoff),]
            }
      })
      
      output$Preprocessexcludeui <- renderUI(selectInput("Preprocessexclude","Selected unwanted cells",multiple = T,choices = colnames(Maindata$rawdata)[-1]))
      
      output$Preprocesshowproctable <- renderDataTable({
            if (!is.null(Maindata$fullprocdata)) {
                  tmpdata <- Maindata$fullprocdata[row.names(Maindata$fullprocdata) %in% input$Preprocesschoosegene,,drop=F]
                  tmpdata <- cbind(row.names(tmpdata),tmpdata)
                  colnames(tmpdata)[1] <- "GENE"
                  tmpdata
            }
      })
      
      output$Preprocesshowgenesummary <- renderPrint(summary(t(Maindata$fullprocdata[row.names(Maindata$fullprocdata) %in% input$Preprocesschoosegene,,drop=F])))
      
      output$Preprocesshowcellsummary <- renderPrint(summary(Maindata$fullprocdata))
      
      output$Preprocessstatusui <- renderUI({
            if (!is.null(Maindata$procdata)) {
                  tagList(
                        p(paste(nrow(Maindata$fullprocdata),"rows out of",nrow(Maindata$rawdata),"rows are preserved, which is",round(nrow(Maindata$fullprocdata)/nrow(Maindata$rawdata)*100,3),"percent.")),
                        p(h5("Choose genes to display")),
                        selectInput("Preprocesschoosegene","Choose gene",choices = row.names(Maindata$fullprocdata),selected = head(row.names(Maindata$fullprocdata),n=2),multiple = T),
                        p(h5("Data Table:")),
                        dataTableOutput("Preprocesshowproctable"),
                        p(h5("Gene level summary table:")),
                        verbatimTextOutput("Preprocesshowgenesummary"),
                        p(h5("Cell level summary table (for all retained genes after filtering):")),
                        verbatimTextOutput("Preprocesshowcellsummary")
                  )
            }
      })
      
      ### Cell Ordering ###
      
      observe({
            if (input$MainMenu != "Ordering")
                  updateRadioButtons(session,"Orderingchoosestep","",list("Step 1: Reduce dimension"="reduction","Step 2: Calculate pseudotime"="ptime","Step 3: Manually adjust starting point (optional)"="start","Save results (optional)"="save"),selected = "reduction")
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
                        p("Pseudotemporal cell ordering data should contain exactly three columns. First column: cell name (character value). Second column: cell state (numeric value). Third column: Pseudotime (numeric value)"),
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
                  plot(Maindata$fullreducepc$sdev[1:20]/sum(Maindata$fullreducepc$sdev),xlab="PC",ylab="Variance proportion",main="PC variance proportion")
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
                                         selectInput("Orderingptimetrimbranchselect","Select branch id on original plot",choices = sort(unique(Maindata$scapdata$State)),multiple = T)
                        ),
                        conditionalPanel(condition="input.Orderingptimetrimmethod=='cell'",
                                         selectInput("Orderingptimetrimcellselect","Select cell name",choices = colnames(Maindata$procdata),multiple = T)
                        ),
                        conditionalPanel(condition="input.Orderingptimetrimmethod=='expression'",
                                          helpText("Cells meeting following criterion simultaneously will be trimmed. Refer to 'trim expression' tab on the right."),
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
                  Maindata$reduceres <- Maindata$fullreduceres[,!colnames(Maindata$fullreduceres) %in% Maindata$trimlist]
                  Maindata$procdata <- Maindata$fullprocdata[!row.names(Maindata$fullprocdata) %in% Maindata$trimlist,] 
                  Maindata$rawlogdata <- Maindata$fullrawlogdata[!colnames(Maindata$fullrawlogdata) %in% Maindata$trimlist,] 
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
            if (input$Orderingptimechoosemethod=="MST") {
                  tagList(
                        sliderInput("OrderingMSTpathnum","Choose number of paths",min=1,max=20,step=1,value=3),
                        checkboxInput("OrderingMSTreversetf","Reverse the ordering",value=F),
                        checkboxInput("OrderingMSTshow_tree","Show tree",value = T),
                        checkboxInput("OrderingMSTshow_backbone","Show diameter path (backbone)",value = T),
                        selectInput("OrderingMSTbackbone_color","Select backbone color",choices = c("black","red","blue","green","yellow")),
                        checkboxInput("OrderingMSTshow_cell_names","Show cell names",value = T),
                        conditionalPanel(condition="input.OrderingMSTshow_cell_names=='1'",textInput("OrderingMSTcell_name_size","Choose the size of cell name labels",value = 3)),
                        checkboxInput("OrderingMSTmarkertf","Use marker gene to define node size",value=F),
                        conditionalPanel(condition="input.OrderingMSTmarkertf=='1'",uiOutput("OrderingMSTmarkerui")),
                        checkboxInput("OrderingMSTrootcelltf","Choose root cell",value=F),
                        conditionalPanel(condition="input.OrderingMSTrootcelltf=='1'",uiOutput("OrderingMSTrootcellui")))
                  
            } else if (input$Orderingptimechoosemethod=="TSP") {
                  tagList(
                        sliderInput("OrderingTSPclunum","Choose number of cell clusters",min=1,max=20,step=1,value=3),
                        checkboxInput("OrderingTSPshow_tree","Show tree",value = T),
                        checkboxInput("OrderingTSPshow_cell_names","Show cell name",value=T),
                        conditionalPanel(condition="input.OrderingTSPshow_cell_names=='1'",textInput("OrderingTSPcell_name_size","Choose the size of cell name labels",value = 3)),
                        checkboxInput("OrderingTSPmarkertf","Use marker gene to define node size",value=F),
                        conditionalPanel(condition="input.OrderingTSPmarkertf=='1'",uiOutput("OrderingTSPmarkerui"))
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
                  tagList(
                        sliderInput("Orderingptimezoominxaxis","X-axis range",min=min(Maindata$fullreduceres[1,]),max=max(Maindata$fullreduceres[1,]),value=c(min(Maindata$fullreduceres[1,]),max(Maindata$fullreduceres[1,]))),
                        sliderInput("Orderingptimezoominyaxis","Y-axis range",min=min(Maindata$fullreduceres[2,]),max=max(Maindata$fullreduceres[2,]),value=c(min(Maindata$fullreduceres[2,]),max(Maindata$fullreduceres[2,])))
                  )
            }
      })
      
      observe({
            if (!is.null(Maindata$reduceres)) {
                  if (input$Orderingptimechoosemethod=="MST" && !is.null(input$OrderingMSTpathnum) && !is.null(Maindata$reduceres)) { 
                        if (!is.null(input$OrderingMSTrootcell) && input$OrderingMSTrootcelltf && nchar(input$OrderingMSTrootcell) != 0) {
                              root_cell <- input$OrderingMSTrootcell 
                        } else {
                              root_cell <- NULL      
                        }
                        num_paths <- as.numeric(input$OrderingMSTpathnum)
                        reverse <- input$OrderingMSTreversetf                            
                        dp <- as.matrix(dist(t(Maindata$reduceres)))
                        gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
                        dp_mst <- minimum.spanning.tree(gp)
                        Maindata$dp_mst <- dp_mst
                        next_node <<- 0
                        res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
                        cc_ordering <- extract_good_branched_ordering(res$subtree, res$root, dp, num_paths, reverse)
                        colnames(cc_ordering) <- c("sample_name","State","Pseudotime")
                        cc_ordering$Pseudotime <- cc_ordering$Pseudotime/max(cc_ordering$Pseudotime) * as.numeric(input$Orderingptimescale)
                        Maindata$scapdata <- cc_ordering
                  } else if (input$Orderingptimechoosemethod=="TSP" && !is.null(input$OrderingTSPclunum)) {
                        set.seed(12345)
                        datadist <- dist(t(Maindata$reduceres))
                        dataTSP <- TSP(datadist)
                        datasolve <- solve_TSP(dataTSP)
                        datalab <- labels(datasolve)
                        distmat <- as.matrix(datadist)
                        alldist <- sapply(1:(length(datalab)-1), function(x) {
                              distmat[datalab[x],datalab[x+1]]
                        })
                        
                        state <- kmeans(t(Maindata$reduceres),centers=as.numeric(input$OrderingTSPclunum))$cluster
                        ptime <- c(0,cumsum(alldist))
                        ptime <- ptime/max(ptime) * as.numeric(input$Orderingptimescale)
                        Maindata$scapdata <- Maindata$pdata <- data.frame(sample_name=datalab,State=state[datalab],Pseudotime=ptime,stringsAsFactors = F)
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
      
      output$OrderingMSTmarkerui <- renderUI({
            selectInput("OrderingMSTmarker","select marker gene used for node size",choices = row.names(Maindata$fullrawlogdata))
      })
      
      output$OrderingTSPmarkerui <- renderUI({
            selectInput("OrderingTSPmarker","select marker gene used for node size",choices = row.names(Maindata$fullrawlogdata))
      })
      
      output$OrderingMSTrootcellui <- renderUI({
            selectInput("OrderingMSTrootcell","select root cell",choices = colnames(Maindata$procdata))
      })
      
      MSTdrawplot <- function(xlabtext="Component 2",ylabtext="Component 1",titletext="Pseudotime ordering plot") {
            x = 1
            y = 2
            color_by = "State"
            show_tree = input$OrderingMSTshow_tree
            show_backbone = input$OrderingMSTshow_backbone
            backbone_color = input$OrderingMSTbackbone_color
            if (!is.null(input$OrderingMSTmarkertf) && input$OrderingMSTmarkertf) {
                  markers <- input$OrderingMSTmarker
            } else {
                  markers <- NULL
            }
            show_cell_names = input$OrderingMSTshow_cell_names
            cell_name_size = as.numeric(input$OrderingMSTcell_name_size)
            
            lib_info_with_pseudo <- Maindata$scapdata
            lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
            S_matrix <- Maindata$reduceres
            ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
            colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
            ica_space_df$sample_name <- row.names(ica_space_df)
            ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
            dp_mst <- Maindata$dp_mst
            edge_list <- as.data.frame(get.edgelist(dp_mst))
            colnames(edge_list) <- c("source", "target")
            edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", ICA_dim_2 = "source_ICA_dim_2"))
            edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", ICA_dim_2 = "target_ICA_dim_2"))
            diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, weights = NA)]$name))
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
            g <- g + theme_minimal(base_size = as.numeric(input$Orderingsaveplotfontsize)) + theme(panel.border = element_blank(), axis.line = element_line()) + 
                  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
                  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
                  labs(title=titletext) + ylab(ylabtext) + xlab(xlabtext) + theme(legend.position = "top", 
                                                                                  legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
                  theme(panel.background = element_rect(fill = "white"))
            
            if (input$Orderingptimezoomintf) {
                  g <- g + coord_cartesian(xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
            } else {
                  sc1 <- 0.05 * (max(Maindata$fullreduceres[1,])-min(Maindata$fullreduceres[1,]))
                  sc2 <- 0.05 * (max(Maindata$fullreduceres[2,])-min(Maindata$fullreduceres[2,]))
                  g <- g + coord_cartesian(xlim = c(min(Maindata$fullreduceres[1,]) - sc1, max(Maindata$fullreduceres[1,])+sc1),ylim=c(min(Maindata$fullreduceres[2,])-sc2, max(Maindata$fullreduceres[2,])+sc2))                  
            }                  
            g 
      }
      
      TSPdrawplot <- function(xlabtext="Component 2",ylabtext="Component 1",titletext="Pseudotime ordering plot") {
            x = 1
            y = 2
            color_by = "State"
            
            show_tree = input$OrderingTSPshow_tree
            show_cell_names = input$OrderingTSPshow_cell_names
            cell_name_size = as.numeric(input$OrderingTSPcell_name_size)
            if (!is.null(input$OrderingTSPmarkertf) && input$OrderingTSPmarkertf) {
                  markers <- input$OrderingTSPmarker
            } else {
                  markers <- NULL
            }
            
            lib_info_with_pseudo <- Maindata$scapdata
            lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
            S_matrix <- Maindata$reduceres
            ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
            colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
            ica_space_df$sample_name <- row.names(ica_space_df)
            ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name")
            
            TSPorder <- Maindata$scapdata[,1]
            edge_list <- data.frame(TSPorder,c(TSPorder[2:length(TSPorder)],TSPorder[1]),stringsAsFactors = F)
            colnames(edge_list) <- c("source", "target")
            edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", ICA_dim_2 = "source_ICA_dim_2"))
            edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", all = TRUE)
            edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", ICA_dim_2 = "target_ICA_dim_2"))
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
            if (show_cell_names) {
                  g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
            }
            g <- g + theme_minimal(base_size = as.numeric(input$Orderingsaveplotfontsize))+theme(panel.border = element_blank(), axis.line = element_line()) + 
                  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
                  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
                  labs(title=titletext) + ylab(ylabtext) + xlab(xlabtext) + theme(legend.position = "top", 
                                                                                  legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
                  theme(panel.background = element_rect(fill = "white"))
            
            if (input$Orderingptimezoomintf) {
                  g <- g + coord_cartesian(xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
            } else {
                  sc1 <- 0.05 * (max(Maindata$fullreduceres[1,])-min(Maindata$fullreduceres[1,]))
                  sc2 <- 0.05 * (max(Maindata$fullreduceres[2,])-min(Maindata$fullreduceres[2,]))
                  g <- g + coord_cartesian(xlim = c(min(Maindata$fullreduceres[1,]) - sc1, max(Maindata$fullreduceres[1,])+sc1),ylim=c(min(Maindata$fullreduceres[2,])-sc2, max(Maindata$fullreduceres[2,])+sc2))
            }  
            g       
      }
      
      output$Orderingptimeshowplot <- renderPlot({
            if (!is.null(Maindata$reduceres)) {
                  if (input$Orderingptimechoosemethod=="MST" && !is.null(Maindata$dp_mst)) {
                        MSTdrawplot()
                  } else if (input$Orderingptimechoosemethod=="TSP" && !is.null(input$OrderingTSPshow_cell_names)) {
                        TSPdrawplot()
                  } else if (input$Orderingptimechoosemethod=="PC") {
                        coord <- t(Maindata$reduceres)
                        if (input$Orderingptimezoomintf) {
                              plot(Maindata$lpcobj,datcol=Maindata$pdata$State,datpch=19,xlim = c(as.numeric(input$Orderingptimezoominxaxis[1]), as.numeric(input$Orderingptimezoominxaxis[2])),ylim=c(as.numeric(input$Orderingptimezoominyaxis[1]), as.numeric(input$Orderingptimezoominyaxis[2])))
                        } else {
                              plot(Maindata$lpcobj,datcol=Maindata$pdata$State,datpch=19)
                        }
                        if (input$OrderingPCshowcellname)
                              text(coord,row.names(coord))     
                  }
            }
      })
      
      output$Orderingptimeshowptime <- renderDataTable(Maindata$scapdata)
      
      #Step 3: starting point
      
      output$Orderingstartmainui <- renderUI({
            if (input$Orderingptimechoosemethod == "TSP") {
                  tabsetPanel(
                        tabPanel("Ordering", uiOutput("Orderingstartshowheatui")
                        ),
                        tabPanel("Genesets",dataTableOutput("Orderingstartsuggestshowgeneset")),
                        tabPanel("Pseudotime",dataTableOutput("Orderingstartsuggestshowpseudotime"))
                  )
            } else {
                  h5("This function is not available for minimum spanning tree approach.")
            }
      })
      
      output$Orderingstartchoosemarkerui <- renderUI({
            if (!is.null(Maindata$procdata))
                  selectInput("Orderingstartchoosemarker","select genes",choices = row.names(Maindata$rawdata),multiple = T)
      })
      
      observe({
            if (!is.null(input$Orderingstartaddbutton) && input$Orderingstartaddbutton > 0) {
                  isolate({
                        if (is.null(Maindata$fullstartgeneset)) {
                              Maindata$fullstartgeneset <- data.frame(gene=input$Orderingstartchoosemarker,trend=input$Orderingstartchoosegenetrend,genesetname="Geneset 1",stringsAsFactors = F)
                        } else {
                              genesetid <- max(as.numeric(sub("Geneset ","",unique(Maindata$fullstartgeneset$genesetname))))+1
                              Maindata$fullstartgeneset <- rbind(Maindata$fullstartgeneset,data.frame(gene=input$Orderingstartchoosemarker,trend=input$Orderingstartchoosegenetrend,genesetname=paste("Geneset",genesetid),stringsAsFactors = F))
                        }
                  })                        
            }          
      })
      
      output$Orderingstartincludegenesetui <- renderUI({
            selectInput("Orderingstartincludegeneset","Select geneset to include",choices = unique(Maindata$fullstartgeneset$genesetname),multiple = T,selected = unique(Maindata$fullstartgeneset$genesetname))
      })
      
      observe({
            if (!is.null(Maindata$fullstartgeneset) && !is.null(input$Orderingstartincludegeneset)) {
                  Maindata$startgeneset <- Maindata$fullstartgeneset[Maindata$fullstartgeneset$genesetname %in% input$Orderingstartincludegeneset,]      
                  increaseexpr <- NULL
                  increasename <- NULL
                  decreaseexpr <- NULL
                  decreasename <- NULL
                  for (i in unique(Maindata$startgeneset$genesetname)) {
                        tmp <- Maindata$startgeneset[Maindata$startgeneset$genesetname == i,]             
                        if (tmp[1,2] != "No") {
                              tmpexpr <- Maindata$rawlogdata[tmp[,1],,drop=F]
                              if (input$Orderingstartscalegeneset)
                                    tmpexpr <- t(scale(t(tmpexpr)))
                              tmpexpr <- colMeans(tmpexpr)
                              
                              if (tmp[1,2] == "increasing") {
                                    increaseexpr <- rbind(increaseexpr,tmpexpr)
                                    increasename <- c(increasename, i)
                              } else if (tmp[1,2] == "decreasing") {
                                    decreaseexpr <- rbind(decreaseexpr,tmpexpr)
                                    decreasename <- c(decreasename, i)
                              }
                        }                        
                  }
                  row.names(increaseexpr) <- increasename
                  row.names(decreaseexpr) <- decreasename
                  Maindata$increaseexpr <- increaseexpr
                  Maindata$decreaseexpr <- decreaseexpr                                          
            }
      })
      
      observe({
            if (input$Orderingchoosestep=='start' && !is.null(input$Orderingstartslider) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
                  pdata <- Maindata$pdata
                  pdata <- pdata[order(pdata$Pseudotime),]
                  start <- as.numeric(input$Orderingstartslider)
                  if (start == 1) {
                        order <- 1:nrow(pdata)
                  } else {
                        order <- c(start:(nrow(pdata)),1:(start-1))
                  }
                  if (input$Orderingstartflip)
                        order <- rev(order)
                  pdata <- pdata[order,]
                  datadist <- dist(t(Maindata$reduceres))
                  distmat <- as.matrix(datadist)
                  alldist <- sapply(1:(nrow(pdata)-1), function(x) {
                        distmat[pdata[x,1],pdata[x+1,1]]
                  })
                  ptime <- c(0,cumsum(alldist))
                  ptime <- ptime/max(ptime) * as.numeric(input$Orderingptimescale)
                  pdata[,3] <- ptime
                  Maindata$scapdata <- pdata
            }
      })
      
      observe({
            if (input$Orderingchoosestep=='start' && !is.null(Maindata$pdata) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
                  pdata <- Maindata$pdata[,-2]
                  pdata <- pdata[order(pdata$Pseudotime),]
                  if (is.null(Maindata$decreaseexpr)) {
                        allexpr <- Maindata$increaseexpr
                  } else {
                        allexpr <- rbind(Maindata$increaseexpr,-Maindata$decreaseexpr)     
                  }
                  allres <- NULL
                  for (flip in c(TRUE,FALSE)) {
                        for (start in 1:nrow(pdata)) {
                              if (start == 1) {
                                    order <- 1:nrow(pdata)
                              } else {
                                    order <- c(start:(nrow(pdata)),1:(start-1))
                              }
                              if (flip)
                                    order <- rev(order)            
                              cellorder <- pdata[order,1]                              
                              pcos <- sum(apply(allexpr[,cellorder,drop=F],1,function(expr) {
                                    sum(sapply(1:(length(expr)-1),function(x) {
                                          sum(expr[(x+1):length(expr)] - expr[x])
                                    }))                                    
                              }))
                              allres <- rbind(allres,c(start,flip,pcos))
                        }
                  }
                  Maindata$Orderingsuggestres <- allres
            }
      })
      
      output$Orderingstartshowheatmap <- renderPlot({
            if (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr)) {
                  allexpr <- rbind(Maindata$increaseexpr,Maindata$decreaseexpr)
                  par(mar=c(0,0,0,0))
                  if (nrow(allexpr) == 1) {
                        plot(allexpr[,Maindata$scapdata[,1]])
                  } else {
                        image(t(allexpr[,Maindata$scapdata[,1]]),axes=F,col=bluered(100))      
                        names <- rev(row.names(allexpr))
                        dim <- length(names)
                        pos <- seq(0,1,length.out=dim)
                        for (i in 1:dim) {
                              text(0.8,pos[i],names[i],cex=1.3)
                        }
                  }
            }            
      })
      
      output$Orderingstartshowheatui <- renderUI({
            if (!is.null(Maindata$procdata) && (!is.null(Maindata$increaseexpr) || !is.null(Maindata$decreaseexpr))) {
                  tagList(                  
                        p(actionButton("Orderingstartsetoptimalbutton","Set optimal value")),
                        checkboxInput("Orderingstartflip","Flip the ordering"),
                        sliderInput("Orderingstartslider","Slide to select the starting point",min=1,max=ncol(Maindata$procdata),step=1,value=1,width='800px'),
                        plotOutput("Orderingstartshowheatmap",width='800px'),
                        p(textOutput("Orderingstartsuggestshowres"))
                  )      
            }
      })
      
      observe({
            if (!is.null(input$Orderingstartsetoptimalbutton) && input$Orderingstartsetoptimalbutton) {
                  isolate({
                        maxpos <- Maindata$Orderingsuggestres[which.max(Maindata$Orderingsuggestres[,3]),]
                        updateCheckboxInput(session = session,"Orderingstartflip",value=as.logical(maxpos[2]))
                        updateSliderInput(session = session, "Orderingstartslider",value=as.integer(maxpos[1]))
                  })
            }
      })
      
      output$Orderingstartsuggestshowres <- renderText({
            if (!is.null(input$Orderingstartsetoptimalbutton)) {
                  pcos <- Maindata$Orderingsuggestres[Maindata$Orderingsuggestres[,1]== as.numeric(input$Orderingstartslider) & Maindata$Orderingsuggestres[,2] == input$Orderingstartflip,3]
                  paste("Pseudotemporal cell ordering score:",pcos)      
            }
      })
      
      output$Orderingstartsuggestshowgeneset <- renderDataTable(Maindata$fullstartgeneset)
      
      output$Orderingstartsuggestshowpseudotime <- renderDataTable(Maindata$scapdata)
      
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
            if (input$Orderingptimechoosemethod=="MST" && !is.null(Maindata$dp_mst)) {
                  MSTdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle)
            } else if (input$Orderingptimechoosemethod=="TSP" && !is.null(input$OrderingTSPshow_cell_names)) {
                  TSPdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle)
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
                  if (input$Orderingptimechoosemethod=="MST") {
                        print(MSTdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle))
                  } else if (input$Orderingptimechoosemethod=="TSP") {
                        print(TSPdrawplot(xlabtext=input$Orderingsaveplotchangexlab,ylabtext=input$Orderingsaveplotchangeylab,titletext=input$Orderingsaveplotchangetitle))
                  }
                  dev.off()
            }
      )
      
      ### Changpoint ###
      
      #choose uploaded pdata or sca pdata
      
      output$Changpointpdataselectui <- renderUI({
            if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata)) {
                  radioButtons("Changepointpdataselect","Select the cell ordering you want to use",choices = list("Uploaded pseudotime"="Uploaded","SCA generated pseudotime"="SCA"))     
            }            
      })
      
      observe({
            if (input$MainMenu=="Changepoint") {
                  if (!is.null(Maindata$uploadpdata) && !is.null(Maindata$scapdata) && !is.null(input$Changepointpdataselect)) {
                        if (input$Changepointpdataselect == 'Uploaded') {
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
            if (!is.null(input$Changepointfilterbutton) && input$Changepointfilterbutton > 0) {
                  isolate({
                        tmpdata <- Maindata$procdata[,Maindata$finalpdata[,1]]
                        ptime <- Maindata$finalpdata[,3]
                        i <- 0
                        filterpval <- apply(tmpdata,1,function(x) {                                                            
                              i <<- i + 1
                              Maindata$filtercalculatecount <<- i
                              anova(gam(x~s(ptime,k=3)),lm(x~1))$s.table[4]
                        })
                        Maindata$filteradjpval <- p.adjust(filterpval,method="fdr")                        
                  })
            }
      })
      
      observe({
            if (!is.null(Maindata$filteradjpval)) {
                  Maindata$filterdata <- Maindata$procdata[Maindata$filteradjpval < as.numeric(input$Changpointfilterfdrval),Maindata$finalpdata[,1]]
                  Maindata$filterresalltable <- data.frame(GENE=row.names(Maindata$procdata),adjusted.p.value=Maindata$filteradjpval)      
                  Maindata$filterrestable <- data.frame(GENE=row.names(Maindata$filterdata),adjusted.p.value=Maindata$filteradjpval[Maindata$filteradjpval < as.numeric(input$Changpointfilterfdrval)])      
            }
      })
      
      output$Changpointfiltersummaryui <- renderUI({
            if (!is.null(Maindata$filterdata) && nrow(Maindata$filterdata)!=0)
                  p(paste(nrow(Maindata$filterdata),"genes out of",nrow(Maindata$procdata),"genes are differentially expressed, which is",round(nrow(Maindata$filterdata)/nrow(Maindata$procdata)*100,3),"percent."))
      })      
      
      output$Changpointfiltercalculatecomplete <- renderText({
            if (!is.null(Maindata$filtercalculatecount) && Maindata$filtercalculatecount == nrow(Maindata$procdata))
                  "Calculation Completed!"
      })
      
      output$Changpointmainui <- renderUI({
            if (input$Changepointchoosemethod == "Filter") {
                  tagList(
                        textOutput("Changpointfiltercalculatecomplete"),
                        uiOutput("Changpointfiltersummaryui"),
                        dataTableOutput("ChangepointFiltershowresult")
                  )
            } else {
                  plotOutput("Changepointviewshowplot",width = "800px",height = "800px")
            }
            
      })
      
      output$ChangepointFiltershowresult <- renderDataTable({
            if (!is.null(Maindata$filterdata) && nrow(Maindata$filterdata)!=0)
                  if (input$Changpointfiltershowresultopt == "all") {
                        Maindata$filterresalltable   
                  } else {
                        Maindata$filterrestable
                  }                  
      })
      
      output$Changpointsavepvaltable <- downloadHandler(
            filename = function() { paste0("Test_for_differentially_expressed.",ifelse(input$Changepointsavepvaltabletype == 'txt','txt','csv')) },
            content = function(file) {
                  if (input$Changepointsavepvaltabletype == 'txt') {
                        if (input$Changpointfiltershowresultopt == "all") {
                              write.table(Maindata$filterresalltable,file,row.names=F,quote=F)                              
                        } else {
                              write.table(Maindata$filterrestable,file,row.names=F,quote=F)     
                        }                          
                  } else {
                        if (input$Changpointfiltershowresultopt == "all") {
                              write.csv(Maindata$filterresalltable,file,row.names=F,quote=F)                              
                        } else {
                              write.csv(Maindata$filterrestable,file,row.names=F,quote=F)     
                        }
                  }                  
            }
      )
      
      output$Changepointviewgeneselectui <- renderUI({            
            selectInput("ChangepointViewgeneselect","Select gene",choices = row.names(Maindata$procdata))
      })
      
      output$Changepointviewshowplot <- renderPlot({
            if (nchar(input$ChangepointViewgeneselect)!=0) {
                  tmpdata <- Maindata$procdata[input$ChangepointViewgeneselect,Maindata$finalpdata[,1]]
                  ptime <- Maindata$finalpdata[,3]
                  if (input$Changepointviewshowstatus) {
                        plot(ptime,tmpdata,xlab="Pseudotime",ylab="Expression value",pch=19,col=Maindata$finalpdata[,2])      
                  } else {
                        plot(ptime,tmpdata,xlab="Pseudotime",ylab="Expression value",pch=19)      
                  }
            }            
      })
      
      ###  Miscellaneous ###
            
      observe({
            if (input$MainMenu != "Ordering")
                  updateRadioButtons(session,"Orderingchoosestep","",list("Step 1: Reduce dimension"="reduction","Step 2: Calculate pseudotime"="ptime","Step 3: Manually adjust starting point (optional)"="start","Save results (optional)"="save"),selected = "reduction")
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
                        p("Sub-population id could only be 0 and 1. 0 stands for subpopulation collected at an early time point and 1 stands for subpopulaiton collected at a latter time point"),
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
                  sum(sapply(1:(length(scoreorder)-1),function(x) {
                        sum(scoreorder[(x+1):length(scoreorder)] - scoreorder[x])
                  })) / (sum(scoreorder==1)*sum(scoreorder==0))
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
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      ###monocle internal functions
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