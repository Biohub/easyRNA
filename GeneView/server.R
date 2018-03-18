
library(shiny)
library(pheatmap)

# Define server logic
shinyServer(function(input, output, session) {
  
  # load dataset
  dataset <- reactive({
    # validate(
    #   need(input$rData$datapath != "" | length(input$sData) != 0, "")
    # )
    if(is.null(input$rData) && isTRUE(input$useData)){
      return(mget(load(input$sData)))
    } else {
      if (!is.null(input$rData) && file.exists(paste0(input$rData$datapath))){
        return(mget(load(input$rData$datapath)))
      } else {
        return(NULL)
      }
    }
    # if(file.exists(paste0(input$rData$datapath))){
    #   mget(load(input$rData$datapath))
    # } else {
    #   if (isTRUE(input$useData)){
    #     mget(load(input$sData))
    #   } else {
    #     return(NULL)
    #   }
    # }
  })
  
  exp <- reactive({
    if(is.null(dataset())){
      return(NULL)
    }
    dataset()$exp
  })
  meta <- reactive({
    if(is.null(dataset())){
      return(NULL)
    }
    dataset()$meta
  })
  
  meta.attr <- reactive({
    if(is.null(meta())){
      return(NULL)
    }
    names(meta())
  })
  
  genes <- reactive({
    if(is.null(input$gene)){
      return(NULL)
    }
    this.gene <- unlist(strsplit(x = input$gene, split = "\\n|\\,|\\,\\s|\\s", perl = T))
    this.gene[ this.gene %in% row.names(exp())]
  })
  
  observe({
    updateSelectInput(session = session, inputId = "attr", choices = meta.attr())
  })
  
  observe({
    this.values <- NULL
    if(!is.null(input$attr)){
      this.values <- unique(meta()[, which(names(meta()) %in% input$attr)])
    }
    updateSelectInput(session = session, inputId = "attr.values", choices = this.values)
  })
  
  output$attrUI <- renderUI({
    if(!input$addAttr){
      return(NULL)
    }
    tagList(
      selectInput(inputId = "attr2", label = "Meta attributes 2", choices = meta.attr())
    )
  })

  output$attrValuesUI <- renderUI({
    if(!input$addAttr){
      return(NULL)
    }
    tagList(
      selectInput(inputId = "attr.values2", label = "Attribute values 2", choices = unique(meta()[, input$attr2]), multiple = T)
    )
  })
  
  attrSelect <- reactive({
    if(is.null(input$attr.values)){
      return(NULL)
    }
    input$attr.values
  })
  
  attrSelect2 <- reactive({
    if(is.null(input$attr.values2)){
      return(NULL)
    }
    input$attr.values2
  })
  
  attrSelectIndex <- reactive({
    attr1.index <- unlist(sapply(attrSelect(), function(x){
      which(meta()[, input$attr] %in% x)
    }, simplify = T))
    if(input$addAttr){
      attr2.index <- unlist(sapply(attrSelect2(), function(x){
        which(meta()[, input$attr2] %in% x)
      }, simplify = T))
      return(intersect(attr1.index, attr2.index))
    }
    return(attr1.index)
  })
  
  thisData <- reactive({
    data.exp <- exp()[genes(), attrSelectIndex()]
    # data.var <- unlist(apply(data.exp, 1, var))
    # data.exp[data.var > 0,]
    return(data.exp)
  })
  
  # clusterRow <- reactive({
  #   input.dist.row <- as.dist(1-cor(t(thisData()), method = "pearson"))
  #   input.clust.row <- hclust(input.dist.row, method = "ward.D2")
  # })
  # 
  # clusterCol <- reactive({
  #   input.dist.col <- as.dist(1-cor(thisData(), method = "pearson"))
  #   input.clust.col <- hclust(input.dist.col, method = "ward.D2")
  # })
  
  heatmap.plot <- function(){
    this.data <- thisData()
    if(is.null(this.data)){
      return(NULL)
    }
    # remove genes with sd = 0
    # data.var <- unlist(apply(this.data, 1, var))
    # this.data <- this.data[data.var > 0,]
    # color
    col.panel <- colorRampPalette(colors = c("blue", "white","red"))
    # cluster
    input.clust.row <- "row" %in% input$heatmap.clust.args
    input.clust.col <- "col" %in% input$heatmap.clust.args
    # annotation
    col.annotation <- data.frame( "Cell type" = meta()[attrSelectIndex(), input$attr], row.names= row.names(meta()[attrSelectIndex(),]))
    
    pheatmap(this.data, color=col.panel(100),
             cluster_rows = input.clust.row,
             cluster_cols = input.clust.col,
             clustering_distance_rows = input$heatmap.dist.args,
             clustering_distance_cols = input$heatmap.dist.args,
             clustering_method = input$heatmap.clust.method,
             legend = T, show_colnames = F, show_rownames = T, 
             annotation_col = col.annotation
    )
  }
  
  bar.plot <- function(){
    this.data <- thisData()
    if(is.null(this.data)){
      return(NULL)
    }
    
    # colors
    library(RColorBrewer)
    cols <- colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(input$attr.values))
    celltypes <- meta()[attrSelectIndex(), input$attr]
    cols.cells <- sapply(celltypes, function(x){
      cols[match(x, input$attr.values)]
    })
    
    # plot(x=0,y=0, type = "n", axes = F, xlab = "", ylab = "",
    #      xlim = c(0,ncol(this.data) +1), ylim = c(0, nrow(this.data) +1))
    # for(i in 1:nrow(this.data)){
    #   segments(x0 = 0, y0 = i, x1 = ncol(this.data), y1 = i, col = "grey")
    #   segments(x0 = 0, y0 = i, x1 = 0, y1 = ncol(this.data) +1, col = "grey")
    #   segments(x0 = 1:ncol(this.data), y0 = i, )
    # }
    if(is.null(nrow(this.data))){
      par(las=1)
      barplot(height = this.data, names.arg = F, col = cols.cells, border = cols.cells)
    } else {
      par(mar=c(0,4,2,2), oma=c(2,2,2,2), las=1)
      layout(matrix(1:nrow(this.data), ncol = 1))
      for(i in 1:nrow(this.data)){
        barplot(height = this.data[i,], names.arg = F, col = cols.cells, border = cols.cells)
        mtext(x=3, text = row.names(this.data)[i])
      }
    }
  }
  
  heatmap.width <- reactive({
    as.numeric(input$heatmap.width)
  })
  
  heatmap.height <- reactive({
    as.numeric(input$heatmap.height)
  })
  
  output$plot <- renderPlot(
    expr = {
      heatmap.plot()
    },
    width = function(){heatmap.width() * 75},
    height = function(){heatmap.height() * 75}
  )
  
  output$heatmap.download <- downloadHandler(
    filename = function(){
      paste("heatmap", Sys.Date(), "pdf", sep = ".")
    },
    content = function(file){
      pdf(file = file, width = heatmap.width(), height = heatmap.height(), onefile=FALSE)
      heatmap.plot()
      dev.off()
    }
  )
  
  output$barplot <- renderPlot(
    expr = {
      bar.plot()
    },
    width = function(){heatmap.width() * 75},
    height = function(){heatmap.height() * 75}
  )
  
  output$text <- renderText({
    input$attr.values
  })
  
})
