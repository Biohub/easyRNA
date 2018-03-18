
library(shiny)
options(shiny.maxRequestSize=100*1024^2) # maximum upload size 100M

shinyUI(fluidPage(
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      # upload Rdata
      fileInput(inputId = "rData", label = "Upload Rdata (expression and meta data)"),
      selectInput(inputId = "sData", label = "Or select available dataset", choices = list("mHSPC_transplantation" = "data/background_transplant.tpm_and_meta.RData")),
      checkboxInput(inputId = "useData", label = "Use avalaible dataset", value = F),
      # --- line ---
      tags$hr(),
      
      # meta data specification
      selectInput(inputId = "attr", label = "Meta attributes", choices = NULL),
      selectInput(inputId = "attr.values", label = "Attribute values", choices = NULL, multiple = T),
      
      ## Add second attribute
      checkboxInput(inputId = "addAttr", label = "add another attribute", value = F),
      # dynamic generation of meta attributes 
      uiOutput(outputId = "attrUI"),
      # dynamic generation of values for attributes
      uiOutput(outputId = "attrValuesUI"),
      # genes
      textAreaInput(inputId = "gene", label = "Genes (by row or comma delimited)", resize = "vertical"),
      fluidRow(
        column(
          width = 6,
          numericInput(inputId = "heatmap.width", label = "Figure width", value = 7)
        ),
        column(
          width = 6,
          numericInput(inputId = "heatmap.height", label = "Figure Height", value = 4)
        )
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Heatmap",
          checkboxGroupInput(inputId = "heatmap.clust.args", label = "Clustering",inline = T, selected = "row",
                             choices = list("By row" = "row", "By column" = "col")),
          radioButtons(inputId = "heatmap.dist.args", label = "Distance", 
                       choices = list("Correlation" = "correlation", "Euclidean" = "euclidean", "Manhattan" = "manhattan", "Minkowski" = "minkowski"), inline = T),
          radioButtons(inputId = "heatmap.clust.method", label = "Clustering method",  selected = "ward.D2",
                       choices = list("Ward.D2" = "ward.D2", "Complete" = "complete", "Single" = "single", "Average" = "average", "Centroid" = "centroid"), inline = T),
          icon("download"),
          downloadLink(outputId = "heatmap.download", label = "Download plot (pdf)"),
          plotOutput("plot")
        ),
        tabPanel(
          "Barplot",
          plotOutput("barplot")
        ),
        tabPanel(
          "Text",
          textOutput("text")
        )
      )
    )
  )
))
