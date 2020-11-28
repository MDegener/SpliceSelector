library(shiny)
library(tidyverse)
library(corrplot)
library(psych)
library(RColorBrewer)

# Define UI ----
ui <- fluidPage(
  titlePanel("Explorartory Analysis of PSI Data"),

  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),

  sidebarLayout(
    sidebarPanel(
      h3("Data Input"),
      fileInput(inputId = "psiData",
                label = "Upload a .csv file with your PSI data"),

      hr(),
      h3("Analysis Options"),

      selectInput(inputId = "analysis",
                  label = "Please select your desired analysis",
                  choices = c("Datatable", "Boxplot", "Pairwise Correlations"),
                  selected = "Datatable"),

      uiOutput(outputId = "exonSelection"),

      width = 2
    ),

    mainPanel(

      conditionalPanel(condition = "input.analysis == 'Datatable'",
                       DT::dataTableOutput(outputId = "psiTable")),

      conditionalPanel(condition = "input.analysis == 'Boxplot'",
                       plotOutput(outputId = "boxplot")),

      conditionalPanel(condition = "input.analysis == 'Pairwise Correlations'",
                       plotOutput(outputId = "corrplot"),
                       selectInput(inputId = "corrMethod",
                                   label = "Correlation method",
                                   choices = c("pearson", "spearman", "kendall")),
                       selectInput(inputId = "pAdjustMethod",
                                   label = "P adjustment method",
                                   choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),

      width = 10
    )
  )
)

# Define server logic ----
server <- function(input, output) {

  loadData <- reactive({
    file <- input$psiData
    ext <- tools::file_ext(file$datapath)

    req(file)
    validate(need(ext == "csv", "Please upload .csv of your PSI data"))

    read_csv(file$datapath)

    #read_csv(system.file("extdata", "miso_psi_data.csv", package = "SpliceSelector"))
  })

  output$psiTable <- DT::renderDataTable({
    loadData()
  },
  options = list(
    pageLength = 10
  ))

  output$exonSelection <- renderUI({
    if(is.null(input$psiData)) return(NULL)
    if(input$analysis == "Datatable") return(NULL)

    psiTable <- loadData()
    uniqueExons <- n_distinct(psiTable$exonID)

    tagList(
      hr(),
      h3("Exon Selection"),

      numericInput(inputId = "randomExons",
                   label = "Select a number of random exons",
                   value = 2, min = 2, max = min(50, uniqueExons)),

      fileInput(inputId = "selectedExons",
                label = "Upload .txt with exonIDs of interest (max. 50)")
    )
  })

  output$boxplot <- renderPlot({
    if(is.null(input$randomExons)) return(NULL)

    psiTable <- loadData()
    exons <- sample(unique(psiTable$exonID), input$randomExons)
    psiTable <- filter(psiTable, exonID %in% exons)

    ggplot(psiTable,
           aes(x = eventName, y = exonInclusion, fill = eventName)) +
      geom_boxplot() +
      theme_bw() +
      theme(text = element_text(size = 20),
            strip.background =element_rect(color="white", fill="white"),
            legend.position = "none",
            axis.text.x = element_text(angle = 60, hjust = 1)) +
      xlab("Exon-Skipping Event") + ylab("Percent Spliced-In (PSI)")
  })

  output$corrplot <- renderPlot({
    if(is.null(input$randomExons)) return(NULL)

    psiTable <- loadData()
    exons <- sample(unique(psiTable$exonID), input$randomExons)
    psiTable <- filter(psiTable, exonID %in% exons)

    psiData <- psiTable %>%
      select(eventName, sampleID, exonInclusion) %>%
      spread(key = "eventName", value = "exonInclusion") %>%
      select(-sampleID) %>% as.matrix()

    makeCorrplot(psiData, method=input$corrMethod, adjust=input$pAdjustMethod)

  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
