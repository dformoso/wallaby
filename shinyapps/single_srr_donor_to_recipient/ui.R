library(shiny)
library(shinyFiles)

shinyUI(fluidPage(
  titlePanel("Output Analysis - Single SRR - Donor to Recipient Workflow"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Directory"),
      textOutput("dir"),
      br(),
      shinyDirButton(
        id = "dir",
        label = "Change Directory",
        title = "Select directory",
        icon = icon('folder')
      ),
      br(),
      
    ),
    
    mainPanel(
      h4("SRRs"),
      htmlOutput("srrs"),
      br(),
      h4("Reference Genomes"),
      htmlOutput("ref_genomes"),
      br(),
      h4("Files"),
      verbatimTextOutput("files")
    )
    
  )
)) 