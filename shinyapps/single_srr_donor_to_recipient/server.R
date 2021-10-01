library(shiny)
library(shinyFiles)
library(stringr)

shinyServer(function(input, output, session) {
  # Directory Object Definition
  shinyDirChoose(input,
                 'dir',
                 roots = c(home = '~/wallaby/workflows'))
  
  # Path for Default Folder
  global <-
    reactiveValues(datapath = '~/wallaby/workflows/cromwell-final-outputs')
  
  # Path for Selected Folder
  dir <- reactive(input$dir)
  
  output$dir <- renderText(global$datapath)
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir()))
                   return()
                 home <- normalizePath('~/wallaby/workflows')
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })
  
  # path
  path <- reactive({
    home <- normalizePath("~/wallaby/workflows/cromwell-final-outputs")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })
  
  # srrs
  output$srrs <- renderPrint({
    HTML(paste(unique(na.omit(
      str_extract(list.files(global$datapath), "SRR+[0-9]+")
    )), sep = "<br/>"))
  })
  
  # ref genomes
  output$ref_genomes <- renderPrint({
    HTML(paste(unique(na.omit(
      str_extract(list.files(global$datapath), "(?<=-to-)[^_]+")
    )), sep = "<br/>"))
  })
  
  # files
  output$files <- renderPrint(list.files(global$datapath))
}) 