


library(shiny)

ui <- shinyUI(pageWithSidebar(
  headerPanel("Model Terms"),
  sidebarPanel(width=4,
               fluidRow(column(12,
                               h2('term'),
                               uiOutput('uiOutpt')
               )), # END fluidRow
               fluidRow(
                 column(4,div()),
                 column(4,actionButton("add", "Add")),
                 column(4,actionButton("remove", "Remove"))
               ) # END fluidRow
  ), # END sidebarPanel
  mainPanel(
    textOutput("text2"),
    tableOutput('tbl')
  )
))

server <- shinyServer(function(input, output) {
  features <- reactiveValues(renderd=c(1),
                             conv=c(50),
                             inlabels=c('A'),
                             outlabels=c('B'))
  
  df <- eventReactive(input$goButton, {
    out <- lapply(features$renderd,function(i){
      fv <- paste0('numInp_',i)
      vn <- paste0('InLabel',i)
      data.frame(Variable=input[[vn]], Value=input[[fv]] )
    })
    do.call(rbind,out)
  })
  
  output$nText <- renderText({
    ntext()
  })
  output$text2 <- renderText({ 
    paste(sprintf("You have selected feature: %s", paste(features$renderd,collapse=", ")))
  })
  
  output$tbl <- renderTable({
    df()
  })
  
  # Increment reactive values array used to store how may rows we have rendered
  observeEvent(input$add,{
    out <- lapply(features$renderd,function(i){
      fv <- paste0('numInp_',i)
      vn <- paste0('InLabel',i)
      vo <- paste0('OutLabel',i)
      data.frame(inlabels=input[[vn]],outlabels=input[[vo]], conv=input[[fv]] )
    })
    df<-do.call(rbind,out)
    print(df)
    features$inlabels <- c(as.character(df$inlabels),' ')
    features$outlabels <- c(as.character(df$outlabels),' ')
    print(c(features$inlabels,features$outlabels))
    
    features$renderd <- c(features$renderd, length(features$renderd)+1)
    print(features$renderd)
    print(names(features))
    features$conv<-c(df$conv,51-length(features$renderd))
  })
  
  observeEvent(input$remove,{
    features$renderd <- features$renderd[-length(features$renderd)]
  })
  
  # If reactive vector updated we render the UI again
  observe({
    output$uiOutpt <- renderUI({
      # Create rows
      rows <- lapply(features$renderd,function(i){
        fluidRow(
          # duplicate choices make selectize poop the bed, use unique():
          column(4,  selectizeInput(paste0('vars',i), 
                                    label = 'Input Name',selected=features$inlabels[i],
                                    choices=unique(c(features$inlabels[i],features$outlabels[!features$outlabels %in% features$inlabels])),
                                    options = list(create = TRUE))),
          column(4,  sliderInput(paste0('numInp_',i), label="Conversion",min = 0, max = 100, value = features$conv[i])),
          column(4, selectizeInput(paste0('OutLabel',i), 
                                   label = "Output Name", selected=features$outlabels[i],
                                   choices=unique(c(features$inlabels,features$outlabels)),
                                   options = list(create = TRUE)))
        )
      })
      do.call(shiny::tagList,rows)
    })
  })
})

shinyApp(ui=ui,server=server)  