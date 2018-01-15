

design_editor <- function(design, formula="", sframe) {
  ui <- function() {
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("run_variable", "Run Variable", names(design)),
          selectInput("show_run", "Show Run", 1),
          textInput("formula", "Event Formula", value=formula),
          textInput("tr", "Repetition Time", value=2),
          
          textInput("blocklens", "Run Lengths", value=426+25)
        ),
        mainPanel(
          plotOutput("dplot")
        )
      )
    )
  }
  
  server <- function(input, output, session) {
    
    output$dplot <- renderPlot({
      if (input$formula == "") {
        plot()
      } else {
        srun <- as.integer(input$show_run)
        keep <- which(design[[input$run_variable]] == srun)
        des <- design[keep,]
        print(nrow(des))
        blocklens <- as.numeric(rep(input$blocklens, length.out=length(unique(design[[input$run_variable]]))))
        bl <- blocklens[srun]
       
        sframe <- sampling_frame(bl, as.numeric(input$tr))
        ev <- event_model(as.formula(input$formula), block=as.formula(paste("~ ", input$run_variable)), data=des, sampling_frame=sframe)
        dmat <- design_matrix(ev)
        matplot(dmat, type='l')
      }
    })
  }
  
  
  shinyApp(ui = ui, server = server)
}