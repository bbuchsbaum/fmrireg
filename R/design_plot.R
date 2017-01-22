
#' @importFrom ggplot2 ggplot aes_string aes
#' @import shiny
design_plot <- function(fmrimod) {
  et <- terms(fmrimod$event_model)
  bt <- terms(fmrimod$baseline_model)
  
  all_terms <- terms(fmrimod)
  term_names <- names(all_terms)
  
  sframe <- fmrimod$event_model$sampling_frame
  
  dflist <- lapply(all_terms, function(term) {
    dm1 <- tibble::as_tibble(design_matrix(term))
    dm1$.block <- sframe$blockids
    dm1$.time <- sframe$time
    cnames <- conditions(term)
    gather(dm1, condition, value, -.time, -.block)
  })
  
  names(dflist) <- term_names
  
  
  
  ui <- function() {
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("term", "Term", term_names),
          selectInput("block", "Blocks", c("All", sort(unique(df2$.block)))),
          sliderInput("range", "Time Range:", min = min(df2$.time), 
                                              max = max(df2$.time), 
                                              value = c(0, 100))
        ),
        mainPanel(
          plotOutput("dplot")
        )
      )
    )
  }
  
  server <- function(input, output, session) {
      output$dplot <- renderPlot({
        dfx <- dflist[[input$term]]
        if (input$block == "All") {
          df3 <- dfx %>% dplyr::filter(.time > input$range[1] & .time < input$range[2])
          ggplot(df3, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
            xlab("Time") + theme_bw(14)
          
        } else {
          df3 <- dfx %>% dplyr::filter(.block == input$block & (.time > input$range[1]) & (.time < input$range[2]))
          ggplot(df3, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
            xlab("Time") + theme_bw(24)
         
        }
      })
  }
  
  shinyApp(ui = ui, server = server)
  
}


#design_plot(fmrimod)

