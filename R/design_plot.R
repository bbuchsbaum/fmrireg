

#' design_plot
#' 
#' @importFrom ggplot2 ggplot aes_string aes
#' @import shiny
design_plot <- function(fmrimod, longnames=FALSE) {
  et <- terms(fmrimod$event_model)
  bt <- terms(fmrimod$baseline_model)
  
  all_terms <- terms(fmrimod)
  term_names <- names(all_terms)
  
  sframe <- fmrimod$event_model$sampling_frame
  
  dflist <- lapply(all_terms, function(term) {
    dm1 <- tibble::as_tibble(design_matrix(term))
    dm1$.block <- sframe$blockids
    dm1$.time <- sframe$time
    
    if (longnames) {
      cnames <- conditions(term)
    } else {
      cnames <- shortnames(term)
    }
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

lookup_hrf <- function(label) {
  switch(label,
         Gamma=HRF_GAMMA,
         Gaussian=HRF_GAUSSIAN,
         SPMG1=HRF_SPMG1,
         SPMG2=HRF_SPMG2,
         SPMG3=HRF_SPMG3)
  
}

hrf_plot <- function() {
  ui <- function() {
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("hrf", "Hemodynamic Response Function", c("Gamma", "Gaussian", "SPMG1", "SPMG2", "SPMG3")),
          sliderInput("duration", "Event Duration", min = 0, 
                      max = 10, 
                      value = 0)
        ),
        mainPanel(
          plotOutput("dplot")
        )
      )
    )
  }
  
  server <- function(input, output, session) {
    output$dplot <- renderPlot({
      hrf <- lookup_hrf(input$hrf)
      duration <- input$duration
      reg <- regressor(0, hrf, duration)
      sgrid <- seq(0,24,by=.1)
      
      if (nbasis(hrf) == 1) {
        df1 <- data.frame(Y=evaluate(reg, sgrid), time=sgrid)
        ggplot2::ggplot(df1, ggplot2::aes(time, Y)) + geom_line() + xlab("Time") + theme_bw(14)
      } else {
        Y <- evaluate(reg, sgrid)
        df1 <- data.frame(Y =Y, time=sgrid)
        p <- ggplot2::ggplot(df1, ggplot2::aes(time, Y.1)) + geom_line() + xlab("Time") + theme_bw(14)
        
        for (i in 2:nbasis(hrf)) {
          df1$Yi <- Y[,i]
          cat("i = ", i, "\n")
          cat("y: ", df1$Yi, "\n")
          p <- p + geom_line(aes_string("time", paste0("Y.", i)), colour=i)
        }
        
        print(p)
        
      }
    })
  }
  shinyApp(ui = ui, server = server)
  
}
#design_plot(fmrimod)
hrf_plot()
