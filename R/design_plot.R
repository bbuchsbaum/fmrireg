

#' Design Plot for fMRI Model
#' 
#' @description
#' Generate an interactive plot of the design matrix for a given fMRI model.
#'
#' @param fmrimod The `fmri_model` object.
#' @param longnames Use long names in the legend (default: FALSE).
#'
#' @importFrom ggplot2 ggplot aes_string aes
#' 
#' @return An interactive plot using Shiny.
#' @export
design_plot <- function(fmrimod, longnames=FALSE) {
  with_package("shiny")
  stopifnot(inherits(fmrimod, "fmri_model"))
  
  et <- terms(fmrimod$event_model)
  bt <- terms(fmrimod$baseline_model)
  
  all_terms <- terms(fmrimod)
  term_names <- names(all_terms)
  
  sframe <- fmrimod$event_model$sampling_frame
  
  condition = value = .time = .block = NULL
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
  dfx <- dflist[[1]]
  
  
  ui <- function() {
    shiny::fluidPage(
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("term", "Term", term_names),
          shiny::selectInput("block", "Blocks", c("All", sort(unique(dfx$.block)))),
          shiny::sliderInput("range", "Time Range:", min = min(dfx$.time), 
                                              max = max(dfx$.time), 
                                              value = c(0, 100))
        ),
        shiny::mainPanel(
          shiny::plotOutput("dplot")
        )
      )
    )
  }
  
  server <- function(input, output, session) {
      output$dplot <- shiny::renderPlot({
        dfx <- dflist[[input$term]]
        p <- if (input$block == "All") {
          df3 <- dfx %>% dplyr::filter(.time > input$range[1] & .time < input$range[2])
          ggplot(df3, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
            xlab("Time") + theme_bw(14)
          
        } else {
          df3 <- dfx %>% dplyr::filter(.block == input$block & (.time > input$range[1]) & (.time < input$range[2]))
          ggplot(df3, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
            xlab("Time") + theme_bw(24)
         
        }
        
        if (length(unique(dfx$condition)) > 25) {
          p <- p + ggplot2::guides(colour=FALSE)
        }
        
        p
      })
  }
  
  shiny::shinyApp(ui = ui, server = server)
  
}

#' @keywords internal
#' @noRd
lookup_hrf <- function(label, lag) {
  switch(label,
         Gamma=getHRF("gamma", lag=lag),
         Gaussian=getHRF("gaussian", lag=lag),
         SPMG1=getHRF("spmg1", lag=lag),
         SPMG2=getHRF("spmg2", lag=lag),
         SPMG3=getHRF("spmg3", lag=lag))
  
}


#' @noRd
hrf_plot <- function() {
  with_package("shiny")
  ui <- function() {
    shiny::fluidPage(
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("hrf", "Hemodynamic Response Function", c("Gamma", "Gaussian", "SPMG1", "SPMG2", "SPMG3")),
          shiny::sliderInput("offset", "Lag", min=0, max=10, value=0),
          shiny::sliderInput("duration", "Event Duration", min = 0, 
                      max = 10, 
                      value = 0)
        ),
        shiny::mainPanel(
          shiny::plotOutput("dplot")
        )
      )
    )
  }
  
  server <- function(input, output, session) {
    output$dplot <- shiny::renderPlot({
      hrf <- lookup_hrf(input$hrf, as.numeric(input$offset))
      duration <- input$duration
      reg <- regressor(0, hrf, duration)
      sgrid <- seq(0,24,by=.1)
      
      if (nbasis(hrf) == 1) {
        df1 <- data.frame(Y=evaluate(reg, sgrid), time=sgrid)
        ggplot2::ggplot(df1, ggplot2::aes(time, Y)) + geom_line() + xlab("Time") + theme_bw(14)
      } else {
        Y.1 <- NULL
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
  shiny::shinyApp(ui = ui, server = server)
  
}
