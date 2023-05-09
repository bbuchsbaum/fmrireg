

#' Design Plot for fMRI Model
#' 
#' @description
#' Generate an interactive plot of the design matrix for a given fMRI model.
#'
#' @param fmrimod The `fmri_model` object.
#' @param longnames Use long names in the legend (default: FALSE).
#'
#' @importFrom ggplot2 ggplot aes_string aes
#' @import shiny
#' 
#' @return An interactive plot using Shiny.
design_plot <- function(fmrimod, longnames=FALSE) {
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
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("term", "Term", term_names),
          selectInput("block", "Blocks", c("All", sort(unique(dfx$.block)))),
          sliderInput("range", "Time Range:", min = min(dfx$.time), 
                                              max = max(dfx$.time), 
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
  
  shinyApp(ui = ui, server = server)
  
}

#' @keywords internal
lookup_hrf <- function(label, lag) {
  switch(label,
         Gamma=getHRF("gamma", lag=lag),
         Gaussian=getHRF("gaussian", lag=lag),
         SPMG1=getHRF("spmg1", lag=lag),
         SPMG2=getHRF("spmg2", lag=lag),
         SPMG3=getHRF("spmg3", lag=lag))
  
}


#' @importFrom shiny fluidPage sidebarLayout sliderInput mainPanel plotOutput
hrf_plot <- function() {
  ui <- function() {
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("hrf", "Hemodynamic Response Function", c("Gamma", "Gaussian", "SPMG1", "SPMG2", "SPMG3")),
          sliderInput("offset", "Lag", min=0, max=10, value=0),
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
  shinyApp(ui = ui, server = server)
  
}
#design_plot(fmrimod)
#hrf_plot()


# @keywords internal
# term_form <- function() {
#   div(class="ui sidebar inverted vertical visible menu",
#       div(class="item",
#           div(class="active header", "Regression Term"),
#           div(class="menu",
#               shiny::uidropdown(),
#               a(class="item", href="#divider", "HRF"),
#               a(class="item", href="#input", "nbasis"),
#           )
#       )
#   )
# }


# design_editor <- function(design, sframe) {
#   requireNamespace("shiny.semantic")
#   term_component <- function() {
#     fluidRow(
#       column(4, selectInput("onset_var", "Onsets", names(design))),
#       column(4, selectInput("factors", "Variable(s)", names(design), multiple=TRUE)),
#       column(4, selectInput("hrf", "HRF", choices=c("Gamma", "Gaussian",
#                                                 "SPMG1", "SPMG2", "SPMG3"),
#                                                 selected="SPMG1"))
#     )
#   }
# 
#   ui <- function() {
#     shinyUI(
#       semanticPage(
#         suppressDependencies("bootstrap"),
#         sidebar()
#       )
#     )
#   }
# 
#   server <- function(input, output, session) {
# 
#     output$dplot <- renderPlot({
#       if (TRUE) {
#         plot(1:100)
#       } else {
#         plot(1:100)
#         # srun <- as.integer(input$show_run)
#         # keep <- which(design[[input$run_variable]] == srun)
#         # des <- design[keep,]
#         # print(nrow(des))
#         # blocklens <- sframe$blocklens
#         # bl <- blocklens[srun]
#         # sframe2 <- sampling_frame(bl, TR=sframe$TR)
#         # print(input$formula)
#         # ev <- event_model(as.formula(input$formula), block=as.formula(paste("~ ", input$run_variable)), data=des, sampling_frame=sframe2)
#         # dmat <- design_matrix(ev)
#         # matplot(dmat, type='l')
#       }
#     })
#   }
# 
# 
#   shinyApp(ui = ui, server = server)
# }
