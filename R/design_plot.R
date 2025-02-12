#' Design Plot for fMRI Model
#'
#' @description
#' Generates an interactive Shiny app that plots the design matrix for a given
#' fMRI model. The design matrix is first converted into a long-format tibble and
#' then plotted over time, faceted by block. Several customization options allow
#' you to adjust the title, axis labels, line thickness, color palette, and more.
#'
#' @param fmrimod An \code{fmri_model} object.
#' @param term_name Optional: Name of the term to plot. If \code{NULL} (the default),
#'   the first term is used.
#' @param longnames Logical; if TRUE, use long condition names in the legend. Default is FALSE.
#' @param plot_title Optional plot title. If \code{NULL}, a default title is generated.
#' @param x_label Label for the x-axis. Default is "Time".
#' @param y_label Label for the y-axis. Default is "Value".
#' @param line_size Numeric; line thickness for the plot. Default is 1.
#' @param color_palette Character; name of a ColorBrewer palette to use (e.g., "Set1"). Default is "Set1".
#' @param facet_ncol Number of columns for facet_wrap. Default is 1.
#' @param theme_custom A ggplot2 theme to apply. Default is \code{theme_bw(base_size = 14)}.
#' @param legend_threshold Numeric; if the number of unique conditions exceeds this value,
#'   the legend is hidden. Default is 25.
#' @param ... Additional arguments passed to ggplot2::geom_line().
#'
#' @return A Shiny app that displays the design plot.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap labs theme_bw scale_color_brewer guides theme_minimal
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @examples
#' if (interactive()) {
#'   ## --- Construct a sampling frame ---
#'   sframe <- sampling_frame(blocklens = c(100, 100), TR = 2, precision = 0.5)
#'
#'   ## --- Create a dummy event table ---
#'   set.seed(123)
#'   event_table <- data.frame(
#'     onset = seq(10, 190, length.out = 20),
#'     x = rnorm(20),
#'     y = rnorm(20),
#'     run = rep(1:2, each = 10)
#'   )
#'
#'   ## --- Construct a baseline model ---
#'   base_mod <- baseline_model(basis = "bs", degree = 3, sframe = sframe, intercept = "runwise")
#'
#'   ## --- Construct an event model using a formula ---
#'   ev_mod <- event_model(x = onset ~ hrf(x) + hrf(y), data = event_table,
#'                         block = ~ run, sampling_frame = sframe,
#'                         drop_empty = TRUE, durations = rep(0, nrow(event_table)))
#'
#'   ## --- Combine into an fMRI model ---
#'   fmri_mod <- fmri_model(ev_mod, base_mod)
#'
#'   ## --- Launch the interactive design plot ---
#'   design_plot(fmrimod = fmri_mod,
#'               term_name = NULL,
#'               longnames = TRUE,
#'               plot_title = "fMRI Design Matrix",
#'               x_label = "Time (s)",
#'               y_label = "Signal",
#'               line_size = 1.5,
#'               color_palette = "Set2",
#'               facet_ncol = 1)
#' }
design_plot <- function(fmrimod, term_name = NULL, longnames = FALSE,
                        plot_title = NULL,
                        x_label = "Time",
                        y_label = "Value",
                        line_size = 1,
                        color_palette = "Set1",
                        facet_ncol = 1,
                        theme_custom = ggplot2::theme_bw(base_size = 14),
                        legend_threshold = 25,
                        ...) {
  with_package("shiny")
  stopifnot(inherits(fmrimod, "fmri_model"))
  
  # Extract event terms from the fmri_model.
  all_terms <- terms(fmrimod)
  term_names <- sapply(all_terms, "[[", "varname")
  if (is.null(term_names) || length(term_names) == 0) {
    stop("No terms found in the fMRI model.")
  }
  
  sframe <- fmrimod$event_model$sampling_frame
  if (!all(c("time", "blockids") %in% names(sframe))) {
    stop("The sampling_frame must contain 'time' and 'blockids' components.")
  }
  
  # Convert each term's design matrix into long format.
  dflist <- lapply(all_terms, function(term) {
    dm <- suppressMessages(tibble::as_tibble(design_matrix(term), .name_repair = "check_unique"))
    dm$.block <- sframe$blockids
    dm$.time <- sframe$time
    cnames <- if (longnames) conditions(term) else shortnames(term)
    if (!is.null(cnames) && length(cnames) == (ncol(dm) - 2)) {
      names(dm)[1:(ncol(dm) - 2)] <- cnames
    }
    tidyr::pivot_longer(dm, cols = -c(.time, .block), names_to = "condition", values_to = "value")
  })
  names(dflist) <- term_names
  
  # Define UI for Shiny app.
  ui <- shiny::fluidPage(
    shiny::titlePanel("Design Plot for fMRI Model"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("term", "Select Term", choices = term_names, selected = term_names[1]),
        shiny::selectInput("block", "Select Block", choices = c("All", sort(unique(dflist[[1]]$.block)))),
        shiny::sliderInput("time_range", "Time Range:",
                           min = min(sframe$time),
                           max = max(sframe$time),
                           value = c(min(sframe$time), max(sframe$time)))
      ),
      shiny::mainPanel(
        shiny::plotOutput("dplot", height = "600px")
      )
    )
  )
  
  # Define server for Shiny app.
  server <- function(input, output, session) {
    output$dplot <- shiny::renderPlot({
      dfx <- dflist[[input$term]]
      df_filtered <- dfx %>% dplyr::filter(dplyr::between(.time, input$time_range[1], input$time_range[2]))
      if (input$block != "All") {
        df_filtered <- dplyr::filter(df_filtered, .block == input$block)
      }
      p <- ggplot2::ggplot(df_filtered, ggplot2::aes(x = .time, y = value, colour = condition)) +
        ggplot2::geom_line(size = line_size, ...) +
        ggplot2::facet_wrap(~ .block, ncol = facet_ncol) +
        ggplot2::labs(title = if (!is.null(plot_title)) plot_title else paste("Design Plot:", input$term),
                      x = x_label, y = y_label, colour = "Condition") +
        ggplot2::scale_color_brewer(palette = color_palette) +
        theme_custom
      if (length(unique(df_filtered$condition)) > legend_threshold) {
        p <- p + ggplot2::guides(colour = "none")
      }
      p
    })
  }
  
  # Launch the Shiny app.
  shiny::shinyApp(ui = ui, server = server)
}
