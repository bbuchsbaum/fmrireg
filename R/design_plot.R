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
#'   sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2, precision = 0.5)
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
#' }-------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#' @importFrom bslib bs_theme card card_header layout_sidebar
design_plot <- function(fmrimod, term_name = NULL, longnames = FALSE,
                         plot_title = NULL,
                         x_label = "Time (s)", y_label = "Amplitude",
                         line_size = 1,
                         color_palette = "viridis",      # <- colour-blind safe
                         facet_ncol   = 2,               # <- sensible default
                         theme_custom = ggplot2::theme_minimal(base_size = 15) +
                                        ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines")),
                         legend_threshold = 30, ...){

  with_package(c("shiny", "plotly", "bslib", "thematic"))
  stopifnot(inherits(fmrimod, "fmri_model"))

  # -- prep ------------------------------------------------------------------
  terms_all  <- terms(fmrimod)
  term_names <- vapply(terms_all, `[[`, character(1), "varname")

  if (is.null(term_name)) term_name <- term_names[1]
  if (!(term_name %in% term_names))
      stop("term_name not found in model: ", term_name)

  sframe <- fmrimod$event_model$sampling_frame
  df_time <- sframe$time
  df_block<- sframe$blockids

  longify <- function(term){
    dm   <- tibble::as_tibble(design_matrix(term), .name_repair = "unique")
    dm$.block <- df_block
    dm$.time  <- df_time

    # pretty column names
    cn <- if (longnames) conditions(term) else shortnames(term)
    if (!is.null(cn) && length(cn)==ncol(dm)-2) names(dm)[1:(ncol(dm)-2)] <- cn

    tidyr::pivot_longer(dm, -c(.time,.block),
                        names_to = "condition", values_to = "value")
  }
  df_long <- lapply(terms_all, longify)
  names(df_long) <- term_names

  # -- shiny UI --------------------------------------------------------------
  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bg = "#fafafa", fg = "#222", primary = "#4c72b0"),
    bslib::card(
      bslib::card_header("fmrireg design viewer"),
      bslib::layout_sidebar(
        sidebar = list(
          shiny::selectInput("term",   "Term",      term_names, term_name),
          shiny::selectInput("block",  "Block",     c("all", sort(unique(df_block)))),
          shiny::sliderInput("timer",  "Time-window",
                             min(df_time), max(df_time),
                             value = range(df_time), step = diff(range(df_time))/200),
          shiny::checkboxInput("zero", "Y‑axis starts at zero", TRUE),
          shiny::hr(),
          shiny::helpText("Drag to zoom, double‑click to reset.")
        ),
        shiny::mainPanel(
          plotly::plotlyOutput("plot", height = "650px", inline = TRUE)
        )
      )
    )
  )

  # ── server ----------------------------------------------------------------
  server <- function(input, output, session){

    reactive_df <- shiny::reactive({
      d <- df_long[[input$term]]
      d <- d[d$.time >= input$timer[1] & d$.time <= input$timer[2], ]
      if (input$block != "all") d <- d[d$.block == input$block, ]
      d
    })

    output$plot <- plotly::renderPlotly({
      d <- reactive_df()
      gg <- ggplot2::ggplot(
              d, ggplot2::aes(.time, value,
                              colour = condition,
                              text = paste0("t = ", round(.time,2),
                                            "<br>cond = ", condition,
                                            "<br>val = ", signif(value,3)))
            ) +
            ggplot2::geom_line(size = line_size, ...) +
            ggplot2::facet_wrap(~ .block, ncol = facet_ncol) +
            ggplot2::labs(title   = plot_title %||% paste("Design:", input$term),
                          x = x_label, y = y_label, colour = "Condition") +
            theme_custom +
            { if (tolower(color_palette) == "viridis")
                  ggplot2::scale_colour_viridis_d(option = "C")
              else ggplot2::scale_colour_brewer(palette = color_palette) } +
            { if (input$zero) ggplot2::expand_limits(y = 0) } +
            { if (length(unique(d$condition)) > legend_threshold)
                  ggplot2:: guides(colour = "none") }

      plotly::ggplotly(gg, tooltip = "text") |>
        plotly::layout(hovermode = "closest") |>
        plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d"))
    })
  }

  shiny::shinyApp(ui, server)
}
