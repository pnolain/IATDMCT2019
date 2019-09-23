library(shiny)
library(shinydashboard)
library(rhandsontable)
library(shinyjs)

tagList(
  tags$head(useShinyjs(),
            tags$link(href = "style.css", rel = "stylesheet")),
  dashboardPage(
    dashboardHeader(title = "TDM Alirocumab Demo"),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      div(id = "loading-content",
          img(src = "balls.svg"),
          h5("Loading...")),
      fluidRow(
        box(
          title = "TDM information",
          width = 12,
          
          fluidRow(
            column(
              4,
              h4("Previous Doses"),
              actionLink("add_admin", "Add"),
              rHandsontableOutput("administration_table"),
              actionLink("clear_admin", "Clear")
            ),
            column(
              4,
              h4("Covariates"),
              actionLink("reset", "Reset"),
              rHandsontableOutput("covariates_table")
            ),
            column(
              4,
              h4("Samples"),
              actionLink("add_sample", "Add"),
              rHandsontableOutput("observations_table"),
              actionLink("clear_samples", "Clear")
            )
          ),
          br(),
          fluidRow(
            column(4,
                   h4("")),
            column(
              4,
              (actionButton("go", "Estimate", class = "btn-primary btn-lg")),
              conditionalPanel(
                condition = "$('#dosing_plot').hasClass('recalculating') |$('#parameters_table').hasClass('recalculating')
                       |$('#distributions_plot').hasClass('recalculating') |$('#forecast_plot').hasClass('recalculating')"  ,
                tags$div(id = "plotSpinner",
                         img(src = "ajax-loader-bar.gif"))
              ),
              uiOutput("model_problem"),
              textOutput("obj_Fn_value")
            ),
            column(4,
                   rHandsontableOutput("convergence_table"))
          )
        )
      ),
      fluidRow(box(
        title = "Profiles",
        width = 12,
        conditionalPanel(
          condition = "$('#dosing_plot').hasClass('recalculating') |$('#parameters_table').hasClass('recalculating')
                       |$('#distributions_plot').hasClass('recalculating') |$('#forecast_plot').hasClass('recalculating')"  ,
          tags$div(id = "plotSpinner",
                   img(src = "ajax-loader.gif"))
        ),
        div(id = "mainPlotContainer",
            plotOutput(
              "dosing_plot",
              height = 500,
              click = clickOpts(id = "dosing_plot_click")
            ))
      )),
      
      
      
      
      fluidRow(box(
        title = "Statistics",
        width = 12,
        fluidRow(column(
          width = 12,
          h4("Parameters"),
          rHandsontableOutput("parameters_table")
        )),
        fluidRow(column(
          width = 12,
          h4("Distributions"),
          plotOutput("distributions_plot", height = 500)
        ))
      )),
      fluidRow(box(
        title = "Forecast",
        width = 12,
        plotOutput("forecast_plot")
      ))
    )
  )
  
)