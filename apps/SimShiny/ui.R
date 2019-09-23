library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)

appCSS <- "
#loading-content {
  position: absolute;
  background: #FFFFFF;
  opacity: 0.9;
  z-index: 1000000;
  left: 0;
  right: 0;
  top:0;
  height: 100%;
  text-align: center;
  color: #00000;
  padding-top: 200px;
}
.sidebar{margin-left:5px;margin-right:5px;}
.sidebar h4{ padding:10px; border-bottom: solid #357ca5;font-weight:bold;}
.sidebar h4:hover{cursor:default; }
#design-container { height: 505px; overflow: auto; overflow-x: hidden; }
#subject-container { height: 505px; overflow: auto; overflow-x: hidden; padding: 5px 0px 0px 0px;}
#subject-container .shiny-input-container { min-width:200px; }
#results_box + div { min-height : 600px}
.shiny-output-error-validation {  color: red; font-weight: bold}
"


dashboardPage(title = "SimShiny - Alirocumab",
  dashboardHeader(title = "SimShiny - Alirocumab", titleWidth = 400),
  dashboardSidebar(
    div(align = "center",
        img(src = "iatdmct2019.jpg", width = "200px")
        ),
    h4(icon("cogs"), span("Model")),
    radioButtons(
      "selected_model", "Selection", selected = "PKPD",
      choices = c(
        "PK/PD" = "PKPD",
        "TMDD" = "TMDD"
      )
    ),
    em(textOutput("model_problem")),
    tags$label("References:"),
    tags$ul(
      tags$li(tags$strong(
        tags$a(target = "_blank", 
               href = "https://dx.doi.org/10.1007%2Fs40262-018-0670-5", 
               "doi: 10.1007/s40262-018-0670-5"))),
      tags$li(tags$strong(
        tags$a(target = "_blank", 
               href = "https://dx.doi.org/10.1007%2Fs40262-016-0505-1", 
               "doi: 10.1007/s40262-016-0505-1")))),
    h4(icon("clock-o"), span("Time settings")),
    uiOutput("design_time_unit"),
    uiOutput("design_time_range"),
    uiOutput("precision"),
    h4(icon("history"), span("Previous simulations")),
    uiOutput("previous_simulations_selection"),
    div(align = "center",
        actionButton("clear_history", "Clear history"))
  ),
  dashboardBody(
    id = "main_tab",
    useShinyjs(),
    shinytoastr::useToastr(),
    inlineCSS(appCSS),
    # Loading message
    div(
      id = "loading-content",
      img(src = "balls.svg"),
      h5("Loading...")
    ),
      div(
      id = "app-content",
      
    fluidRow(column(
      width = 6,
      box(
        id = "design-container",
        title = list(icon("dashboard"), span("Design")),
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        textInput("run_name", label = "Simulation name"),
        fluidRow(
          column(
            width = 5,
            fluidRow(column(6, selectInput("dose", "Dose (mg)", choices = c(75, 150, 300))),
                     column(6, uiOutput("route"))),
            numericInput("dosing_interval",
                                                             sprintf("Dosing interval (%s)", 14),
                                                             value = 14),
            sliderInput(
              "n_admins",
              "Number of administrations",
              value = 6,
              min = 1 ,
              max = 10,
              step = 1
            ),
            div(align = "center", actionButton("add_administrations", "Add"))
          ),
          column(
            width = 7,
            strong("Administrations"),
            dataTableOutput("administrations_table", height = 310),
            div(
              align = "center",
              style = "padding:5px;",
              actionButton("remove_administrations", "Remove selected"),
              actionButton("clear_administrations", "Clear")
            )
          )
        ),
        div(align = "center",
            actionButton("simulate_button", "Simulate"))
      )
    ), column(
      width = 6,
      box(
        id = "subject-container",
        title = list(icon("female"), icon("male"), span("Subject")),
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        tabBox(
          width = 12,
          tabPanel(
            value = "covariates_tab",
            title = list(icon("sliders"), span("Covariates")),
            uiOutput("covariates_ui_output")
          ),
          tabPanel(
            value = "parameters_tab",
            title = list(icon("wrench"), span("Parameters")),
            uiOutput("parameters_ui_output"),
            div(align = "center", actionButton("reset_parameters", "Reset"))
          )
        ),
        column(6, checkboxInput("simulate_iiv", "Simulate with variability", value = FALSE)),
                 column(6, sliderInput("iiv_n", label = "N individuals", value = 500, min = 100, max = 1000, step = 100))
      )
    )),
    fluidRow(column(
      width = 12,
      div(
        id = "results-container",
        tabBox(
          id = "results_box",
          title = list(icon("table"), span("Results")),
          width = 12,
          side = "right",
          tabPanel(
            value = "plots_tab",
            title = list(icon("line-chart"), span("Plots")),
            fluidRow(
              column(width = 9, div(
                plotOutput(
                  "simulation_plot",
                  height = 700,
                  dblclick = "simulation_plot_dblclick",
                  brush = brushOpts(id = "simulation_plot_brush",
                                    resetOnNew = TRUE)
                )
              )),
              column(
                width = 3,
                h4("Displayed output"),
                uiOutput("selected_output"),
                h5("Reference lines"),
                fluidRow(column(
                  width = 6,
                  textInput("reference_line_label", "Label")
                ),
                column(
                  width = 6,
                  numericInput("reference_line_value", "Value", value = 0)
                )),
                div(
                  align = "center",
                  actionButton("add_reference_line", "Add"),
                  actionButton("clear_reference_lines", "Clear")
                ),
                hr(),
                h4("Plot settings"),
                fluidRow(column(
                  width = 6,
                  selectInput("x_scale", "X-Scale", choices = c("Linear", "Logarithmic"))
                ),
                column(
                  width = 6,
                  selectInput("y_scale", "Y-Scale", choices = c("Linear", "Logarithmic"))
                )),
                h6("*Brush and double-click to zoom"),
                checkboxInput("show_iiv", "Show IIV 5%-95% percentile range", value = TRUE)
              )
            )
          ),
          tabPanel(
            title = list(icon("database"), span("Data")),
            value = "data_tab",
            fluidRow(column(
              width = 12,
              h4("Exposition time range"),
              flowLayout(
                numericInput("exposition_start", "From", value = NA),
                numericInput("exposition_end", "To", value = NA)
              )
            )),
            fluidRow(tabBox(
              width = 12,
              tabPanel(
                title = list(icon("area-chart"), span("Exposure")),
                value = "exposure_tab",
                fluidRow(column(
                  width = 12,
                  fluidRow(column(width =
                                    6, h4(
                                      "Individual exposure indices"
                                    )),
                           column(
                             width = 6, downloadButton("download_exposure_parameters")
                           )),
                  uiOutput("individual_exposure_table")
                ))
              ),
              tabPanel(
                title = list(icon("file"), span("Outputs")),
                value = "outputsTab",
                fluidRow(column(width =
                                  6, h4("Outputs")),
                         column(
                           width = 6, downloadButton("download_outputs")
                         )),
                dataTableOutput("outputs_table")
              )
            ))
          )
        )
      )
    )),
    div(align = "right",
        h6(textOutput(
          "computation_message"
        )))
  ))
)
