library(shiny)

ui <- fluidPage(
  titlePanel("Theophylline PK data visualization"),
  sidebarLayout(
    sidebarPanel(
      # Create a select drop-down
      selectInput("subject_id", "Select an individual:", choices = 1:12),
      img(src = "iatdmct2019.jpg")),
    mainPanel(
      # Place-holder for a plot to
      plotOutput("pk_profile")
    )
  )
)

server <- function(input, output) {
  # Filter ID data
  individual_data <- reactive({
    # Subset the Theophylline dataset
    subset(Theoph, Subject == input$subject_id)
  })
  
  # Generate a plot named "pk_profile"
  output$pk_profile <- renderPlot({
    cc_data <- individual_data()
    
    op <- par(mar = c(5, 4, 0, 2))
    plot(conc ~ Time, data = cc_data,
         xlab = "Time since drug administration (hr)",
         ylab = "Theophylline concentration (mg/L)",
         sub  = paste0("Theophylline data - ID: ", input$subject_id),
         type = "o")
    par(op)
  })
}

shinyApp(ui, server)
