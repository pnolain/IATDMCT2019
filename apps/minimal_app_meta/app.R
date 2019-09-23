library(shiny)
library(shinymeta)

ui <- fluidPage(
  titlePanel("Theophylline PK data visualization - Reproducible"),
  sidebarLayout(
    sidebarPanel(
      selectInput("subject_id", "Select an individual:", choices = 1:12),
      img(src = "iatdmct2019.jpg")),
    mainPanel(
      outputCodeButton(plotOutput("pk_profile"))
    )
  )
)

server <- function(input, output) {
  individual_data <- metaReactive({
    subset(Theoph, Subject == ..(input$subject_id))
  })
  
  output$pk_profile <- metaRender(renderPlot, {
    cc_data <- ..(individual_data())
    
    plot(conc ~ Time, data = cc_data,
         xlab = "Time since drug administration (hr)",
         ylab = "Theophylline concentration (mg/L)",
         sub  = paste0("Theophylline data - ID: ", ..(input$subject_id)),
         type = "o")
  })
  
  observeEvent(input$pk_profile_output_code, {
    code <- expandChain(output$pk_profile())
    displayCodeModal(code, wordWrap = TRUE)
  })
}

shinyApp(ui, server)
