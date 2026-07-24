functions_module_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Funktioner der skal dækkes"),
    fluidRow(
      column(4, textInput(ns("name"), "Funktionsnavn", placeholder = "Fx stuegang, ambulatorium, akut")),
      column(4, selectInput(ns("coverage_type"), "Dækkes", choices = c("Dagligt", "Månedligt", "Efter behov"))),
      column(2, numericInput(ns("required_per_day"), "Pr. dag", value = 1, min = 0, step = 1)),
      column(2, numericInput(ns("required_per_month"), "Pr. måned", value = 0, min = 0, step = 1))
    ),
    checkboxInput(ns("requires_specialist"), "Kræver speciallæge", value = FALSE),
    actionButton(ns("save"), "Gem funktion", class = "btn-primary"),
    br(), br(),
    tableOutput(ns("functions"))
  )
}

functions_module_server <- function(id, conn) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$save, {
      req(nzchar(input$name))
      add_clinical_function(
        conn,
        input$name,
        input$coverage_type,
        input$required_per_day,
        input$required_per_month,
        input$requires_specialist
      )
      showNotification("Funktionen er gemt", type = "message")
    })

    output$functions <- renderTable({
      functions <- load_clinical_functions(conn)
      if (nrow(functions) == 0) {
        return(data.frame(Besked = "Der er endnu ikke oprettet funktioner."))
      }
      functions$requires_specialist <- ifelse(functions$requires_specialist == 1, "Ja", "Nej")
      names(functions) <- c("ID", "Navn", "Dækkes", "Pr. dag", "Pr. måned", "Kræver speciallæge", "Aktiv")
      functions
    })
  })
}
