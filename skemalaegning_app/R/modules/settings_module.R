settings_module_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Grunddata for normperioden"),
    dateInput(ns("period_start"), "Normperiode starter", value = Sys.Date()),
    numericInput(ns("period_length_weeks"), "Normperiodens længde i uger", value = 12, min = 1, step = 1),
    textAreaInput(ns("notes"), "Noter til skemalægger", rows = 4),
    actionButton(ns("save"), "Gem grunddata", class = "btn-primary"),
    br(), br(),
    tableOutput(ns("current_settings"))
  )
}

settings_module_server <- function(id, conn) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$save, {
      save_schedule_settings(conn, input$period_start, input$period_length_weeks, input$notes)
      showNotification("Grunddata er gemt", type = "message")
    })

    output$current_settings <- renderTable({
      settings <- load_schedule_settings(conn)
      if (nrow(settings) == 0) {
        return(data.frame(Besked = "Der er endnu ikke gemt grunddata."))
      }
      names(settings) <- c("Startdato", "Længde i uger", "Noter")
      settings
    })
  })
}
