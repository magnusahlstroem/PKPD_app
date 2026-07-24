doctors_module_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Læger på afdelingen"),
    fluidRow(
      column(4, textInput(ns("name"), "Navn")),
      column(
        4,
        selectInput(
          ns("employment_type"),
          "Ansættelse",
          choices = c("Speciallæge", "Uddannelseslæge", "Introduktionslæge", "Reservelæge", "Vikar", "Andet")
        )
      ),
      column(4, numericInput(ns("weekly_hours"), "Timetal pr. uge", value = 37, min = 0, step = 0.5))
    ),
    uiOutput(ns("function_picker")),
    textAreaInput(ns("notes"), "Noter", rows = 3, placeholder = "Fx særlige hensyn, forskning, undervisning eller faste fridage"),
    actionButton(ns("save"), "Gem læge", class = "btn-primary"),
    br(), br(),
    tableOutput(ns("doctors"))
  )
}

doctors_module_server <- function(id, conn) {
  moduleServer(id, function(input, output, session) {
    output$function_picker <- renderUI({
      functions <- load_clinical_functions(conn)
      choices <- stats::setNames(functions$id, functions$name)
      checkboxGroupInput(session$ns("function_ids"), "Funktioner lægen kan tage", choices = choices)
    })

    observeEvent(input$save, {
      req(nzchar(input$name))
      add_doctor(
        conn,
        input$name,
        input$employment_type,
        input$weekly_hours,
        input$function_ids %||% integer(),
        input$notes
      )
      showNotification("Lægen er gemt", type = "message")
    })

    output$doctors <- renderTable({
      doctors <- load_doctors(conn)
      if (nrow(doctors) == 0) {
        return(data.frame(Besked = "Der er endnu ikke oprettet læger."))
      }
      names(doctors) <- c("ID", "Navn", "Ansættelse", "Timer/uge", "Aktiv", "Noter", "Funktioner")
      doctors
    })
  })
}
