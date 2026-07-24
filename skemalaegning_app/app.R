library(shiny)
library(DBI)
library(RSQLite)

source("R/db.R", local = TRUE)
source("R/scheduler.R", local = TRUE)
source("R/modules/settings_module.R", local = TRUE)
source("R/modules/functions_module.R", local = TRUE)
source("R/modules/doctors_module.R", local = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

conn <- connect_db()
init_db(conn)
onStop(function() DBI::dbDisconnect(conn))

ui <- fluidPage(
  titlePanel("Dermatologisk skemalægning"),
  p("Grunddata til skemalægning for læger på en dermatologisk afdeling."),
  tabsetPanel(
    tabPanel("Normperiode", settings_module_ui("settings")),
    tabPanel("Funktioner", functions_module_ui("functions")),
    tabPanel("Læger", doctors_module_ui("doctors")),
    tabPanel(
      "Første skemaforslag",
      h2("Simpel test af skemamotor"),
      p("Dette faneblad er midlertidigt og kan senere kobles direkte til de gemte databaseoplysninger."),
      textAreaInput(
        "people",
        "Personer",
        value = "Anna: man-08, tir-08, ons-08\nBo: man-08, ons-08\nClara: tir-08, ons-08",
        rows = 6,
        placeholder = "Format: Navn: tidsrum, tidsrum"
      ),
      textAreaInput(
        "shifts",
        "Opgaver/vagter/lektioner",
        value = "Ambulatorium: man-08\nStuegang: tir-08\nAkut: ons-08",
        rows = 6,
        placeholder = "Format: Opgave: tidsrum"
      ),
      numericInput("max_assignments", "Maks. tildelinger pr. person", value = 999, min = 1, step = 1),
      actionButton("create", "Lav skema", class = "btn-primary"),
      h3("Forslag til skema"),
      tableOutput("schedule"),
      h3("Ikke planlagt"),
      tableOutput("unassigned")
    )
  )
)

server <- function(input, output, session) {
  settings_module_server("settings", conn)
  functions_module_server("functions", conn)
  doctors_module_server("doctors", conn)

  result <- eventReactive(input$create, {
    people <- parse_people(input$people)
    shifts <- parse_shifts(input$shifts)
    schedule <- create_schedule(people, shifts, max_assignments = input$max_assignments)

    list(
      schedule = schedule,
      unassigned = find_unassigned_shifts(shifts, schedule)
    )
  }, ignoreInit = FALSE)

  output$schedule <- renderTable({
    schedule <- result()$schedule
    if (nrow(schedule) == 0) {
      return(data.frame(Besked = "Der kunne ikke laves nogen tildelinger med de angivne data."))
    }

    names(schedule) <- c("Opgave", "Tidsrum", "Person")
    schedule
  })

  output$unassigned <- renderTable({
    unassigned <- result()$unassigned
    if (nrow(unassigned) == 0) {
      return(data.frame(Besked = "Alle opgaver blev planlagt."))
    }

    names(unassigned) <- c("Opgave", "Tidsrum")
    unassigned
  })
}

shinyApp(ui = ui, server = server)
