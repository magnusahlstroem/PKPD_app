library(shiny)

source("R/scheduler.R", local = TRUE)

ui <- fluidPage(
  titlePanel("Automatisk skemalægning"),
  sidebarLayout(
    sidebarPanel(
      helpText("Indtast personer og opgaver. Appen laver et simpelt skema uden dobbeltbooking i samme tidsrum."),
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
        value = "Matematik: man-08\nBiologi: tir-08\nFysik: ons-08",
        rows = 6,
        placeholder = "Format: Opgave: tidsrum"
      ),
      numericInput("max_assignments", "Maks. tildelinger pr. person", value = 999, min = 1, step = 1),
      actionButton("create", "Lav skema", class = "btn-primary")
    ),
    mainPanel(
      h3("Forslag til skema"),
      tableOutput("schedule"),
      h3("Ikke planlagt"),
      tableOutput("unassigned")
    )
  )
)

server <- function(input, output, session) {
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
