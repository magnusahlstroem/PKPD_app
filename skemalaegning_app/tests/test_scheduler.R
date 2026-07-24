source("skemalaegning_app/R/scheduler.R")

people <- data.frame(
  name = c("Anna", "Bo"),
  available_slots = c("man-08, tir-08", "tir-08"),
  stringsAsFactors = FALSE
)
shifts <- data.frame(
  task = c("Matematik", "Biologi"),
  slot = c("man-08", "tir-08"),
  stringsAsFactors = FALSE
)

schedule <- create_schedule(people, shifts)
stopifnot(identical(schedule$person, c("Anna", "Bo")))

people <- data.frame(name = "Anna", available_slots = "man-08", stringsAsFactors = FALSE)
shifts <- data.frame(task = c("Matematik", "Biologi"), slot = c("man-08", "man-08"), stringsAsFactors = FALSE)
schedule <- create_schedule(people, shifts)
stopifnot(nrow(schedule) == 1)
stopifnot(schedule$task[[1]] == "Matematik")

people <- data.frame(name = "Anna", available_slots = "man-08, tir-08", stringsAsFactors = FALSE)
shifts <- data.frame(task = c("Matematik", "Biologi"), slot = c("man-08", "tir-08"), stringsAsFactors = FALSE)
schedule <- create_schedule(people, shifts, max_assignments = 1)
stopifnot(nrow(schedule) == 1)

parsed_people <- parse_people("Anna: man-08, tir-08\nBo: ons-08")
stopifnot(identical(parsed_people$name, c("Anna", "Bo")))

parsed_shifts <- parse_shifts("Matematik: man-08\nBiologi: tir-08")
stopifnot(identical(parsed_shifts$task, c("Matematik", "Biologi")))

cat("All scheduler tests passed\n")
