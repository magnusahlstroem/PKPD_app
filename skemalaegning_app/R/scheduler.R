# Simple scheduling logic for the Shiny app.

create_schedule <- function(people, shifts, max_assignments = Inf) {
  if (nrow(people) == 0 || nrow(shifts) == 0) {
    return(empty_schedule())
  }

  assignment_counts <- stats::setNames(rep(0, nrow(people)), people$name)
  booked_slots <- stats::setNames(vector("list", nrow(people)), people$name)
  booked_slots <- lapply(booked_slots, function(x) character())

  assignments <- list()

  for (shift_index in seq_len(nrow(shifts))) {
    shift <- shifts[shift_index, , drop = FALSE]
    candidates <- people[
      vapply(
        seq_len(nrow(people)),
        function(person_index) {
          person <- people[person_index, , drop = FALSE]
          available_slots <- split_slots(person$available_slots)
          person_name <- person$name

          shift$slot %in% available_slots &&
            !(shift$slot %in% booked_slots[[person_name]]) &&
            assignment_counts[[person_name]] < max_assignments
        },
        logical(1)
      ),
      ,
      drop = FALSE
    ]

    if (nrow(candidates) == 0) {
      next
    }

    candidate_counts <- assignment_counts[candidates$name]
    selected <- candidates[order(candidate_counts, candidates$name), , drop = FALSE][1, ]
    selected_name <- selected$name

    assignments[[length(assignments) + 1]] <- data.frame(
      task = shift$task,
      slot = shift$slot,
      person = selected_name,
      stringsAsFactors = FALSE
    )

    assignment_counts[[selected_name]] <- assignment_counts[[selected_name]] + 1
    booked_slots[[selected_name]] <- c(booked_slots[[selected_name]], shift$slot)
  }

  if (length(assignments) == 0) {
    return(empty_schedule())
  }

  do.call(rbind, assignments)
}

parse_people <- function(text) {
  lines <- clean_lines(text)
  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      return(NULL)
    }

    data.frame(
      name = trimws(parts[[1]]),
      available_slots = trimws(paste(parts[-1], collapse = ":")),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(name = character(), available_slots = character()))
  }

  do.call(rbind, rows)
}

parse_shifts <- function(text) {
  lines <- clean_lines(text)
  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      return(NULL)
    }

    data.frame(
      task = trimws(parts[[1]]),
      slot = trimws(paste(parts[-1], collapse = ":")),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame(task = character(), slot = character()))
  }

  do.call(rbind, rows)
}

find_unassigned_shifts <- function(shifts, schedule) {
  if (nrow(shifts) == 0) {
    return(shifts)
  }

  if (nrow(schedule) == 0) {
    return(shifts)
  }

  shifts[!(shifts$task %in% schedule$task), , drop = FALSE]
}

split_slots <- function(slots) {
  trimws(unlist(strsplit(slots, ",", fixed = TRUE)))
}

clean_lines <- function(text) {
  lines <- trimws(unlist(strsplit(text, "\n", fixed = TRUE)))
  lines[nzchar(lines)]
}

empty_schedule <- function() {
  data.frame(task = character(), slot = character(), person = character())
}
