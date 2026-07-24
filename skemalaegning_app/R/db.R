# Database helpers for the dermatology scheduling app.

get_db_path <- function() {
  path <- Sys.getenv("SKEMA_DB_PATH", unset = file.path("data", "skemalaegning.sqlite"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

connect_db <- function(path = get_db_path()) {
  DBI::dbConnect(RSQLite::SQLite(), path)
}

init_db <- function(conn) {
  DBI::dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS schedule_settings (
      id INTEGER PRIMARY KEY CHECK (id = 1),
      period_start TEXT NOT NULL,
      period_length_weeks INTEGER NOT NULL,
      notes TEXT
    )
  ")

  DBI::dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS clinical_functions (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      name TEXT NOT NULL UNIQUE,
      coverage_type TEXT NOT NULL,
      required_per_day INTEGER NOT NULL DEFAULT 1,
      required_per_month INTEGER NOT NULL DEFAULT 0,
      requires_specialist INTEGER NOT NULL DEFAULT 0,
      active INTEGER NOT NULL DEFAULT 1
    )
  ")

  DBI::dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS doctors (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      name TEXT NOT NULL UNIQUE,
      employment_type TEXT NOT NULL,
      weekly_hours REAL NOT NULL,
      active INTEGER NOT NULL DEFAULT 1,
      notes TEXT
    )
  ")

  DBI::dbExecute(conn, "
    CREATE TABLE IF NOT EXISTS doctor_functions (
      doctor_id INTEGER NOT NULL,
      function_id INTEGER NOT NULL,
      PRIMARY KEY (doctor_id, function_id),
      FOREIGN KEY (doctor_id) REFERENCES doctors(id) ON DELETE CASCADE,
      FOREIGN KEY (function_id) REFERENCES clinical_functions(id) ON DELETE CASCADE
    )
  ")
}

save_schedule_settings <- function(conn, period_start, period_length_weeks, notes = "") {
  DBI::dbExecute(
    conn,
    "INSERT INTO schedule_settings (id, period_start, period_length_weeks, notes)
     VALUES (1, ?, ?, ?)
     ON CONFLICT(id) DO UPDATE SET
       period_start = excluded.period_start,
       period_length_weeks = excluded.period_length_weeks,
       notes = excluded.notes",
    params = list(as.character(period_start), as.integer(period_length_weeks), notes)
  )
}

load_schedule_settings <- function(conn) {
  DBI::dbGetQuery(conn, "SELECT period_start, period_length_weeks, notes FROM schedule_settings WHERE id = 1")
}

add_clinical_function <- function(conn, name, coverage_type, required_per_day, required_per_month, requires_specialist) {
  DBI::dbExecute(
    conn,
    "INSERT INTO clinical_functions
      (name, coverage_type, required_per_day, required_per_month, requires_specialist)
     VALUES (?, ?, ?, ?, ?)
     ON CONFLICT(name) DO UPDATE SET
       coverage_type = excluded.coverage_type,
       required_per_day = excluded.required_per_day,
       required_per_month = excluded.required_per_month,
       requires_specialist = excluded.requires_specialist,
       active = 1",
    params = list(
      name,
      coverage_type,
      as.integer(required_per_day),
      as.integer(required_per_month),
      as.integer(requires_specialist)
    )
  )
}

load_clinical_functions <- function(conn, active_only = TRUE) {
  sql <- "SELECT id, name, coverage_type, required_per_day, required_per_month, requires_specialist, active FROM clinical_functions"
  if (active_only) {
    sql <- paste(sql, "WHERE active = 1")
  }
  DBI::dbGetQuery(conn, paste(sql, "ORDER BY name"))
}

add_doctor <- function(conn, name, employment_type, weekly_hours, function_ids = integer(), notes = "") {
  DBI::dbWithTransaction(conn, {
    DBI::dbExecute(
      conn,
      "INSERT INTO doctors (name, employment_type, weekly_hours, notes)
       VALUES (?, ?, ?, ?)
       ON CONFLICT(name) DO UPDATE SET
         employment_type = excluded.employment_type,
         weekly_hours = excluded.weekly_hours,
         notes = excluded.notes,
         active = 1",
      params = list(name, employment_type, as.numeric(weekly_hours), notes)
    )

    doctor_id <- DBI::dbGetQuery(conn, "SELECT id FROM doctors WHERE name = ?", params = list(name))$id[[1]]
    DBI::dbExecute(conn, "DELETE FROM doctor_functions WHERE doctor_id = ?", params = list(doctor_id))

    for (function_id in as.integer(function_ids)) {
      DBI::dbExecute(
        conn,
        "INSERT OR IGNORE INTO doctor_functions (doctor_id, function_id) VALUES (?, ?)",
        params = list(doctor_id, function_id)
      )
    }
  })
}

load_doctors <- function(conn, active_only = TRUE) {
  sql <- "
    SELECT
      doctors.id,
      doctors.name,
      doctors.employment_type,
      doctors.weekly_hours,
      doctors.active,
      doctors.notes,
      COALESCE(GROUP_CONCAT(clinical_functions.name, ', '), '') AS functions
    FROM doctors
    LEFT JOIN doctor_functions ON doctor_functions.doctor_id = doctors.id
    LEFT JOIN clinical_functions ON clinical_functions.id = doctor_functions.function_id
  "
  if (active_only) {
    sql <- paste(sql, "WHERE doctors.active = 1")
  }
  sql <- paste(sql, "GROUP BY doctors.id ORDER BY doctors.name")
  DBI::dbGetQuery(conn, sql)
}
