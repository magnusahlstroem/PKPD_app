from pathlib import Path


def test_database_schema_contains_core_tables():
    db_file = Path("skemalaegning_app/R/db.R").read_text()

    for table in [
        "schedule_settings",
        "clinical_functions",
        "doctors",
        "doctor_functions",
    ]:
        assert table in db_file


def test_app_sources_expected_modules():
    app_file = Path("skemalaegning_app/app.R").read_text()

    for module in [
        "settings_module.R",
        "functions_module.R",
        "doctors_module.R",
    ]:
        assert module in app_file
