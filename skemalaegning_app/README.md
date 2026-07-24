# Dermatologisk skemalægning

Koden ligger i mappen `skemalaegning_app/` i roden af dette repository. Projektet er en selvstændig R/Shiny-app til at opbygge grunddata for skemalægning på en dermatologisk afdeling.

## Formål

Appen skal på sigt kunne hjælpe skemalæggeren med at lave lægeskemaer ud fra:

- normperiode, startdato og længde
- funktioner der skal dækkes dagligt eller månedligt
- læger på afdelingen
- timetal og ansættelsestype, fx speciallæge eller uddannelseslæge
- hvilke funktioner den enkelte læge kan varetage
- gemte grunddata i en lokal database

## Hvor ligger koden?

```text
skemalaegning_app/
├── app.R                         # Shiny-appens samlede indgang
├── DESCRIPTION                   # R-afhængigheder: shiny, DBI og RSQLite
├── R/
│   ├── db.R                      # SQLite database og gem/load-funktioner
│   ├── scheduler.R               # Første simple skemalægningsmotor
│   └── modules/
│       ├── doctors_module.R      # Modul til læger og deres kompetencer
│       ├── functions_module.R    # Modul til funktioner der skal dækkes
│       └── settings_module.R     # Modul til normperiode og grunddata
├── data/                         # Lokal database oprettes her ved kørsel
└── tests/
    └── test_scheduler.R          # Simple tests af skemalægningslogikken
```

## Kør appen

Åbn R/RStudio i repoets rodmappe og kør:

```r
shiny::runApp("skemalaegning_app")
```

Hvis du mangler pakkerne, kan du først køre:

```r
install.packages(c("shiny", "DBI", "RSQLite"))
```

## Database

Appen bruger SQLite via `DBI` og `RSQLite`. Som standard gemmes databasen lokalt i:

```text
skemalaegning_app/data/skemalaegning.sqlite
```

Du kan ændre placeringen ved at sætte miljøvariablen `SKEMA_DB_PATH`, før appen startes.

## Moduler i første version

### Normperiode

Her kan skemalæggeren gemme startdato, længde i uger og noter til normperioden.

### Funktioner

Her kan skemalæggeren oprette funktioner som fx ambulatorium, stuegang, akut, supervision eller operationsdag. For hver funktion kan man angive, om den dækkes dagligt/månedligt, antal pr. dag/måned og om den kræver speciallæge.

### Læger

Her kan skemalæggeren oprette læger, angive ansættelsestype, timetal og vælge hvilke funktioner lægen kan tage.

### Første skemaforslag

Dette faneblad indeholder fortsat en simpel test af skemamotoren. Næste naturlige skridt er at koble den direkte til de gemte læger, funktioner og normperioden fra databasen.

## Kør tests

Fra repoets rodmappe kan du køre:

```bash
Rscript skemalaegning_app/tests/test_scheduler.R
```
