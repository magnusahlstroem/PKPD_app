# Automatisk skemalægning

Koden ligger i den nye mappe `skemalaegning_app/` i roden af dette repository. Projektet er lavet som en selvstændig R/Shiny-app, så det passer bedre til en R-baseret arbejdsgang.

## Hvor ligger koden?

```text
skemalaegning_app/
├── app.R                  # Shiny-appens brugerflade og server
├── R/
│   └── scheduler.R        # Selve skemalægningslogikken
└── tests/
    └── test_scheduler.R   # Simple tests af skemalægningslogikken
```

## Kør appen

Åbn R/RStudio i repoets rodmappe og kør:

```r
shiny::runApp("skemalaegning_app")
```

Hvis du ikke har Shiny installeret endnu, kan du først køre:

```r
install.packages("shiny")
```

## Kør tests

Fra repoets rodmappe kan du køre:

```bash
Rscript skemalaegning_app/tests/test_scheduler.R
```

## Hvad kan første version?

Appen kan lave et simpelt automatisk skema ud fra:

- personer, der kan planlægges
- opgaver, vagter eller lektioner, der skal dækkes
- tidsrum, hvor personer er tilgængelige
- en regel om at undgå dobbeltbooking i samme tidsrum
- en simpel balancering, hvor personen med færrest tildelinger vælges først

Det er bevidst en lille startversion, som er nem at forstå og bygge videre på i R. Senere kan vi udvide med flere regler, prioriteringer, import fra Excel/CSV og mere avanceret optimering.
