# PKPD_app

En Shiny-app til illustration af simple én-kompartment PK/PD-profiler.

## Faglig validering

Modellen bør valideres i små, dokumenterede trin, før resultater bruges i en
klinisk sammenhæng:

1. Aftal modellens antagelser med en klinisk farmakolog (én-kompartmentmodel,
   førsteordens-elimination, lineær kinetik og hvordan steady state defineres).
2. Udvælg anonymiserede referencecases med kendte koncentrationer og PK/PD-mål.
3. Registrér forventet resultat, enhed, kilde og accepteret tolerance for hver case.
4. Tilføj casene som automatiserede tests, og få afvigelser fagligt vurderet.
5. Versionsstyr lægemiddelparametre og dokumentér dato og kilde for ændringer.

Appen beregner fri koncentration. `protein_binding` angives internt som en
fraktion fra 0 til 1, mens brugerfladen viser procent. fAUC/MIC beregnes med
trapezmetoden på den faktisk simulerede tidsakse og skaleres til 24 timer.

> **Bemærk:** Appen er et beregnings- og visualiseringsværktøj og er ikke i sig
> selv valideret som medicinsk udstyr eller som grundlag for patientbehandling.

## Automatiske tests

Testene dækker en analytisk trapez-case, uafhængighed af tidsopløsning, simple
grænsetilfælde for tid over MIC samt afvisning af ugyldige input.

```sh
Rscript tests/testthat.R
```

Testmiljøet kræver R-pakkerne `testthat`, `tidyverse`, `rlang` og `zoo`.
