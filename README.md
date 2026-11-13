# El Niño and the annual probability of reaching alert and evacuation levels on the middle Uruguay River

Code and data for the article of the same title (submitted to *Atmósfera*). For seven gauges of the Uruguay River in
Misiones and Corrientes (Argentina), the analysis fits a Poisson process for flood frequency, with a rate that depends
on the Niño 3.4 index, and a generalised Pareto distribution for flood height. It then gives the annual probability of
reaching each official alert and evacuation level for several index values.

Every number in the article comes from a CSV file in `tables/`. Comments in the scripts are in Spanish.

## Requirements

R 4.5.1 with base R only for the fit. The figures also need ggplot2 4.0.3, sf 1.0.24, ggrepel 0.9.8 and patchwork 1.3.2.

## Running

Run the scripts from the repository root, in this order. Each one writes to `tables/` (and `validation/`, `figures/`).

```
Rscript scripts/frecuencia.R            # fit per gauge and pooled, bootstrap (seed 2026), full model -> Table S1
Rscript scripts/frecuencia_sola.R       # probabilities with d = 0 -> Table II, Table S7
Rscript scripts/robustez.R              # flood definition, dams, trend, constancy of the relationship
Rscript scripts/robustez_indice.R       # index windows, ERSST, lagged index, 90% profile intervals
Rscript scripts/calibracion.R           # calibration and Poisson dispersion (Tables S2, S3)
Rscript scripts/calibracion_garruchos.R # diagnosis of the Garruchos alert level (Table S4)
Rscript scripts/figuras.R               # Figures 1-4 and the graphical abstract
```

`scripts/modelo.R` holds the data loading and model functions and is sourced by the others; `scripts/tema.R` holds the
figure theme. `Rscript scripts/frecuencia.R --sintetica` checks the estimator on simulated data, and
`--semilla 1` / `--semilla 2` give the control runs with other bootstrap seeds (`validation/*_semilla*.csv`).

Figure 1 also needs three Natural Earth 1:10m layers, which are not redistributed here. Download
`ne_10m_admin_0_countries`, `ne_10m_admin_1_states_provinces` and `ne_10m_rivers_lake_centerlines` from
<https://www.naturalearthdata.com/downloads/10m-cultural-vectors/> (and the physical vectors page for rivers) and unzip
each into `data/geo/<layer>/`.

## Data

The files come from a regional database of the author that also holds six Paraná and Iguazú gauges (series 8–14).
The scripts use only the seven Uruguay gauges (series 60–68); the daily readings are given for those seven only.

| File | Content |
|---|---|
| `data/crecidas_diario.csv` | Daily stage readings (m) of the seven Uruguay gauges, from the start of each record to 7 November 2026 |
| `data/crecidas_mes.csv` | Days with readings per gauge and month, and whether the month is valid (readings on at least 80% of days) |
| `data/crecidas_enso.csv` | Monthly Niño 3.4 anomaly: HadISST, OISST v2.1 and ERSST v6, the series used (`nino34`) and its six-month mean (`x6`) |
| `data/crecidas_evento.csv` | Floods identified by the threshold rule (98th percentile, 10-day separation), with peak date and height |
| `data/crecidas_max_anual.csv` | Annual maximum per gauge and hydrological year |
| `data/crecidas_homogeneidad.csv` | Pettitt tests on annual maxima and on differences between neighbouring gauges (jumps used in Section 2.1) |
| `data/crecidas_serie.csv` | Gauges, threshold *u*, and official alert and evacuation levels with their source |
| `data/crecidas_ajuste.csv`, `data/crecidas_prob_condicional.csv` | Output of an earlier Python implementation, used only for the cross-validation in `validation/crossval_report.csv` |
| `data/geo/hidrometros.geojson` | Gauge locations and levels |
| `data/geo/rio_uruguay_ign.geojson` | Course of the Uruguay River in Argentine territory |

Sources: daily stage readings of the Prefectura Naval Argentina scales, distributed by the Information and Hydrological
Alert System of the Instituto Nacional del Agua (INA), downloaded on 9 November 2026; alert and evacuation levels from
the INA and from the hydrological monitor of the Government of Misiones; HadISST (Met Office Hadley Centre), OISST v2.1
and ERSST v6 (NOAA); river course from the Instituto Geográfico Nacional (Argentina), perennial watercourses layer.
Users of the stage readings should cite the INA.

## Licence

Code: MIT (`LICENSE`). Derived tables and figures: CC BY 4.0. The input series remain subject to the terms of their
original providers.
