# probabilidad anual de alcanzar las cotas con el modelo de «frecuencia sola» (d = 0).
#
# Regla en el registro de decisiones del autor (10-nov-2026, noche), escrita antes de correr. La tabla principal
# (tables/prob_condicional.csv) usa la b y la d de cada hidrómetro. Como el artículo afirma la frecuencia y no la
# altura, acá se repite con la misma tasa (a, b) y con la GPD estacionaria del mismo ajuste (c0, ξ0; d = 0).
# Mismo remuestreo por años hidrológicos, 500 réplicas, semilla 2026.
#
# Uso (desde la raíz del repositorio):  Rscript scripts/frecuencia_sola.R

source("scripts/modelo.R", encoding = "UTF-8")
N_BOOT <- 500L; set.seed(2026)
tabla <- read.csv("tables/prob_condicional.csv", fileEncoding = "UTF-8", check.names = FALSE)
q5 <- function(v) quantile(v, 0.05, na.rm = TRUE); q95 <- function(v) quantile(v, 0.95, na.rm = TRUE)

filas <- list()
for (sid in unique(tabla$serie)) {
  m0 <- serie[as.character(sid), "mes_inicio_hidro"]
  D <- datos(sid, "x6"); f <- ajustar(D); t <- tabla[tabla$serie == sid, ]
  ye <- anio_h(D$ye, m0); ym <- anio_h(D$mm, m0); anios <- sort(unique(ym))
  idx_e <- split(seq_along(ye), ye); idx_m <- split(seq_along(ym), ym)
  P <- matrix(NA_real_, nrow(t), N_BOOT)
  for (r in seq_len(N_BOOT)) {
    s <- remuestrear(idx_e, idx_m, sample(anios, length(anios), replace = TRUE))
    if (length(s$ie) < 10) next
    fr <- ajustar(list(xe = D$xe[s$ie], y = D$y[s$ie], xm = D$xm[s$im]))
    P[, r] <- prob(t$z_m, D$u, fr$a, fr$b, fr$c0, 0, fr$xi0, t$x)
  }
  p0 <- prob(t$z_m, D$u, f$a, f$b, f$c0, 0, f$xi0, t$x)
  stopifnot(max(abs(prob(t$z_m, D$u, f$a, f$b, f$c, f$d, f$xi, t$x) - t$p)) < 1e-8)   # reproduce la principal
  filas[[length(filas) + 1]] <- data.frame(serie = sid, nombre = t$nombre, nivel = t$nivel, z_m = t$z_m,
    clase = t$clase, x = t$x, p_frec = p0, p_frec_lo = apply(P, 1, q5), p_frec_hi = apply(P, 1, q95),
    p_principal = t$p, dif_puntos = 100 * (p0 - t$p), d_hidrometro = f$d)
  cat(sid, t$nombre[1], ": listo\n")
}
out <- do.call(rbind, filas)
write.csv(out, "tables/prob_frecuencia_sola.csv", row.names = FALSE, fileEncoding = "UTF-8")
ev <- out[out$nivel == "evacuación", ]
cat(sprintf("evacuación: diferencia (frecuencia sola − principal) entre %+.1f y %+.1f puntos; a +1,5 °C entre %+.1f y %+.1f\n",
            min(ev$dif_puntos), max(ev$dif_puntos), min(ev$dif_puntos[abs(ev$x - 1.5) < 1e-9]),
            max(ev$dif_puntos[abs(ev$x - 1.5) < 1e-9])))
