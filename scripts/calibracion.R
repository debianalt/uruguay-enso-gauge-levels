# calibración de las probabilidades anuales y supuesto de Poisson.
#
# Regla en el registro de decisiones del autor (11-nov-2026), escrita y aprobada antes de correr.
#   A. Calibración: en años hidrológicos con los doce meses válidos, años observados que alcanzan cada cota (O)
#      contra la suma de las probabilidades del modelo de frecuencia sola con la trayectoria mensual del índice
#      (E). Calibrado si |O − E|/√V < 1,645, con V = ΣP(1 − P). Dentro de la muestra.
#   B. Poisson: dispersión de Pearson de los conteos anuales, con valor p por simulación (2.000, semilla 2026).
#      Se rechaza si p < 0,05.
#   Si B rechaza: probabilidad con conteo binomial negativo, P = 1 − (1 + λS/k)^(−k), k por momentos.
#
# Uso (desde la raíz del repositorio):  Rscript scripts/calibracion.R

source("scripts/modelo.R", encoding = "UTF-8")
set.seed(2026); N_SIM <- 2000L; N_MIN <- 10L
tabla <- read.csv("tables/prob_frecuencia_sola.csv", fileEncoding = "UTF-8", check.names = FALSE)

# probabilidad de que un exceso supere z − u con la GPD estacionaria (d = 0)
S_gpd <- function(z, u, c0, xi0) {
  s <- exp(c0); t <- 1 + xi0 * (z - u) / s
  if (abs(xi0) > 1e-6) ifelse(t > 0, pmax(t, 0)^(-1 / xi0), 0) else exp(-(z - u) / s)
}

cal <- list(); disp <- list(); nb <- list()
for (sid in URUGUAY) {
  m0 <- serie[as.character(sid), "mes_inicio_hidro"]; nombre <- serie[as.character(sid), "nombre"]
  D <- datos(sid, "x6"); f <- ajustar(D)
  am <- anio_h(D$mm, m0); ae <- anio_h(D$ye, m0)
  completos <- as.integer(names(which(table(am) == 12)))
  lam_mes <- exp(f$a + f$b * D$xm) / 12
  L <- tapply(lam_mes, am, sum)[as.character(completos)]                 # tasa esperada de cada año
  xbar <- tapply(D$xm, am, mean)[as.character(completos)]
  O_n <- sapply(completos, function(a) sum(ae == a))                     # crecidas observadas en el año
  n <- length(completos)
  if (n < N_MIN) {     # posterior a la regla: con menos de diez años completos no se calcula (Pepirí Miní, Panambí)
    disp[[length(disp) + 1]] <- data.frame(serie = sid, nombre = nombre, anios = n, crecidas = sum(O_n),
      esperadas = sum(L), media = NA, varianza = NA, phi = NA, p = NA, rechaza = NA, k_nb = NA)
    next
  }

  # --- B. dispersión
  phi <- sum((O_n - L)^2 / L) / (n - 2)
  phi_sim <- replicate(N_SIM, { o <- rpois(n, L); sum((o - L)^2 / L) / (n - 2) })
  p_phi <- mean(phi_sim >= phi)
  exceso <- sum((O_n - L)^2) - sum(L)
  k <- if (exceso > 0) sum(L^2) / exceso else Inf
  disp[[length(disp) + 1]] <- data.frame(serie = sid, nombre = nombre, anios = n, crecidas = sum(O_n),
    esperadas = sum(L), media = mean(O_n), varianza = var(O_n), phi = phi, p = p_phi, rechaza = p_phi < 0.05, k_nb = k)

  # --- A. calibración, por cota
  t <- tabla[tabla$serie == sid, ]
  for (niv in unique(t$nivel)) {
    z <- t$z_m[t$nivel == niv][1]
    S <- S_gpd(z, D$u, f$c0, f$xi0)
    P <- 1 - exp(-L * S)
    obs <- sapply(completos, function(a) any(D$y[ae == a] + D$u >= z))
    O <- sum(obs); E <- sum(P); V <- sum(P * (1 - P)); zeta <- (O - E) / sqrt(V)
    ter <- cut(xbar, quantile(xbar, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("bajo", "medio", "alto"))
    fila <- data.frame(serie = sid, nombre = nombre, nivel = niv, z_m = z, anios = n, observados = O, esperados = E,
      desvio = sqrt(V), zeta = zeta, calibrado = abs(zeta) < 1.645)
    for (g in levels(ter)) {
      fila[[paste0("n_", g)]] <- sum(ter == g); fila[[paste0("x_", g)]] <- mean(xbar[ter == g])
      fila[[paste0("obs_", g)]] <- mean(obs[ter == g]); fila[[paste0("pred_", g)]] <- mean(P[ter == g])
    }
    cal[[length(cal) + 1]] <- fila

    # binomial negativa en los valores de la tabla (sólo se usa si B rechaza)
    tt <- t[t$nivel == niv, ]; lam <- exp(f$a + f$b * tt$x)
    p_nb <- if (is.finite(k)) 1 - (1 + lam * S / k)^(-k) else 1 - exp(-lam * S)
    stopifnot(max(abs((1 - exp(-lam * S)) - tt$p_frec)) < 1e-8)          # reproduce la tabla principal
    nb[[length(nb) + 1]] <- data.frame(serie = sid, nombre = nombre, nivel = niv, z_m = z, x = tt$x,
      p_frec = tt$p_frec, p_nb = p_nb, dif_puntos = 100 * (p_nb - tt$p_frec), k_nb = k)
  }
}
cal <- do.call(rbind, cal); disp <- do.call(rbind, disp); nb <- do.call(rbind, nb)
write.csv(cal, "tables/calibracion.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(disp, "tables/poisson_dispersion.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(nb, "tables/calibracion_binomial_negativa.csv", row.names = FALSE, fileEncoding = "UTF-8")

options(width = 200)
cat("--- B. Dispersión de los conteos anuales\n"); print(format(disp, digits = 3), row.names = FALSE)
cat("\n--- A. Calibración\n")
print(format(cal[, c("nombre", "nivel", "z_m", "anios", "observados", "esperados", "desvio", "zeta", "calibrado")],
             digits = 3), row.names = FALSE)
cat("\n--- A. Por tercil del índice medio del año (frecuencia observada / predicha)\n")
print(format(cal[, c("nombre", "nivel", "obs_bajo", "pred_bajo", "obs_medio", "pred_medio", "obs_alto", "pred_alto")],
             digits = 2), row.names = FALSE)
ev15 <- nb[nb$nivel == "evacuación" & abs(nb$x - 1.5) < 1e-9, ]
cat("\n--- Binomial negativa, evacuación a +1,5 °C (puntos contra la tabla principal)\n")
print(format(ev15[, c("nombre", "p_frec", "p_nb", "dif_puntos", "k_nb")], digits = 3), row.names = FALSE)
